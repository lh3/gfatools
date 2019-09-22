#include <stdio.h>
#include <assert.h>
#include "gfa-priv.h"
#include "khash.h"
KHASH_INIT(gt, uint32_t, char, 0, __ac_Wang_hash, kh_int_hash_equal)

typedef struct {
	int32_t j, i, a;
	double sc;
} gt_sc_t;

#define GT_MAX_SC 3

typedef struct {
	int32_t n;
	float vsc;
	gt_sc_t s[GT_MAX_SC]; // top 2 scores
} gt_max_t;

typedef struct {
	int32_t is_arc;
	int32_t x;
	float w;
} gt_elem_t;

typedef struct {
	int32_t l, flt;
	gt_elem_t *w;
	double s;
} gt_walk_t;

typedef struct {
	int32_t n_al;
	int32_t al[2];
	double s_alt, s_het;
} gt_call_t;

static inline float gt_get_dc(const gfa_aux_t *aux)
{
	const uint8_t *dc;
	dc = gfa_aux_get(aux->l_aux, aux->aux, "dc");
	return dc && dc[0] == 'f'? *(float*)(dc + 1) : 0.0f;
}

static void gt_all_weights(const gfa_t *g, const gfa_sub_t *sub, float *wv, float *wa)
{
	int32_t j, k;
	for (j = 0; j < sub->n_v; ++j) {
		const gfa_subv_t *t = &sub->v[j];
		wv[j] = gt_get_dc(&g->seg[sub->v[j].v>>1].aux);
		for (k = 0; k < t->n; ++k)
			wa[t->off + k] = gt_get_dc(&g->link_aux[g->arc[(uint32_t)sub->a[t->off + k]].link_id]);
	}
}

static int32_t gt_get_ref_walk(const gfa_t *g, const gfa_sub_t *sub, int32_t jst, int32_t jen, gt_elem_t *walk)
{
	int32_t i, j, k;
	j = jst, k = 0;
	while (1) {
		const gfa_subv_t *t = &sub->v[j];
		uint64_t a;
		for (i = 0; i < t->n; ++i) {
			a = sub->a[t->off + i];
			if (g->arc[(uint32_t)a].rank == 0)
				break;
		}
		assert(i < t->n);
		j = a>>32;
		walk[k].is_arc = 1, walk[k++].x = t->off + i;
		if (j == jen) break;
		walk[k].is_arc = 0, walk[k++].x = j;
	}
	return k;
}

static double gt_cal_weight(const gfa_t *g, const gfa_sub_t *sub, const float *wv, const float *wa, int32_t len, gt_elem_t *walk)
{
	double s = 0.0;
	int32_t j;
	for (j = 0; j < len; ++j) {
		gt_elem_t *w = &walk[j];
		w->w = w->is_arc? wa[w->x] : wv[w->x];
		s += w->w;
	}
	return s;
}

static int32_t gt_best_walks(const gfa_t *g, const gfa_sub_t *sub, const float *wv, const float *wa, int32_t jst, int32_t jen, gt_walk_t *walk)
{
	int32_t j, k, n_walk;
	gt_max_t *sc;

	// fill sc[]
	GFA_CALLOC(sc, jen - jst + 1);
	sc[jen - jst].n = 1;
	sc[jen - jst].s[0].j = sc[jen - jst].s[0].i = -1;
	for (j = jen - 1; j >= jst; --j) {
		const gfa_subv_t *t = &sub->v[j];
		gt_max_t *s0 = &sc[j - jst];
		s0->vsc = j == jst? 0.0f : wv[j];
		for (k = 0; k < t->n; ++k) { // iterate over neighbors
			uint64_t a = sub->a[t->off + k];
			uint32_t jp = (uint32_t)(a>>32);
			int32_t i;
			gt_max_t *s1;
			float dc;
			assert(t->off + k < sub->n_a);
			if (jp <= j || jp > jen) continue; // ignore cycles or arcs outside this subgraph
			dc = wa[t->off + k];
			s1 = &sc[jp - jst];
			for (i = 0; i < s1->n; ++i) { // iterate over path ends
				double score;
				int32_t ins = -1;
				score = s1->s[i].sc + dc + s0->vsc;
				if (s0->n == GT_MAX_SC) {
					int32_t x, y;
					for (x = 0; x < s0->n; ++x)
						if (s0->s[x].sc < score)
							break;
					if (x < s0->n) { // make room for the new partial path
						for (y = s0->n - 1; y > x; --y)
							s0->s[y] = s0->s[y - 1];
						ins = x;
					}
				} else ins = s0->n;
				if (ins >= 0) {
					s0->s[ins].j = jp;
					s0->s[ins].i = i;
					s0->s[ins].a = t->off + k;
					s0->s[ins].sc = score;
					if (s0->n < GT_MAX_SC) ++s0->n;
				}
			}
		}
	}

	// allocate path[]
	n_walk = sc->n + 1;
	for (k = 0; k < n_walk; ++k) // k==sc->n for the reference path
		GFA_MALLOC(walk[k].w, (jen - jst + 1) * 2); // this is over-allocating, but it should not be an issue

	// backtrack
	for (k = 0; k < sc->n; ++k) {
		int32_t l = 0, j = jst, i = k;
		gt_elem_t *w = walk[k + 1].w;
		while (1) {
			gt_sc_t *s = &sc[j - jst].s[i];
			j = s->j, i = s->i;
			w[l].is_arc = 1, w[l].x = s->a, w[l].w = 0.0, ++l;
			if (j < 0 || j == jen) break;
			w[l].is_arc = 0, w[l].x = j, w[l].w = 0.0, ++l;
		}
		walk[k+1].l = l;
	}
	free(sc);

	// the reference walk
	walk[0].l = gt_get_ref_walk(g, sub, jst, jen, walk[0].w);
	return n_walk;
}

static double gt_relative(int32_t l0, const gt_elem_t *w0, int32_t l1, const gt_elem_t *w1)
{
	int32_t i;
	khash_t(gt) *h;
	double s = 0.0;
	h = kh_init(gt);
	kh_resize(gt, h, l0<<1);
	for (i = 0; i < l0; ++i) {
		uint32_t x = w0[i].is_arc? 1U<<31 | w0[i].x : w0[i].x;
		int absent;
		kh_put(gt, h, x, &absent);
	}
	for (i = 0; i < l1; ++i) {
		uint32_t x = w1[i].is_arc? 1U<<31 | w1[i].x : w1[i].x;
		if (kh_get(gt, h, x) == kh_end(h))
			s += w1[i].w;
	}
	kh_destroy(gt, h);
	return s;
}

static int32_t gt_filter_walk(int32_t n_walk, gt_walk_t *walk)
{
	int32_t j, k;
	for (k = 1, j = 1; k < n_walk; ++k) {
		double s;
		s = gt_relative(walk[0].l, walk[0].w, walk[k].l, walk[k].w);
		if ((float)s + 1.0f > 1.0f) walk[j++] = walk[k];
		else free(walk[k].w);
	}
	return j;
}

static void gt_genotype(int32_t n_walk, const gt_walk_t *walk, gt_call_t *c)
{
	int32_t j, k, max_j = -1, max_k = -1;
	double max = 0.0;
	assert(n_walk >= 1);
	for (j = 0; j < n_walk; ++j) { // choose the pair of walks that maximize the weight
		int32_t i;
		double sj = 0.0;
		khash_t(gt) *h;
		const gt_walk_t *wj = &walk[j];
		h = kh_init(gt);
		kh_resize(gt, h, wj->l);
		for (i = 0; i < wj->l; ++i) {
			uint32_t x = wj->w[i].is_arc? 1U<<31|wj->w[i].x : wj->w[i].x;
			int absent;
			kh_put(gt, h, x, &absent);
			sj += wj->w[i].w;
		}
		for (k = j; k < n_walk; ++k) {
			const gt_walk_t *wk = &walk[k];
			double sk = sj;
			for (i = 0; i < wk->l; ++i) {
				uint32_t x = wk->w[i].is_arc? 1U<<31|wk->w[i].x : wk->w[i].x;
				if (kh_get(gt, h, x) == kh_end(h))
					sk += wk->w[i].w;
			}
			if (max < sk) max = sk, max_j = j, max_k = k;
		}
		kh_destroy(gt, h);
	}
	if (max_j != max_k) {
		double skj, sjk;
		sjk = gt_relative(walk[max_j].l, walk[max_j].w, walk[max_k].l, walk[max_k].w);
		skj = gt_relative(walk[max_k].l, walk[max_k].w, walk[max_j].l, walk[max_j].w);
		if (sjk < skj) { // max_j is the better allele
			if ((float)sjk + 1.0f == 1.0f) c->n_al = 1, c->al[0] = max_j, c->s_het = 0.0;
			else c->n_al = 2, c->al[0] = max_j, c->al[1] = max_k, c->s_het = sjk;
		} else {
			if ((float)skj + 1.0f == 1.0f) c->n_al = 1, c->al[0] = max_k, c->s_het = 0.0;
			else c->n_al = 2, c->al[0] = max_k, c->al[1] = max_j, c->s_het = skj;
		}
	} else c->n_al = 1, c->al[0] = max_j, c->s_het = 0.0;
}

static void gt_call(int32_t n_walk, gt_walk_t *walk, gt_call_t *c)
{
	double max_s;
	int32_t k, max_k;
	assert(n_walk > 0);
	memset(c, 0, sizeof(gt_call_t));
	c->n_al = 1, c->al[0] = 0, c->s_alt = c->s_het = 0.0; // reference allele only
	if (n_walk == 1) return;
	max_s = -1.0, max_k = -1;
	for (k = 1; k < n_walk; ++k) {
		double s;
		s = gt_relative(walk[0].l, walk[0].w, walk[k].l, walk[k].w);
		if (max_s < s) max_s = s, max_k = k;
	}
	c->s_alt = max_s;
	gt_genotype(n_walk, walk, c);
}

static int32_t gt_walk_compact(int32_t n_walk, gt_walk_t *walk, gt_call_t *c)
{
	int32_t n, i, k, n2o[GT_MAX_SC+1], o2n[GT_MAX_SC+1];
	walk[0].flt = 0;
	for (k = 1; k < n_walk; ++k) {
		walk[k].flt = 1;
		for (i = 0; i < c->n_al; ++i)
			if (k == c->al[i]) break;
		walk[k].flt = i == c->n_al? 1 : 0;
	}
	for (k = 0, n = 0; k < n_walk; ++k) {
		if (!walk[k].flt)
			walk[n] = walk[k], n2o[n++] = k;
		else free(walk[k].w);
		o2n[k] = -1;
	}
	for (k = 0; k < n; ++k) o2n[n2o[k]] = k;
	for (k = 0; k < c->n_al; ++k) c->al[k] = o2n[c->al[k]];
	return n;
}

static void gt_print(const gfa_t *g, const gfa_sub_t *sub, int32_t jst, int32_t jen, int32_t n_walk, const gt_walk_t *walk, const gt_call_t *call, int32_t is_path)
{
	FILE *fp = stdout;
	int32_t i, k;
	const gfa_seg_t *seg_st, *seg_en;
	seg_st = &g->seg[sub->v[jst].v>>1];
	seg_en = &g->seg[sub->v[jen].v>>1];
	fprintf(fp, "%s\t%d\t%d", g->sseq[seg_st->snid].name, seg_st->soff + seg_st->len, seg_en->soff);
	if (is_path) fprintf(fp, "\t%c%s\t%c%s", "><"[sub->v[jst].v&1], g->seg[sub->v[jst].v>>1].name, "><"[sub->v[jen].v&1], g->seg[sub->v[jen].v>>1].name);
	fprintf(fp, "\t%.2f\t%.2f", call->s_alt, call->s_het);
	if (n_walk == 1) { // only the reference path is present
		fputs("\t0/0\n", stdout);
		return;
	}
	if (call->n_al == 1) fprintf(fp, "\t%d/%d", call->al[0], call->al[0]);
	else fprintf(fp, "\t%d/%d", call->al[0], call->al[1]);
	fprintf(fp, "\t%d", n_walk);
	for (k = 0; k < n_walk; ++k) {
		const gt_walk_t *w = &walk[k];
		int32_t c = 0;
		fputc('\t', fp);
		if (is_path) {
			for (i = 0; i < w->l; ++i) {
				if (!w->w[i].is_arc) {
					uint32_t v = sub->v[w->w[i].x].v;
					fprintf(fp, "%c%s", "><"[v&1], g->seg[v>>1].name);
					++c;
				}
			}
		} else {
			int32_t l;
			char *buf;
			for (i = 0, l = 0; i < w->l; ++i)
				if (!w->w[i].is_arc)
					l += g->seg[sub->v[w->w[i].x].v>>1].len;
			GFA_CALLOC(buf, l + 1);
			for (i = 0, l = 0; i < w->l; ++i) {
				if (!w->w[i].is_arc) {
					uint32_t v = sub->v[w->w[i].x].v;
					const gfa_seg_t *seg = &g->seg[v>>1];
					if (v&1) {
						for (k = seg->len - 1; k >= 0; --k)
							buf[l++] = gfa_comp_table[(uint8_t)seg->seq[k]];
					} else {
						memcpy(&buf[l], seg->seq, seg->len);
						l += seg->len;
					}
					++c;
				}
			}
			fputs(buf, fp);
			free(buf);
		}
		if (c == 0) printf("*");
	}
	fputc('\n', fp);
}

static void gfa_gt_simple_interval(const gfa_t *g, const gfa_sub_t *sub, const float *wv, const float *wa, int32_t jst, int32_t jen, float min_dc, int32_t is_path)
{
	int32_t k, n_walk;
	gt_walk_t walk[GT_MAX_SC+1];
	gt_call_t call;
	const gfa_seg_t *seg_st, *seg_en;

	memset(walk, 0, sizeof(gt_walk_t) * (GT_MAX_SC + 1));
	seg_st = &g->seg[sub->v[jst].v>>1];
	seg_en = &g->seg[sub->v[jen].v>>1];
	assert(seg_st->rank == 0 && seg_en->rank == 0 && seg_st->snid == seg_en->snid);
	//fprintf(stderr, "XX\t%s\t%s\n", g->seg[sub->v[jst].v>>1].name, g->seg[sub->v[jen].v>>1].name);

	n_walk = gt_best_walks(g, sub, wv, wa, jst, jen, walk);
	for (k = 0; k < n_walk; ++k)
		walk[k].s = gt_cal_weight(g, sub, wv, wa, walk[k].l, walk[k].w);
	n_walk = gt_filter_walk(n_walk, walk);
	gt_call(n_walk, walk, &call);
	n_walk = gt_walk_compact(n_walk, walk, &call);
	gt_print(g, sub, jst, jen, n_walk, walk, &call, is_path);
	for (k = 0; k < n_walk; ++k)
		free(walk[k].w);
}

void gfa_gt_simple_print(const gfa_t *g, float min_dc, int32_t is_path) // FIXME: doesn't work with translocations
{
	uint32_t i, *vs, *vmin;
	GFA_MALLOC(vs, g->n_sseq);
	GFA_MALLOC(vmin, g->n_sseq);
	for (i = 0; i < g->n_sseq; ++i)
		vs[i] = (uint32_t)-1, vmin[i] = UINT32_MAX;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		if (s->rank != 0 || s->snid < 0) continue;
		if ((uint32_t)s->soff < vmin[s->snid])
			vmin[s->snid] = s->soff, vs[s->snid] = i<<1;
	}
	free(vmin);
	for (i = 0; i < g->n_sseq; ++i) {
		gfa_sub_t *sub;
		int32_t j, jst, max_a;
		float *wa, *wv;
		if (vs[i] == (uint32_t)-1) continue;
		sub = gfa_sub_from(0, g, vs[i], 0);
		GFA_MALLOC(wv, sub->n_v);
		GFA_MALLOC(wa, sub->n_a);
		gt_all_weights(g, sub, wv, wa);
		for (j = 0, jst = 0, max_a = -1; j < sub->n_v; ++j) {
			gfa_subv_t *t = &sub->v[j];
			int32_t k;
			if (j == max_a) {
				const gfa_seg_t *sst = &g->seg[sub->v[jst].v>>1];
				const gfa_seg_t *sen = &g->seg[t->v>>1];
				if (sst->snid == i && sen->snid == i)
					gfa_gt_simple_interval(g, sub, wv, wa, jst, j, min_dc, is_path);
				max_a = -1, jst = j;
			}
			for (k = 0; k < t->n; ++k)
				if ((int32_t)(sub->a[t->off + k]>>32) > max_a)
					max_a = sub->a[t->off + k]>>32;
		}
		free(wa); free(wv);
		gfa_sub_destroy(sub);
	}
	free(vs);
}
