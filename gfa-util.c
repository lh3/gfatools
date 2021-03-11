#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include "gfa-priv.h"
#include "kvec.h"
#include "ksort.h"
#include "kdq.h"
KDQ_INIT(uint64_t)

void gfa_sub(gfa_t *g, int n, char *const* seg, int step)
{
	int32_t i;
	int8_t *flag;
	kdq_t(uint64_t) *q;
	if (n == 0) return;
	q = kdq_init(uint64_t);
	GFA_CALLOC(flag, g->n_seg * 2);
	for (i = 0; i < n; ++i) {
		int32_t s;
		s = gfa_name2id(g, seg[i]);
		if (s >= 0) {
			kdq_push(uint64_t, q, ((uint64_t)s<<1|0)<<32);
			kdq_push(uint64_t, q, ((uint64_t)s<<1|1)<<32);
		}
	}
	for (i = 0; i < g->n_seg; ++i) // mark all segments to be deleted
		g->seg[i].del = 1;
	while (kdq_size(q) > 0) {
		uint64_t x = *kdq_shift(uint64_t, q);
		uint32_t v = x>>32;
		int r = (int32_t)x;
		if (flag[v]) continue; // already visited
		flag[v] = 1;
		g->seg[v>>1].del = 0;
		if (r < step) {
			uint32_t nv = gfa_arc_n(g, v);
			gfa_arc_t *av = gfa_arc_a(g, v);
			for (i = 0; i < nv; ++i) {
				if (flag[av[i].w] == 0)
					kdq_push(uint64_t, q, (uint64_t)av[i].w<<32 | (r + 1));
				if (flag[av[i].w^1] == 0)
					kdq_push(uint64_t, q, (uint64_t)(av[i].w^1)<<32 | (r + 1));
			}
		}
	}
	kdq_destroy(uint64_t, q);
	free(flag);
	gfa_arc_rm(g);
}

static uint64_t find_join(const gfa_t *g, uint32_t v)
{
	gfa_seg_t *t, *s = &g->seg[v>>1];
	int32_t i, nv, n_low, n_r;
	uint32_t w;
	gfa_arc_t *av;
	if (s->rank == 0) return (uint64_t)-1;
	nv = gfa_arc_n(g, v);
	av = gfa_arc_a(g, v);
	for (i = 0, n_low = n_r = 0, w = 0; i < nv; ++i) {
		gfa_arc_t *q = &av[i];
		if (q->rank >= 0 && q->rank == s->rank) {
			++n_r, w = q->w;
		} else {
			t = &g->seg[q->w>>1];
			if (t->rank >= 0 && t->rank < s->rank)
				++n_low, w = q->w;
		}
	}
	if (n_r != 1 && gfa_verbose >= 2)
		fprintf(stderr, "[W] failed to find the associated arc for vertex %c%s[%d]: %d,%d\n", "><"[v&1], g->seg[v>>1].name, v, n_r, n_low);
	if (n_r != 1 && n_low != 1) return (uint64_t)-1;
	t = &g->seg[w>>1];
	return (uint64_t)t->snid<<32 | (uint32_t)(w&1? t->soff + t->len : t->soff) << 1 | (w&1);
}

gfa_sfa_t *gfa_gfa2sfa(const gfa_t *g, int32_t *n_sfa_, int32_t write_seq)
{
	int32_t i, j, k, *scnt, *soff, n_sfa;
	gfa_sfa_t *sfa = 0;
	uint64_t *a;

	*n_sfa_ = 0;
	if (g->n_sseq == 0) return 0;

	// precount
	GFA_CALLOC(scnt, g->n_sseq);
	for (i = 0; i < g->n_seg; ++i)
		if (g->seg[i].snid >= 0)
			++scnt[g->seg[i].snid];
	GFA_MALLOC(soff, g->n_sseq + 1);
	for (soff[0] = 0, i = 1; i <= g->n_sseq; ++i)
		soff[i] = soff[i - 1] + scnt[i - 1];

	// fill a[]
	GFA_BZERO(scnt, g->n_sseq);
	GFA_MALLOC(a, g->n_seg);
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		if (s->snid < 0) continue;
		a[soff[s->snid] + scnt[s->snid]] = (uint64_t)s->soff<<32 | i;
		++scnt[s->snid];
	}
	for (i = 0; i < g->n_sseq; ++i)
		if (scnt[i] > 1)
			radix_sort_gfa64(&a[soff[i]], &a[soff[i+1]]);
	free(scnt);

	// check
	n_sfa = g->n_sseq;
	for (i = 0; i < g->n_sseq; ++i) {
		const gfa_seg_t *s;
		if (soff[i] == soff[i+1]) --n_sfa;
		if (soff[i] == soff[i+1]) continue;
		s = &g->seg[(int32_t)a[soff[i]]];
		if (s->rank == 0 && s->soff != 0) {
			if (gfa_verbose >= 2)
				fprintf(stderr, "[W] rank-0 stable sequence \"%s\" not started with 0\n", g->sseq[s->snid].name);
			goto end_check;
		}
		for (j = soff[i] + 1; j < soff[i+1]; ++j) {
			const gfa_seg_t *s = &g->seg[(int32_t)a[j-1]];
			const gfa_seg_t *t = &g->seg[(int32_t)a[j]];
			if (s->soff + s->len > t->soff) {
				if (gfa_verbose >= 2)
					fprintf(stderr, "[W] overlap on stable sequence \"%s\"\n", g->sseq[s->snid].name);
				goto end_check;
			}
			if (s->rank == 0 && s->soff + s->len != t->soff) {
				if (gfa_verbose >= 2)
					fprintf(stderr, "[W] rank-0 stable sequence \"%s\" is not contiguous\n", g->sseq[s->snid].name);
				goto end_check;
			}
			if (s->rank != t->rank) {
				if (gfa_verbose >= 2)
					fprintf(stderr, "[W] stable sequence \"%s\" associated with different ranks\n", g->sseq[s->snid].name);
				goto end_check;
			}
			if (s->soff + s->len == t->soff) {
				int32_t k, nv;
				const gfa_arc_t *av;
				nv = gfa_arc_n(g, (uint32_t)a[j-1]<<1);
				av = gfa_arc_a(g, (uint32_t)a[j-1]<<1);
				for (k = 0; k < nv; ++k)
					if (av[k].w == (uint32_t)a[j]<<1)
						break;
				if (s->rank == 0 && k == nv) {
					if (gfa_verbose >= 2)
						fprintf(stderr, "[W] nearby segments on rank-0 stable sequence \"%s\" are not connected\n", g->sseq[s->snid].name);
					goto end_check;
				}
				if (k == nv) ++n_sfa;
			} else ++n_sfa;
		}
	}

	// fill sfa[]
	*n_sfa_ = n_sfa;
	GFA_CALLOC(sfa, n_sfa);
	for (i = 0, k = 0; i < g->n_sseq; ++i) {
		int32_t jst;
		if (soff[i] == soff[i+1]) continue;
		for (j = soff[i] + 1, jst = j - 1; j <= soff[i+1]; ++j) {
			int32_t is_cont = 0;
			if (j < soff[i+1]) {
				const gfa_seg_t *s = &g->seg[(int32_t)a[j-1]];
				const gfa_seg_t *t = &g->seg[(int32_t)a[j]];
				if (s->soff + s->len == t->soff) {
					int32_t k, nv;
					const gfa_arc_t *av;
					nv = gfa_arc_n(g, (uint32_t)a[j-1]<<1);
					av = gfa_arc_a(g, (uint32_t)a[j-1]<<1);
					for (k = 0; k < nv; ++k)
						if (av[k].w == (uint32_t)a[j]<<1)
							break;
					if (k < nv) is_cont = 1;
				}
			}
			if (!is_cont) {
				int32_t l;
				const gfa_seg_t *s = &g->seg[(int32_t)a[jst]];
				gfa_sfa_t *p = &sfa[k++];
				assert(jst < j);
				p->snid = s->snid, p->soff = s->soff, p->rank = s->rank;
				p->end[0] = find_join(g, (uint32_t)a[jst]<<1|1);
				if (p->end[0] != (uint64_t)-1) p->end[0] ^= 1;
				p->end[1] = find_join(g, (uint32_t)a[j-1]<<1);
				for (l = jst, p->len = 0; l < j; ++l)
					p->len += g->seg[(int32_t)a[l]].len;
				if (write_seq) {
					GFA_MALLOC(p->seq, p->len + 1);
					for (l = jst, p->len = 0; l < j; ++l) {
						s = &g->seg[(int32_t)a[l]];
						memcpy(&p->seq[p->len], s->seq, s->len);
						p->len += s->len;
					}
					p->seq[p->len] = 0;
				}
				jst = j;
			}
		}
	}
	assert(k == n_sfa);

end_check:
	free(soff);
	free(a);
	return sfa;
}

void gfa_sort_ref_arc(gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	for (v = 0; v < n_vtx; ++v) {
		gfa_seg_t *s = &g->seg[v>>1];
		int32_t i, nv;
		gfa_arc_t *av, b;
		if (s->rank != 0) continue;
		nv = gfa_arc_n(g, v);
		av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].w;
			gfa_seg_t *t = &g->seg[w>>1];
			if (t->rank == 0 && t->snid == s->snid && (v&1) == (w&1)) {
				if (((v&1) == 0 && s->soff + s->len == t->soff) || ((v&1) == 1 && t->soff + t->len == s->soff))
					break;
			}
		}
		if (nv > 0 && i == nv) fprintf(stderr, "X\t%c%s\t%d\t%s\t%d\n", "><"[v&1], s->name, i, g->sseq[s->snid].name, s->soff);
		assert(nv == 0 || i < nv);
		if (i > 0 && i < nv) b = av[i], av[i] = av[0], av[0] = b;
	}
}

typedef struct {
	int32_t ld, sd, rd;
	int32_t lp, sp;
} bb_aux_t;

static void bb_write_seq(const gfa_t *g, int32_t n, const uint32_t *v, int32_t l_seq, char *seq)
{
	int32_t k, l;
	for (k = n - 1, l = 0; k >= 0; --k) {
		const gfa_seg_t *s = &g->seg[v[k]>>1];
		if (v[k]&1) {
			int32_t p;
			for (p = s->len - 1; p >= 0; --p)
				seq[l++] = gfa_comp_table[(uint8_t)s->seq[p]];
		} else {
			memcpy(&seq[l], s->seq, s->len);
			l += s->len;
		}
	}
	assert(l == l_seq);
	seq[l] = 0;
}

static int32_t bb_n_paths(const gfa_t *g, const gfa_sub_t *sub, int32_t js, int32_t je)
{
	int32_t j, k;
	int64_t *cnt, c;
	GFA_CALLOC(cnt, je - js + 1);
	cnt[0] = 1;
	for (j = js; j < je; ++j) {
		const gfa_subv_t *t = &sub->v[j];
		for (k = 0; k < t->n; ++k) {
			uint64_t a = sub->a[t->off + k];
			int32_t jv = (int32_t)(a>>32);
			if (jv <= j || jv > je) continue;
			if (cnt[jv - js] + cnt[j - js] > INT32_MAX)
				cnt[jv - js] = INT32_MAX;
			else cnt[jv - js] += cnt[j - js];
		}
	}
	c = cnt[je - js];
	free(cnt);
	return c < INT32_MAX? c : INT32_MAX;
}

gfa_bubble_t *gfa_bubble(const gfa_t *g, int32_t *n_bb_)
{
	extern void radix_sort_gfa32(uint32_t*, uint32_t*);
	uint32_t i, *vs, *vmin, *vtmp = 0;
	int32_t n_bb = 0, m_bb = 0, m_vtmp = 0;
	gfa_bubble_t *bb = 0;
	gfa_scbuf_t *scbuf;

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

	scbuf = gfa_scbuf_init(g);
	for (i = 0; i < g->n_sseq; ++i) {
		gfa_sub_t *sub;
		int32_t j, jst, max_a;
		bb_aux_t *ba;

		if (vs[i] == (uint32_t)-1) continue;
		#if 0
		sub = gfa_sub_from(0, g, vs[i], 0);
		#else
		sub = gfa_scc1(0, g, scbuf, vs[i]);
		#endif
		//gfa_sub_print(stderr, g, sub);
		GFA_CALLOC(ba, sub->n_v);
		for (j = 0; j < sub->n_v; ++j)
			ba[j].sd = INT32_MAX, ba[j].lp = ba[j].sp = -1;
		ba[0].sd = 0;
		for (j = 0; j < sub->n_v; ++j) {
			gfa_subv_t *t = &sub->v[j];
			int32_t k;
			for (k = 0; k < t->n; ++k) {
				uint64_t a = sub->a[t->off + k];
				int32_t jv = (int32_t)(a>>32);
				int32_t l = (int32_t)g->arc[(uint32_t)a].v_lv;
				if (jv <= j) continue; // skip loop or cycle
				if (ba[jv].sd >= ba[j].sd + l)
					ba[jv].sd = ba[j].sd + l, ba[jv].sp = j;
				if (ba[jv].ld < ba[j].ld + l)
					ba[jv].ld = ba[j].ld + l, ba[jv].lp = j;
			}
		}
		for (j = 0, jst = 0, max_a = -1; j < sub->n_v; ++j) {
			gfa_subv_t *t = &sub->v[j];
			int32_t k;
			if (j == max_a) {
				const gfa_seg_t *sst = &g->seg[sub->v[jst].v>>1];
				const gfa_seg_t *sen = &g->seg[t->v>>1];
				if (sst->snid == i && sen->snid == i) {
					int32_t n, l;
					uint32_t *v;
					float f;
					gfa_bubble_t *b;

					// basic information
					if (j - jst <= 1) continue;
					if (n_bb == m_bb) GFA_EXPAND(bb, m_bb);
					b = &bb[n_bb++];
					b->snid = i;
					b->vs = sub->v[jst].v;
					b->ve = t->v;
					b->ss = sst->soff + sst->len;
					b->se = sen->soff;
					b->len_min = ba[j].sd - ba[jst].sd - sst->len;
					b->len_max = ba[j].ld - ba[jst].ld - sst->len;
					b->n_paths = bb_n_paths(g, sub, jst, j);
					assert(b->len_min >= 0);
					assert(b->len_max >= 0 && b->len_max >= b->len_min);
					b->n_seg = j - jst + 1;
					l = (b->len_min + 1) + (b->len_max + 1);
					l = (l + 3) / 4 + b->n_seg;
					GFA_CALLOC(b->v, l);
					b->seq_min = (char*)(b->v + b->n_seg);
					b->seq_max = b->seq_min + b->len_min + 1;
					for (k = jst; k <= j; ++k)
						b->v[k - jst] = sub->v[k].v;

					// test bubble involving both strands (mostly inversions)
					if (b->n_seg > m_vtmp) {
						m_vtmp = b->n_seg;
						kroundup32(m_vtmp);
						GFA_REALLOC(vtmp, m_vtmp);
					}
					for (k = 0; k < b->n_seg; ++k) vtmp[k] = b->v[k]>>1;
					radix_sort_gfa32(vtmp, vtmp + b->n_seg);
					for (k = 1; k < b->n_seg; ++k)
						if (vtmp[k] == vtmp[k-1]) break;
					b->is_bidir = (k < b->n_seg);

					// generate sequences and cf_min/cf_max
					GFA_MALLOC(v, j - jst);
					k = j, n = 0, f = 1.0f;
					while (k > jst) {
						if (k < j) v[n++] = sub->v[k].v;
						k = ba[k].sp;
					}
					bb_write_seq(g, n, v, b->len_min, b->seq_min);
					k = j, n = 0, f = 1.0f;
					while (k > jst) {
						if (k < j) v[n++] = sub->v[k].v;
						k = ba[k].lp;
					}
					bb_write_seq(g, n, v, b->len_max, b->seq_max);
					free(v);
				}
				max_a = -1, jst = j;
			}
			for (k = 0; k < t->n; ++k)
				if ((int32_t)(sub->a[t->off + k]>>32) > max_a)
					max_a = sub->a[t->off + k]>>32;
		}
		free(ba);
		gfa_sub_destroy(sub);
	}
	free(vtmp);
	gfa_scbuf_destroy(scbuf);
	free(vs);
	*n_bb_ = n_bb;
	return bb;
}

const char *gfa_parse_reg(const char *s, int32_t *beg, int32_t *end)
{
	int32_t i, k, l, name_end;
	*beg = *end = -1;
	name_end = l = strlen(s);
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
	}
	// parse the interval
	if (name_end < l) {
		char *tmp, *tmp0;
		tmp0 = tmp = (char*)malloc(l - name_end + 1);
		for (i = name_end + 1, k = 0; i < l; ++i)
			if (s[i] != ',') tmp[k++] = s[i];
		tmp[k] = 0;
		if ((*beg = strtol(tmp, &tmp, 10) - 1) < 0) *beg = 0;
		*end = *tmp? strtol(tmp + 1, &tmp, 10) : 1<<29;
		if (*beg > *end) name_end = l;
		free(tmp0);
	}
	if (name_end == l) *beg = 0, *end = 1<<29;
	return s + name_end;
}

static char **gfa_append_list(char **a, uint32_t *n, uint32_t *m, const char *p)
{
	if (*n == *m) GFA_EXPAND(a, *m);
	a[(*n)++] = gfa_strdup(p);
	return a;
}

char **gfa_query_by_id(const gfa_t *g, int32_t n_bb, const gfa_bubble_t *bb, int32_t snid, int32_t start, int32_t end, int *n_seg_)
{ // TODO: This is an inefficient implementationg. Faster query requires to index the bubble intervals first.
	int32_t i, j, last = 0, bb_st = -1, bb_st_on = -1, bb_en = -1, bb_en_on = -1, bb_last = -1;
	uint32_t n_seg = 0, m_seg = 0;
	char **seg = 0;
	assert(start <= end && start >= 0);
	*n_seg_ = 0;
	for (i = 0; i < n_bb; ++i) {
		const gfa_bubble_t *b = &bb[i];
		if (i == 0 || bb[i].snid != bb[i-1].snid) last = 0;
		if (b->snid != snid) continue;
		assert(b->n_seg > 0);
		bb_last = i;
		if (last <= start && start < b->ss) {
			assert(bb_st < 0);
			bb_st = i, bb_st_on = 1;
		} else if (b->ss <= start && start < b->se) {
			bb_st = i, bb_st_on = 0;
		}
		if (last < end && end <= b->ss) {
			assert(bb_st >= 0);
			bb_en = i, bb_en_on = 1;
		} else if (b->ss < end && end <= b->se) {
			bb_en = i, bb_en_on = 0;
		}
		last = b->se;
	}
	if (bb_last < 0) return 0; // snid not found
	if (bb_st < 0) { // on the last stem
		const gfa_seg_t *s = &g->seg[bb[bb_last].v[bb[bb_last].n_seg - 1]>>1];
		assert(s->snid == snid && start >= s->soff);
		if (start < s->soff + s->len)
			seg = gfa_append_list(seg, &n_seg, &m_seg, s->name);
	} else if (bb_st_on && bb_st == bb_en && bb_en_on) { // on one stem
		seg = gfa_append_list(seg, &n_seg, &m_seg, g->seg[bb[bb_st].v[0]>>1].name);
	} else { // extract bubbles
		if (bb_en < 0) bb_en = bb_last;
		for (i = bb_st; i <= bb_en; ++i) {
			int32_t s = i == bb_st? 0 : 1;
			for (j = s; j < bb[i].n_seg; ++j)
				seg = gfa_append_list(seg, &n_seg, &m_seg, g->seg[bb[i].v[j]>>1].name);
		}
	}
	*n_seg_ = n_seg;
	return seg;
}

char **gfa_query_by_reg(const gfa_t *g, int32_t n_bb, const gfa_bubble_t *bb, const char *reg, int *n_seg)
{
	int32_t snid, start, end;
	const char *p;
	char *tmp;
	*n_seg = 0;
	p = gfa_parse_reg(reg, &start, &end);
	if (p == 0) return 0;
	tmp = gfa_strndup(reg, p - reg);
	snid = gfa_sseq_get(g, tmp);
	free(tmp);
	if (snid < 0) return 0;
	return gfa_query_by_id(g, n_bb, bb, snid, start, end, n_seg);
}
