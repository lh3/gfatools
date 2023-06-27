#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include "gfa-priv.h"
#include "kvec.h"
#include "ksort.h"
#include "kdq.h"
KDQ_INIT(uint64_t)

#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, gfa_map64_t, gfa_map64, uint64_t, int32_t, kh_hash_uint64, kh_eq_generic)

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

static int32_t *gfa_append_list(int32_t *a, uint32_t *n, uint32_t *m, int32_t seg)
{
	if (*n == *m) GFA_EXPAND(a, *m);
	a[(*n)++] = seg;
	return a;
}

int32_t *gfa_query_by_id(const gfa_t *g, int32_t n_bb, const gfa_bubble_t *bb, int32_t snid, int32_t start, int32_t end, int *n_seg_)
{ // TODO: This is an inefficient implementationg. Faster query requires to index the bubble intervals first.
	uint32_t n_seg = 0, m_seg = 0;
	int32_t i, j, last = 0, bb_st = -1, bb_st_on = -1, bb_en = -1, bb_en_on = -1, bb_last = -1, *seg = 0;
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
		uint32_t v = bb[bb_last].v[bb[bb_last].n_seg - 1];
		const gfa_seg_t *s = &g->seg[v>>1];
		assert(s->snid == snid && start >= s->soff);
		if (start < s->soff + s->len)
			seg = gfa_append_list(seg, &n_seg, &m_seg, v>>1);
	} else if (bb_st_on && bb_st == bb_en && bb_en_on) { // on one stem
		seg = gfa_append_list(seg, &n_seg, &m_seg, bb[bb_st].v[0]>>1);
	} else { // extract bubbles
		if (bb_en < 0) bb_en = bb_last;
		for (i = bb_st; i <= bb_en; ++i) {
			int32_t s = i == bb_st? 0 : 1;
			for (j = s; j < bb[i].n_seg; ++j)
				seg = gfa_append_list(seg, &n_seg, &m_seg, bb[i].v[j]>>1);
		}
	}
	*n_seg_ = n_seg;
	return seg;
}

int32_t *gfa_query_by_reg(const gfa_t *g, int32_t n_bb, const gfa_bubble_t *bb, const char *reg, int *n_seg)
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

gfa_t *gfa_subview(gfa_t *g, int32_t n_seg, const int32_t *seg)
{
	gfa_map64_t *h;
	gfa_t *f;
	int32_t i, j, absent;
	uint32_t k, l;

	h = gfa_map64_init();
	for (i = j = 0; i < n_seg; ++i)
		if (seg[i] < g->n_seg) {
			k = gfa_map64_put(h, seg[i], &absent);
			if (absent) kh_val(h, k) = j++;
		}

	GFA_CALLOC(f, 1);
	f->n_seg = kh_size(h);
	GFA_CALLOC(f->seg, f->n_seg);
	for (k = 0; k < kh_end(h); ++k) {
		int32_t nv, s, i, j, t;
		const gfa_arc_t *av;
		if (!kh_exist(h, k)) continue;
		s = kh_key(h, k), t = kh_val(h, k);
		f->seg[t] = g->seg[s];
		for (j = 0; j < 2; ++j) {
			uint32_t v = (uint32_t)s<<1 | j;
			nv = gfa_arc_n(g, v);
			av = gfa_arc_a(g, v);
			for (i = 0; i < nv; ++i) {
				gfa_arc_t *a;
				l = gfa_map64_get(h, av[i].w>>1);
				if (l != kh_end(h)) {
					if (f->n_arc == f->m_arc) {
						f->m_arc += (f->m_arc>>1) + 16;
						GFA_REALLOC(f->arc, f->m_arc);
					}
					a = &f->arc[f->n_arc++];
					*a = av[i];
					a->v_lv = (uint64_t)t<<33 | (uint64_t)j<<32 | a->v_lv<<32>>32;
					a->w = (uint32_t)kh_val(h, l)<<1 | (a->w&1);
				}
			}
		}
	}

	gfa_map64_destroy(h);
	gfa_arc_sort(f);
	gfa_arc_index(f);
	f->sseq = g->sseq;
	f->link_aux = g->link_aux;
	return f;
}

void gfa_subview_destroy(gfa_t *f)
{
	free(f->idx); free(f->arc); free(f->seg); free(f);
}

int32_t *gfa_sub_extend(const gfa_t *g, int n_seg, const int32_t *seg, int step, int32_t *n_ret)
{
	uint32_t v;
	int32_t i, k, *ret;
	int8_t *flag, *sflag;
	kdq_t(uint64_t) *q;
	if (n_seg == 0) return 0;
	q = kdq_init(uint64_t, 0);
	GFA_CALLOC(flag, g->n_seg * 2);
	for (i = 0; i < n_seg; ++i) {
		kdq_push(uint64_t, q, ((uint64_t)seg[i]<<1|0)<<32);
		kdq_push(uint64_t, q, ((uint64_t)seg[i]<<1|1)<<32);
	}
	while (kdq_size(q) > 0) {
		uint64_t x = *kdq_shift(uint64_t, q);
		uint32_t v = x>>32;
		int r = (int32_t)x;
		if (flag[v]) continue; // already visited
		flag[v] = 1;
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
	GFA_CALLOC(sflag, g->n_seg);
	for (v = 0; v < gfa_n_vtx(g); ++v)
		if (flag[v]) sflag[v>>1] = 1;
	free(flag);
	for (i = 0, k = 0; i < g->n_seg; ++i)
		if (sflag[i]) ++k;
	GFA_CALLOC(ret, k);
	for (i = 0, k = 0; i < g->n_seg; ++i)
		if (sflag[i]) ret[k++] = i;
	free(sflag);
	*n_ret = k;
	return ret;
}

int32_t *gfa_list2seg(const gfa_t *g, int32_t n_seg, char *const* seg, int32_t *n_ret)
{
	int32_t i, k, *ret;
	GFA_MALLOC(ret, n_seg);
	for (i = k = 0; i < n_seg; ++i) {
		int32_t s;
		s = gfa_name2id(g, seg[i]);
		if (s >= 0) ret[k++] = s;
	}
	*n_ret = k;
	return ret;
}

void gfa_walk_flip(gfa_t *g)
{
	int32_t i, j;
	int8_t *strand;
	if (g->n_walk == 0) return;
	GFA_CALLOC(strand, g->n_seg);
	for (i = 0; i < g->n_walk; ++i) {
		gfa_walk_t *w = &g->walk[i];
		for (j = 0; j < w->n_v; ++j)
			if (strand[w->v[j]>>1] == 0)
				strand[w->v[j]>>1] = w->v[j]&1? -1 : 1;
	}
	for (i = 0; i < g->n_walk; ++i) {
		gfa_walk_t *w = &g->walk[i];
		int32_t n[2];
		n[0] = n[1] = 0;
		for (j = 0; j < w->n_v; ++j) {
			int8_t s;
			assert(strand[w->v[j]>>1] != 0);
			s = w->v[j]&1? -1 : 1;
			if (s == strand[w->v[j]>>1]) ++n[0];
			else ++n[1];
		}
		if (n[0] >= n[1]) continue;
		for (j = 0; j < w->n_v>>1; ++j) {
			uint32_t t = w->v[j]^1;
			w->v[j] = w->v[w->n_v - 1 - j]^1;
			w->v[w->n_v - 1 - j] = t;
		}
		if (w->n_v&1) w->v[w->n_v>>1] ^= 1;
	}
	free(strand);
}
