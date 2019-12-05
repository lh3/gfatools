#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "gfa-priv.h"
#include "kvec.h"

/********************************************
 * Preprocessing up to transitive reduction *
 ********************************************/

typedef struct { uint32_t n, m; uint32_t *a; } gfa32_v;

// delete short arcs
int gfa_arc_del_short(gfa_t *g, int min_ovlp_len, float drop_ratio)
{
	uint32_t v, n_vtx = gfa_n_vtx(g), n_short = 0;
	for (v = 0; v < n_vtx; ++v) {
		gfa_arc_t *av = gfa_arc_a(g, v);
		uint32_t i, thres, nv = gfa_arc_n(g, v);
		if (nv < 2) continue;
		thres = (uint32_t)(av[0].ov * drop_ratio + .499);
		if (thres < min_ovlp_len) thres = min_ovlp_len;
		for (i = nv - 1; i >= 1 && av[i].ov < thres; --i);
		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
	}
	if (gfa_verbose >= 3) fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
	if (n_short) {
		gfa_cleanup(g);
		gfa_fix_symm_del(g);
	}
	return n_short;
}

// delete multi-arcs
int gfa_arc_del_multi_risky(gfa_t *g)
{
	uint32_t *cnt, n_vtx = gfa_n_vtx(g), n_multi = 0, v;
	GFA_CALLOC(cnt, n_vtx);
	for (v = 0; v < n_vtx; ++v) {
		gfa_arc_t *av = gfa_arc_a(g, v);
		int32_t i, nv = gfa_arc_n(g, v);
		if (nv < 2) continue;
		for (i = nv - 1; i >= 0; --i) ++cnt[av[i].w];
		for (i = nv - 1; i >= 0; --i)
			if (--cnt[av[i].w] != 0)
				av[i].del = 1, ++n_multi;
	}
	free(cnt);
	if (gfa_verbose >= 3 && n_multi > 0) fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
	if (n_multi) gfa_cleanup(g);
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
int gfa_arc_del_asymm_risky(gfa_t *g)
{
	uint32_t e, n_asymm = 0;
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].w^1, u = g->arc[e].v_lv>>32^1;
		uint32_t i, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].w == u) break;
		if (i == nv) g->arc[e].del = 1, ++n_asymm;
	}
	if (gfa_verbose >= 3 && n_asymm > 0) fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
	if (n_asymm) gfa_cleanup(g);
	return n_asymm;
}

void gfa_fix_symm_del(gfa_t *g)
{
	gfa_arc_del_multi_risky(g);
	gfa_arc_del_asymm_risky(g);
}

// transitive reduction; see Myers, 2005
int gfa_arc_del_trans(gfa_t *g, int fuzz)
{
	uint8_t *mark;
	uint32_t v, n_vtx = gfa_n_vtx(g), n_reduced = 0;

	mark = (uint8_t*)calloc(n_vtx, 1);
	for (v = 0; v < n_vtx; ++v) {
		uint32_t L, i, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		if (nv == 0) continue; // no hits
		if (g->seg[v>>1].del) {
			for (i = 0; i < nv; ++i) av[i].del = 1, ++n_reduced;
			continue;
		}
		for (i = 0; i < nv; ++i)
			mark[av[i].w] = g->seg[av[i].w>>1].del? 2 : 1;
		L = gfa_arc_len(av[nv-1]) + fuzz;
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].w;
			uint32_t j, nw = gfa_arc_n(g, w);
			gfa_arc_t *aw = gfa_arc_a(g, w);
			if (mark[av[i].w] != 1) continue;
			for (j = 0; j < nw && gfa_arc_len(aw[j]) + gfa_arc_len(av[i]) <= L; ++j)
				if (mark[aw[j].w]) mark[aw[j].w] = 2;
		}
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].w;
			uint32_t j, nw = gfa_arc_n(g, w);
			gfa_arc_t *aw = gfa_arc_a(g, w);
			for (j = 0; j < nw && (j == 0 || gfa_arc_len(aw[j]) < fuzz); ++j)
				if (mark[aw[j].w]) mark[aw[j].w] = 2;
		}
		for (i = 0; i < nv; ++i) {
			if (mark[av[i].w] == 2) av[i].del = 1, ++n_reduced;
			mark[av[i].w] = 0;
		}
	}
	free(mark);
	if (gfa_verbose >= 3) fprintf(stderr, "[M::%s] transitively reduced %d arcs\n", __func__, n_reduced);
	if (n_reduced) {
		gfa_cleanup(g);
		gfa_fix_symm_del(g);
	}
	return n_reduced;
}

int gfa_arc_del_weak(gfa_t *g)
{
	uint32_t n_vtx = gfa_n_vtx(g), v, n_abnormal = 0, n_del = 0;
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, nv = gfa_arc_n(g, v), n_strong = 0, top_strong = 0;
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			if (!av[i].strong) continue;
			if (n_strong == i) ++top_strong;
			++n_strong;
		}
		if (n_strong == 0) continue;
		if (top_strong == 0) {
			++n_abnormal;
		} else {
			for (i = 0; i < nv; ++i)
				if (!av[i].strong)
					av[i].del = 1, ++n_del;
		}
	}
	if (gfa_verbose >= 3 && n_del > 0)
		fprintf(stderr, "[M::%s] removed %d weak arcs; %d abnormal arcs\n", __func__, n_del, n_abnormal);
	if (n_del) {
		gfa_cleanup(g);
		gfa_fix_symm_del(g);
	}
	return n_del;
}

int gfa_arc_pair_strong(gfa_t *g)
{
	uint32_t e, n_flip = 0;
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v, u, i, nv;
		gfa_arc_t *av;
		if (g->arc[e].strong == 0) continue;
		v = g->arc[e].w^1;
		u = g->arc[e].v_lv>>32^1;
		nv = gfa_arc_n(g, v);
		av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].w == u && av[i].strong == 0)
				av[i].strong = 1, ++n_flip;
	}
	if (gfa_verbose >= 3 && n_flip > 0) fprintf(stderr, "[M::%s] added %d strong arcs\n", __func__, n_flip);
	return n_flip;
}

/********************************
 * Probe unitig ends from reads *
 ********************************/

#define GFA_VT_MERGEABLE 0
#define GFA_VT_TIP       1
#define GFA_VT_MULTI_OUT 2
#define GFA_VT_MULTI_IN  3

static inline int32_t gfa_deg(const gfa_t *g, uint32_t v) // compute out-degree, excluding deleted edges
{
	uint32_t i, n, nv;
	const gfa_arc_t *av;
	if (g->seg[v>>1].del) return 0;
	nv = gfa_arc_n(g, v);
	av = gfa_arc_a(g, v);
	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
}

static inline int32_t gfa_deg2(const gfa_t *g, uint32_t v, uint32_t *w, int32_t *l)
{
	uint32_t i, nv, nv0, k;
	int32_t min_l;
	const gfa_arc_t *av;
	*l = 0, *w = (uint32_t)-1;
	if (g->seg[v>>1].del) return 0;
	min_l = g->seg[v>>1].len;
	nv0 = k = gfa_arc_n(g, v);
	av = gfa_arc_a(g, v);
	for (i = nv = 0; i < nv0; ++i)
		if (!av[i].del)
			++nv, k = i, min_l = gfa_arc_len(av[i]) < min_l? gfa_arc_len(av[i]) : min_l;
	*l = min_l;
	*w = nv == 1? av[k].w : (uint32_t)-1;
	return nv;
}

static inline int32_t gfa_vtype(const gfa_t *g, uint32_t v, uint32_t *w_, int32_t *l_)
{
	int32_t nv, nw, l;
	uint32_t w;
	nv = gfa_deg2(g, v, &w, &l);
	if (l_) *l_ = l;
	if (w_) *w_ = w;
	if (nv == 0) return GFA_VT_TIP;
	if (nv > 1) return GFA_VT_MULTI_OUT;
	nw = gfa_deg(g, w^1);
	return nw == 1? GFA_VT_MERGEABLE : GFA_VT_MULTI_IN;
}

static inline int32_t gfa_uext(const gfa_t *g, uint32_t v, int32_t max_ext, int32_t *ne, int32_t *le, uint32_t *end_v, gfa32_v *a)
{
	int32_t vt, n_ext = 0, l_ext = 0;
	if (a) a->n = 0;
	if (a) kv_push(uint32_t, *a, v);
	if (end_v) *end_v = (uint32_t)-1;
	do {
		uint32_t w;
		int32_t l;
		vt = gfa_vtype(g, v, &w, &l);
		l_ext += l;
		if (end_v && (vt == GFA_VT_MERGEABLE || vt == GFA_VT_MULTI_IN)) *end_v = w;
		if (vt != GFA_VT_MERGEABLE) break;
		++n_ext;
		if (a) kv_push(uint32_t, *a, w);
		v = w;
	} while (--max_ext > 0);
	if (ne) *ne = n_ext;
	if (le) *le = l_ext;
	return vt;
}

/*************************
 * Topological extension *
 *************************/

typedef struct {
	uint32_t p; // the optimal parent vertex
	int32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of reads
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: whether the vertex has been visited before
} gfa_tinfo_t;

typedef struct {
	gfa_tinfo_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
	int32_t n_short_tip, n_bb, dist; // n_short_tip: number of short tips; n_bb: number of bubbles; dist: min distance from one of the input vertices
	uint32_t v_end;     // end vertex; dist and v_end only set when n_bb>0
} gfa_tbuf_t;

#define GFA_TE_THRU_SHORT_TIP  0x1
#define GFA_TE_THRU_BUBBLE     0x2

#define GFA_TE_TYPE_NONE       0x0
#define GFA_TE_TYPE_BUBBLE     0x1
#define GFA_TE_TYPE_SHORT_TIP  0x2
#define GFA_TE_TYPE_SELF_CYC   0x4
#define GFA_TE_TYPE_SELF_BICYC 0x8

static gfa_tbuf_t *gfa_tbuf_init(const gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	gfa_tbuf_t *b;
	GFA_CALLOC(b, 1);
	GFA_CALLOC(b->a, gfa_n_vtx(g));
	for (v = 0; v < n_vtx; ++v)
		b->a[v].p = (uint32_t)-1, b->a[v].d = INT32_MIN;
	return b;
}

static void gfa_tbuf_destroy(gfa_tbuf_t *b)
{
	if (b == 0) return;
	free(b->a); free(b->S.a); free(b->b.a); free(b->e.a);
	free(b);
}

static void gfa_tbuf_reset(gfa_tbuf_t *b)
{
	uint32_t i;
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		gfa_tinfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->r = 0, t->p = (uint32_t)-1, t->d = INT32_MIN;
	}
}

static int32_t gfa_topo_ext(const gfa_t *g, uint32_t n_v0, const uint32_t *v0, int32_t max_dist, int32_t thru_flag, gfa_tbuf_t *b)
{
	uint32_t i, type = 0, n_pending = 0; // n_pending: number of visited vertices that are not sorted
	int32_t max_d = INT32_MIN; // max_d: max gfa_tinfo_t::d of visited vertices

	b->S.n = b->b.n = b->e.n = 0;
	b->n_short_tip = b->n_bb = 0, b->dist = -1, b->v_end = (uint32_t)-1;
	if (n_v0 == 0) return 0;
	for (i = 0; i < n_v0; ++i) {
		uint32_t v = v0[i];
		gfa_tinfo_t *t = &b->a[v];
		if (g->seg[v>>1].del) continue;
		t->p = (uint32_t)-1, t->r = 0, t->c = t->s = 0; // ->s has to be 0, as gfa_tbuf_reset() doesn't reset the initial vertices
		t->d = -(int32_t)g->seg[v>>1].len;
		kv_push(uint32_t, b->S, v);
	}

	while (b->S.n > 0 && max_d <= max_dist) {
		uint32_t v = kv_pop(b->S), nv = gfa_arc_n(g, v);
		int32_t d = b->a[v].d, c = b->a[v].c, end_dist = d + (int32_t)g->seg[v>>1].len;
		gfa_arc_t *av = gfa_arc_a(g, v);
		if (n_pending == 0 && end_dist > b->dist)
			b->dist = end_dist, b->v_end = v;
		if (gfa_deg(g, v) == 0) { // a tip
			if (end_dist < max_dist) { // a tip shorter than max_dist
				++b->n_short_tip;
				if (thru_flag & GFA_TE_THRU_SHORT_TIP) continue;
				else break;
			} else break; // if we come here, we have a tip beyond max_dist; we stop
		}
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			uint32_t j, w = av[i].w;
			int32_t l = (int32_t)av[i].v_lv; // u->w with length l
			gfa_tinfo_t *t = &b->a[w];
			if (av[i].del) continue;
			for (j = 0; j < n_v0; ++j) // test cycles or bidirected cycles involving the starting vertices
				if (w>>1 == v0[j]>>1 && !g->seg[v0[j]>>1].del) break;
			if (j < n_v0) {
				type |= w == v0[j]? GFA_TE_TYPE_SELF_CYC : GFA_TE_TYPE_SELF_BICYC; // cycle or bidirected cycle
				break;
			}
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i); // save the edge
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for gfa_tbuf_reset()
				t->p = v, t->s = 1, t->d = d + l, t->c = c + 1;
				t->r = gfa_deg(g, w^1);
				++n_pending;
			} else { // visited before
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				if (d + l < t->d) t->d = d + l; // update dist
			}
			max_d = max_d > t->d? max_d : t->d;
			assert(t->r > 0);
			if (--(t->r) == 0) {
				kv_push(uint32_t, b->S, w);
				--n_pending;
			}
		}
		if (i < nv || b->S.n == 0) break;
		if (b->S.n == 1 && n_pending == 0) { // end of a bubble
			uint32_t v = b->S.a[0];
			b->dist = b->a[v].d, b->v_end = v;
			++b->n_bb;
			if (!(thru_flag & GFA_TE_THRU_BUBBLE)) break;
		}
	}
	if (b->n_bb)  type |= GFA_TE_TYPE_BUBBLE;
	if (b->n_short_tip) type |= GFA_TE_TYPE_SHORT_TIP;
	return type;
}

/************************
 * Basic graph cleaning *
 ************************/

int gfa_drop_tip(gfa_t *g, int tip_cnt, int tip_len)
{
	gfa32_v a = {0,0,0};
	uint32_t n_vtx = gfa_n_vtx(g), v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		int32_t l_ext, vt;
		if (g->seg[v>>1].del) continue;
		if (gfa_deg(g, v^1) != 0) continue; // not a tip
		vt = gfa_uext(g, v, tip_cnt, 0, &l_ext, 0, &a);
		if (vt == GFA_VT_MERGEABLE) continue; // not a short unitig
		if (l_ext > tip_len) continue; // tip too long
		for (i = 0; i < a.n; ++i)
			gfa_seg_del(g, a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M::%s] drop %d tips\n", __func__, cnt);
	return cnt;
}

int gfa_drop_internal(gfa_t *g, int max_ext)
{
	gfa32_v a = {0,0,0};
	uint32_t n_vtx = gfa_n_vtx(g), v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		int32_t l_ext;
		if (g->seg[v>>1].del) continue;
		if (gfa_vtype(g, v^1, 0, 0) != GFA_VT_MULTI_IN) continue;
		if (gfa_uext(g, v, max_ext, 0, &l_ext, 0, &a) != GFA_VT_MULTI_IN) continue;
		for (i = 0; i < a.n; ++i)
			gfa_seg_del(g, a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	fprintf(stderr, "[M::%s] cut %d internal sequences\n", __func__, cnt);
	if (cnt > 0) gfa_cleanup(g);
	return cnt;
}

static inline int32_t gfa_is_extended(const gfa_t *g, uint32_t v, int32_t min_dist, int32_t max_dist, gfa_tbuf_t *b)
{
	if (g->seg[v>>1].len >= min_dist) return 1;
	gfa_topo_ext(g, 1, &v, max_dist - g->seg[v>>1].len, GFA_TE_THRU_BUBBLE, b);
	gfa_tbuf_reset(b);
	return b->dist + g->seg[v>>1].len > min_dist? 1 : 0;
}

int gfa_cut_z(gfa_t *g, int32_t min_dist, int32_t max_dist)
{
	uint32_t n_cut = 0;
	int64_t e;
	gfa_tbuf_t *b;
	b = gfa_tbuf_init(g);
	for (e = 0; e < g->n_arc; ++e) {
		gfa_arc_t *a = &g->arc[e], *av, *aw;
		uint32_t v, w, u[2];
		if (a->del) continue;
		v = a->v_lv>>32, w = a->w^1;
		if (gfa_arc_n(g, v) != 2 || gfa_arc_n(g, w) != 2) continue;
		av = gfa_arc_a(g, v);
		aw = gfa_arc_a(g, w);
		if (av[0].del || av[1].del || aw[0].del || aw[1].del) continue;
		assert(av[0].w == (w^1) || av[1].w == (w^1));
		assert(aw[0].w == (v^1) || aw[1].w == (v^1));
		u[0] = av[0].w != (w^1)? av[0].w : av[1].w;
		u[1] = aw[0].w != (v^1)? aw[0].w : aw[1].w;
		if (gfa_deg(g, u[0]^1) != 1) continue;
		if (gfa_deg(g, u[1]^1) != 1) continue;
		if (!gfa_is_extended(g, v^1, min_dist, max_dist, b)) continue;
		if (!gfa_is_extended(g, w^1, min_dist, max_dist, b)) continue;
		if (!gfa_is_extended(g, u[0], min_dist, max_dist, b)) continue;
		if (!gfa_is_extended(g, u[1], min_dist, max_dist, b)) continue;
		a->del = 1;
		gfa_arc_del(g, w, v^1, 1);
		++n_cut;
	}
	gfa_tbuf_destroy(b);
	fprintf(stderr, "[M::%s] cut %d Z-junctions\n", __func__, n_cut);
	if (n_cut > 0) gfa_cleanup(g);
	return n_cut;
}

int gfa_topocut(gfa_t *g, float drop_ratio, int32_t tip_cnt, int32_t tip_len)
{
	uint32_t n_vtx = gfa_n_vtx(g), v, n_cut = 0;
	uint64_t k, n_b = 0, *b;

	// collect all overlaps at bifurcations
	assert(g->n_arc < UINT32_MAX);
	GFA_MALLOC(b, g->n_arc); // FIXME: we don't need an array this large, but probably it doesn't matter.
	for (v = 0; v < n_vtx; ++v) {
		gfa_arc_t *av = gfa_arc_a(g, v);
		uint32_t i, nv = gfa_arc_n(g, v);
		if (nv < 2) continue;
		for (i = 0; i < nv; ++i)
			b[n_b++] = (uint64_t)av[i].ov << 32 | (av - g->arc + i);
	}
	radix_sort_gfa64(b, b + n_b);

	// test each edge
	for (k = 0; k < n_b; ++k) {
		gfa_arc_t *a = &g->arc[(uint32_t)b[k]];
		uint32_t i, iv, iw, v = a->v_lv>>32, w = a->w^1, to_del = 0;
		uint32_t nv = gfa_arc_n(g, v), nw = gfa_arc_n(g, w), kv, kw;
		uint32_t ov_max = 0, ow_max = 0;
		int32_t vt, l_ext;
		gfa_arc_t *av, *aw;

		if (nv == 1 && nw == 1) continue;
		av = gfa_arc_a(g, v);
		aw = gfa_arc_a(g, w);

		for (i = 0, kv = 0; i < nv; ++i) {
			if (av[i].del) continue;
			if (ov_max < av[i].ov) ov_max = av[i].ov;
			++kv;
		}
		if (kv >= 2 && a->ov > ov_max * drop_ratio) continue;
		for (i = 0, kw = 0; i < nw; ++i) {
			if (aw[i].del) continue;
			if (ow_max < aw[i].ov) ow_max = aw[i].ov;
			++kw;
		}
		if (kw >= 2 && a->ow > ow_max * drop_ratio) continue;
		if (kv == 1 && kw == 1) continue;

		for (iv = 0; iv < nv; ++iv)
			if (av[iv].w == (w^1)) break;
		assert(iv < nv);
		for (iw = 0; iw < nw; ++iw)
			if (aw[iw].w == (v^1)) break;
		assert(iw < nw);
		assert(av[iv].del == aw[iw].del);
		if (av[iv].del && aw[iw].del) continue;

		if (kv > 1 && kw > 1) {
			if (a->ov < ov_max * drop_ratio && a->ow < ow_max * drop_ratio)
				to_del = 1;
		} else if (kw == 1) {
			vt = gfa_uext(g, w^1, tip_cnt, 0, &l_ext, 0, 0);
			//if (vt != GFA_VT_MERGEABLE && l_ext < tip_len) to_del = 1;
			if (vt == GFA_VT_MULTI_IN && l_ext < tip_len) to_del = 1;
		} else if (kv == 1) {
			vt = gfa_uext(g, v^1, tip_cnt, 0, &l_ext, 0, 0);
			//if (vt != GFA_VT_MERGEABLE && l_ext < tip_len) to_del = 1;
			if (vt == GFA_VT_MULTI_IN && l_ext < tip_len) to_del = 1;
		}
		if (to_del)
			av[iv].del = aw[iw].del = 1, ++n_cut;
	}

	free(b);
	if (gfa_verbose >= 3) fprintf(stderr, "[M::%s] %d topology-aware cuts\n", __func__, n_cut);
	if (n_cut) {
		gfa_cleanup(g);
		gfa_fix_symm_del(g);
	}
	return n_cut;
}

/******************
 * Bubble popping *
 ******************/

// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void gfa_bub_backtrack(gfa_t *g, uint32_t v0, int max_del, gfa_tbuf_t *b)
{
	uint32_t i, v;
	assert(b->S.n == 1);
	if (max_del > 0) {
		int32_t n_kept = 0;
		v = b->v_end;
		do { ++n_kept, v = b->a[v].p; } while (v != v0);
		if (b->b.n - n_kept > max_del) return;
	}
	for (i = 0; i < b->b.n; ++i)
		g->seg[b->b.a[i]>>1].del = 1;
	for (i = 0; i < b->e.n; ++i) {
		gfa_arc_t *a = &g->arc[b->e.a[i]];
		a->del = 1;
		gfa_arc_del(g, a->w^1, a->v_lv>>32^1, 1);
	}
	v = b->v_end;
	do {
		uint32_t u = b->a[v].p; // u->v
		g->seg[v>>1].del = 0;
		gfa_arc_del(g, u, v, 0);
		gfa_arc_del(g, v^1, u^1, 0);
		v = u;
	} while (v != v0);
}

// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static int32_t gfa_bub_pop1(gfa_t *g, uint32_t v0, int max_dist, int max_del, int protect_tip, gfa_tbuf_t *b)
{
	uint64_t ret = 0;
	if (gfa_deg(g, v0) < 2) return 0; // no bubbles
	gfa_topo_ext(g, 1, &v0, max_dist, protect_tip? 0 : GFA_TE_THRU_SHORT_TIP, b);
	if (b->n_bb) {
		gfa_bub_backtrack(g, v0, max_del, b);
		ret = 1 | (uint64_t)b->n_short_tip << 32;
	}
	gfa_tbuf_reset(b);
	return ret;
}

// pop bubbles
int gfa_pop_bubble(gfa_t *g, int max_dist, int max_del, int protect_tip)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	uint64_t n_pop = 0;
	gfa_tbuf_t *b;
	b = gfa_tbuf_init(g);
	for (v = 0; v < n_vtx; ++v)
		if (!g->seg[v>>1].del && gfa_deg(g, v) >= 2)
			n_pop += gfa_bub_pop1(g, v, max_dist, max_del, protect_tip, b);
	gfa_tbuf_destroy(b);
	if (n_pop) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M::%s] popped %d bubbles and trimmed short %d tips\n", __func__, (uint32_t)n_pop, (uint32_t)(n_pop>>32));
	return n_pop;
}

/****************
 * Unitig graph *
 ****************/

#include "kdq.h"
KDQ_INIT(uint64_t)

#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])

gfa_t *gfa_ug_gen(const gfa_t *g)
{
	int32_t *mark;
	uint32_t i, v, n_vtx = gfa_n_vtx(g);
	kdq_t(uint64_t) *q;
	gfa_t *ug;

	ug = gfa_init();
	mark = (int32_t*)calloc(n_vtx, 4);

	q = kdq_init(uint64_t);
	for (v = 0; v < n_vtx; ++v) {
		uint32_t w, x, l, start, end, len, tmp, len_r;
		char utg_name[11];
		gfa_seg_t *u;
		gfa_arc_t *a;

		if (g->seg[v>>1].del || mark[v]) continue;
		mark[v] = 1;
		q->count = 0, start = v, end = v^1, len = len_r = 0;
		// forward
		w = v;
		while (1) {
			if (arc_cnt(g, w) != 1) break;
			x = arc_first(g, w).w; // w->x
			if (arc_cnt(g, x^1) != 1) break;
			mark[x] = mark[w^1] = 1;
			a = &arc_first(g, w);
			l = gfa_arc_len(*a);
			kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
			end = x^1, len += l, len_r += gfa_arc_lw(g, *a);
			w = x;
			if (x == v) break;
		}
		if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
			l = g->seg[end>>1].len;
			kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
			len += l;
		} else { // circular unitig
			start = end = UINT32_MAX;
			goto add_unitig; // then it is not necessary to do the backward
		}
		// backward
		x = v;
		while (1) { // similar to forward but not the same
			if (arc_cnt(g, x^1) != 1) break;
			w = arc_first(g, x^1).w ^ 1; // w->x
			if (arc_cnt(g, w) != 1) break;
			mark[x] = mark[w^1] = 1;
			a = &arc_first(g, w);
			l = gfa_arc_len(*a);
			kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
			start = w, len += l, len_r += gfa_arc_lw(g, *a);
			x = w;
		}
		len_r += g->seg[start>>1].len;
add_unitig:
		if (start != UINT32_MAX) mark[start] = mark[end] = 1;
		sprintf(utg_name, "utg%.7d", ug->n_seg + 1);
		tmp = gfa_add_seg(ug, utg_name);
		u = &ug->seg[tmp];
		u->seq = 0, u->len = len;
		GFA_MALLOC(u->utg, 1);
		u->utg->start = start, u->utg->end = end, u->utg->n = kdq_size(q), u->circ = (start == UINT32_MAX);
		u->utg->m = u->utg->n;
		kroundup32(u->utg->m);
		u->utg->a = (uint64_t*)malloc(8 * u->utg->m);
		u->utg->name = (char**)malloc(sizeof(char*) * u->utg->m);
		u->utg->len2 = len_r;
		for (i = 0; i < kdq_size(q); ++i) {
			u->utg->a[i] = kdq_at(q, i);
			u->utg->name[i] = gfa_strdup(g->seg[u->utg->a[i]>>33].name);
		}
	}
	kdq_destroy(uint64_t, q);

	// add arcs between unitigs; reusing mark for a different purpose
	for (v = 0; v < n_vtx; ++v) mark[v] = -1;
	for (i = 0; i < ug->n_seg; ++i) {
		if (ug->seg[i].circ) continue;
		mark[ug->seg[i].utg->start] = i<<1 | 0;
		mark[ug->seg[i].utg->end] = i<<1 | 1;
	}
	for (i = 0; i < g->n_arc; ++i) {
		gfa_arc_t *p = &g->arc[i];
		if (p->del) continue;
		if (mark[p->v_lv>>32^1] >= 0 && mark[p->w] >= 0) {
			gfa_seg_t *s1 = &ug->seg[mark[p->v_lv>>32^1]>>1];
			gfa_seg_t *s2 = &ug->seg[mark[p->w]>>1];
			int ov = p->ov, ow = p->ow;
			if (ov >= s1->len) ov = s1->len - 1;
			if (ow >= s2->len) ow = s2->len - 1;
			gfa_add_arc1(ug, mark[p->v_lv>>32^1]^1, mark[p->w], ov, ow, -1, 0);
		}
	}
	free(mark);
	gfa_finalize(ug);
	return ug;
}
