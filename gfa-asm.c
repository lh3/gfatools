#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "gfa.h"
#include "kvec.h"

/**********************
 * Graph manipulation *
 **********************/

typedef struct { size_t n, m; uint64_t *a; } gfa64_v;

// delete short arcs
int gfa_arc_del_short(gfa_t *g, float drop_ratio)
{
	uint32_t v, n_vtx = gfa_n_vtx(g), n_short = 0;
	for (v = 0; v < n_vtx; ++v) {
		gfa_arc_t *av = gfa_arc_a(g, v);
		uint32_t i, thres, nv = gfa_arc_n(g, v);
		if (nv < 2) continue;
		thres = (uint32_t)(av[0].ov * drop_ratio + .499);
		for (i = nv - 1; i >= 1 && av[i].ov < thres; --i);
		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
	}
	if (n_short) {
		gfa_cleanup(g);
		gfa_symm(g);
	}
	if (gfa_verbose >= 3) fprintf(stderr, "[M] removed %d short overlaps\n", n_short);
	return n_short;
}

// delete multi-arcs
static int gfa_arc_del_multi(gfa_t *g)
{
	uint32_t *cnt, n_vtx = gfa_n_vtx(g), n_multi = 0, v;
	cnt = (uint32_t*)calloc(n_vtx, 4);
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
	if (n_multi) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] removed %d multi-arcs\n", n_multi);
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
static int gfa_arc_del_asymm(gfa_t *g)
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
	if (n_asymm) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] removed %d asymmetric arcs\n", n_asymm);
	return n_asymm;
}

void gfa_symm(gfa_t *g)
{
	gfa_arc_del_multi(g);
	gfa_arc_del_asymm(g);
	g->is_symm = 1;
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
		for (i = 0; i < nv; ++i) mark[av[i].w] = 1;
		L = gfa_arc_len(av[nv-1]) + fuzz;
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].w;
			uint32_t j, nw = gfa_arc_n(g, w);
			gfa_arc_t *aw = gfa_arc_a(g, w);
			if (mark[av[i].w] != 1) continue;
			for (j = 0; j < nw && gfa_arc_len(aw[j]) + gfa_arc_len(av[i]) <= L; ++j)
				if (mark[aw[j].w]) mark[aw[j].w] = 2;
		}
		#if 0
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].w;
			uint32_t j, nw = gfa_arc_n(g, w);
			gfa_arc_t *aw = gfa_arc_a(g, w);
			for (j = 0; j < nw && (j == 0 || gfa_arc_len(aw[j]) < fuzz); ++j)
				if (mark[aw[j].w]) mark[aw[j].v] = 2;
		}
		#endif
		for (i = 0; i < nv; ++i) {
			if (mark[av[i].w] == 2) av[i].del = 1, ++n_reduced;
			mark[av[i].w] = 0;
		}
	}
	free(mark);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] transitively reduced %d arcs\n", n_reduced);
	if (n_reduced) {
		gfa_cleanup(g);
		gfa_symm(g);
	}
	return n_reduced;
}

/**********************************
 * Filter short potential unitigs *
 **********************************/

#define GFA_ET_MERGEABLE 0
#define GFA_ET_TIP       1
#define GFA_ET_MULTI_OUT 2
#define GFA_ET_MULTI_NEI 3

static inline int gfa_is_utg_end(const gfa_t *g, uint32_t v, uint64_t *lw)
{
	uint32_t w, nv, nw, nw0, nv0 = gfa_arc_n(g, v^1);
	int i, i0 = -1;
	gfa_arc_t *aw, *av = gfa_arc_a(g, v^1);
	for (i = nv = 0; i < nv0; ++i)
		if (!av[i].del) i0 = i, ++nv;
	if (nv == 0) return GFA_ET_TIP; // tip
	if (nv > 1) return GFA_ET_MULTI_OUT; // multiple outgoing arcs
	if (lw) *lw = av[i0].v_lv<<32 | av[i0].w;
	w = av[i0].w ^ 1;
	nw0 = gfa_arc_n(g, w);
	aw = gfa_arc_a(g, w);
	for (i = nw = 0; i < nw0; ++i)
		if (!aw[i].del) ++nw;
	if (nw != 1) return GFA_ET_MULTI_NEI;
	return GFA_ET_MERGEABLE;
}

int gfa_extend(const gfa_t *g, uint32_t v, int max_ext, gfa64_v *a)
{
	int ret;
	uint64_t lw;
	a->n = 0;
	kv_push(uint64_t, *a, v);
	do {
		ret = gfa_is_utg_end(g, v^1, &lw);
		if (ret != 0) break;
		kv_push(uint64_t, *a, lw);
		v = (uint32_t)lw;
	} while (--max_ext > 0);
	return ret;
}

int gfa_cut_tip(gfa_t *g, int max_ext)
{
	gfa64_v a = {0,0,0};
	uint32_t n_vtx = gfa_n_vtx(g), v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		if (g->seg[v>>1].del) continue;
		if (gfa_is_utg_end(g, v, 0) != GFA_ET_TIP) continue; // not a tip
		if (gfa_extend(g, v, max_ext, &a) == GFA_ET_MERGEABLE) continue; // not a short unitig
		for (i = 0; i < a.n; ++i)
			gfa_seg_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] cut %d tips\n", cnt);
	return cnt;
}

int gfa_cut_internal(gfa_t *g, int max_ext)
{
	gfa64_v a = {0,0,0};
	uint32_t n_vtx = gfa_n_vtx(g), v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		if (g->seg[v>>1].del) continue;
		if (gfa_is_utg_end(g, v, 0) != GFA_ET_MULTI_NEI) continue;
		if (gfa_extend(g, v, max_ext, &a) != GFA_ET_MULTI_NEI) continue;
		for (i = 0; i < a.n; ++i)
			gfa_seg_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] cut %d internal sequences\n", cnt);
	return cnt;
}

int gfa_cut_biloop(gfa_t *g, int max_ext)
{
	gfa64_v a = {0,0,0};
	uint32_t n_vtx = gfa_n_vtx(g), v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		uint32_t nv, nw, w = UINT32_MAX, x, ov = 0, ox = 0;
		gfa_arc_t *av, *aw;
		if (g->seg[v>>1].del) continue;
		if (gfa_is_utg_end(g, v, 0) != GFA_ET_MULTI_NEI) continue;
		if (gfa_extend(g, v, max_ext, &a) != GFA_ET_MULTI_OUT) continue;
		x = (uint32_t)a.a[a.n - 1] ^ 1;
		nv = gfa_arc_n(g, v ^ 1), av = gfa_arc_a(g, v ^ 1);
		for (i = 0; i < nv; ++i)
			if (!av[i].del) w = av[i].w ^ 1;
		assert(w != UINT32_MAX);
		nw = gfa_arc_n(g, w), aw = gfa_arc_a(g, w);
		for (i = 0; i < nw; ++i) { // we are looking for: v->...->x', w->v and w->x
			if (aw[i].del) continue;
			if (aw[i].w == x) ox = aw[i].ov;
			if (aw[i].w == v) ov = aw[i].ov;
		}
		if (ov == 0 && ox == 0) continue;
		if (ov > ox) {
			gfa_arc_del(g, w, x, 1);
			gfa_arc_del(g, x^1, w^1, 1);
			++cnt;
		}
	}
	free(a.a);
	if (cnt > 0) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] cut %d small bi-loops\n", cnt);
	return cnt;
}

/******************
 * Bubble popping *
 ******************/

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of reads
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
} binfo_t;

typedef struct {
	binfo_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents
	kvec_t(uint32_t) T; // set of tips
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;

// count the number of outgoing arcs, excluding reduced arcs
static inline int count_out(const gfa_t *g, uint32_t v)
{
	uint32_t i, n, nv = gfa_arc_n(g, v);
	const gfa_arc_t *av = gfa_arc_a(g, v);
	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
}

// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void gfa_bub_backtrack(gfa_t *g, uint32_t v0, buf_t *b)
{
	uint32_t i, v;
	assert(b->S.n == 1);
	for (i = 0; i < b->b.n; ++i)
		g->seg[b->b.a[i]>>1].del = 1;
	for (i = 0; i < b->e.n; ++i) {
		gfa_arc_t *a = &g->arc[b->e.a[i]];
		a->del = 1;
		gfa_arc_del(g, a->w^1, a->v_lv>>32^1, 1);
	}
	v = b->S.a[0];
	do {
		uint32_t u = b->a[v].p; // u->v
		g->seg[v>>1].del = 0;
		gfa_arc_del(g, u, v, 0);
		gfa_arc_del(g, v^1, u^1, 0);
		v = u;
	} while (v != v0);
}

// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static uint64_t gfa_bub_pop1(gfa_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	if (g->seg[v0>>1].del) return 0; // already deleted
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	b->a[v0].c = b->a[v0].d = 0;
	kv_push(uint32_t, b->S, v0);
	do {
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		assert(nv > 0);
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			uint32_t w = av[i].w, l = (uint32_t)av[i].v_lv; // u->w with length l
			binfo_t *t = &b->a[w];
			if (w == v0) goto pop_reset;
			if (av[i].del) continue;
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			if (d + l > max_dist) break; // too far
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				t->p = v, t->s = 1, t->d = d + l;
				t->r = count_out(g, w^1);
				++n_pending;
			} else { // visited before
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				if (d + l < t->d) t->d = d + l; // update dist
			}
			assert(t->r > 0);
			if (--(t->r) == 0) {
				uint32_t x = gfa_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);
	gfa_bub_backtrack(g, v0, b);
	n_pop = 1 | (uint64_t)b->T.n<<32;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}

// pop bubbles
int gfa_pop_bubble(gfa_t *g, int max_dist)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	uint64_t n_pop = 0;
	buf_t b;
	if (!g->is_symm) gfa_symm(g);
	memset(&b, 0, sizeof(buf_t));
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		if (nv < 2 || g->seg[v>>1].del) continue;
		for (i = 0; i < nv; ++i) // gfa_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += gfa_bub_pop1(g, v, max_dist, &b);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) gfa_cleanup(g);
	if (gfa_verbose >= 3) fprintf(stderr, "[M] popped %d bubbles and trimmed %d tips\n", (uint32_t)n_pop, (uint32_t)(n_pop>>32));
	return n_pop;
}

/****************
 * Unitig graph *
 ****************/

#include "kdq.h"
KDQ_INIT(uint64_t)

void gfa_arc_sort(gfa_t *g);
void gfa_arc_index(gfa_t *g);
uint32_t gfa_fix_semi_arc(gfa_t *g);
void gfa_fix_arc_len(gfa_t *g);
uint64_t gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp);

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

		if (g->seg[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
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
			end = x^1, len += l, len_r += a->lw;
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
			start = w, len += l, len_r += a->lw;
			x = w;
		}
		len_r += g->seg[start>>1].len;
add_unitig:
		if (start != UINT32_MAX) mark[start] = mark[end] = 1;
		sprintf(utg_name, "utg%.7d", ug->n_seg + 1);
		tmp = gfa_add_seg(ug, utg_name);
		u = &ug->seg[tmp];
		u->seq = 0, u->len = len;
		u->utg.start = start, u->utg.end = end, u->utg.n = kdq_size(q), u->circ = (start == UINT32_MAX);
		u->utg.m = u->utg.n;
		kroundup32(u->utg.m);
		u->utg.a = (uint64_t*)malloc(8 * u->utg.m);
		u->utg.name = (char**)malloc(sizeof(char*) * u->utg.m);
		u->utg.len2 = len_r;
		for (i = 0; i < kdq_size(q); ++i) {
			u->utg.a[i] = kdq_at(q, i);
			u->utg.name[i] = strdup(g->seg[u->utg.a[i]>>33].name);
		}
	}
	kdq_destroy(uint64_t, q);

	// add arcs between unitigs; reusing mark for a different purpose
	for (v = 0; v < n_vtx; ++v) mark[v] = -1;
	for (i = 0; i < ug->n_seg; ++i) {
		if (ug->seg[i].circ) continue;
		mark[ug->seg[i].utg.start] = i<<1 | 0;
		mark[ug->seg[i].utg.end] = i<<1 | 1;
	}
	for (i = 0; i < g->n_arc; ++i) {
		gfa_arc_t *p = &g->arc[i];
		if (p->del) continue;
		if (mark[p->v_lv>>32^1] >= 0 && mark[p->w] >= 0) {
			gfa_seg_t *s1 = &ug->seg[mark[p->v_lv>>32^1]>>1];
			gfa_seg_t *s2 = &ug->seg[mark[p->w]>>1];
			int ov = p->ov;
			if (ov >= s1->len || ov >= s2->len)
				ov = (s1->len < s2->len? s1->len : s2->len) - 1;
			gfa_add_arc1(ug, mark[p->v_lv>>32^1]^1, mark[p->w], ov, INT32_MAX, -1, 0);
		}
	}
	free(mark);
	gfa_arc_sort(ug);
	gfa_arc_index(ug);
	gfa_fix_semi_arc(ug);
	gfa_fix_arc_len(ug);
	gfa_cleanup(ug);
	return ug;
}
