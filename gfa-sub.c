#include <assert.h>
#include <stdio.h>
#include "gfa-priv.h"
#include "kalloc.h"
#include "kavl.h"
#include "khash.h"
#include "ksort.h"

/*********************************************
 * Extract a subgraph starting from a vertex *
 *********************************************/

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa32, int32_t, generic_key, 4)

typedef struct tnode_s {
	uint64_t nd;
	uint32_t v, in_tree:31, forced:1;
	KAVL_HEAD(struct tnode_s) head;
} tnode_t;

typedef tnode_t *tnode_p;

#define tn_lt(a, b) ((a)->nd < (b)->nd || ((a)->nd == (b)->nd && (a)->v < (b)->v))
#define tn_cmp(a, b) (tn_lt(b, a) - tn_lt(a, b))

KAVL_INIT(v, tnode_t, head, tn_cmp)
KHASH_MAP_INIT_INT(v, tnode_p)

static inline tnode_t *gen_tnode(void *km, const gfa_t *g, uint32_t v, int32_t d)
{
	tnode_t *p;
	KMALLOC(km, p, 1);
	p->v = v, p->in_tree = 1, p->forced = 0;
	p->nd = (uint64_t)gfa_arc_n(g, v^1) << 32 | d;
	return p;
}

/* Extract a subgraph extended from a vertex within a radius. If the subgraph
 * is DAG, vertices are in the topological sorting order. The algorithm is
 * modified from Kahn's algorithm.
 */
gfa_sub_t *gfa_sub_from(void *km0, const gfa_t *g, uint32_t v0, int32_t max_dist)
{
	void *km;
	tnode_t *p, *root = 0, **L = 0;
	khash_t(v) *h;
	khint_t k;
	int32_t j, n_L = 0, m_L = 0, n_arc = 0, m_arc = 0, off, n_bidir = 0, orphan_inv = 0;
	int absent;
	gfa_sub_t *sub = 0;

	km = km_init2(km0, 0x10000);
	h = kh_init2(v, km);

	k = kh_put(v, h, v0, &absent);
	p = kh_val(h, k) = gen_tnode(km, g, v0, 0);
	kavl_insert(v, &root, p, 0);

	while (kavl_size(head, root) > 0) {
		tnode_t *q = 0;
		int32_t i, nv, d;
		gfa_arc_t *av;
		const tnode_t *r;
		kavl_itr_t(v) itr;

		kavl_itr_first(v, root, &itr);
		r = kavl_at(&itr);
		if (r->nd>>32 > 0 || orphan_inv) {
			while ((r = kavl_at(&itr)) != 0) {
				if (kh_get(v, h, r->v^1) != kh_end(h)) {
					q = kavl_erase(v, &root, r, 0);
					break;
				}
				if (kavl_itr_next(v, &itr) == 0) break;
			}
			orphan_inv = 0;
		}
		if (q == 0) q = kavl_erase_first(v, &root); // take out the "smallest" vertex
		fprintf(stderr, "OUT vertex:%c%s[%u], remained:%d, orphan_inv:%d\n", "><"[q->v&1], g->seg[q->v>>1].name, q->v, kavl_size(head, root), orphan_inv);
		q->forced = (q->nd >> 32 > 0);
		q->in_tree = 0;
		if (n_L == m_L) KEXPAND(km, L, m_L);
		L[n_L++] = q;

		k = kh_get(v, h, q->v^1);
		if (k != kh_end(h) && kh_val(h, k)->in_tree)
			orphan_inv = 1;

		d = (uint32_t)q->nd;
		nv = gfa_arc_n(g, q->v);
		av = gfa_arc_a(g, q->v);
		for (i = 0; i < nv; ++i) {
			gfa_arc_t *avi = &av[i];
			int32_t dt = d + (uint32_t)avi->v_lv;
			if (max_dist > 0 && dt > max_dist) continue;
			k = kh_get(v, h, avi->w^1);
			if (k != kh_end(h) && !kh_val(h, k)->in_tree && !kh_val(h, k)->forced) {
				++n_bidir;
				continue;
			}
			++n_arc;
			k = kh_put(v, h, avi->w, &absent);
			if (absent) { // a vertex that hasn't been visited before
				p = kh_val(h, k) = gen_tnode(km, g, avi->w, dt);
			} else { // visited before; then update the info
				p = kh_val(h, k);
				if (!p->in_tree) continue; // when there is a cycle, a vertex may be added to L[] already
				kavl_erase(v, &root, p, 0);
				if (dt < (uint32_t)p->nd)
					p->nd = p->nd>>32<<32 | dt;
			}
			assert(p->nd>>32 > 0);
			p->nd -= 1ULL<<32;
			kavl_insert(v, &root, p, 0); // insert/re-insert to the tree
		}
	}
	assert(kh_size(h) == n_L);

	KCALLOC(km0, sub, 1);
	sub->km = km0;
	sub->n_v = n_L;
	sub->n_a = n_arc;
	KCALLOC(sub->km, sub->v, n_L);
	KCALLOC(sub->km, sub->a, n_arc);
	m_arc = n_arc;
	sub->is_dag = 1;

	for (j = 0; j < n_L; ++j) L[j]->in_tree = j; // reuse ->in_tree for a different purpose
	for (j = 0, off = 0; j < sub->n_v; ++j) {
		int32_t i, nv, o0 = off;
		gfa_arc_t *av;
		nv = gfa_arc_n(g, L[j]->v);
		av = gfa_arc_a(g, L[j]->v);
		for (i = 0; i < nv; ++i) {
			gfa_arc_t *avi = &av[i];
			k = kh_get(v, h, avi->w);
			if (k == kh_end(h)) continue;
			if (off == m_arc) KEXPAND(sub->km, sub->a, m_arc);
			sub->a[off++] = (uint64_t)kh_val(h, k)->in_tree << 32 | (avi - g->arc);
		}
		sub->v[j].v = L[j]->v;
		sub->v[j].d = (uint32_t)L[j]->nd;
		sub->v[j].off = o0;
		sub->v[j].n = off - o0;
		if (o0 < off) {
			radix_sort_gfa64(&sub->a[o0], &sub->a[off]);
			if (sub->a[o0]>>32 <= j) sub->is_dag = 0;
		}
	}
	if (off != n_arc) {
		assert(n_bidir > 0); // off != n_arc should only happen when n_bidir>0
		fprintf(stderr, "[W::%s] unusual bubble chain starting at %c%s: off=%d, n_arc=%d, n_bidir=%d\n", __func__, "><"[v0&1], g->seg[v0>>1].name, off, n_arc, n_bidir);
	}

	km_destroy(km);
	//gfa_sub_print(stderr, g, sub);
	return sub;
}

void gfa_sub_destroy(gfa_sub_t *sub)
{
	void *km;
	if (sub == 0) return;
	km = sub->km;
	kfree(km, sub->v); kfree(km, sub->a); kfree(km, sub);
}

void gfa_sub_print(FILE *fp, const gfa_t *g, const gfa_sub_t *sub)
{
	int32_t i, j;
	for (i = 0; i < sub->n_v; ++i) {
		gfa_subv_t *p = &sub->v[i];
		fprintf(fp, "[%d]\t%d\t%c%s\t%d\t%d", i, p->v, "><"[p->v&1], g->seg[p->v>>1].name, p->d, p->n);
		if (p->n > 0) {
			fputc('\t', fp);
			for (j = 0; j < p->n; ++j) {
				if (j) fputc(',', fp);
				fprintf(fp, "%d", (uint32_t)(sub->a[p->off + j]>>32));
			}
		}
		fputc('\n', fp);
	}
}
