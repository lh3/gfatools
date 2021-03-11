#include <assert.h>
#include <stdio.h>
#include "gfa-priv.h"
#include "kalloc.h"
#include "kavl.h"
#include "khash.h"
#include "ksort.h"
#include "kvec.h"

/*********************************************
 * Extract a subgraph starting from a vertex *
 *********************************************/

typedef struct tnode_s {
	uint64_t xnd;
	uint32_t v, in_tree:31, forced:1;
	KAVL_HEAD(struct tnode_s) head;
} tnode_t, *tnode_p;

#define tn_n(p) ((uint32_t)((p)->xnd<<1>>33))
#define tn_lt(a, b) ((a)->xnd < (b)->xnd || ((a)->xnd == (b)->xnd && (a)->v < (b)->v))
#define tn_cmp(a, b) (tn_lt(b, a) - tn_lt(a, b))

KAVL_INIT(v, tnode_t, head, tn_cmp)
KHASH_MAP_INIT_INT(v, tnode_p)

static inline tnode_t *gen_tnode(void *km, const gfa_t *g, uint32_t v, int32_t d)
{
	tnode_t *p;
	KMALLOC(km, p, 1);
	p->v = v, p->in_tree = 1, p->forced = 0;
	p->xnd = 1ULL<<63 | (uint64_t)gfa_arc_n(g, v^1) << 32 | d;
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

		#if 0
		kavl_itr_first(v, root, &itr);
		fprintf(stderr, "PEEK");
		while ((r = kavl_at(&itr)) != 0) {
			fprintf(stderr, " %c%s:x=%d:n=%d:d=%d", "><"[r->v&1], g->seg[r->v>>1].name, (int)(r->xnd>>63), tn_n(r), (uint32_t)r->xnd);
			if (kavl_itr_next(v, &itr) == 0) break;
		}
		fputc('\n', stderr);
		#endif
		kavl_itr_first(v, root, &itr);
		r = kavl_at(&itr);
		if (orphan_inv) { // then prioritize on vertices whose complements have been moved out of the tree
			while ((r = kavl_at(&itr)) != 0) {
				k = kh_get(v, h, r->v^1);
				if (k != kh_end(h) && !kh_val(h, k)->in_tree) {
					--orphan_inv;
					q = kavl_erase(v, &root, r, 0);
					break;
				}
				if (kavl_itr_next(v, &itr) == 0) break;
			}
		} else if (tn_n(r) > 0) { // FIXME: be careful of the worst-case time complexity!
			int n = 0;
			nv = gfa_arc_n(g, r->v^1);
			av = gfa_arc_a(g, r->v^1);
			for (i = 0; i < nv; ++i) {
				gfa_arc_t *avi = &av[i];
				khint_t k1, k2;
				k1 = kh_get(v, h, avi->w^1);
				k2 = kh_get(v, h, avi->w);
				if ((k1 == kh_end(h) && k2 != kh_end(h) && !kh_val(h, k2)->in_tree) || (k2 == kh_end(h) && k1 != kh_end(h) && !kh_val(h, k1)->in_tree))
					++n;
				else break;
			}
			if (i < nv) {
				tnode_p *a;
				KMALLOC(km, a, kavl_size(head, root));
				n = 0;
				while ((r = kavl_at(&itr)) != 0) {
					a[n++] = p = (tnode_t*)r;
					p->xnd &= ~(1ULL<<63);
					if (kavl_itr_next(v, &itr) == 0) break;
				}
				root = 0;
				for (i = 0; i < n; ++i)
					kavl_insert(v, &root, a[i], 0);
				kfree(km, a);
			}
		}
		if (q == 0) q = kavl_erase_first(v, &root); // take out the "smallest" vertex
		q->forced = (tn_n(q) > 0 || q->xnd>>63 == 0);
		q->in_tree = 0;
		if (n_L == m_L) KEXPAND(km, L, m_L);
		L[n_L++] = q;

		k = kh_get(v, h, q->v^1);
		if (k != kh_end(h) && kh_val(h, k)->in_tree)
			++orphan_inv;
		//fprintf(stderr, "OUT vertex:%c%s[%u], remained:%d, orphan_inv:%d\n", "><"[q->v&1], g->seg[q->v>>1].name, q->v, kavl_size(head, root), orphan_inv);

		d = (uint32_t)q->xnd;
		nv = gfa_arc_n(g, q->v);
		av = gfa_arc_a(g, q->v);
		for (i = 0; i < nv; ++i) {
			gfa_arc_t *avi = &av[i];
			int32_t dt = d + g->seg[avi->w>>1].len;
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
				if (dt < (uint32_t)p->xnd)
					p->xnd = p->xnd>>32<<32 | dt;
			}
			assert(tn_n(p) > 0);
			p->xnd -= 1ULL<<32;
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
		sub->v[j].d = (uint32_t)L[j]->xnd;
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
