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

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa32, int32_t, generic_key, 4)

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

/****************
 * Tarjan's SCC *
 ****************/

typedef struct {
	uint32_t index, low;
	uint32_t i;     // index in gfa_sub_t::v[]; a temporary field
	uint32_t start; // starting vertex
	uint16_t stack, out;
} gfa_scinfo_t;

typedef struct gfa_scbuf_t {
	uint32_t index;
	gfa_scinfo_t *a;     // node information
	kvec_t(uint32_t) ts; // Tarjan's stack
	kvec_t(uint64_t) ds; // DFS stack
} gfa_scbuf_t;

gfa_scbuf_t *gfa_scbuf_init(const gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	gfa_scbuf_t *b;
	GFA_CALLOC(b, 1);
	GFA_CALLOC(b->a, n_vtx);
	for (v = 0; v < n_vtx; ++v)
		b->a[v].index = b->a[v].start = (uint32_t)-1;
	return b;
}

void gfa_scbuf_destroy(gfa_scbuf_t *b)
{
	free(b->a); free(b->ts.a); free(b->ds.a); free(b);
}

gfa_sub_t *gfa_scc1(void *km0, const gfa_t *g, gfa_scbuf_t *b, uint32_t v0)
{
	gfa_sub_t *sub;
	uint32_t k, off, m_v = 0;
	void *km;

	km = km_init2(km0, 0x10000);
	KCALLOC(km0, sub, 1);
	sub->km = km0;

	kv_push(uint64_t, b->ds, (uint64_t)v0<<32);
	while (b->ds.n > 0) {
		uint64_t x = kv_pop(b->ds);
		uint32_t i = (uint32_t)x, v = x>>32, nv;
		if (i == 0) { // i is the number of outgoing edges already visited
			b->a[v].low = b->a[v].index = b->index++;
			b->a[v].stack = 1;
			kv_push(uint32_t, b->ts, v);
		}
		nv = gfa_arc_n(g, v);
		if (i == nv) { // done with v
			if (b->a[v].low == b->a[v].index) {
				int32_t i, j = b->ts.n - 1;
				while (b->ts.a[j] != v) --j;
				for (i = b->ts.n - 1; i >= j; --i) {
					uint32_t w = b->ts.a[i];
					if (b->a[w^1].stack == 0 && !b->a[w^1].out) {
						gfa_subv_t *p;
						if (sub->n_v == m_v) KEXPAND(sub->km, sub->v, m_v);
						p = &sub->v[sub->n_v++];
						p->v = w;
						b->a[w].out = 1;
					}
					b->a[w].stack = 0;
				}
				b->ts.n = j;
			}
			if (b->ds.n > 0) { // if call stack is not empty, update the top element
				uint32_t w = v;
				v = b->ds.a[b->ds.n - 1] >> 32;
				b->a[v].low = b->a[v].low < b->a[w].low? b->a[v].low : b->a[w].low;
			}
		} else { // process v's neighbor av[i].w
			gfa_arc_t *av = gfa_arc_a(g, v);
			uint32_t w = av[i].w;
			kv_push(uint64_t, b->ds, (uint64_t)v<<32 | (i+1)); // update the old top of the stack
			if (b->a[w].index == (uint32_t)-1)
				kv_push(uint64_t, b->ds, (uint64_t)w<<32);
			else if (b->a[w].stack)
				b->a[v].low = b->a[v].low < b->a[w].index? b->a[v].low : b->a[w].index;
		}
	}

	// reverse the vertex array
	for (k = 0; k < sub->n_v>>1; ++k) {
		gfa_subv_t x;
		x = sub->v[k], sub->v[k] = sub->v[sub->n_v - k - 1], sub->v[sub->n_v - k - 1] = x;
	}

	// fill other fields in sub
	for (k = 0; k < sub->n_v; ++k)
		b->a[sub->v[k].v].i = k, b->a[sub->v[k].v].start = v0;
	for (k = 0, off = 0; k < sub->n_v; ++k) { // precompute the length of gfa_sub_t::a[]
		uint32_t v = sub->v[k].v;
		int32_t i, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (b->a[av[i].w].start == v0)
				++off;
	}
	sub->n_a = off;
	KCALLOC(sub->km, sub->a, sub->n_a);
	for (k = 0, off = 0; k < sub->n_v; ++k) {
		uint32_t o0, v = sub->v[k].v;
		int32_t i, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0, o0 = off; i < nv; ++i)
			if (b->a[av[i].w].start == v0)
				sub->a[off++] = (uint64_t)b->a[av[i].w].i << 32 | (&av[i] - g->arc);
		sub->v[k].d = 0;
		sub->v[k].off = o0;
		sub->v[k].n = off - o0;
		if (o0 < off) {
			radix_sort_gfa64(&sub->a[o0], &sub->a[off]);
			if (sub->a[o0]>>32 <= k) sub->is_dag = 0;
		}
	}
	return sub;
}

void gfa_scc_all(const gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	gfa_scbuf_t *b;
	b = gfa_scbuf_init(g);
	for (v = 0; v < n_vtx; ++v)
		if (b->a[v].index == (uint32_t)-1 && b->a[v^1].index == (uint32_t)-1) {
			gfa_sub_t *sub;
			sub = gfa_scc1(0, g, b, v);
			gfa_sub_print(stderr, g, sub);
			gfa_sub_destroy(sub);
		}
	gfa_scbuf_destroy(b);
}
