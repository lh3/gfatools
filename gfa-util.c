#include "gfa.h"
#include "kvec.h"

void gfa_arc_rm(gfa_t *g);

void gfa_sub(gfa_t *g, int n, char *const* seg, int step)
{
	int32_t i;
	kvec_t(uint64_t) stack = {0,0,0};
	for (i = 0; i < n; ++i) {
		int32_t s;
		s = gfa_name2id(g, seg[i]);
		if (s >= 0) {
			kv_push(uint64_t, stack, (uint64_t)(s<<1|0)<<32);
			kv_push(uint64_t, stack, (uint64_t)(s<<1|1)<<32);
		}
	}
	for (i = 0; i < g->n_seg; ++i) // mark all segments to be deleted
		g->seg[i].del = 1;
	while (stack.n) {
		uint64_t x = kv_pop(stack);
		uint32_t v = x>>32, r = (uint32_t)x;
		if (g->seg[v>>1].del == 0) continue;
		g->seg[v>>1].del = 0;
		if (r < step) {
			uint32_t nv = gfa_arc_n(g, v);
			gfa_arc_t *av = gfa_arc_a(g, v);
			for (i = 0; i < nv; ++i)
				if (g->seg[av[i].w>>1].del)
					kv_push(uint64_t, stack, (uint64_t)av[i].w<<32 | (r + 1));
		}
	}
	free(stack.a);
	gfa_arc_rm(g);
}
