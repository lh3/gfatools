#include "gfa-priv.h"

typedef struct {
	int32_t j, i;
	double sc;
} gt_sc_t;

#define GT_MAX_SC 2

typedef struct {
	int32_t n;
	gt_sc_t s[GT_MAX_SC]; // top 2 scores
} gt_max_t;

static void gfa_genotype_simple_interval(const gfa_t *g, const gfa_sub_t *sub, int32_t jst, int32_t jen)
{
	int32_t j, k, path_len[GT_MAX_SC];
	gt_max_t *sc;
	uint32_t *path[GT_MAX_SC];

	GFA_CALLOC(sc, jen - jst + 1);
	sc[jen - jst].n = 1;
	sc[jen - jst].s[0].j = sc[jen - jst].s[0].i = -1;

	// fill sc[]
	for (j = jen - 1; j >= jst; ++j) {
		const gfa_subv_t *t = &sub->v[j];
		const gfa_aux_t *aux = &g->seg[t->v].aux;
		const uint8_t *dc;
		double svj;
		dc = gfa_aux_get(aux->l_aux, aux->aux, "dc");
		svj = dc && dc[0] == 'f'? *(float*)(dc + 1) : 0.0f;
		gt_max_t *s0 = &sc[j - jst];
		for (k = 0; k < t->n; ++k) {
			uint64_t a = sub->a[t->off + k];
			uint32_t jp = (uint32_t)a;
			const gfa_arc_t *arc;
			int32_t i;
			gt_max_t *s1;
			if (jp <= j) continue;
			arc = &g->arc[a>>32];
			s1 = &sc[jp - jst];
			for (i = 0; i < s1->n; ++i) {
				double score, dc_val;
				aux = &g->link_aux[arc->link_id];
				dc = gfa_aux_get(aux->l_aux, aux->aux, "dc");
				dc_val = dc && dc[0] == 'f'? *(float*)(dc + 1) : 0.0f;
				score = s1->s[i].sc + dc_val + svj;
				if (s0->n < GT_MAX_SC) {
					s0->s[s0->n].j = jp;
					s0->s[s0->n].i = i;
					s0->s[s0->n].sc = score;
					++s0->n;
				} else {
					int32_t x, y;
					for (x = 0; x < s0->n; ++x)
						if (s0->s[x].sc < score)
							break;
					if (x < s0->n) {
						for (y = s0->n - 1; y > x; --y)
							s0->s[y] = s0->s[y - 1];
						s0->s[x].j = jp;
						s0->s[x].i = i;
						s0->s[x].sc = score;
					}
				}
			}
		}
	}

	// backtrack
	for (k = 0; k < sc->n; ++k) {
		int32_t l = 0, j = 0, i = k;
		// count
		while (1) {
			gt_sc_t *s = &sc[j].s[i];
			j = s->j, i = s->i;
			if (j < 0) break;
			++l;
		}
		// fill path[k][]
		path_len[k] = l;
		GFA_MALLOC(path[k], l);
		l = 0, j = 0, i = k;
		while (1) {
			gt_sc_t *s = &sc[j].s[i];
			j = s->j, i = s->i;
			if (j < 0) break;
			path[k][l++] = sub->v[j].v;
		}
	}

	// free
	for (k = 0; k < sc->n; ++k)
		free(path[k]);
	free(sc);
}

void gfa_genotype_simple(const gfa_t *g) // FIXME: doesn't work with translocations
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
		if (vs[i] == (uint32_t)-1) continue;
		sub = gfa_sub_from(0, g, vs[i], 0);
		for (j = 0, jst = 0, max_a = -1; j < sub->n_v; ++j) {
			gfa_subv_t *t = &sub->v[j];
			int32_t k;
			if (j == max_a) {
				const gfa_seg_t *sst = &g->seg[sub->v[jst].v>>1];
				const gfa_seg_t *sen = &g->seg[t->v>>1];
				if (sst->snid == i && sen->snid == i)
					gfa_genotype_simple_interval(g, sub, jst, j);
				max_a = -1, jst = j;
			}
			for (k = 0; k < t->n; ++k)
				if ((int32_t)(sub->a[t->off + k]>>32) > max_a)
					max_a = sub->a[t->off + k]>>32;
		}
		gfa_sub_destroy(sub);
	}
	free(vs);
}
