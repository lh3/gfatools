#include <stdio.h>
#include <assert.h>
#include "gfa-priv.h"
#include "khash.h"

typedef struct {
	int32_t j, i, a;
	double sc;
} gt_sc_t;

#define GT_MAX_SC 2

typedef struct {
	int32_t n;
	float vsc;
	gt_sc_t s[GT_MAX_SC]; // top 2 scores
} gt_max_t;

typedef struct {
	int32_t is_arc;
	int32_t x;
	float w;
} gt_node_t;

static inline float get_dc(const gfa_aux_t *aux)
{
	const uint8_t *dc;
	dc = gfa_aux_get(aux->l_aux, aux->aux, "dc");
	return dc && dc[0] == 'f'? *(float*)(dc + 1) : 0.0f;
}

static void gfa_gt_simple_interval(const gfa_t *g, const gfa_sub_t *sub, int32_t jst, int32_t jen)
{
	int32_t j, k, n_path, path_len[GT_MAX_SC+1];
	gt_max_t *sc;
	gt_node_t *path[GT_MAX_SC+1], *pa;
	double score[GT_MAX_SC+1];

	assert(g->seg[sub->v[jst].v>>1].rank == 0);
	assert(g->seg[sub->v[jen].v>>1].rank == 0);
	//fprintf(stderr, "XX\t%s\t%s\n", g->seg[sub->v[jst].v>>1].name, g->seg[sub->v[jen].v>>1].name);

	GFA_CALLOC(sc, jen - jst + 1);
	sc[jen - jst].n = 1;
	sc[jen - jst].s[0].j = sc[jen - jst].s[0].i = -1;
	for (k = 0; k <= GT_MAX_SC; ++k)
		path[k] = 0, path_len[k] = 0, score[k] = 0.0;

	// fill sc[]
	//gfa_sub_print(stderr, g, sub);
	//fprintf(stderr, "n_v=%d,st=%d,en=%d\n", sub->n_v, jst, jen);
	for (j = jen - 1; j >= jst; --j) {
		const gfa_subv_t *t = &sub->v[j];
		gt_max_t *s0 = &sc[j - jst];
		s0->vsc = j == jst? 0.0f : get_dc(&g->seg[t->v>>1].aux);
		for (k = 0; k < t->n; ++k) { // iterate over neighbors
			uint64_t a = sub->a[t->off + k];
			uint32_t jp = (uint32_t)(a>>32);
			int32_t i;
			gt_max_t *s1;
			float dc;
			assert(t->off + k < sub->n_a);
			//fprintf(stderr, "j=%d off=%d k=%d jp=%d\n", j, t->off, k, jp);
			if (jp <= j || jp > jen) continue; // ignore cycles or arcs outside this subgraph
			dc = get_dc(&g->link_aux[g->arc[(uint32_t)a].link_id]);
			s1 = &sc[jp - jst];
			for (i = 0; i < s1->n; ++i) { // iterate over path ends
				double score;
				score = s1->s[i].sc + dc + s0->vsc;
				if (s0->n < GT_MAX_SC) {
					s0->s[s0->n].j = jp;
					s0->s[s0->n].i = i;
					s0->s[s0->n].a = t->off + k;
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
						s0->s[x].a = t->off + k;
						s0->s[x].sc = score;
					}
				}
			}
		}
		//fprintf(stderr, "j=%d v=%c%s[%d] n_best=%d: ", j, "><"[t->v&1], g->seg[t->v>>1].name, t->v, s0->n); for (k = 0; k < s0->n; ++k) fprintf(stderr, "%f,", s0->s[k].sc); fprintf(stderr, "\n");
	}

	// allocate path[]
	n_path = sc->n + 1;
	for (k = 0; k < n_path; ++k) // k==sc->n for the reference path
		GFA_MALLOC(path[k], (jen - jst + 1) * 2); // this is over-allocating, but it should not be an issue

	// find the reference path
	j = jst, k = 0, pa = path[0];
	while (1) {
		const gfa_subv_t *t = &sub->v[j];
		uint64_t a;
		int32_t i;
		float w;
		for (i = 0; i < t->n; ++i) {
			a = sub->a[t->off + i];
			if (g->arc[(uint32_t)a].rank == 0)
				break;
		}
		assert(i < t->n);
		j = a>>32;
		w = get_dc(&g->link_aux[g->arc[(uint32_t)a].link_id]);
		pa[k].is_arc = 1, pa[k].x = t->off + i, pa[k].w = w, ++k;
		if (j == jen) break;
		w = get_dc(&g->seg[sub->v[j].v>>1].aux);
		pa[k].is_arc = 0, pa[k].x = j, pa[k].w = w, ++k;
	}
	path_len[0] = k;

	// backtrack
	for (k = 0; k < sc->n; ++k) {
		int32_t l = 0, j = jst, i = k;
		pa = path[k + 1];
		while (1) {
			gt_sc_t *s = &sc[j - jst].s[i];
			j = s->j, i = s->i;
			pa[l].is_arc = 1, pa[l].x = s->a, pa[l].w = 0.0, ++l;
			if (j < 0 || j == jen) break;
			pa[l].is_arc = 0, pa[l].x = j, pa[l].w = 0.0, ++l;
		}
		path_len[k+1] = l;
	}

	// update weight and score[]
	for (k = 0; k < n_path; ++k) {
		double s = 0.0;
		for (j = 0; j < path_len[k]; ++j) {
			gt_node_t *p = &path[k][j];
			if (p->is_arc) p->w = get_dc(&g->link_aux[g->arc[(uint32_t)sub->a[p->x]].link_id]);
			else p->w = get_dc(&g->seg[sub->v[p->x].v>>1].aux);
			s += p->w;
		}
		score[k] = s;
	}

	// squeeze out ALT ref paths
	for (k = 1, j = 1; k < n_path; ++k) {
		int32_t is_drop = 0;
		if (path_len[k] == path_len[0]) {
			int32_t i;
			for (i = 0; i < path_len[0]; ++i)
				if (path[k][i].x != path[0][i].x)
					break;
			is_drop = (i == path_len[0]);
		}
		if ((float)score[k] + 1.0f == 1.0f)
			is_drop = 1;
		if (is_drop == 0)
			path_len[j] = path_len[k], path[j] = path[k], score[j++] = score[k];
		else free(path[k]);
	}
	n_path = j;

	// print
	printf("VP\t%c%s\t%c%s\t%d", "><"[sub->v[jst].v&1], g->seg[sub->v[jst].v>>1].name, "><"[sub->v[jen].v&1], g->seg[sub->v[jen].v>>1].name, n_path);
	for (k = 0; k < n_path; ++k) {
		int32_t i;
		printf("\t%.2f:", score[k]);
		for (i = 0; i < path_len[k]; ++i) {
			if (!path[k][i].is_arc) {
				uint32_t v = sub->v[path[k][i].x].v;
				printf("%c%s", "><"[v&1], g->seg[v>>1].name);
			}
		}
	}
	printf("\n");

	// free
	for (k = 0; k < n_path; ++k)
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
					gfa_gt_simple_interval(g, sub, jst, j);
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
