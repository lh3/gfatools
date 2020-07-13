#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ketopt.h"
#include "kvec.h"
#include "paf.h"
#include "miniasm.h"
#include "gfa-priv.h"

#define MA_VERSION "0.3-r166-dirty"

extern double gfa_realtime0;
extern double gfa_realtime(void);
extern double gfa_cputime(void);

static inline int64_t gfa_str2num(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int i, c;
	sdict_t *d;
	ma_hit_t *hit;
	ma_sub_t *sub = 0;
	size_t n_hits;
	float cov = 40.0;
	double max_cut_ratio = 0.9;
	char *fn_reads = 0;

	double int_frac = 0.8, min_iden = 0.0;
	int add_dual = 1, flt = 0, gen_ug = 0, clean = 0, aggre = 0, keep_uni_edge = 0, small_n = 3;
	int min_dp = 2, min_ovlp = 500, min_match = 0, max_hang = 100;

	gfa_t *sg = 0;
	ma_ug_t *ug = 0;

	while ((c = ketopt(&o, argc, argv, 1, "bfh:o:ucUi:an:r:", 0)) >= 0) {
		if (c == 'b') add_dual = 0;
		else if (c == 'f') ++flt;
		else if (c == 'h') max_hang = gfa_str2num(o.arg, 0);
		else if (c == 'o') min_ovlp = gfa_str2num(o.arg, 0);
		else if (c == 'u') gen_ug = 1;
		else if (c == 'n') small_n = atoi(o.arg);
		else if (c == 'U') keep_uni_edge = 1;
		else if (c == 'c') ++clean;
		else if (c == 'a') ++aggre;
		else if (c == 'i') fn_reads = o.arg;
		else if (c == 'r') max_cut_ratio = atof(o.arg);
		else if (c == 'V') {
			printf("%s\n", MA_VERSION);
			return 0;
		}
	}
	if (add_dual) keep_uni_edge = 1;
	gfa_verbose = ma_verbose;
	if (argc == o.ind) {
		fprintf(stderr, "Usage: paf2gfa [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -n INT      threshold for tips and small bubbles [%d]\n", small_n);
		fprintf(stderr, "  -b          both directions of an arc are present in input\n");
		fprintf(stderr, "  -U          keep unidirectional edges (effective with -b)\n");
		fprintf(stderr, "  -f          cut and filter initial hits\n");
		fprintf(stderr, "  -h NUM      max overhang length [%d]\n", max_hang);
		fprintf(stderr, "  -o NUM      min overlap length [%d]\n", min_ovlp);
		fprintf(stderr, "  -c          apply graph cleaning (up to 3)\n");
		fprintf(stderr, "  -r FLOAT    max edge cut ratio (between 0.5 and 1) [%g]\n", max_cut_ratio);
		fprintf(stderr, "  -u          generate unitigs\n");
		fprintf(stderr, "  -i FILE     input reads []\n");
		return 1;
	}

	ma_sys_init();
	d = sd_init();

	hit = ma_hit_read(argv[o.ind], min_ovlp, min_match, d, &n_hits, add_dual, 0);
	if (flt >= 1) {
		sub = ma_hit_sub(min_dp, min_iden, 0, n_hits, hit, d->n_seq);
		n_hits = ma_hit_cut(sub, min_ovlp, n_hits, hit);
		n_hits = ma_hit_flt(sub, (int)(max_hang * 1.5 + .499), (int)(min_ovlp * .5 + .499), n_hits, hit, &cov);
	}
	if (flt >= 1) GFA_REALLOC(hit, n_hits);

	sg = ma_sg_gen(max_hang, int_frac, min_ovlp, d, sub, n_hits, hit);
	if (!keep_uni_edge) gfa_fix_symm_del(sg);
	else gfa_fix_symm_add(sg);
	gfa_arc_del_multi_risky(sg);
	if (clean >= 1) {
		gfa_arc_pair_strong(sg);
		gfa_arc_del_weak(sg);
		gfa_arc_del_trans(sg, 100);
		gfa_drop_tip(sg, 1, INT32_MAX);
	}
	if (clean >= 2) {
		gfa_topocut(sg, 0.3, small_n < 3? small_n : 3, INT32_MAX);
		gfa_drop_tip(sg, 2, INT32_MAX);
		gfa_topocut(sg, 0.5, small_n, INT32_MAX);
		gfa_drop_tip(sg, small_n, INT32_MAX);
		if (max_cut_ratio > 0.5) {
			gfa_topocut(sg, 0.7 < max_cut_ratio? 0.7 : max_cut_ratio, small_n, INT32_MAX);
			gfa_drop_tip(sg, small_n, INT32_MAX);
			gfa_arc_del_short(sg, 2000, 0.5);
			gfa_drop_tip(sg, small_n, INT32_MAX);
			if (max_cut_ratio > 0.7) {
				gfa_topocut(sg, max_cut_ratio, small_n, INT32_MAX);
				gfa_drop_tip(sg, small_n, INT32_MAX);
			}
		}
		gfa_pop_bubble(sg, 1000, small_n, 1);
	}
	if (clean >= 3) {
		gfa_drop_internal(sg, 1);
		gfa_drop_tip(sg, small_n, INT32_MAX);
		gfa_cut_z(sg, 50000, 100000); // this is not well tested
	}
	if (aggre >= 1) {
		gfa_pop_bubble(sg, 100000, INT32_MAX, 1);
		gfa_drop_tip(sg, 100, 100000);
	}

	if (gen_ug) {
		ug = ma_ug_gen(sg);
		if (fn_reads) ma_ug_seq(ug, d, sub, fn_reads);
		ma_ug_print(ug, d, sub, stdout);
	} else {
		ma_sg_print(sg, d, sub, stdout);
	}

	free(sub); free(hit);
	sd_destroy(d);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, gfa_realtime() - gfa_realtime0, gfa_cputime());
	return 0;
}
