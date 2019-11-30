#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ketopt.h"
#include "kvec.h"
#include "paf.h"
#include "miniasm.h"
#include "gfa-priv.h"

#define MA_VERSION "0.3-r179"

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
	char *fn_reads = 0;

	double int_frac = 0.8, min_iden = 0.0;
	int add_dual = 1, flt = 0, gen_ug = 0, clean = 0;
	int min_dp = 2, min_ovlp = 500, min_match = 0, max_hang = 100;

	gfa_t *sg = 0;
	ma_ug_t *ug = 0;

	while ((c = ketopt(&o, argc, argv, 1, "bfh:o:uc", 0)) >= 0) {
		if (c == 'b') add_dual = 0;
		else if (c == 'f') ++flt;
		else if (c == 'h') max_hang = gfa_str2num(o.arg, 0);
		else if (c == 'o') min_ovlp = gfa_str2num(o.arg, 0);
		else if (c == 'u') gen_ug = 1;
		else if (c == 'c') ++clean;
		else if (c == 'V') {
			printf("%s\n", MA_VERSION);
			return 0;
		}
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: paf2gfa [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b          both directions of an arc are present in input\n");
		fprintf(stderr, "  -f          cut and filter initial hits\n");
		fprintf(stderr, "  -h NUM      max overhang length [%d]\n", max_hang);
		fprintf(stderr, "  -o NUM      min overlap length [%d]\n", min_ovlp);
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
	if (clean >= 1) {
		gfa_arc_del_trans(sg, 100);
		gfa_cut_tip(sg, 1, INT32_MAX);
		gfa_topocut(sg, 0.3, 3, INT32_MAX);
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
