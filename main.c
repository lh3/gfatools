#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include "ketopt.h"
#include "gfa-priv.h"

#include <zlib.h>
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

char **gv_read_list(const char *o, int *n_)
{
	int n = 0, m = 0;
	char **s = 0;
	*n_ = 0;
	if (*o != '@') {
		const char *q = o, *p;
		for (p = q;; ++p) {
			if (*p == ',' || *p == 0) {
				if (n == m) {
					m = m? m<<1 : 16;
					s = (char**)realloc(s, m * sizeof(char*));
				}
				s[n++] = gfa_strndup(q, p - q);
				if (*p == 0) break;
				q = p + 1;
			}
		}
	} else {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		int dret;

		fp = gzopen(o + 1, "r");
		if (fp == 0) return 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			if (n == m) {
				m = m? m<<1 : 16;
				s = (char**)realloc(s, m * sizeof(char*));
			}
			s[n++] = gfa_strndup(str.s, p - str.s);
		}
		ks_destroy(ks);
		gzclose(fp);
	}
	if (s) s = (char**)realloc(s, n * sizeof(char*));
	*n_ = n;
	return s;
}

int main_view(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, M_only = 1, step = 0, is_del = 0;
	char *list_arg = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "v:dr:l:", 0)) >= 0) {
		if (c == 'v') gfa_verbose = atoi(o.arg);
		else if (c == 'd') is_del = 1;
		else if (c == 'r') step = atoi(o.arg);
		else if (c == 'l') list_arg = o.arg;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools view [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -v INT        verbose level [%d]\n", gfa_verbose);
		fprintf(stderr, "  -l STR/@FILE  segment list to subset []\n");
		fprintf(stderr, "  -r INT        subset radius (effective with -l) [%d]\n", step);
		fprintf(stderr, "  -d            delete the list of segments (requiring -l; ignoring -r)\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	if (list_arg) {
		int i, n;
		char **list;
		list = gv_read_list(list_arg, &n);
		if (!is_del) {
			gfa_sub(g, n, list, step);
		} else {
			for (i = 0; i < n; ++i) {
				int32_t seg;
				seg = gfa_name2id(g, list[i]);
				if (seg >= 0) gfa_seg_del(g, seg);
			}
		}
		for (i = 0; i < n; ++i) free(list[i]);
		free(list);
	}
	gfa_print(g, stdout, M_only);
	gfa_destroy(g);
	return 0;
}

int main_gfa2bed(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c, merged = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "m", 0)) >= 0)
		if (c == 'm') merged = 1;
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools gfa2bed [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m     merge adjacent intervals on stable sequences\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	if (merged == 0) {
		for (i = 0; i < g->n_seg; ++i) {
			gfa_seg_t *s = &g->seg[i];
			if (s->snid >= 0 && s->soff >= 0)
				printf("%s\t%d\t%d\t%s\n", g->sseq[s->snid].name, s->soff, s->soff + s->len, s->name);
		}
	} else {
		int32_t n_sfa;
		gfa_seg_t *r;
		r = gfa_gfa2sfa(g, &n_sfa, 0);
		for (i = 0; i < n_sfa; ++i) {
			gfa_seg_t *s = &r[i];
			printf("%s\t%d\t%d\n", g->sseq[s->snid].name, s->soff, s->soff + s->len);
		}
		free(r);
	}
	gfa_destroy(g);
	return 0;
}

int main_gfa2fa(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c, is_stable = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "s", 0)) >= 0)
		if (c == 's') is_stable = 1;
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools gfa2fa [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s     output stable sequences (rGFA only)\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	if (is_stable == 0) {
		for (i = 0; i < g->n_seg; ++i) {
			gfa_seg_t *s = &g->seg[i];
			printf(">%s\n%s\n", s->name, s->seq);
		}
	} else {
		int32_t n_sfa;
		gfa_seg_t *r;
		r = gfa_gfa2sfa(g, &n_sfa, 1);
		for (i = 0; i < n_sfa; ++i) {
			gfa_seg_t *s = &r[i];
			if (s->rank == 0)
				printf(">%s\n", g->sseq[s->snid].name);
			else
				printf(">%s_%d_%d\n", g->sseq[s->snid].name, s->soff, s->soff + s->len);
			printf("%s\n", s->seq);
		}
		free(r);
	}
	gfa_destroy(g);
	return 0;
}

int main_asm(int argc, char *argv[])
{
	const char *tr_opts = "v:R:T:B:O:rtbomu";
	ketopt_t o = KETOPT_INIT;
	int c, gap_fuzz = 1000, max_ext = 4, bub_dist = 50000, M_only = 1;
	float ovlp_drop_ratio = .7f;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, tr_opts, 0)) >= 0);
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools asm [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -v INT      verbose level [%d]\n", gfa_verbose);
		fprintf(stderr, "  -u          generate unitig graph (unambiguous merge)\n");
		fprintf(stderr, "  -r          transitive reduction\n");
		fprintf(stderr, "  -R INT      fuzzy length for -r [%d]\n", gap_fuzz);
		fprintf(stderr, "  -t          trim tips\n");
		fprintf(stderr, "  -T INT      tip length for -t [%d]\n", max_ext);
		fprintf(stderr, "  -b          pop bubbles\n");
		fprintf(stderr, "  -B INT      max bubble dist for -b [%d]\n", bub_dist);
		fprintf(stderr, "  -o          drop shorter overlaps\n");
		fprintf(stderr, "  -O FLOAT    dropped/longest<FLOAT, for -o [%g]\n", ovlp_drop_ratio);
		fprintf(stderr, "  -m          misc trimming\n");
		fprintf(stderr, "Note: the order of options matters; one option may be applied >1 times.\n");
		return 1;
	}

	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}

	o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, tr_opts, 0)) >= 0) {
		if (c == 'v') gfa_verbose = atoi(o.arg);
		else if (c == 'R') gap_fuzz = atoi(o.arg);
		else if (c == 'r') gfa_arc_del_trans(g, gap_fuzz);
		else if (c == 'T') max_ext = atoi(o.arg);
		else if (c == 't') gfa_cut_tip(g, max_ext);
		else if (c == 'B') bub_dist = atoi(o.arg);
		else if (c == 'b') gfa_pop_bubble(g, bub_dist);
		else if (c == 'O') ovlp_drop_ratio = atof(o.arg);
		else if (c == 'o') {
			if (gfa_arc_del_short(g, ovlp_drop_ratio) != 0) {
				gfa_cut_tip(g, max_ext);
				gfa_pop_bubble(g, bub_dist);
			}
		} else if (c == 'm') {
			gfa_cut_internal(g, 1);
			gfa_cut_biloop(g, max_ext);
			gfa_cut_tip(g, max_ext);
			gfa_pop_bubble(g, bub_dist);
		} else if (c == 'u') {
			gfa_t *ug;
			ug = gfa_ug_gen(g);
			gfa_destroy(g);
			g = ug;
		}
	}

	gfa_print(g, stdout, M_only);
	gfa_destroy(g);
	return 0;
}

int main(int argc, char *argv[])
{
	extern double realtime(void);
	extern double cputime(void);
	double t_start;
	int ret = 0, i;

	if (argc == 1) {
		fprintf(stderr, "Usage: gfatools <command> <arguments>\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  view        read a GFA file\n");
		fprintf(stderr, "  gfa2fa      convert GFA to FASTA\n");
		fprintf(stderr, "  gfa2bed     convert GFA to BED (requiring rGFA)\n");
		fprintf(stderr, "  asm         miniasm-like graph transformation\n");
		fprintf(stderr, "  version     print version number\n");
		return 1;
	}

	t_start = realtime();
	if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "gfa2bed") == 0) ret = main_gfa2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "gfa2fa") == 0) ret = main_gfa2fa(argc-1, argv+1);
	else if (strcmp(argv[1], "asm") == 0) ret = main_asm(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(GFA_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, GFA_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_start, cputime());
	}
	return ret;
}
