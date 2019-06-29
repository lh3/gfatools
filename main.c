#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include "ketopt.h"
#include "gfa.h"

#include <zlib.h>
#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

char **gv_read_list(const char *o, int *n_)
{
	int n = 0, m = 0;
	char **s = 0;
	*n_ = 0;
	if (*o == ',') {
		const char *q = o + 1, *p;
		for (p = q;; ++p) {
			if (*p == ',' || *p == 0) {
				if (n == m) {
					m = m? m<<1 : 16;
					s = (char**)realloc(s, m * sizeof(char*));
				}
				s[n++] = strndup(q, p - q);
				if (*p == 0) break;
				q = p + 1;
			}
		}
	} else {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		int dret;

		fp = gzopen(o, "r");
		if (fp == 0) return 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			if (n == m) {
				m = m? m<<1 : 16;
				s = (char**)realloc(s, m * sizeof(char*));
			}
			s[n++] = strndup(str.s, p - str.s);
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
	int c;
	gfa_t *g;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "v:", 0)) >= 0) {
		if (c == 'v') gfa_verbose = atoi(optarg);
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools view [-v verbose] <in.gfa>\n");
		return 1;
	}
	g = gfa_read(argv[optind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	gfa_print(g, stdout, 1);
	gfa_destroy(g);
	return 0;
}

int main_asm(int argc, char *argv[])
{
	const char *tr_opts = "v:R:T:B:O:rtbom1s:S:d:u";
	int c;
	int gap_fuzz = 1000, max_ext = 4, bub_dist = 50000, M_only = 0, sub_step = 0;
	float ovlp_drop_ratio = .7f;
	gfa_t *g;

	while ((c = getopt(argc, argv, tr_opts)) >= 0);
	if (optind == argc) {
		fprintf(stderr, "Usage: gfaview [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  General:\n");
		fprintf(stderr, "    -v INT      verbose level [%d]\n", gfa_verbose);
		fprintf(stderr, "    -1          only output CIGAR-M operators (for compatibility)\n");
		fprintf(stderr, "    -u          generate unitig graph (unambiguous merge)\n");
		fprintf(stderr, "  Subgraph:\n");
		fprintf(stderr, "    -s EXPR     list of segment names to extract []\n");
		fprintf(stderr, "    -S INT      include neighbors in a radius [%d]\n", sub_step);
		fprintf(stderr, "    -d EXPR     list of segment names to delete []\n");
		fprintf(stderr, "  Graph simplification:\n");
		fprintf(stderr, "    -r          transitive reduction\n");
		fprintf(stderr, "    -R INT      fuzzy length for -r [%d]\n", gap_fuzz);
		fprintf(stderr, "    -t          trim tips\n");
		fprintf(stderr, "    -T INT      tip length for -t [%d]\n", max_ext);
		fprintf(stderr, "    -b          pop bubbles\n");
		fprintf(stderr, "    -B INT      max bubble dist for -b [%d]\n", bub_dist);
		fprintf(stderr, "    -o          drop shorter overlaps\n");
		fprintf(stderr, "    -O FLOAT    dropped/longest<FLOAT, for -o [%g]\n", ovlp_drop_ratio);
		fprintf(stderr, "    -m          misc trimming\n");
		fprintf(stderr, "Note: the order of options matters; one option may be applied >1 times.\n");
		return 1;
	}

	g = gfa_read(argv[optind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	
	optind = 1;
	while ((c = getopt(argc, argv, tr_opts)) >= 0) {
		if (c == 'v') gfa_verbose = atoi(optarg);
		else if (c == '1') M_only = 1;
		else if (c == 'R') gap_fuzz = atoi(optarg);
		else if (c == 'r') gfa_arc_del_trans(g, gap_fuzz);
		else if (c == 'T') max_ext = atoi(optarg);
		else if (c == 't') gfa_cut_tip(g, max_ext);
		else if (c == 'B') bub_dist = atoi(optarg);
		else if (c == 'b') gfa_pop_bubble(g, bub_dist);
		else if (c == 'O') ovlp_drop_ratio = atof(optarg);
		else if (c == 'S') sub_step = atoi(optarg);
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
		} else if (c == 's') {
			int i, n;
			char **s;
			s = gv_read_list(optarg, &n);
			gfa_sub(g, n, s, sub_step);
			for (i = 0; i < n; ++i) free(s[i]);
			free(s);
		} else if (c == 'd') {
			int i, n;
			char **s;
			s = gv_read_list(optarg, &n);
			for (i = 0; i < n; ++i) {
				int32_t seg;
				seg = gfa_name2id(g, s[i]);
				if (seg >= 0) gfa_seg_del(g, seg);
				free(s[i]);
			}
			free(s);
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
		fprintf(stderr, "  view      read a GFA file\n");
		return 1;
	}

	t_start = realtime();
	if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else {
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
