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

static inline int64_t gfa_str2num(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	*q = p;
	return (int64_t)(x + .499);
}

int main_view(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, out_flag = 0, step = 0, is_del = 0, fix_multi = 0;
	char *list_arg = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "v:dr:l:SM", 0)) >= 0) {
		if (c == 'v') gfa_verbose = atoi(o.arg);
		else if (c == 'd') is_del = 1;
		else if (c == 'r') step = atoi(o.arg);
		else if (c == 'l') list_arg = o.arg;
		else if (c == 'S') out_flag |= GFA_O_NO_SEQ;
		else if (c == 'M') fix_multi = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools view [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -v INT        verbose level [%d]\n", gfa_verbose);
		fprintf(stderr, "  -l STR/@FILE  segment list to subset []\n");
		fprintf(stderr, "  -r INT        subset radius (effective with -l) [%d]\n", step);
		fprintf(stderr, "  -d            delete the list of segments (requiring -l; ignoring -r)\n");
		fprintf(stderr, "  -M            remove multiple edges\n");
		fprintf(stderr, "  -S            don't print sequences\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	if (fix_multi) gfa_fix_multi(g);
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
	gfa_print(g, stdout, out_flag);
	gfa_destroy(g);
	return 0;
}

int main_stat(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, max_deg = 0;
	int64_t i, n_vtx;
	uint64_t tot_seg_len = 0, seg0_len = 0, n_link = 0, tot_deg = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools stat <in.gfa>\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	printf("Number of segments: %d\n", g->n_seg);
	for (i = 0; i < g->n_arc; ++i)
		if (!g->arc[i].comp) ++n_link;
	printf("Number of links: %lld\n", (long long)n_link);
	printf("Number of arcs: %lld\n", (long long)g->n_arc);
	printf("Max rank: %d\n", g->max_rank);
	for (i = 0; i < g->n_seg; ++i) {
		tot_seg_len += g->seg[i].len;
		if (g->seg[i].rank == 0) seg0_len += g->seg[i].len;
	}
	printf("Total segment length: %lld\n", (long long)tot_seg_len);
	if (g->n_seg)
		printf("Average segment length: %.3f\n", (double)tot_seg_len / g->n_seg);
	printf("Sum of rank-0 segment lengths: %lld\n", (long long)seg0_len);
	n_vtx = gfa_n_vtx(g);
	for (i = 0; i < n_vtx; ++i) {
		int32_t nv = gfa_arc_n(g, i);
		if (nv > max_deg) max_deg = nv;
		tot_deg += nv;
	}
	printf("Max degree: %d\n", max_deg);
	if (n_vtx > 0)
		printf("Average degree: %.3f\n", (double)tot_deg / n_vtx);
	gfa_destroy(g);
	return 0;
}

int main_gfa2bed(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c, merged = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "s", 0)) >= 0)
		if (c == 's') merged = 1;
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools gfa2bed [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s     merge adjacent intervals on stable sequences\n");
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
		int32_t j, n_sfa;
		gfa_sfa_t *r;
		r = gfa_gfa2sfa(g, &n_sfa, 0);
		for (i = 0; i < n_sfa; ++i) {
			gfa_sfa_t *s = &r[i];
			printf("%s\t%d\t%d", g->sseq[s->snid].name, s->soff, s->soff + s->len);
			if (s->rank > 0) {
				for (j = 0; j < 2; ++j) {
					if (s->end[j] == (uint64_t)-1) printf("\t*\t*\t*");
					else printf("\t%c\t%s\t%d", "><"[s->end[j]&1], g->sseq[s->end[j]>>32].name, (uint32_t)s->end[j]>>1);
				}
			} else printf("\t*\t*\t*\t*\t*\t*");
			putchar('\n');
		}
		free(r);
	}
	gfa_destroy(g);
	return 0;
}

static void print_seq(FILE *fp, const char *seq, int32_t line_len)
{
	if (line_len <= 0) {
		fputs(seq, fp);
	} else {
		int32_t i, l;
		l = strlen(seq);
		for (i = 0; i < l; i += line_len) {
			if (i) fputc('\n', fp);
			if (i + line_len < l) fwrite(&seq[i], 1, line_len, fp);
			else fputs(&seq[i], fp);
		}
	}
	fputc('\n', fp);
}

int main_gfa2fa(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c, is_stable = 0, skip_rank0 = 0, line_len = 60, ref_only = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "sPl:0", 0)) >= 0) {
		if (c == 's') is_stable = 1;
		else if (c == 'P') skip_rank0 = is_stable = 1;
		else if (c == 'l') line_len = atoi(o.arg);
		else if (c == '0') ref_only = is_stable = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools gfa2fa [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT   line length [%d]\n", line_len);
		fprintf(stderr, "  -s       output stable sequences (rGFA only)\n");
		fprintf(stderr, "  -P       skip rank-0 sequences (rGFA only; force -s)\n");
		fprintf(stderr, "  -0       only output rank-0 sequences (rGFA only; force -s)\n");
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
			printf(">%s\n", s->name);
			print_seq(stdout, s->seq, line_len);
		}
	} else {
		int32_t j, n_sfa;
		gfa_sfa_t *r;
		r = gfa_gfa2sfa(g, &n_sfa, 1);
		for (i = 0; i < n_sfa; ++i) {
			gfa_sfa_t *s = &r[i];
			if (s->rank == 0) {
				if (!skip_rank0) {
					printf(">%s\n", g->sseq[s->snid].name);
					print_seq(stdout, s->seq, line_len);
				}
			} else if (!ref_only) {
				printf(">%s_%d_%d", g->sseq[s->snid].name, s->soff, s->soff + s->len);
				for (j = 0; j < 2; ++j) {
					if (s->end[j] == (uint64_t)-1) printf("\t*");
					else printf("\t%c%s:%d", "><"[s->end[j]&1], g->sseq[s->end[j]>>32].name, (uint32_t)s->end[j]>>1);
				}
				putchar('\n');
				print_seq(stdout, s->seq, line_len);
			}
			free(s->seq);
		}
		free(r);
	}
	gfa_destroy(g);
	return 0;
}

int main_blacklist(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, j, c, n_bb, min_len = 100;
	gfa_t *g;
	gfa_bubble_t *bb;

	while ((c = ketopt(&o, argc, argv, 1, "l:", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools blacklist [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    min region length [%d]\n", min_len);
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	bb = gfa_bubble(g, &n_bb);
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		int32_t rst = b->ss, ren = b->se;
		if (ren - rst < min_len) {
			int32_t ext = (min_len - (ren - rst) + 1) / 2;
			rst -= ext, ren += ext;
			if (rst < 0) rst = 0;
			if (ren > g->sseq[b->snid].max)
				ren = g->sseq[b->snid].max;
		}
		printf("%s\t%d\t%d\t%d\t", g->sseq[b->snid].name, rst, ren, b->n_seg);
		for (j = 0; j < b->n_seg; ++j) {
			if (j) fputc(',', stdout);
			printf("%s", g->seg[b->v[j]>>1].name);
		}
		putchar('\n');
		free(b->v);
	}
	free(bb);
	gfa_destroy(g);
	return 0;
}

int main_bubble(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, j, c, n_bb;
	gfa_t *g;
	gfa_bubble_t *bb;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools bubble <in.gfa>\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	bb = gfa_bubble(g, &n_bb);
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%.3g\t%.3g\t", g->sseq[b->snid].name, b->ss, b->se, b->n_seg,
			   b->len_min, b->len_max, b->cf_min, b->cf_max);
		for (j = 0; j < b->n_seg; ++j) {
			if (j) fputc(',', stdout);
			printf("%s", g->seg[b->v[j]>>1].name);
		}
		if (b->len_min == 0) printf("\t*");
		else printf("\t%s", b->seq_min);
		printf("\t%s\n", b->seq_max);
		free(b->v);
	}
	free(bb);
	gfa_destroy(g);
	return 0;
}

int main_gt(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, is_path = 0;
	float min_dc = 5.0f;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "d:p", 0)) >= 0) {
		if (c == 'd') min_dc = atof(o.arg);
		else if (c == 'p') is_path = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools gt [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -d FLOAT      min depth [%g]\n", min_dc);
		fprintf(stderr, "  -p            print path instead of allele sequences (for debugging)\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	gfa_gt_simple_print(g, min_dc, is_path);
	gfa_destroy(g);
	return 0;
}

int main_asm(int argc, char *argv[])
{
	const char *tr_opts = "v:ur:t:b:B:s:o:c:";
	ketopt_t o = KETOPT_INIT;
	int c, oflag = 0;
	gfa_t *g;
	char *p;

	while ((c = ketopt(&o, argc, argv, 1, tr_opts, 0)) >= 0);
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools asm [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -r INT          transitive reduction (fuzzy length)\n");
		fprintf(stderr, "  -t INT1[,INT2]  cut tips (tip seg count, tip length [inf])\n");
		fprintf(stderr, "  -b INT          pop bubbles along with small tips (max dist)\n");
		fprintf(stderr, "  -B INT          pop bubbles but protect small tips (max dist)\n");
		fprintf(stderr, "  -s INT          pop simple bubbles (max seg count)\n");
		fprintf(stderr, "  -o FLOAT[,INT]  cut short overlaps (ratio to the longest overlap, overlap length [0])\n");
		fprintf(stderr, "  -c FLOAT[,INT1[,INT2]]\n");
		fprintf(stderr, "                  cut overlaps, topology aware (ratio, tip seg count [3], tip length [inf])\n");
		fprintf(stderr, "  -u              generate unitigs\n");
		fprintf(stderr, "  -v INT          verbose level [%d]\n", gfa_verbose);
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
		if (c == 'v') {
			gfa_verbose = atoi(o.arg);
		} else if (c == 'u') {
			gfa_t *ug;
			ug = gfa_ug_gen(g);
			gfa_destroy(g);
			g = ug;
		} else if (c == 'r') {
			int32_t fuzz;
			fuzz = gfa_str2num(o.arg, &p);
			gfa_arc_del_trans(g, fuzz);
		} else if (c == 't') {
			int32_t max_ext, max_len = INT32_MAX;
			max_ext = gfa_str2num(o.arg, &p);
			if (*p == ',') max_len = gfa_str2num(p + 1, &p);
			gfa_cut_tip(g, max_ext, max_len);
		} else if (c == 'b') {
			int32_t dist;
			dist = gfa_str2num(o.arg, &p);
			gfa_pop_bubble(g, dist, 0);
		} else if (c == 'B') {
			int32_t dist;
			dist = gfa_str2num(o.arg, &p);
			gfa_pop_bubble(g, dist, 1);
		} else if (c == 's') {
			int32_t min_side, max_side = 20;
			min_side = gfa_str2num(o.arg, &p);
			if (*p == ',') max_side = gfa_str2num(p + 1, &p);
			gfa_bub_simple(g, min_side, max_side);
		} else if (c == 'o') {
			double ratio;
			int32_t min_len = 0;
			ratio = strtod(o.arg, &p);
			if (*p == ',') min_len = gfa_str2num(p + 1, &p);
			gfa_arc_del_short(g, min_len, ratio);
		} else if (c == 'c') {
			double ratio;
			int32_t tip_cnt = 3, tip_len = INT32_MAX;
			ratio = strtod(o.arg, &p);
			if (*p == ',') tip_cnt = gfa_str2num(p + 1, &p);
			if (*p == ',') tip_len = gfa_str2num(p + 1, &p);
			gfa_topocut(g, ratio, tip_cnt, tip_len);
		}
	}

	gfa_print(g, stdout, oflag);
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
		fprintf(stderr, "  stat        statistics about a GFA file\n");
		fprintf(stderr, "  gfa2fa      convert GFA to FASTA\n");
		fprintf(stderr, "  gfa2bed     convert GFA to BED (requiring rGFA)\n");
		fprintf(stderr, "  blacklist   blacklist regions\n");
		fprintf(stderr, "  gt          genotype from the \"dc\" tag (requring rGFA)\n");
		fprintf(stderr, "  asm         miniasm-like graph transformation\n");
		fprintf(stderr, "  version     print version number\n");
		return 1;
	}

	t_start = realtime();
	if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "stat") == 0) ret = main_stat(argc-1, argv+1);
	else if (strcmp(argv[1], "gfa2bed") == 0) ret = main_gfa2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "gfa2fa") == 0) ret = main_gfa2fa(argc-1, argv+1);
	else if (strcmp(argv[1], "blacklist") == 0) ret = main_blacklist(argc-1, argv+1);
	else if (strcmp(argv[1], "bubble") == 0) ret = main_bubble(argc-1, argv+1);
	else if (strcmp(argv[1], "gt") == 0) ret = main_gt(argc-1, argv+1);
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
