#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include "ketopt.h"
#include "gfa-priv.h"
#include "kalloc.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define GFATOOLS_VERSION "0.5-r275-dirty"

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
	int c, out_flag = 0, step = 0, is_del = 0, fix_multi = 0, flip_walk = 0;
	char *list_arg = 0, *reg_arg = 0;
	gfa_t *g, *f = 0;

	while ((c = ketopt(&o, argc, argv, 1, "v:dr:l:SMR:w", 0)) >= 0) {
		if (c == 'v') gfa_verbose = atoi(o.arg);
		else if (c == 'd') is_del = 1;
		else if (c == 'r') step = atoi(o.arg);
		else if (c == 'l') list_arg = o.arg;
		else if (c == 'S') out_flag |= GFA_O_NO_SEQ;
		else if (c == 'R') reg_arg = o.arg;
		else if (c == 'M') fix_multi = 1;
		else if (c == 'w') flip_walk = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools view [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -v INT        verbose level [%d]\n", gfa_verbose);
		fprintf(stderr, "  -l STR/@FILE  segment list to subset []\n");
		fprintf(stderr, "  -R STR        a region like chr1:101-200 (a 1-based closed region) []\n");
		fprintf(stderr, "  -r INT        subset radius (effective with -l) [%d]\n", step);
		fprintf(stderr, "  -d            delete the list of segments (requiring -l; ignoring -r)\n");
		fprintf(stderr, "  -M            remove multiple edges\n");
		fprintf(stderr, "  -S            don't print sequences\n");
		fprintf(stderr, "  -w            flip walk\n");
		return 1;
	}
	if (list_arg && reg_arg) {
		fprintf(stderr, "ERROR: -R and -l can't be used at the same time\n");
		return 3;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	if (fix_multi) gfa_fix_multi(g);
	if (list_arg || reg_arg) {
		int32_t i, n_seg, *seg;
		if (list_arg) {
			char **list;
			int32_t n;
			list = gv_read_list(list_arg, &n);
			seg = gfa_list2seg(g, n, list, &n_seg);
			for (i = 0; i < n; ++i) free(list[i]);
			free(list);
		} else {
			gfa_bubble_t *bb;
			int32_t n_bb;
			gfa_sort_ref_arc(g);
			bb = gfa_bubble(g, &n_bb);
			seg = gfa_query_by_reg(g, n_bb, bb, reg_arg, &n_seg);
			for (i = 0; i < n_bb; ++i) free(bb[i].v);
			free(bb);
		}
		if (n_seg == 0) goto end_view; // nothing to extract
		if (step > 0) {
			int32_t n_seg0 = n_seg, *seg0 = seg;
			seg = gfa_sub_extend(g, n_seg0, seg0, step, &n_seg);
			free(seg0);
		}
		if (is_del) {
			int8_t *flag;
			int32_t n, *seg_rev;
			GFA_CALLOC(flag, g->n_seg);
			for (i = 0; i < n_seg; ++i)
				if (seg[i] < g->n_seg)
					flag[seg[i]] = 1;
			for (i = 0, n = 0; i < g->n_seg; ++i)
				if (flag[i]) ++n;
			GFA_CALLOC(seg_rev, n);
			for (i = 0, n = 0; i < g->n_seg; ++i)
				if (flag[i]) seg_rev[n++] = i;
			free(flag);
			free(seg);
			seg = seg_rev, n_seg = n;
		}
		f = gfa_subview2(g, n_seg, seg, 1);
		free(seg);
	}
	if (flip_walk) gfa_walk_flip(f);
	gfa_print(f, stdout, out_flag);
	gfa_subview_destroy(f);
end_view:
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
				printf("%s\t%d\t%d\t%s\t%d\n", g->sseq[s->snid].name, s->soff, s->soff + s->len, s->name, s->rank);
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
	int32_t i, j, c, n_bb, min_len = 100, with_bidir = 0;
	gfa_t *g;
	gfa_bubble_t *bb;

	while ((c = ketopt(&o, argc, argv, 1, "l:b", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
		else if (c == 'b') with_bidir = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools blacklist [options] <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    min region length [%d]\n", min_len);
		fprintf(stderr, "  -b        include regions involving both strands (mostly inversions)\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	gfa_sort_ref_arc(g);
	bb = gfa_bubble(g, &n_bb);
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		int32_t rst = b->ss, ren = b->se;
		if (!with_bidir && b->is_bidir) continue;
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
	int32_t i, j, c, n_bb, sub_gfa = 0;
	gfa_t *g;
	gfa_bubble_t *bb;

	while ((c = ketopt(&o, argc, argv, 1, "g", 0)) >= 0) {
		if (c == 'g') sub_gfa = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools bubble [-g] <in.gfa>\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	gfa_sort_ref_arc(g);
	bb = gfa_bubble(g, &n_bb);
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		if (sub_gfa) printf("bb\t");
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t-1\t-1\t-1\t", g->sseq[b->snid].name, b->ss, b->se, b->n_seg, b->n_paths,
			   b->is_bidir, b->len_min, b->len_max);
		for (j = 0; j < b->n_seg; ++j) {
			if (j) fputc(',', stdout);
			printf("%s", g->seg[b->v[j]>>1].name);
		}
		if (b->len_min == 0) fputs("\t*", stdout);
		else printf("\t%s", b->seq_min);
		if (b->len_max == 0) fputs("\t*", stdout);
		else printf("\t%s", b->seq_max);
		putchar('\n');
		if (sub_gfa) {
			int32_t i, *seg;
			gfa_t *f;
			GFA_CALLOC(seg, b->n_seg);
			for (i = 0; i < b->n_seg; ++i)
				seg[i] = b->v[i]>>1;
			f = gfa_subview(g, b->n_seg, seg);
			gfa_print(f, stdout, GFA_O_NO_SEQ);
			gfa_subview_destroy(f);
			printf("be\n");
		}
		free(b->v);
	}
	free(bb);
	gfa_destroy(g);
	return 0;
}

int main_sql(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, write_seq = 0;
	gfa_t *g;

	while ((c = ketopt(&o, argc, argv, 1, "s", 0)) >= 0) {
		if (c == 's') write_seq = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gfatools sql <in.gfa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s      write sequence\n");
		return 1;
	}
	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "ERROR: failed to read the graph\n");
		return 2;
	}
	gfa_sort_ref_arc(g);
	gfa_sql_write(stdout, g, write_seq);
	gfa_destroy(g);
	return 0;
}

int main_asm(int argc, char *argv[])
{
	const char *tr_opts = "v:ur:t:b:B:o:c:z:y";
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
		fprintf(stderr, "  -b INT1[,INT2]  pop bubbles (max radius, max deletions [inf])\n");
		fprintf(stderr, "  -B INT1[,INT2]  pop bubbles along with small tips (max radius, max del [inf])\n");
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
	gfa_arc_del_multi_risky(g);

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
			gfa_drop_tip(g, max_ext, max_len);
		} else if (c == 'b') {
			int32_t dist, max_del = INT32_MAX;
			dist = gfa_str2num(o.arg, &p);
			if (*p == ',') max_del = gfa_str2num(p + 1, &p);
			gfa_pop_bubble(g, dist, max_del, 1);
		} else if (c == 'B') {
			int32_t dist, max_del = INT32_MAX;
			dist = gfa_str2num(o.arg, &p);
			if (*p == ',') max_del = gfa_str2num(p + 1, &p);
			gfa_pop_bubble(g, dist, max_del, 0);
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
		} else if (c == 'z') {
			int32_t min_dist, max_dist = -1;
			min_dist = gfa_str2num(o.arg, &p);
			if (*p == ',') max_dist = gfa_str2num(p + 1, &p);
			else max_dist = min_dist * 2;
			gfa_cut_z(g, min_dist, max_dist);
		} else if (c == 'y') {
			gfa_scc_all(g);
		}
	}

	gfa_print(g, stdout, oflag);
	gfa_destroy(g);
	return 0;
}

int main_ed(int argc, char *argv[])
{
	extern int gfa_ed_dbg;
	gzFile fp;
	kseq_t *ks;
	ketopt_t o = KETOPT_INIT;
	gfa_t *g;
	gfa_edrst_t rst;
	gfa_edseq_t *es;
	gfa_edopt_t opt;
	int c, no_kalloc = 0;
	uint32_t v0 = 0<<0|0;
	void *km = 0;
	char *sname = 0;

	gfa_edopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "cl:s:d:w:m:n:K", 0)) >= 0) {
		if (c == 'l') opt.max_lag = atoi(o.arg);
		else if (c == 'w') opt.bw_dyn = atoi(o.arg);
		else if (c == 'n') opt.max_chk = atoi(o.arg);
		else if (c == 'm') opt.s_term = atoi(o.arg);
		else if (c == 'c') opt.traceback = 1;
		else if (c == 's') sname = o.arg;
		else if (c == 'd') gfa_ed_dbg = atoi(o.arg);
		else if (c == 'K') no_kalloc = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: gfa ed [options] <target.gfa|fa> <query.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    max lag behind the furthest wavefront; 0 to disable [0]\n");
		fprintf(stderr, "  -n INT    check max lag if there are more than INT diagnoals [%d]\n", opt.max_chk);
		fprintf(stderr, "  -w INT    band width [%d]\n", opt.bw_dyn);
		fprintf(stderr, "  -m INT    max edit distance; -1 to disable [-1]\n");
		fprintf(stderr, "  -s STR    starting segment name [first]\n");
		fprintf(stderr, "  -c        report the alignment path\n");
		return 1;
	}

	//km = km_init();

	g = gfa_read(argv[o.ind]);
	assert(g);
	if (sname) {
		int32_t sid;
		sid = gfa_name2id(g, sname);
		if (sid < 0) fprintf(stderr, "ERROR: failed to find segment '%s'\n", sname);
		else v0 = sid<<1 | 0; // TODO: also allow to change the orientation
	}
	es = gfa_edseq_init(g);

	km = no_kalloc? 0 : km_init();
	fp = gzopen(argv[o.ind+1], "r");
	assert(fp);
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int32_t s;
		s = gfa_edit_dist(km, &opt, g, es, ks->seq.l, ks->seq.s, v0, 0, &rst);
		if (opt.traceback) {
			int32_t i, last_len = -1, len = 0;
			printf("%s\t%d\t0\t%d\t+\t", ks->name.s, ks->seq.l, ks->seq.l);
			for (i = 0; i < rst.nv; ++i) {
				uint32_t v = rst.v[i];
				printf("%c%s", "><"[v&1], g->seg[v>>1].name);
				last_len = g->seg[v>>1].len;
				len += last_len;
			}
			printf("\t%d\t0\t%d\t%d\n", len, len - (last_len - rst.end_off) + 1, rst.s);
			kfree(km, rst.v);
		} else printf("%s\t%d\n", ks->name.s, s);
	}
	kseq_destroy(ks);
	gzclose(fp);

	gfa_edseq_destroy(g->n_seg, es);
	gfa_destroy(g);
	km_destroy(km);
	return 0;
}

int main(int argc, char *argv[])
{
	extern double gfa_realtime(void);
	extern double gfa_cputime(void);
	double t_start;
	int ret = 0, i;

	if (argc == 1) {
		fprintf(stderr, "Usage: gfatools <command> <arguments>\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  view        read a GFA file\n");
		fprintf(stderr, "  stat        statistics about a GFA file\n");
		fprintf(stderr, "  gfa2fa      convert GFA to FASTA\n");
		fprintf(stderr, "  gfa2bed     convert rGFA to BED (requiring rGFA)\n");
		fprintf(stderr, "  blacklist   blacklist regions\n");
		fprintf(stderr, "  bubble      print bubble-like regions (EXPERIMENTAL)\n");
		fprintf(stderr, "  asm         miniasm-like graph transformation\n");
		fprintf(stderr, "  sql         export rGFA to SQLite (requiring rGFA)\n");
		fprintf(stderr, "  ed          GWFA prefix alignment (for evaluation only)\n");
		fprintf(stderr, "  version     print version number\n");
		return 1;
	}

	t_start = gfa_realtime();
	if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "stat") == 0) ret = main_stat(argc-1, argv+1);
	else if (strcmp(argv[1], "gfa2bed") == 0) ret = main_gfa2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "gfa2fa") == 0) ret = main_gfa2fa(argc-1, argv+1);
	else if (strcmp(argv[1], "blacklist") == 0) ret = main_blacklist(argc-1, argv+1);
	else if (strcmp(argv[1], "bubble") == 0) ret = main_bubble(argc-1, argv+1);
	else if (strcmp(argv[1], "asm") == 0) ret = main_asm(argc-1, argv+1);
	else if (strcmp(argv[1], "sql") == 0) ret = main_sql(argc-1, argv+1);
	else if (strcmp(argv[1], "ed") == 0) ret = main_ed(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("gfa.h: %s\ngfatools: %s\n", GFA_VERSION, GFATOOLS_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, GFATOOLS_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, gfa_realtime() - t_start, gfa_cputime());
	}
	return ret;
}
