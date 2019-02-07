#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "kstring.h"
#include "gfa.h"

#include "kseq.h"
KSTREAM_INIT2(, gzFile, gzread, 65536)

#include "khash.h"
KHASH_INIT2(seg,, kh_cstr_t, uint32_t, 1, kh_str_hash_func, kh_str_hash_equal)
typedef khash_t(seg) seghash_t;

#include "ksort.h"
#define gfa_arc_key(a) ((a).v_lv)
KRADIX_SORT_INIT(arc, gfa_arc_t, gfa_arc_key, 8)

int gfa_verbose = 2;

/********************
 * Tag manipulation *
 ********************/

int gfa_aux_parse(char *s, uint8_t **data, int *max)
{
	char *q, *p;
	kstring_t str;
	if (s == 0) return 0;
	str.l = 0, str.m = *max, str.s = (char*)*data;
	if (*s == '\t') ++s;
	for (p = q = s;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
				int type = q[3];
				kputsn_(q, 2, &str);
				q += 5;
				if (type == 'A') {
					kputc_('A', &str);
					kputc_(*q, &str);
				} else if (type == 'i') {
					int32_t x;
					x = strtol(q, &q, 10);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'f') {
					float x;
					x = strtod(q, &q);
					kputc_('f', &str); kputsn_(&x, 4, &str);
				} else if (type == 'Z') {
					kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
				} else if (type == 'B') {
					type = *q++; // q points to the first ',' following the typing byte
					if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
						int32_t n;
						char *r;
						for (r = q, n = 0; *r; ++r)
							if (*r == ',') ++n;
						kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
						// TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
						if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'I') while (q + 1 < p) { uint32_t x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
					}
				} // should not be here, as we have tested all types
			}
			q = p + 1;
			if (c == 0) break;
		}
	}
	if (str.s) str.s[str.l] = 0;
	*max = str.m, *data = (uint8_t*)str.s;
	return str.l;
}

int gfa_aux_format(int l_aux, const uint8_t *aux, char **t, int *max)
{
	kstring_t str;
	const uint8_t *s = aux;
	str.l = 0, str.s = *t, str.m = *max;
	while (s < aux + l_aux) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s++;
		kputc('\t', &str); kputsn((char*)key, 2, &str); kputc(':', &str);
		if (type == 'A') { kputsn("A:", 2, &str); kputc(*s, &str); ++s; }
		else if (type == 'i') { kputsn("i:", 2, &str); kputw(*(int32_t*)s, &str); s += 4; }
		else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
		else if (type == 'Z') { kputc(type, &str); kputc(':', &str); while (*s) kputc(*s++, &str); ++s; }
		else if (type == 'B') {
			uint8_t sub_type = *(s++);
			int32_t i, n;
			memcpy(&n, s, 4);
			s += 4; // no point to the start of the array
			kputsn("B:", 2, &str); kputc(sub_type, &str); // write the typing
			for (i = 0; i < n; ++i) { // FIXME: for better performance, put the loop after "if"
				kputc(',', &str);
				if ('c' == sub_type)      { kputw(*(int8_t*)s, &str); ++s; }
				else if ('C' == sub_type) { kputw(*(uint8_t*)s, &str); ++s; }
				else if ('s' == sub_type) { kputw(*(int16_t*)s, &str); s += 2; }
				else if ('S' == sub_type) { kputw(*(uint16_t*)s, &str); s += 2; }
				else if ('i' == sub_type) { kputw(*(int32_t*)s, &str); s += 4; }
				else if ('I' == sub_type) { kputuw(*(uint32_t*)s, &str); s += 4; }
				else if ('f' == sub_type) { ksprintf(&str, "%g", *(float*)s); s += 4; }
			}
		}
	}
	*t = str.s, *max = str.m;
	return str.l;
}

static inline int gfa_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + gfa_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += gfa_aux_type2size(type); \
	} while(0)

uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2])
{
	const uint8_t *s = data;
	int y = tag[0]<<8 | tag[1];
	while (s < data + l_data) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return (uint8_t*)s;
		__skip_tag(s);
	}
	return 0;
}

// s MUST BE returned by gfa_aux_get()
int gfa_aux_del(int l_data, uint8_t *data, uint8_t *s)
{
	uint8_t *p;
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, l_data - (s - data));
	return l_data - (s - p);
}

/******************
 * Basic routines *
 ******************/

gfa_t *gfa_init(void)
{
	gfa_t *g;
	g = (gfa_t*)calloc(1, sizeof(gfa_t));
	g->h_names = kh_init(seg);
	return g;
}

void gfa_destroy(gfa_t *g)
{
	uint32_t i, j;
	uint64_t k;
	if (g == 0) return;
	kh_destroy(seg, (seghash_t*)g->h_names);
	for (i = 0; i < g->n_seg; ++i) {
		gfa_seg_t *s = &g->seg[i];
		free(s->name);
		free(s->aux.aux);
		for (j = 0; j < s->utg.n; ++j)
			free(s->utg.name[j]);
		free(s->utg.name);
		free(s->utg.a);
	}
	for (k = 0; k < g->n_arc; ++k)
		free(g->arc_aux[k].aux);
	free(g->idx); free(g->seg); free(g->arc); free(g->arc_aux);
	free(g);
}

int32_t gfa_add_seg(gfa_t *g, const char *name)
{
	khint_t k;
	int absent;
	seghash_t *h = (seghash_t*)g->h_names;
	k = kh_put(seg, h, name, &absent);
	if (absent) {
		gfa_seg_t *s;
		if (g->n_seg == g->m_seg) {
			uint32_t old_m = g->m_seg;
			g->m_seg = g->m_seg? g->m_seg<<1 : 16;
			g->seg = (gfa_seg_t*)realloc(g->seg, g->m_seg * sizeof(gfa_seg_t));
			memset(&g->seg[old_m], 0, (g->m_seg - old_m) * sizeof(gfa_seg_t));
		}
		s = &g->seg[g->n_seg++];
		kh_key(h, k) = s->name = strdup(name);
		s->del = s->len = 0;
		kh_val(h, k) = g->n_seg - 1;
	}
	return kh_val(h, k);
}

uint64_t gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp)
{
	gfa_arc_t *a;
	if (g->m_arc == g->n_arc) {
		uint64_t old_m = g->m_arc;
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (gfa_arc_t*)realloc(g->arc, g->m_arc * sizeof(gfa_arc_t));
		memset(&g->arc[old_m], 0, (g->m_arc - old_m) * sizeof(gfa_arc_t));
		g->arc_aux = (gfa_aux_t*)realloc(g->arc_aux, g->m_arc * sizeof(gfa_aux_t));
		memset(&g->arc_aux[old_m], 0, (g->m_arc - old_m) * sizeof(gfa_aux_t));
	}
	a = &g->arc[g->n_arc++];
	a->v_lv = (uint64_t)v << 32;
	a->w = w, a->ov = ov, a->ow = ow, a->lw = 0;
	a->link_id = link_id >= 0? link_id : g->n_arc - 1;
	a->del = 0;
	a->comp = comp;
	return a->link_id;
}

void gfa_arc_sort(gfa_t *g)
{
	radix_sort_arc(g->arc, g->arc + g->n_arc);
}

uint64_t *gfa_arc_index_core(size_t max_seq, size_t n, const gfa_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || gfa_arc_head(a[i-1]) != gfa_arc_head(a[i]))
			idx[gfa_arc_head(a[i-1])] = (uint64_t)last<<32 | (i - last), last = i;
	return idx;
}

void gfa_arc_index(gfa_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = gfa_arc_index_core(g->n_seg, g->n_arc, g->arc);
}

/****************
 * Line parsers *
 ****************/

int gfa_parse_S(gfa_t *g, char *s)
{
	int i, is_ok = 0;
	char *p, *q, *seg = 0, *seq = 0, *rest = 0;
	uint32_t sid, len = 0;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) seg = q;
			else if (i == 1) {
				seq = q[0] == '*'? 0 : strdup(q);
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			++i, q = p + 1;
			if (c == 0) break;
		}
	}
	if (is_ok) { // all mandatory fields read
		int l_aux, m_aux = 0;
		uint8_t *aux = 0;
		gfa_seg_t *s;
		l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
		if (seq == 0) {
			uint8_t *s;
			s = gfa_aux_get(l_aux, aux, "LN");
			if (s && s[0] == 'i')
				len = *(int32_t*)(s+1);
		} else len = strlen(seq);
		sid = gfa_add_seg(g, seg);
		s = &g->seg[sid];
		s->len = len, s->seq = seq;
		s->aux.m_aux = m_aux, s->aux.l_aux = l_aux, s->aux.aux = aux;
	} else return -1;
	return 0;
}

int gfa_parse_L(gfa_t *g, char *s)
{
	int i, oriv = -1, oriw = -1, is_ok = 0;
	char *p, *q, *segv = 0, *segw = 0, *rest = 0;
	int32_t ov = INT32_MAX, ow = INT32_MAX;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) {
				segv = q;
			} else if (i == 1) {
				if (*q != '+' && *q != '-') return -2;
				oriv = (*q != '+');
			} else if (i == 2) {
				segw = q;
			} else if (i == 3) {
				if (*q != '+' && *q != '-') return -2;
				oriw = (*q != '+');
			} else if (i == 4) {
				if (*q == ':') {
					ov = INT32_MAX;
					ow = isdigit(*(q+1))? strtol(q+1, &q, 10) : INT32_MAX;
				} else if (isdigit(*q)) {
					char *r;
					ov = strtol(q, &r, 10);
					if (isupper(*r)) { // CIGAR
						ov = ow = 0;
						do {
							long l;
							l = strtol(q, &q, 10);
							if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
							if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
							++q;
						} while (isdigit(*q));
					} else if (*r == ':') { // overlap lengths
						ow = isdigit(*(r+1))? strtol(r+1, &r, 10) : INT32_MAX;
					} else break;
				} else break;
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			++i, q = p + 1;
			if (c == 0) break;
		}
	}
	if (is_ok) {
		uint32_t v, w;
		uint64_t link_id;
		int l_aux, m_aux = 0;
		uint8_t *aux = 0;
		v = gfa_add_seg(g, segv) << 1 | oriv;
		w = gfa_add_seg(g, segw) << 1 | oriw;
		link_id = gfa_add_arc1(g, v, w, ov, ow, -1, 0);
		l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
		if (l_aux) {
			gfa_aux_t *a = &g->arc_aux[link_id];
			uint8_t *s_L1, *s_L2;
			a->l_aux = l_aux, a->m_aux = m_aux, a->aux = aux;
			s_L1 = gfa_aux_get(a->l_aux, a->aux, "L1");
			if (s_L1) {
				if (ov != INT32_MAX && s_L1[0] == 'i')
					g->seg[v>>1].len = g->seg[v>>1].len > ov + *(int32_t*)(s_L1+1)? g->seg[v>>1].len : ov + *(int32_t*)(s_L1+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L1);
			}
			s_L2 = gfa_aux_get(a->l_aux, a->aux, "L2");
			if (s_L2) {
				if (ow != INT32_MAX && s_L2[0] == 'i')
					g->seg[w>>1].len = g->seg[w>>1].len > ow + *(int32_t*)(s_L2+1)? g->seg[w>>1].len : ow + *(int32_t*)(s_L2+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L2);
			}
			if (a->l_aux == 0) {
				free(a->aux);
				a->aux = 0;
			}
		}
	} else return -1;
	return 0;
}

/********************
 * Fix graph issues *
 ********************/

uint32_t gfa_fix_no_seg(gfa_t *g)
{
	uint32_t i, n_err = 0;
	for (i = 0; i < g->n_seg; ++i) {
		gfa_seg_t *s = &g->seg[i];
		if (s->len == 0) {
			++n_err, s->del = 1;
			if (gfa_verbose >= 2)
				fprintf(stderr, "[W] segment '%s' is used on an L-line but not defined on an S-line\n", s->name);
		}
	}
	return n_err;
}

void gfa_fix_arc_len(gfa_t *g)
{
	uint64_t k;
	for (k = 0; k < g->n_arc; ++k) {
		gfa_arc_t *a = &g->arc[k];
		uint32_t v = gfa_arc_head(*a), w = gfa_arc_tail(*a);
		if (g->seg[v>>1].del || g->seg[w>>1].del) {
			a->del = 1;
		} else {
			a->v_lv |= g->seg[v>>1].len - a->ov;
			if (a->ow != INT32_MAX)
				a->lw = g->seg[w>>1].len - a->ow;
		}
	}
}

uint32_t gfa_fix_semi_arc(gfa_t *g)
{
	uint32_t n_err = 0, v, n_vtx = gfa_n_vtx(g);
	int i, j;
	for (v = 0; v < n_vtx; ++v) {
		int nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			if (!av[i].del && (av[i].ow == INT32_MAX || av[i].ov == INT32_MAX)) { // overlap length is missing
				uint32_t w = av[i].w^1;
				int is_multi = 0, c, jv = -1, nw = gfa_arc_n(g, w);
				gfa_arc_t *aw = gfa_arc_a(g, w);
				for (j = 0, c = 0; j < nw; ++j)
					if (!aw[j].del && aw[j].w == (v^1)) ++c, jv = j;
				if (c == 1) {
					if (av[i].ov != INT32_MAX && aw[jv].ow != INT32_MAX && av[i].ov != aw[jv].ow) is_multi = 1;
					if (av[i].ow != INT32_MAX && aw[jv].ov != INT32_MAX && av[i].ow != aw[jv].ov) is_multi = 1;
				}
				if (c == 1 && !is_multi) {
					if (aw[jv].ov != INT32_MAX) av[i].ow = aw[jv].ov;
					if (aw[jv].ow != INT32_MAX) av[i].ov = aw[jv].ow;
				} else {
					if (gfa_verbose >= 2)
						fprintf(stderr, "[W] can't infer overlap length for %s%c -> %s%c\n",
								g->seg[v>>1].name, "+-"[v&1], g->seg[w>>1].name, "+-"[(w^1)&1]);
					++n_err;
					av[i].del = 1;
				}
			}
		}
	}
	return n_err;
}

uint32_t gfa_fix_symm(gfa_t *g)
{
	uint32_t n_err = 0, v, n_vtx = gfa_n_vtx(g);
	int i;
	for (v = 0; v < n_vtx; ++v) {
		int nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			int j, nw;
			gfa_arc_t *aw, *avi = &av[i];
			if (avi->del || avi->comp) continue;
			nw = gfa_arc_n(g, avi->w^1);
			aw = gfa_arc_a(g, avi->w^1);
			for (j = 0; j < nw; ++j) {
				gfa_arc_t *awj = &aw[j];
				if (awj->del || awj->comp) continue;
				if (awj->w == (v^1) && awj->ov == avi->ow && awj->ow == avi->ov) { // complement found
					awj->comp = 1;
					awj->link_id = avi->link_id;
					break;
				}
			}
			if (j == nw) {
				gfa_arc_t *arc_old = g->arc;
				gfa_add_arc1(g, avi->w^1, v^1, avi->ow, avi->ov, avi->link_id, 1);
				if (arc_old != g->arc) av = gfa_arc_a(g, v); // g->arc may be reallocated
			}
		}
	}
	if (n_vtx < gfa_n_vtx(g)) {
		gfa_arc_sort(g);
		gfa_arc_index(g);
	}
	return n_err;
}

void gfa_arc_rm(gfa_t *g)
{
	uint32_t e, n;
	for (e = n = 0; e < g->n_arc; ++e) {
		uint32_t u = g->arc[e].v_lv>>32, v = g->arc[e].w;
		if (!g->arc[e].del && !g->seg[u>>1].del && !g->seg[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

void gfa_cleanup(gfa_t *g)
{
	gfa_arc_rm(g);
	if (!g->is_srt) {
		gfa_arc_sort(g);
		g->is_srt = 1;
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	if (g->idx == 0) gfa_arc_index(g);
}

/****************
 * User-end I/O *
 ****************/

gfa_t *gfa_read(const char *fn)
{
	gzFile fp;
	gfa_t *g;
	kstring_t s = {0,0,0};
	kstream_t *ks;
	int dret;
	uint64_t lineno = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	g = gfa_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
		int ret = 0;
		++lineno;
		if (s.l < 3 || s.s[1] != '\t') continue; // empty line
		if (s.s[0] == 'S') ret = gfa_parse_S(g, s.s);
		else if (s.s[0] == 'L') ret = gfa_parse_L(g, s.s);
		if (ret < 0 && gfa_verbose >= 1)
			fprintf(stderr, "[E] invalid %c-line at line %ld (error code %d)\n", s.s[0], (long)lineno, ret);
	}
	free(s.s);
	gfa_fix_no_seg(g);
	gfa_arc_sort(g);
	gfa_arc_index(g);
	gfa_fix_semi_arc(g);
	gfa_fix_symm(g);
	gfa_fix_arc_len(g);
	gfa_cleanup(g);
	ks_destroy(ks);
	gzclose(fp);
	return g;
}

void gfa_print(const gfa_t *g, FILE *fp, int M_only)
{
	uint32_t i;
	uint64_t k;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		if (s->del) continue;
		fprintf(fp, "S\t%s\t", s->name);
		if (s->seq) fputs(s->seq, fp);
		else fputc('*', fp);
		if (s->aux.l_aux == 0 || !gfa_aux_get(s->aux.l_aux, s->aux.aux, "LN"))
			fprintf(fp, "\tLN:i:%d", s->len);
		if (s->aux.l_aux > 0) {
			char *t = 0;
			int max = 0, len;
			len = gfa_aux_format(s->aux.l_aux, s->aux.aux, &t, &max);
			fputs(t, fp);
			free(t);
		}
		fputc('\n', fp);
		if (s->utg.n) {
			uint32_t j, l;
			for (j = l = 0; j < s->utg.n; ++j) {
				const gfa_utg_t *u = &s->utg;
				fprintf(fp, "a\t%s\t%d\t%s\t%c\t%d\n", s->name, l, u->name[j], "+-"[u->a[j]>>32&1], (uint32_t)u->a[j]);
				l += (uint32_t)u->a[j];
			}
		}
	}
	for (k = 0; k < g->n_arc; ++k) {
		const gfa_arc_t *a = &g->arc[k];
		const gfa_aux_t *aux = &g->arc_aux[a->link_id];
		if (a->del || a->comp) continue;
		fprintf(fp, "L\t%s\t%c\t%s\t%c", g->seg[a->v_lv>>33].name, "+-"[a->v_lv>>32&1], g->seg[a->w>>1].name, "+-"[a->w&1]);
		if (M_only) {
			fprintf(fp, "\t%dM", a->ov < a->ow? a->ov : a->ow);
		} else {
			if (a->ov == a->ow) fprintf(fp, "\t%dM", a->ov);
			else fprintf(fp, "\t%d:%d", a->ov, a->ow);
		}
		fprintf(fp, "\tL1:i:%d", gfa_arc_len(*a));
		fprintf(fp, "\tL2:i:%d", a->lw);
		if (aux->l_aux) {
			char *t = 0;
			int max = 0, len;
			len = gfa_aux_format(aux->l_aux, aux->aux, &t, &max);
			if (t) fputs(t, fp);
			free(t);
		}
		fputc('\n', fp);
	}
}

int32_t gfa_name2id(const gfa_t *g, const char *name)
{
	seghash_t *h = (seghash_t*)g->h_names;
	khint_t k;
	k = kh_get(seg, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

