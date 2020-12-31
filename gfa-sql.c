#include <stdarg.h>
#include <stdio.h>
#include <assert.h>
#include "kstring.h"
#include "gfa-priv.h"

const char *gfa_sql_schema = "\
CREATE TABLE seg (           -- segment\n\
  sid     INTEGER,           -- segmend ID\n\
  name    TEXT NOT NULL,     -- segment name\n\
  len     INTEGER NOT NULL,  -- segment length\n\
  sname   TEXT,              -- stable name (tag SN)\n\
  soff    INTEGER,           -- stable offset (tag SO)\n\
  rank    INTEGER,           -- rank (tag SR)\n\
  PRIMARY KEY (sid)\n\
);\n\
CREATE TABLE vtx (           -- vertex (aka oriented segment)\n\
  vid     INTEGER,           -- vertex ID\n\
  sid     INTEGER NOT NULL,  -- segment ID\n\
  strand  INTEGER NOT NULL,  -- orientation\n\
  PRIMARY KEY (vid)\n\
);\n\
CREATE TABLE arc (           -- arc (edge)\n\
  vid1    INTEGER,           -- first vertex ID\n\
  vid2    INTEGER,           -- second vertex ID\n\
  ori     INTEGER NOT NULL,  -- in the input (0 or 1)\n\
  rank    INTEGER,           -- rank (tag SR)\n\
  PRIMARY KEY (vid1, vid2)\n\
);\n\
CREATE TABLE seq (           -- segment sequence\n\
  sid     INTEGRE,           -- segment ID\n\
  seq     TEXT NOT NULL,     -- sequence\n\
  PRIMARY KEY (sid)\n\
);\n\
CREATE TABLE bbl (           -- bubble\n\
  bid     INTEGER,           -- bubble ID\n\
  sname   TEXT NOT NULL,     -- stable name\n\
  start   INTEGER NOT NULL,  -- start\n\
  end     INTEGER NOT NULL,  -- end\n\
  src     INTEGER NOT NULL,  -- source vertex\n\
  sink    INTEGER NOT NULL,  -- sink vertex\n\
  n_vtx   INTEGER NOT NULL,  -- number of vertices in the bubble\n\
  inv     INTEGER NOT NULL,  -- inversion or not (0 or 1)\n\
  n_path  INTEGER NOT NULL,  -- number of possible paths\n\
  min_len INTEGER NOT NULL,  -- min path length\n\
  max_len INTEGER NOT NULL,  -- max path length\n\
  PRIMARY KEY (bid)\n\
);\n\
CREATE TABLE b2v (           -- bubble-vertex relationship\n\
  bid     INTEGER NOT NULL,  -- bubble ID\n\
  vid     INTEGER NOT NULL   -- vertex ID\n\
);\n\
";

const char *gfa_sql_index = "\
CREATE INDEX idx_name  ON seg (name);\n\
CREATE INDEX idx_scoor ON seg (sname, soff);\n\
CREATE INDEX idx_s2v   ON vtx (sid, strand);\n\
CREATE INDEX idx_bst   ON bbl (sname, start);\n\
CREATE INDEX idx_ben   ON bbl (sname, end);\n\
CREATE INDEX idx_b2v   ON b2v (bid, vid);\n\
";

const char *gfa_sql_view = "\
CREATE VIEW named_vtx AS SELECT v.vid, v.sid, s.name, v.strand, s.len, s.rank FROM vtx v, seg s WHERE v.sid = s.sid;\n\
CREATE VIEW named_arc AS\n\
  SELECT a.vid1, v1.sid AS sid1, s1.name AS name1, v1.strand AS strand1, a.vid2, v2.sid AS sid2, s2.name AS name2, v2.strand AS strand2, a.ori, a.rank\n\
    FROM arc a, vtx v1, vtx v2, seg s1, seg s2\n\
    WHERE a.vid1 = v1.vid AND a.vid2 = v2.vid AND v1.sid = s1.sid AND v2.sid = s2.sid;\n\
";

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void mg_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else abort();
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

void gfa_sql_write_seg(FILE *fp, const gfa_t *g, kstring_t *out)
{
	uint32_t i;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		out->l = 0;
		if (s->rank >= 0 && s->snid >= 0 && s->soff >= 0) {
			mg_sprintf_lite(out, "INSERT INTO seg (sid,name,len,sname,soff,rank) VALUES ('%u','%s','%d','%s','%d','%d');\n",
				i, s->name, s->len, g->sseq[s->snid].name, s->soff, s->rank);
		} else {
			mg_sprintf_lite(out, "INSERT INTO seg (sid,name,len) VALUES ('%u','%s','%d')\n", i, s->name, s->len);
		}
		fputs(out->s, fp);
	}
}

void gfa_sql_write_vtx(FILE *fp, const gfa_t *g, kstring_t *out)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	for (v = 0; v < n_vtx; ++v) {
		out->l = 0;
		mg_sprintf_lite(out, "INSERT INTO vtx (vid,sid,strand) VALUES ('%u','%d','%d');\n", v, v>>1, v&1? -1 : 1);
		fputs(out->s, fp);
	}
}

void gfa_sql_write_arc(FILE *fp, const gfa_t *g, kstring_t *out)
{
	uint64_t i;
	for (i = 0; i < g->n_arc; ++i) {
		const gfa_arc_t *a = &g->arc[i];
		out->l = 0;
		if (a->rank >= 0)
			mg_sprintf_lite(out, "INSERT INTO arc (vid1,vid2,ori,rank) VALUES ('%u','%d','%d','%d');\n", (uint32_t)(a->v_lv>>32), a->w, !a->comp, a->rank);
		else
			mg_sprintf_lite(out, "INSERT INTO arc (vid1,vid2,ori) VALUES ('%u','%d');\n", (uint32_t)(a->v_lv>>32), a->w, !a->comp);
		fputs(out->s, fp);
	}
}

void gfa_sql_write_bbl(FILE *fp, const gfa_t *g, kstring_t *out)
{
	int32_t i, n_bb;
	gfa_bubble_t *bb;
	bb = gfa_bubble(g, &n_bb);
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		out->l = 0;
		assert(b->n_seg >= 3);
		mg_sprintf_lite(out, "INSERT INTO bbl (bid,sname,start,end,src,sink,n_vtx,inv,n_path,min_len,max_len) VALUES ('%d','%s','%d','%d','%d','%d','%d','%d','%d','%d','%d');\n",
			i, g->sseq[b->snid].name, b->ss, b->se, b->vs, b->ve, b->n_seg - 2, !!b->is_bidir, b->n_paths, b->len_min, b->len_max);
		fputs(out->s, fp);
	}
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		int32_t j;
		for (j = 0; j < b->n_seg; ++j) {
			out->l = 0;
			mg_sprintf_lite(out, "INSERT INTO b2v (bid,vid) VALUES (%d,%d);\n", i, b->v[j]);
			fputs(out->s, fp);
		}
		free(b->v);
	}
	free(bb);
}

void gfa_sql_write(FILE *fp, const gfa_t *g)
{
	kstring_t out = {0,0,0};
	fputs(gfa_sql_schema, fp);
	fputs("BEGIN TRANSACTION;\n", fp);
	gfa_sql_write_seg(fp, g, &out);
	gfa_sql_write_vtx(fp, g, &out);
	gfa_sql_write_arc(fp, g, &out);
	gfa_sql_write_bbl(fp, g, &out);
	fputs("END TRANSACTION;\n", fp);
	fputs(gfa_sql_index, fp);
	fputs(gfa_sql_view, fp);
	free(out.s);
}
