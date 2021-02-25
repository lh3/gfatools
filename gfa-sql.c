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
CREATE TABLE bwalk (         -- walk in each bubble\n\
  wid     INTEGER,           -- walk ID\n\
  bid     INTEGER NOT NULL,  -- bubble ID\n\
  len     INTEGER NOT NULL,  -- walk length\n\
  walk    TEXT NOT NULL,     -- walk\n\
  PRIMARY KEY (wid)\n\
);\n\
CREATE TABLE call (\n\
  bid     INTEGER,           -- bubble ID\n\
  sample  TEXT,              -- sample name\n\
  wid     INTEGER NOT NULL,  -- walk ID\n\
  ctg     TEXT NOT NULL,     -- sample contig name\n\
  start   INTEGER NOT NULL,  -- contig start (BED-like)\n\
  end     INTEGER NOT NULL,  -- contig end\n\
  strand  INTEGER NOT NULL,  -- contig strand\n\
  PRIMARY KEY (bid, sample)\n\
);\n\
";

const char *gfa_sql_index = "\
CREATE UNIQUE INDEX idx_name  ON seg   (name);\n\
CREATE        INDEX idx_scoor ON seg   (sname, soff);\n\
CREATE UNIQUE INDEX idx_s2v   ON vtx   (sid, strand);\n\
CREATE UNIQUE INDEX idx_bst   ON bbl   (sname, start);\n\
CREATE UNIQUE INDEX idx_ben   ON bbl   (sname, end);\n\
CREATE        INDEX idx_b2v   ON b2v   (bid, vid);\n\
CREATE        INDEX idx_v2b   ON b2v   (vid);\n\
CREATE        INDEX idx_bwalk ON bwalk (bid);\n\
CREATE        INDEX idx_cst   ON call  (ctg, start);\n\
CREATE        INDEX idx_cen   ON call  (ctg, end);\n\
";

const char *gfa_sql_view = "\
CREATE VIEW named_vtx AS SELECT v.vid, v.sid, s.name, v.strand, s.len, s.rank FROM vtx v, seg s WHERE v.sid = s.sid;\n\
CREATE VIEW named_arc AS\n\
  SELECT a.vid1, v1.sid AS sid1, s1.name AS name1, v1.strand AS strand1, a.vid2, v2.sid AS sid2, s2.name AS name2, v2.strand AS strand2, a.ori, a.rank\n\
    FROM arc a, vtx v1, vtx v2, seg s1, seg s2\n\
    WHERE a.vid1 = v1.vid AND a.vid2 = v2.vid AND v1.sid = s1.sid AND v2.sid = s2.sid;\n\
CREATE VIEW var_call AS\n\
  SELECT b.bid, b.sname, b.start, b.end, c.sample, c.wid, w.len, c.strand, c.ctg, c.start AS ctg_start, c.end AS ctg_end\n\
    FROM bbl b, call c, bwalk w WHERE c.bid=b.bid AND c.wid = w.wid AND c.bid = w.bid ORDER BY c.bid;\n\
";

// WITH x AS (SELECT t.vid FROM bbl b, b2v t WHERE b.sname="chr11" AND b.end>=10000000 AND b.start<=11000000 AND b.bid=t.bid) SELECT v.* FROM named_vtx v WHERE v.vid IN x;
// WITH x AS (SELECT t.vid>>1 FROM bbl b, b2v t WHERE b.sname="chr11" AND b.end>=10000000 AND b.start<=11000000 AND b.bid=t.bid) SELECT a.* FROM named_arc a WHERE a.vid1>>1 IN x AND a.vid2>>1 IN x AND a.ori=1;
// SELECT c.bid, c.sname, c.start, c.end, c.len, COUNT(c.wid) FROM var_call c WHERE c.sname="chr11" AND c.end>=10000000 AND c.start<=11000000 group by c.bid, c.wid;

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

void gfa_sql_write_seq(FILE *fp, const gfa_t *g, kstring_t *out)
{
	uint32_t i;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		out->l = 0;
		mg_sprintf_lite(out, "INSERT INTO seq (sid,seq) VALUES ('%d','%s');\n", i, s->seq);
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

void gfa_sql_write(FILE *fp, const gfa_t *g, int ins_seq)
{
	kstring_t out = {0,0,0};
	fputs(gfa_sql_schema, fp);
	fputs("BEGIN TRANSACTION;\n", fp);
	gfa_sql_write_seg(fp, g, &out);
	gfa_sql_write_vtx(fp, g, &out);
	gfa_sql_write_arc(fp, g, &out);
	gfa_sql_write_bbl(fp, g, &out);
	if (ins_seq) gfa_sql_write_seq(fp, g, &out);
	fputs("END TRANSACTION;\n", fp);
	fputs(gfa_sql_index, fp);
	fputs(gfa_sql_view, fp);
	free(out.s);
}
