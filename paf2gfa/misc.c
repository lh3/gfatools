#include "miniasm.h"
#include "gfa-priv.h"

int ma_verbose = 3;

const char *ma_timestamp(void)
{
	extern double gfa_realtime0;
	extern double gfa_realtime(void);
	extern double gfa_cputime(void);
	static char buf[256];
	double rt, ct;
	rt = gfa_realtime() - gfa_realtime0;
	ct = gfa_cputime();
	snprintf(buf, 255, "%.3f*%.2f", rt, ct/rt);
	return buf;
}

void ma_sys_init(void)
{
	extern void gfa_sys_init(void);
	gfa_sys_init();
}

#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, shash_t, h_s2i, kh_cstr_t, uint32_t, kh_hash_str, kh_eq_str)

sdict_t *sd_init(void)
{
	sdict_t *d;
	d = (sdict_t*)calloc(1, sizeof(sdict_t));
	d->h = h_s2i_init();
	return d;
}

void sd_destroy(sdict_t *d)
{
	uint32_t i;
	if (d == 0) return;
	if (d->h) h_s2i_destroy((shash_t*)d->h);
	for (i = 0; i < d->n_seq; ++i)
		free(d->seq[i].name);
	free(d->seq);
	free(d);
}

int32_t sd_put(sdict_t *d, const char *name, uint32_t len)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	int absent;
	k = h_s2i_put(h, name, &absent);
	if (absent) {
		sd_seq_t *s;
		if (d->n_seq == d->m_seq) {
			d->m_seq = d->m_seq? d->m_seq<<1 : 16;
			d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));
		}
		s = &d->seq[d->n_seq];
		s->len = len, s->aux = 0, s->del = 0;
		kh_key(h, k) = s->name = gfa_strdup(name);
		kh_val(h, k) = d->n_seq++;
	} // TODO: test if len is the same;
	return kh_val(h, k);
}

int32_t sd_get(const sdict_t *d, const char *name)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	k = h_s2i_get(h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

void sd_hash(sdict_t *d)
{
	uint32_t i;
	shash_t *h;
	if (d->h) return;
	d->h = h = h_s2i_init();
	for (i = 0; i < d->n_seq; ++i) {
		int absent;
		khint_t k;
		k = h_s2i_put(h, d->seq[i].name, &absent);
		kh_val(h, k) = i;
	}
}

int32_t *sd_squeeze(sdict_t *d)
{
	int32_t *map, i, j;
	if (d->h) {
		h_s2i_destroy((shash_t*)d->h);
		d->h = 0;
	}
	map = (int32_t*)calloc(d->n_seq, 4);
	for (i = j = 0; i < d->n_seq; ++i) {
		if (d->seq[i].del) {
			free(d->seq[i].name);
			map[i] = -1;
		} else d->seq[j] = d->seq[i], map[i] = j++;
	}
	d->n_seq = j;
	sd_hash(d);
	return map;
}
