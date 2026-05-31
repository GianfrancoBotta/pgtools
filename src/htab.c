#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

# include "htab.h"
#include "khashl.h" // hash table
#include "utils.h"

KHASHL_MAP_INIT(, pg_ht_t, pg_ht, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// Operations on hash tables and bloom filter.

pg_mht_t *pg_mht_init(int k, int pre)
{
	pg_mht_t *h;
	int i;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre); // allocate the array of partitions.
	for (i = 0; i < 1<<h->pre; ++i) {
		h->h[i].h = pg_ht_init(); // initialize hash table for each bucket.

		pthread_mutex_init(&h->h[i].lock, NULL); // initialize lock for each partition
	}
	return h;
}

void pg_mht_destroy(pg_mht_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i) {
		pthread_mutex_destroy(&h->h[i].lock);
		pg_ht_destroy(h->h[i].h); // destroy hash table for each bucket.
	}
	free(h->h); free(h);
}

long long pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq, int ff)
{	
	long long n_del = 0;
    int i, n = 1 << h->pre;

    for (i = 0; i < n; ++i) {
		// store entries to delete
        pg_ht1_t *g = &h->h[i];
		uint64_t *del_part = malloc(kh_size(g->h) * sizeof(uint64_t));
		int n_del_part = 0;
        khint_t k;
        // pthread_mutex_lock(&g->lock);
        for (k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;
			uint32_t v = kh_val(g->h, k);
			if ((double)(n_proc - val_count(v)) / n_tot > (1.0 - min_freq)) {
				del_part[n_del_part++] = kh_key(g->h, k);
			}
			else if (ff && (!val_snp1(v) || val_snp2(v))) {
				del_part[n_del_part++] = kh_key(g->h, k);
			}
		}

		// delete entries
		for (int d = 0; d < n_del_part; ++d) {
			k = pg_ht_get(g->h, del_part[d]);
    		if (k != kh_end(g->h)) {
				pg_ht_del(g->h, k);
			}
		}
		n_del += n_del_part;
		free(del_part);
        // pthread_mutex_unlock(&g->lock);
    }

	return n_del;
}

int pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int filt)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	pg_ht1_t *g;
	if (n == 0) return 0;

	g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers.
	pthread_mutex_lock(&g->lock);
	
	for (j = 0; j < n; ++j) {
		int absent;
		uint64_t key;
		uint32_t cb = a[j].cb;
		key = a[j].h_flanks >> h->pre;
		if ((a[j].h_flanks & mask) != (a[0].h_flanks & mask)) continue;

		if (filt) { // after filtering, do not insert k-mers
			khint_t k = pg_ht_get(g->h, key);
			absent = (k == kh_end(g->h));
			if (!absent) { 
				uint32_t v = kh_val(g->h, k);
				uint32_t cnt = val_count(v);
				uint32_t snp1, snp2, cb2 = val_cb2(v);
				if (val_snp1(v) ^ val_snp2(v)) { // already known as SNP, check if it is multi-allelic
					snp1 = 1;
					if (cb != val_cb1(v) && cb != val_cb2(v)) {
						snp2 = 1; // multi-allelic SNP
					} else {
						snp2 = 0; // bi-allelic SNP
					}
				} else {
					if (val_snp1(v) & val_snp2(v)) { // already known as multi-allelic SNP
						snp1 = 1; snp2 = 1; // already known as multi-allelic SNP
					} else if (cb != val_cb1(v)) { // newly identified SNP
						snp1 = 1; snp2 = 0;
						cb2 = cb; // store the second central base
					} else { // still non-SNP
						snp1 = 0; snp2 = 0;
					}
				}
				kh_val(g->h, k) = val_pack(cnt < COUNTER_MAX ? cnt + 1 : cnt, snp2, snp1, cb2, val_cb1(v));
			}
		} else { // before filtering, insert or update k-mers
			khint_t k = pg_ht_put(g->h, key, &absent);
			if (absent) {
				++n_ins;
				kh_val(g->h, k) = val_pack(1, 0, 0, 0, cb);  // first occurrence, SNP unknown
			} else {
				uint32_t v = kh_val(g->h, k);
				uint32_t cnt = val_count(v);
				uint32_t snp1, snp2, cb2 = 0;
				if (val_snp1(v) ^ val_snp2(v)) { // already known as SNP, check if it is multi-allelic
					snp1 = 1;
					if (cb != val_cb1(v) && cb != val_cb2(v)) {
						snp2 = 1; // multi-allelic SNP
					} else {
						snp2 = 0; // bi-allelic SNP
					}
				} else {
					if (val_snp1(v) & val_snp2(v)) { // already known as multi-allelic SNP
						snp1 = 1; snp2 = 1; // already known as multi-allelic SNP
					} else if (cb != val_cb1(v)) { // newly identified SNP
						snp1 = 1; snp2 = 0;
						cb2 = cb; // store the second central base
					} else { // still non-SNP
						snp1 = 0; snp2 = 0;
					}
				}
				kh_val(g->h, k) = val_pack(cnt < COUNTER_MAX ? cnt + 1 : cnt, snp2, snp1, cb2, val_cb1(v));
			}
		}
	}

	pthread_mutex_unlock(&g->lock);
	
	return n_ins;
}


void pg_mht_tighten(pg_mht_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_t *g = h->h[i].h;
		if (kh_size(g) * 3 < kh_capacity(g))
			pg_ht_m_resize(g, kh_size(g) * 3);
	}
}


void pg_mht_dump(const pg_mht_t *h, const char *fn)
{
	FILE *fp;
	char seq[h->k+1];
	int i;
	uint64_t hash_mask = (1ULL<<((h->k-1)*2)) - 1; // to hash only the flanks
	if ((fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout) == 0) return -1;
	fprintf(fp, "kmer\tcount\n");
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_t *g = h->h[i].h;
		khint_t k;
		for (k = 0; k < kh_end(g); ++k)
			if (kh_exist(g, k)) {
				// reverse hash and obtain kmers
				uint64_t flanks = pg_hash64_inv(kh_key(g, k) << h->pre | i, hash_mask); // reconstruct the full flanks
				uint64_t v = kh_val(g, k);

				int mid = h->k >> 1; // 15 for k=31
				// right flank
				for (int j = 0; j < mid; ++j)
					seq[h->k-1-j] = nt4_seq_table[(flanks >> (j * 2)) & 3];
				// central base
				seq[mid] = nt4_seq_table[val_cb1(v)];
				// left flank
				for (int j = 0; j < mid; ++j)
					seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];
				seq[h->k] = '\0';
				
				char cb1 = nt4_seq_table[val_cb1(v)];
				char cb2 = nt4_seq_table[val_cb2(v)];
				uint64_t cnt = val_count(v);

				fprintf(fp, "%s\t%c|%c\t%llu\n", seq, cb1, cb2, (unsigned long long)cnt);
			}
	}
	fprintf(stderr, "[M::%s] Dumped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
}