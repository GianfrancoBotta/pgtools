#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>

#define PG_MAGIC "PGKM"
#define VAL_INFO_BITS 6
#define COUNTER_BITS (32 - VAL_INFO_BITS)
#define COUNTER_MAX ((1U << COUNTER_BITS) - 1)

#define val_count(v) ((v) >> VAL_INFO_BITS)
#define val_cb1(v) ((v) & 0x3U)
#define val_cb2(v) (((v) >> 2) & 0x3U)
#define val_snp1(v) (((v) >> 4) & 0x1U)
#define val_snp2(v) (((v) >> 5) & 0x1U)
#define val_pack(cnt, snp2, snp1, cb2, cb1) (((uint32_t)(cnt) << VAL_INFO_BITS) | ((snp2) << 5) | ((snp1) << 4) | ((cb2) << 2) | (cb1))

typedef struct __attribute__((packed)) {
	uint64_t h_flanks;
	uint8_t cb; // stores central base of k-mer and SNP information
} ch_seq_t;

typedef struct { // terminal options
    int64_t chunk_size;
    int64_t n_threads;
    double min_freq;
	int32_t k;
    int32_t pre; // number of bits for partitioning.
    int verbose;
} pg_opt_t;

struct pg_ht_t; // see khashl.h and htab.c for the definition of pg_ht_t.

typedef struct { 
    struct pg_ht_t *h; // hash table for each bucket.
    pthread_mutex_t lock;
} pg_ht1_t;

typedef struct {
    pg_ht1_t *h; // array of partitions (size = 1 << pre).
    uint64_t counts; // stores the total k-mer counts in pangenome.
    int32_t k; // k-mers length.
    int32_t pre; // stores the k-mer flanking sequences, encoded as 2 bits per base (used for partitioning).
} pg_mht_t;

// typedef struct __attribute__((packed)) {
//     uint64_t kmer; // stores the k-mer flanking sequences, encoded as 2 bits per base, and the central base (biallelic SNPs).
//     uint16_t count; // stores the count of k-mer with central base a and central base b.
// } pg_snpmer_t;

pg_mht_t *pg_mht_init(int k, int pre);
void pg_mht_destroy(pg_mht_t *h);
int pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int filt);
long long pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq, int ff);
void pg_mht_tighten(pg_mht_t *h);
pg_mht_t *pg_count(const char **fns, int n_fns, const pg_opt_t *opt);
void pg_mht_dump(const pg_mht_t *h, const char *fn);

#endif // HTAB_H