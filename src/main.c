#include <stdio.h>

#include "htab.h"
#include "ketopt.h"
#include "sys.h"
#include "utils.h"

int main_count(int argc, char *argv[])
{   
    pg_mht_t *h = 0;
	char *fn_in = 0, *fn_out = 0;
	int c, i, min_cnt = 1, max_cnt = 1, max_out = 0, check_n = 10, pre_resize = 0;
	pg_opt_t opt;
	ketopt_t o = KETOPT_INIT;
	pg_opt_init(&opt);
	opt.chunk_size = mm_parse_num("1.9g");
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:i:o:c:x:e:s:r", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'c') min_cnt = atoi(o.arg);
		else if (c == 'x') max_cnt = atoi(o.arg);
		else if (c == 'e') max_out = atoi(o.arg);
		else if (c == 's') check_n = atoi(o.arg);
		else if (c == 'r') pre_resize = 1;
		else if (c == 'p') opt.pre = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = mm_parse_num(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'i') fn_in = o.arg;
		else if (c == 'o') fn_out = o.arg;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: yak cntasm [options] <in1.fa> [in2.fa [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -c INT     min count [%d]\n", min_cnt);
		fprintf(stderr, "  -x INT     max count [%d]\n", max_cnt);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -r         resize before merging; use if merging is slow\n");
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -e INT     exclude a k-mer if absent from INT samples [%d]\n", max_out);
		fprintf(stderr, "  -s INT     shrink the hash table every INT samples [%d]\n", check_n);
		fprintf(stderr, "  -K INT     chunk size [1.9g]\n");
		fprintf(stderr, "  -i FILE    input k-mer dump []\n");
		fprintf(stderr, "  -o FILE    output k-mer dump []\n");
		fprintf(stderr, "Note: if input and output file names are identical, input is overwritten\n");
		return 1;
	}
	if (opt.pre < COUNTER_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 32) {
		fprintf(stderr, "ERROR: -k must be <=31\n");
		return 1;
	}

    h = pg_count(&argv[o.ind], argc - o.ind, &opt); // count k-mers in the input files argv

    // fprintf(stderr, "[M::%s::%.3f*%.2f] processed all files; %ld distinct k-mers in the hash table\n", __func__,
	// 			pg_realtime(), pg_cputime() / pg_realtime(), (long)h->counts);

    pg_mht_tighten(h);

    if (fn_out) pg_mht_dump(h, fn_out);
	pg_mht_destroy(h);

    return 0;
}


int main(int argc, char *argv[])
{   
    int ret = 0;
    pg_reset_realtime();
    if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);

    if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PG_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, pg_realtime(), pg_cputime(), pg_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}

    return 0;
}