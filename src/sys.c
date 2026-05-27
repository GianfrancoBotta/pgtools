#include <sys/resource.h>
#include <sys/time.h>
#include <stdio.h>
#include "sys.h"

int pg_verbose = 3;

static double pg_realtime0;

double pg_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline double pg_realtime_core(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void pg_reset_realtime(void)
{
	pg_realtime0 = pg_realtime_core();
}

double pg_realtime(void)
{
	return pg_realtime_core() - pg_realtime0;
}

long pg_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}