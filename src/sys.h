#ifndef SYS_H
#define SYS_H

double pg_cputime(void);
void pg_reset_realtime(void);
double pg_realtime(void);
long pg_peakrss(void);

#endif // SYS_H