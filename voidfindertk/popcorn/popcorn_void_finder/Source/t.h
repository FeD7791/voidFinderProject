#ifndef T_H
#define T_H

#include <ctime>

int CountLines(char *filename);

void Time(clock_t ti, int N);

void Progress(int par, int tot);

long long int find_ngrid(long long int ntracers, float radmin, float boxsize);

#endif
