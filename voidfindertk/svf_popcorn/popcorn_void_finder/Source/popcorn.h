#ifndef POPCORN_H_
#define POPCORN_H_

#include "grid.h"

int popcorn_start(float min_radius, double eps, float limit, int rank, int numrank);
void popcorn_work(int cnt, int nshell, float min_radius, double eps, float limit, grid::Grid *G);

#endif