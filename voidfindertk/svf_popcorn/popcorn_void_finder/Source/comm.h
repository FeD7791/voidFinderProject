#ifndef COMM_H_
#define COMM_H_

#include "grid.h"
#include "objects.h"
#include <vector>

void popcorn_rma(int rank, int nshell, float min_radius, double eps, float limit, grid::Grid *G);
void svf_rma(int rank, float limseed, float rhobar, int nstep, grid::Grid *G, float min_radius);
std::vector<grid::popcorn> gather_voids(int rank, int nproc, const std::vector<grid::popcorn> &voids);
std::vector<grid::Sphere> gather_sphvds(int rank, int nproc, const std::vector<grid::Sphere> &voids);

#endif
