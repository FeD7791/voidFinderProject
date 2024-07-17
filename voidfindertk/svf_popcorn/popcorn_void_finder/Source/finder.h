#ifndef FINDER_H
#define FINDER_H

#include <string>
#include <vector>

#include "grid.h"
#include "objects.h"

int raw_finder(int ic, int jc, int kc, float limit, float rhobar, int nstep, grid::Grid *g, float rad_min);

void FindCenters_fast(float limseed, float rhobar, int nstep, grid::Grid *G, float min_radius);

void FindCenters(float limseed);

void FindSphvds(float limit, std::string rawfname);

void CleanSphvds_Tol0(grid::Grid *G, float maxscale);

void CleanSphvds(grid::Grid *G, float maxscale);

void CleanSphvds_Tol0_opt(grid::Grid *G, float maxscale);

void FindSphvds_fast(float limit, std::string fname, int nshell, grid::Grid *G);

void FindSphvds_fastV2(float limit, std::string fname, int nshell, grid::Grid *G, float min_radius);

float xper(float x, float boxsize);

float logfactorial(int n);

void svf_work(int wrk, float limseed, float rhobar, int nstep, grid::Grid *G, float min_radius);

#endif
