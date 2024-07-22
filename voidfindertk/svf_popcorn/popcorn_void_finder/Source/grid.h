#ifndef GRID_H_
#define GRID_H_

#include <vector>

#include "allvars.h"
#include "objects.h"

struct cell_state {
  bool is_in;
  bool is_out;
};

namespace grid {
class Grid {
public:
  int Ngrid;
  unsigned long long int Npart;
  int nsten;
  float3 boxsize;
  float maxscale;
  int nsearch;
  float gridsize;

  float *vert_dist;
  float *edge_dist;

  long long int *first;
  unsigned long long int *cnt;
  long long int *next;

  Grid(int nbin, unsigned long long int nt, float3 box, float maxsc);

  template <typename T> void build(T *rts);
};

int indx_period(int ind, int ngrid);

cell_state **compute_steps(Grid *g, int nstep);

template <typename T> float distance(T p1, T p2, float3 boxsize);

template <typename T> int fast_counter(T centre, int irad, int nstep, Grid *g, std::vector<T> *pnts);

template <typename T> int slow_counter(T centre, float(rad), Grid *g, std::vector<T> *pnts);

template <typename T> int *raw_counter(T centre, int irad, int nstep, Grid *g, std::vector<T> *pnts);

// template<typename T>
// int raw_finder(int ic,int jc, int kc, float limit,float rhobar,int nstep,
// Grid* g, vector<T> *pnts,cell_state **status_rad,float rad_min);

template <typename T>
float_in fast_finder(T centre, float limit, float rhobar, int nstep, Grid *g, std::vector<T> *pnts, int indx);

template <typename T>
float_in ball_finder(T centre, float limit, float rhobar, int irad, int nstep, Grid *g, std::vector<T> *pnts);

template <typename T>
float_in shell_finder(T centre, float limit, float rhobar, int irad, int nstep, Grid *g, std::vector<T> *pnts);

template <typename T> bool is_touched(T centre, Grid *g, std::vector<T> *svptr, float maxscale);

template <typename T> std::vector<int> touch_list(T centre, Grid *g, std::vector<T> *others, float maxscale);

cell_state is_there(int shifti, int shiftj, int shiftk, float rad);
float periodicity_delta(float coord, float boxsize);

template <typename T>
std::vector<float_in> fast_ball(T centre, float rad, int nstep, Grid *g, std::vector<T> *pnts,
                                std::vector<bool> *visit);

float periodicity(float coord, float box_size);

float distance(float x1, float y1, float z1, float x2, float y2, float z2, float3 box_size);

int index3(int ix, int iy, int iz, int nside);

std::vector<Sphere> pochoclea(Sphere root, std::vector<Sphere> *pop, float3 box_size, double eps);

void check_inner_spheres(double **test_pop, int *sizes, double eps);

void check_inner_spheres_v2(double **test_pop, int *sizes);

bool inflate_fast(Sphere *sph, class popcorn *pope, Grid *g, int nstep, std::vector<bool> *visit,
                  std::vector<int> *in_halos, double *volret, int *nhret, float limit, float density,
                  std::vector<halo> &Tracer, float rad_min);

int init_flags(Sphere *centre, Grid *g, std::vector<halo> *pnts, std::vector<int> *in_halos);

std::vector<int> revertindex(int ind, int nside);

template <typename T> int build_vecinos(int nside, T *rts, std::vector<vecv> *vecinos, float box_size);

template <typename T> int build_grid(int nside, T *rts, std::vector<class cell> *grid, float box_size);

void check_boundary_splits(double *test_pop, int sizes, float3 boxsize);

void RSD_data_z();
} // namespace grid

int iper(int, int);
#include "grid.tpp"

#endif
