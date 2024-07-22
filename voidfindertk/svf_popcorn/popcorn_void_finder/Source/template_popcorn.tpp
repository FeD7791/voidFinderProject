#include <cassert>
#include <vector>

#include "objects.h"

namespace grid {

template <typename T> int build_vecinos(int nside, T *rts, std::vector<grid::vecv> *vecinos, float3 box_size) {
  int ix, iy, iz, iii;
  int ivecino, iix, iiy, iiz;
  int indx, indy, indz;
  int i, j, k;
  typename T::iterator it;
  std::vector<int> ind;
  int cnt;
  class cell dumm;
  std::vector<class cell> grid;
  std::vector<int>::iterator nei;

  dumm = cell(-1, -1, -1);

  cnt = 0;
  for (i = 0; i < (nside * nside * nside); i++)
    grid.push_back(cell());

  for (i = 0; i < (*rts).size(); i++)
    (*vecinos).push_back(vecv());

  for (i = 0; i < (*rts).size(); i++) {

    for (j = 0; j < ((*rts)[i]).membs.size(); j++) {
      if (((*rts)[i]).membs[j].x < 0)
        ((*rts)[i]).membs[j].x = +box_size.x;
      if (((*rts)[i]).membs[j].x > box_size.x)
        ((*rts)[i]).membs[j].x = -box_size.x;

      if (((*rts)[i]).membs[j].y < 0)
        ((*rts)[i]).membs[j].y = +box_size.y;
      if (((*rts)[i]).membs[j].y > box_size.y)
        ((*rts)[i]).membs[j].y = -box_size.y;

      if (((*rts)[i]).membs[j].z < 0)
        ((*rts)[i]).membs[j].z = +box_size.z;
      if (((*rts)[i]).membs[j].z > box_size.z)
        ((*rts)[i]).membs[j].z = -box_size.z;
    }
  }
  // cargando particulas

  for (i = 0; i < (*rts).size(); i++) {
    ind = which_cell_pop((*rts)[i], nside, box_size);
    for (j = 0; j < ind.size(); j++)
      grid[ind[j]].in.push_back(i);

    ind.clear();
  }

  for (i = 0; i < grid.size(); i++)
    sort((grid[i].in).begin(), (grid[i].in).end());

  for (i = 0; i < grid.size(); i++)
    for (j = 0; j < grid[i].in.size(); j++)
      for (k = 0; k < grid[i].in.size(); k++)
        (*vecinos)[grid[i].in[k]].v.push_back(grid[i].in[j]);

  for (i = 0; i < (*vecinos).size(); i++) {
    sort((((*vecinos)[i]).v).begin(), (((*vecinos)[i]).v).end());
    nei = std::unique(((*vecinos)[i].v).begin(), ((*vecinos)[i].v).end());
    ((*vecinos)[i].v).resize(std::distance(((*vecinos)[i].v).begin(), nei));
  }
  return (0);
}

template <typename T> int which_cell(T it, int nside, float3 box_size) {
  int ix, iy, iz, ind;
  double facbinx = (double)nside / box_size.x;
  double facbiny = (double)nside / box_size.y;
  double facbinz = (double)nside / box_size.z;
  float xc, yc, zc;
  xc = it.x;
  yc = it.y;
  zc = it.z;
  ix = (int)floor((double)(xc)*facbinx);
  iy = (int)floor((double)(yc)*facbiny);
  iz = (int)floor((double)(zc)*facbinz);

  if (xc == box_size.x)
    ix = 0;
  if (yc == box_size.y)
    iy = 0;
  if (zc == box_size.z)
    iz = 0;

  assert(ix < nside && ix >= 0);
  assert(iy < nside && iy >= 0);
  assert(iz < nside && iz >= 0);
  ind = index3(ix, iy, iz, nside);
  return ind;
}

template <typename T> int build_grid(int nside, T *rts, std::vector<class cell> *grid, float3 box_size) {
  int ind, ix, iy, iz;
  int ivecino, iix, iiy, iiz;
  int indx, indy, indz;
  int i;
  typename T::iterator it;

  // inicializando
  for (ix = 0; ix < nside; ix++)
    for (iy = 0; iy < nside; iy++)
      for (iz = 0; iz < nside; iz++) {
        (*grid).push_back(cell());
      }

  for (ix = 0; ix < nside; ix++)
    for (iy = 0; iy < nside; iy++)
      for (iz = 0; iz < nside; iz++) {
        ind = index3(ix, iy, iz, nside);
        (*grid)[ind] = (cell(ix, iy, iz));
      }

  // cargando particulas
  for (i = 0; i < (*rts).size(); i++) {
    ind = which_cell((*rts)[i], nside, box_size);
    (*grid)[ind].in.push_back(i);
  }

  // apuntando a vecinos
  for (ix = 0; ix < nside; ix++)
    for (iy = 0; iy < nside; iy++)
      for (iz = 0; iz < nside; iz++) {
        i = 0;
        ind = index3(ix, iy, iz, nside);
        for (iix = ix - 1; iix <= ix + 1; iix++)
          for (iiy = iy - 1; iiy <= iy + 1; iiy++)
            for (iiz = iz - 1; iiz <= iz + 1; iiz++) {
              ivecino = index3(iix, iiy, iiz, nside);
              (*grid)[ind].next[i] = ivecino;
              i++;
            }
      }

  return (0);
}

} // namespace grid
