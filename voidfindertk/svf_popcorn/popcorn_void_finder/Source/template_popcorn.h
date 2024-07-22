#ifndef TEMPLATE_POPCORN_H_
#define TEMPLATE_POPCORN_H_

#include <vector>

#include "objects.h"

namespace grid {

template <typename T> int build_vecinos(int nside, T *rts, std::vector<grid::vecv> *vecinos, float3 box_size);

template <typename T> int which_cell(T it, int nside, float3 box_size);

template <typename T> int build_grid(int nside, T *rts, std::vector<class cell> *grid, float3 box_size);
} // namespace grid

#include "template_popcorn.tpp"

#endif