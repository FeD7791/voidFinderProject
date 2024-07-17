#include "comm.h"

#include <mpi.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include "finder.h"
#include "grid.h"
#include "objects.h"
#include "popcorn.h"
#include <omp.h>

#if MPI_VERSION < 3
#error Popcorn requires MPI 3
#endif

#define CHECK_MPI(MPICALL)                                                                                             \
  {                                                                                                                    \
    int status = MPICALL;                                                                                              \
    if (status != MPI_SUCCESS) {                                                                                       \
      fprintf(stderr, "MPI call failed with %d: " #MPICALL "\n", status);                                              \
      MPI_Finalize();                                                                                                  \
      exit(-1);                                                                                                        \
    }                                                                                                                  \
  }

template <class WorkFunction, class SkipPredicate>
static void distribute_work_rma(int rank, int from, int to, WorkFunction work, SkipPredicate skip) {
  MPI_Win window;
  int *window_ptr;
  MPI_Aint window_size = rank == 0 ? sizeof(int) : 0;

  CHECK_MPI(MPI_Win_allocate(window_size, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window_ptr, &window));

  int curr = from;
  int total = to;
  if (rank == 0) {
    int discard = 0;
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, window);
    CHECK_MPI(MPI_Fetch_and_op(&from, &discard, MPI_INT, 0, 0, MPI_REPLACE, window));
    MPI_Win_unlock(0, window);
  }
  CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

  while (curr < total) {
    int one = 1;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, window);
    CHECK_MPI(MPI_Fetch_and_op(&one, &curr, MPI_INT, 0, 0, MPI_SUM, window));
    MPI_Win_unlock(0, window);

    if (curr < total) {
      if (skip(curr)) {
        std::cout << "Rank " << rank << " skipped #" << curr << "\n";
      } else {
        if (curr % 100 == 0) {
          std::cout << "Rank " << rank << " working on #" << curr << "\n";
        }
        work(curr);
      }
    }
  }
  CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
  CHECK_MPI(MPI_Win_free(&window));
}

template <class WorkFunction> static void distribute_work_2(int rank, int from, int to, WorkFunction work) {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> localWorkIndices;
  int lw = (to - from) / size;
  int localWork = lw;
  if (rank == (size - 1)) {
    int remainder = (to - from) % size;
    localWork += remainder;
  }

  for (int i = 0; i < localWork; i++) {
    localWorkIndices.push_back(i + rank * lw);
  }

#pragma omp parallel for default(none) schedule(dynamic) shared(std::cout, to, from, rank, localWorkIndices, work)
  for (int index : localWorkIndices) {
    if (index % (int)((to - from) / 100) == 0) {
      std::cout << "Rank " << rank << " thread " << omp_get_thread_num() << " working on #" << index << " of "
                << (to - from) << std::endl;
    }
    work(index);
  }
}

void popcorn_rma(int rank, int nshell, float min_radius, double eps, float limit, grid::Grid *G) {
  // distribute_work_rma(
  //     rank, 0, Sphvd.size(), [=](int i) { popcorn_work(i, nshell, min_radius, eps, limit, G); },
  //     [&](int i) -> bool { return Sphvd[i].erase; });
  distribute_work_2(rank, 0, Sphvd.size(), [=](int i) { popcorn_work(i, nshell, min_radius, eps, limit, G); });
}

void svf_rma(int rank, float limseed, float rhobar, int nstep, grid::Grid *G, float min_radius) {
  int Ngrid = (*G).Ngrid;
  int Ngrid3 = Ngrid * Ngrid * Ngrid;

  distribute_work_2(rank, 0, Ngrid3, [=](int i) { svf_work(i, limseed, rhobar, nstep, G, min_radius); });
}
struct popcorn_novec {
  int id = -1;
  int nmem = 0;
  int num_in_halos = 0;
  int npart = 0;
  int isCave = -1;
  double vol = 0.0;

  popcorn_novec() : id(-1), nmem(0), num_in_halos(0), npart(0), isCave(-1), vol(0.0) {}

  popcorn_novec(const grid::popcorn &p)
      : id(p.id), nmem(p.nmem), num_in_halos(p.in_halos.size()), npart(p.npart), isCave(p.isCave), vol(p.vol) {}

  bool operator==(const popcorn_novec &other) const {
    return (id == other.id) && (nmem == other.nmem) && (num_in_halos == other.num_in_halos) && (npart == other.npart) &&
           (isCave == other.isCave) && (vol == other.vol);
  }
};

static void flatten_voids(const std::vector<grid::popcorn> &voids, std::vector<popcorn_novec> &voids_novec,
                          std::vector<int> &all_in_halos, std::vector<grid::Sphere> &all_membs) {

  for (auto &p : voids) {
    voids_novec.emplace_back(p);
    all_in_halos.insert(all_in_halos.end(), p.in_halos.begin(), p.in_halos.end());
    all_membs.insert(all_membs.end(), p.membs.begin(), p.membs.end());
  }
}

static std::vector<grid::popcorn> reconstruct_voids(const std::vector<popcorn_novec> &voids_novec,
                                                    const std::vector<int> &all_in_halos,
                                                    const std::vector<grid::Sphere> &all_membs) {
  std::vector<grid::popcorn> result;

  int membs_pos = 0;
  int in_halos_pos = 0;
  for (auto &p : voids_novec) {
    grid::popcorn tmp;
    tmp.id = p.id;
    tmp.vol = p.vol;
    tmp.npart = p.npart;
    tmp.isCave = p.isCave;
    tmp.nmem = p.nmem;
    for (int m = 0; m < p.nmem; ++m) {
      tmp.membs.push_back(all_membs[membs_pos]);
      ++membs_pos;
    }
    for (int h = 0; h < p.num_in_halos; ++h) {
      tmp.in_halos.push_back(all_in_halos[in_halos_pos]);
      ++in_halos_pos;
    }
    result.push_back(tmp);
  }
  return result;
}

std::vector<grid::popcorn> gather_voids(int rank, int nproc, const std::vector<grid::popcorn> &voids) {
  // Flatten popcorns into 3 arrays
  std::vector<popcorn_novec> voids_novec;
  std::vector<int> all_in_halos;
  std::vector<grid::Sphere> all_membs;
  flatten_voids(voids, voids_novec, all_in_halos, all_membs);

  // Send counts to rank 0
  int my_voids = voids_novec.size();
  int my_in_halos = all_in_halos.size();
  int my_membs = all_membs.size();

  std::vector<int> per_rank_num_voids(nproc);
  std::vector<int> per_rank_num_in_halos(nproc);
  std::vector<int> per_rank_num_membs(nproc);
  CHECK_MPI(MPI_Gather(&my_voids, 1, MPI_INT, per_rank_num_voids.data(), 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_MPI(MPI_Gather(&my_in_halos, 1, MPI_INT, per_rank_num_in_halos.data(), 1, MPI_INT, 0, MPI_COMM_WORLD));
  CHECK_MPI(MPI_Gather(&my_membs, 1, MPI_INT, per_rank_num_membs.data(), 1, MPI_INT, 0, MPI_COMM_WORLD));

  // Calculate per-rank starting displacement and total elements
  std::vector<int> per_rank_void_starts;
  std::vector<int> per_rank_in_halos_starts;
  std::vector<int> per_rank_membs_starts;
  if (rank == 0) {
    std::exclusive_scan(per_rank_num_voids.begin(), per_rank_num_voids.end(), std::back_inserter(per_rank_void_starts),
                        0);
    std::exclusive_scan(per_rank_num_in_halos.begin(), per_rank_num_in_halos.end(),
                        std::back_inserter(per_rank_in_halos_starts), 0);
    std::exclusive_scan(per_rank_num_membs.begin(), per_rank_num_membs.end(), std::back_inserter(per_rank_membs_starts),
                        0);

    int total_voids = std::reduce(per_rank_num_voids.begin(), per_rank_num_voids.end());
    int total_in_halos = std::reduce(per_rank_num_in_halos.begin(), per_rank_num_in_halos.end());
    int total_membs = std::reduce(per_rank_num_membs.begin(), per_rank_num_membs.end());
    voids_novec.resize(total_voids);
    all_in_halos.resize(total_in_halos);
    all_membs.resize(total_membs);
  }

  // Gather popcorns in rank 0
  int mpi_popcorn_count = 6;
  int mpi_popcorn_blocklengths[] = {1, 1, 1, 1, 1, 1};
  MPI_Aint mpi_popcorn_displacements[] = {offsetof(popcorn_novec, id),           offsetof(popcorn_novec, nmem),
                                          offsetof(popcorn_novec, num_in_halos), offsetof(popcorn_novec, npart),
                                          offsetof(popcorn_novec, isCave),       offsetof(popcorn_novec, vol)};
  MPI_Datatype mpi_popcorn_types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
  MPI_Datatype mpi_popcorn_novec;
  CHECK_MPI(MPI_Type_create_struct(mpi_popcorn_count, mpi_popcorn_blocklengths, mpi_popcorn_displacements,
                                   mpi_popcorn_types, &mpi_popcorn_novec));
  CHECK_MPI(MPI_Type_commit(&mpi_popcorn_novec));

  CHECK_MPI(MPI_Gatherv(voids_novec.data(), my_voids, mpi_popcorn_novec, voids_novec.data(), per_rank_num_voids.data(),
                        per_rank_void_starts.data(), mpi_popcorn_novec, 0, MPI_COMM_WORLD));

  CHECK_MPI(MPI_Type_free(&mpi_popcorn_novec));

  // Gather in_halos in rank 0
  CHECK_MPI(MPI_Gatherv(all_in_halos.data(), my_in_halos, MPI_INT, all_in_halos.data(), per_rank_num_in_halos.data(),
                        per_rank_in_halos_starts.data(), MPI_INT, 0, MPI_COMM_WORLD));

  // Gather spheres in rank 0
  int mpi_sphere_count = 12;
  int mpi_sphere_blocklengths[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Aint mpi_sphere_displacements[] = {
      offsetof(grid::Sphere, x),      offsetof(grid::Sphere, y),     offsetof(grid::Sphere, z),
      offsetof(grid::Sphere, radius), offsetof(grid::Sphere, delta), offsetof(grid::Sphere, vol),
      offsetof(grid::Sphere, lvl),    offsetof(grid::Sphere, id),    offsetof(grid::Sphere, erase),
      offsetof(grid::Sphere, ToF),    offsetof(grid::Sphere, np),    offsetof(grid::Sphere, irad)};
  MPI_Datatype mpi_sphere_types[] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,    MPI_FLOAT,    MPI_FLOAT, MPI_FLOAT,
                                     MPI_INT,   MPI_INT,   MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_INT,   MPI_INT};
  MPI_Datatype mpi_sphere;
  CHECK_MPI(MPI_Type_create_struct(mpi_sphere_count, mpi_sphere_blocklengths, mpi_sphere_displacements,
                                   mpi_sphere_types, &mpi_sphere));
  CHECK_MPI(MPI_Type_commit(&mpi_sphere));

  CHECK_MPI(MPI_Gatherv(all_membs.data(), my_membs, mpi_sphere, all_membs.data(), per_rank_num_membs.data(),
                        per_rank_membs_starts.data(), mpi_sphere, 0, MPI_COMM_WORLD));

  CHECK_MPI(MPI_Type_free(&mpi_sphere));

  // Reconstruct popcorn vector in rank 0 from gathered arrays
  std::vector<grid::popcorn> result;
  if (rank == 0) {
    result = reconstruct_voids(voids_novec, all_in_halos, all_membs);
  }
  // non-0 ranks return an empty vector
  return result;
}

std::vector<grid::Sphere> gather_sphvds(int rank, int nproc, const std::vector<grid::Sphere> &voids) {

  std::vector<grid::Sphere> voids_novec(voids);
  // Send counts to rank 0
  int my_voids = voids_novec.size();

  std::vector<int> per_rank_num_voids(nproc);
  CHECK_MPI(MPI_Gather(&my_voids, 1, MPI_INT, per_rank_num_voids.data(), 1, MPI_INT, 0, MPI_COMM_WORLD));

  // Calculate per-rank starting displacement and total elements
  std::vector<int> per_rank_void_starts;
  if (rank == 0) {
    std::exclusive_scan(per_rank_num_voids.begin(), per_rank_num_voids.end(), std::back_inserter(per_rank_void_starts),
                        0);

    int total_voids = std::reduce(per_rank_num_voids.begin(), per_rank_num_voids.end());
    voids_novec.resize(total_voids);
  }

  // Gather spheres in rank 0
  int mpi_sphere_count = 12;
  int mpi_sphere_blocklengths[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Aint mpi_sphere_displacements[] = {
      offsetof(grid::Sphere, x),      offsetof(grid::Sphere, y),     offsetof(grid::Sphere, z),
      offsetof(grid::Sphere, radius), offsetof(grid::Sphere, delta), offsetof(grid::Sphere, vol),
      offsetof(grid::Sphere, lvl),    offsetof(grid::Sphere, id),    offsetof(grid::Sphere, erase),
      offsetof(grid::Sphere, ToF),    offsetof(grid::Sphere, np),    offsetof(grid::Sphere, irad)};
  MPI_Datatype mpi_sphere_types[] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,    MPI_FLOAT,    MPI_FLOAT, MPI_FLOAT,
                                     MPI_INT,   MPI_INT,   MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_INT,   MPI_INT};
  MPI_Datatype mpi_sphere;
  CHECK_MPI(MPI_Type_create_struct(mpi_sphere_count, mpi_sphere_blocklengths, mpi_sphere_displacements,
                                   mpi_sphere_types, &mpi_sphere));
  CHECK_MPI(MPI_Type_commit(&mpi_sphere));

  CHECK_MPI(MPI_Gatherv(voids_novec.data(), my_voids, mpi_sphere, voids_novec.data(), per_rank_num_voids.data(),
                        per_rank_void_starts.data(), mpi_sphere, 0, MPI_COMM_WORLD));

  CHECK_MPI(MPI_Type_free(&mpi_sphere));

  return voids_novec;
}
