// Need to import static functions
#include "../comm.cpp"

#include "gtest/gtest.h"

static void initTestVoids(std::vector<grid::popcorn> &voids) {
  voids.clear();
  voids.resize(4);

  voids[0].id = 1;
  voids[0].vol = 1.0;
  voids[0].npart = 1;
  voids[0].isCave = 1;
  voids[0].nmem = 1;
  voids[0].membs.emplace_back(1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
  voids[0].in_halos.emplace_back(1);

  voids[1].id = 2;
  voids[1].vol = 2.0;
  voids[1].npart = 2;
  voids[1].isCave = 2;
  voids[1].nmem = 2;
  voids[1].membs.emplace_back(1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
  voids[1].membs.emplace_back(2.0f, 2.0f, 2.0f, 2.0f, 2.0f);
  voids[1].in_halos.emplace_back(1);
  voids[1].in_halos.emplace_back(2);

  voids[2].id = 3;
  voids[2].vol = 3.0;
  voids[2].npart = 3;
  voids[2].isCave = 3;
  voids[2].nmem = 0;

  voids[3].id = 4;
  voids[3].vol = 4.0;
  voids[3].npart = 4;
  voids[3].isCave = 4;
  voids[3].nmem = 1;
  voids[3].membs.emplace_back(1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
  voids[3].in_halos.emplace_back(1);
  voids[3].in_halos.emplace_back(2);
}

static void initFlattenedTestVoids(std::vector<popcorn_novec> &voids_novec, std::vector<int> &all_in_halos,
                                   std::vector<grid::Sphere> &all_membs) {

  voids_novec.clear();
  all_in_halos.clear();
  all_membs.clear();

  voids_novec.resize(4);

  voids_novec[0].id = 1;
  voids_novec[0].vol = 1.0;
  voids_novec[0].npart = 1;
  voids_novec[0].isCave = 1;
  voids_novec[0].nmem = 1;
  voids_novec[0].num_in_halos = 1;
  all_membs.emplace_back(1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
  all_in_halos.emplace_back(1);

  voids_novec[1].id = 2;
  voids_novec[1].vol = 2.0;
  voids_novec[1].npart = 2;
  voids_novec[1].isCave = 2;
  voids_novec[1].nmem = 2;
  voids_novec[1].num_in_halos = 2;

  all_membs.emplace_back(1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
  all_membs.emplace_back(2.0f, 2.0f, 2.0f, 2.0f, 2.0f);
  all_in_halos.emplace_back(1);
  all_in_halos.emplace_back(2);

  voids_novec[2].id = 3;
  voids_novec[2].vol = 3.0;
  voids_novec[2].npart = 3;
  voids_novec[2].isCave = 3;
  voids_novec[2].nmem = 0;
  voids_novec[2].num_in_halos = 0;

  voids_novec[3].id = 4;
  voids_novec[3].vol = 4.0;
  voids_novec[3].npart = 4;
  voids_novec[3].isCave = 4;
  voids_novec[3].nmem = 1;
  voids_novec[3].num_in_halos = 2;

  all_membs.emplace_back(1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
  all_in_halos.emplace_back(1);
  all_in_halos.emplace_back(2);
}

TEST(PopcornConversion, FlattensOk) {
  std::vector<grid::popcorn> voids;
  initTestVoids(voids);

  std::vector<popcorn_novec> target_voids_novec;
  std::vector<int> target_all_in_halos;
  std::vector<grid::Sphere> target_all_membs;
  initFlattenedTestVoids(target_voids_novec, target_all_in_halos, target_all_membs);

  std::vector<popcorn_novec> voids_novec;
  std::vector<int> all_in_halos;
  std::vector<grid::Sphere> all_membs;
  flatten_voids(voids, voids_novec, all_in_halos, all_membs);

  EXPECT_TRUE(voids_novec == target_voids_novec);
  EXPECT_TRUE(all_in_halos == target_all_in_halos);
  EXPECT_TRUE(all_membs == target_all_membs);
}

TEST(PopcornConversion, ReconstructsOk) {
  std::vector<popcorn_novec> voids_novec;
  std::vector<int> all_in_halos;
  std::vector<grid::Sphere> all_membs;
  initFlattenedTestVoids(voids_novec, all_in_halos, all_membs);

  std::vector<grid::popcorn> target_voids;
  initTestVoids(target_voids);

  std::vector<grid::popcorn> voids = reconstruct_voids(voids_novec, all_in_halos, all_membs);

  EXPECT_TRUE(voids == target_voids);
}

// Sanity check for one rank
TEST(PopcornMPI, GatherOnlyRoot) {
  int rank = -1;
  int nproc = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::vector<grid::popcorn> input_voids;
  if (rank == 0) {
    initTestVoids(input_voids);
  }
  std::vector<grid::popcorn> result = gather_voids(rank, nproc, input_voids);

  EXPECT_TRUE(result == input_voids);
}

// First and last ranks contribute
TEST(PopcornMPI, GatherTwoProcs) {
  int rank = -1;
  int nproc = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (nproc < 2) {
    GTEST_SKIP();
  }

  std::vector<grid::popcorn> target_voids;
  std::vector<grid::popcorn> input_voids;
  if (rank == 0) {
    initTestVoids(target_voids);

    initTestVoids(input_voids);
    input_voids.erase(input_voids.begin() + 1, input_voids.end());
  } else if (rank == nproc - 1) {
    initTestVoids(input_voids);
    input_voids.erase(input_voids.begin());
  }
  std::vector<grid::popcorn> result = gather_voids(rank, nproc, input_voids);

  EXPECT_TRUE(result == target_voids);
}

// First four ranks contribute
TEST(PopcornMPI, GatherFourProcs) {
  int rank = -1;
  int nproc = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (nproc < 4) {
    GTEST_SKIP();
  }

  std::vector<grid::popcorn> target_voids;
  std::vector<grid::popcorn> input_voids;
  if (rank == 0) {
    initTestVoids(target_voids);
  }

  if (rank < 4) {
    std::vector<grid::popcorn> tmp_voids;
    initTestVoids(tmp_voids);
    input_voids.push_back(tmp_voids[rank]);
  }
  std::vector<grid::popcorn> result = gather_voids(rank, nproc, input_voids);

  EXPECT_TRUE(result == target_voids);
}

TEST(PopcornMPI, DistributeOk) {
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int local_sum = 0;
  distribute_work_rma(
      rank, 0, 20, [&](int i) { local_sum += 1 << i; }, [](int i) { return i % 3 == 0; });

  int result = 0;
  MPI_Allreduce(&local_sum, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int target = 0b10110110110110110110;
  EXPECT_TRUE(result == target);
}