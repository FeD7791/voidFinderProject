#include <mpi.h>

#include <sstream>

#include "allvars.h"
#include "colors.h"
#include "comm.h"
#include "constants.h" //  K_OUTFILE_INTERSECS K_THRSLD K_CRTP K_EPS
#include "grid.h"
#include "io_lib/io.h"
#include "objects.h"
#include "popcorn.h"
#include "t.h"

// POPFILE K_POPFILE_COMPUTE_INTERSECS
#include "configuration_lib/converter.h"
#include "configuration_lib/variables_manager.h"
#include "configuration_lib/variables_tools.h"

using namespace std;
using namespace iate_variables;

int popcorn_mpi_start(float min_radius, double eps, float limit, int rank) {
  // int Ngrid=512;
  int Ngrid = find_ngrid(Tracer.size(), min_radius, boxsize.x);
  int nshell = 30;

  int nvds = Sphvd.size();
  // Cleaning small voids
  float radmax = -1.0;
  for (int i = 0; i < nvds; i++) {
    if (radmax < Sphvd[i].radius) {
      radmax = Sphvd[i].radius;
    }
  }

  float maxscale = radmax * 1.1;

  for (long unsigned int ii = 0; ii < Tracer.size(); ii++) {
    int i = int(Tracer[ii].x / boxsize.x * Ngrid);
    if (unlikely(i > (Ngrid - 1))) {
      i = Ngrid - 1;
    }
    int j = int(Tracer[ii].y / boxsize.y * Ngrid);
    if (unlikely(j > (Ngrid - 1))) {
      j = Ngrid - 1;
    }
    int k = int(Tracer[ii].z / boxsize.z * Ngrid);
    if (unlikely(k > (Ngrid - 1))) {
      k = Ngrid - 1;
    }
    int ind = i + Ngrid * (j + Ngrid * k);
    Tracer[ii].gbin = ind;
  }

  cout << "sorting particles by grid..." << endl;
  cout.flush();
  sort(Tracer.begin(), Tracer.end(), grid::small_grid_index_first());
  cout << "done." << endl;

  grid::Grid *G = new grid::Grid(Ngrid, Tracer.size(), boxsize, maxscale);
  cout << "placing particles into the grid... " << Ngrid << endl;
  cout.flush();
  G->build(&Tracer);
  cout << "done." << endl;
  cout.flush();

  for (int i = 0; i < nvds; i++) {
    if (Sphvd[i].radius < G->maxscale / (float)nshell) {
      Sphvd[i].erase = true;
    }
  }

  int cc = 0;
  for (int i = 0; i < nvds; i++) {
    if (!Sphvd[i].erase) {
      cc++;
    }
  }

  cout << LGREEN << "total de candidatos: " << cc << DEFA << endl;
  /// Main computation /////////////////////////////////////////////////////
  cout << CYAN << "procesando.... " << DEFA << endl;

  clock_t t = clock();

  popcorn_rma(rank, nshell, min_radius, eps, limit, G);

  Time(t, NCORES);

  return 0;
}

int main(int argc, char *argv[]) {
  // #ifdef MPI_VERSION
  //
  MPI_Init(&argc, &argv);

  int numprocs = 1;
  int proc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  cout << "MPI: proc " << proc << " of " << numprocs;
  // #endif

  // initialize instance
  VariablesManager *vars = VariablesManager::getInstance();
  // load variables
  loadConfigurationVars(argc, argv);
  // float boxsize = Converter::stringToFloat(vars->valueFromKey(K_BOXSIZE));
  double eps = Converter::stringToDouble(vars->valueFromKey(K_EPS));
  string fname_tracer = vars->valueFromKey(K_TRSFILE);
  string sph_file = vars->valueFromKey(K_SPHFILE);
  string fmt = vars->valueFromKey(K_FILEFMT);
  int num_files = Converter::stringToInt(vars->valueFromKey(K_NUM_FILE));
  if (!fmt.compare("ASCII") || !fmt.compare("Minerva")) {
    float boxlength = Converter::stringToFloat(vars->valueFromKey(K_BOXSIZE));
    boxsize.x = boxlength;
    boxsize.y = boxlength;
    boxsize.z = boxlength;
  }

  if (!fmt.compare("HDF5_SUBFIND_GROUPS")) {
    masslim = Converter::stringToFloat(vars->valueFromKey(K_MASSMIN));
  }

  if (!fmt.compare("HDF5_SUBFIND_SUBHALOS")) {
    masslim = Converter::stringToFloat(vars->valueFromKey(K_MASSMIN));
    cout << "minimum mass= " << masslim << endl;
  }
  float min_radius = Converter::stringToFloat(vars->valueFromKey(K_MINRADIUS));

  // input parameters
  // underdensity parameter (typical values -0.9 o -0.8)
  float densth = Converter::stringToFloat(vars->valueFromKey(K_DENSTH));
  radmax = Converter::stringToFloat(vars->valueFromKey(K_MAXRADIUS));

  cout << "box....." << endl;
  string rawpop_file = vars->valueFromKey(K_RAWPOPFILE);
  // file with the sperical voids to be used as seeds

  // Lectura de datos
  ReadHeaderTracers(fmt, fname_tracer, num_files);
  ReadTracers(fmt, fname_tracer, num_files);
  cout << "x) boxsize: " << boxsize.x << endl;
  cout << "y) boxsize: " << boxsize.y << endl;
  cout << "z) boxsize: " << boxsize.z << endl;

  Read_sv(sph_file, boxsize);
  popcorn_mpi_start(min_radius, eps, densth, proc);

  // gather results on rank 0
  auto all_voids = gather_voids(proc, numprocs, pop_voids);
  pop_voids = all_voids;

  if (proc == 0) {
    Write_raw_popcorn(rawpop_file);
  }
  MPI_Finalize();

  return 0;
}
