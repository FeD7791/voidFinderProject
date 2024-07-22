#include <mpi.h>

#include <omp.h>

#include <string>

#include "allvars.h" //boxsize in ascii case
#include "comm.h"
#include "constants.h" //  K_OUTFILE_INTERSECS K_THRSLD K_CRTP K_EPS
#include "finder.h"
// POPFILE K_POPFILE_COMPUTE_INTERSECS
#include "configuration_lib/converter.h"
#include "configuration_lib/variables_manager.h"
#include "configuration_lib/variables_tools.h"
#include "io_lib/io.h"
#include "svf.h"
#include "t.h"

using namespace std;
using namespace iate_variables;
using namespace grid;

int svf_mpi_start(float limit, float min_radius, float radmax, int rank) {
  Grid *G;
  long long int Ngrid = find_ngrid((long long int)Tracer.size(), min_radius, boxsize.x);
  // int Ngrid=12;
  int nshell = (int)(radmax / min_radius) + 1;
  float rhobar;

  cout << "using a grid of " << Ngrid << "**3 boxcells" << endl;
  cout << "number of tracers: " << NTRAC << " in  a box of " << boxsize.x << endl;
  cout.flush();

  // cargo el grid en el que está cada tracer
  for (unsigned long int ii = 0; ii < Tracer.size(); ii++) {
    long long int i = (long long int)(Tracer[ii].x / boxsize.x * Ngrid);
    if (unlikely(i > (Ngrid - 1)))
      i = Ngrid - 1;
    long long int j = (long long int)(Tracer[ii].y / boxsize.y * Ngrid);
    if (unlikely(j > (Ngrid - 1)))
      j = Ngrid - 1;
    long long int k = (long long int)(Tracer[ii].z / boxsize.z * Ngrid);
    if (unlikely(k > (Ngrid - 1)))
      k = Ngrid - 1;
    long long int ind = i + Ngrid * (j + Ngrid * k);
    Tracer[ii].gbin = ind;
    // Tracer[ii].id=ii;
  }

  if (Tracer.size() < 16777216) {
    cout << "sorting particles by grid..." << endl;
    cout.flush();
    sort(Tracer.begin(), Tracer.end(), small_grid_index_first());
    cout << "done." << endl;
  }

  NMEAN = ((double)Tracer.size() / (boxsize.x * boxsize.y * boxsize.z));

  // Construyo grid para identificar voids y calculo de velocidades
  G = new Grid(Ngrid, Tracer.size(), boxsize, radmax);

  cout << "placing particles into the grid... " << endl;
  cout.flush();
  (*G).build(&Tracer);
  cout << "done." << endl;
  cout.flush();

  rhobar = (float)Tracer.size() / (boxsize.x * boxsize.y * boxsize.z);

  /// trabajo
  svf_rma(rank, limit, rhobar, nshell, G, min_radius);

  // FindCenters_fast(limit, rhobar, nshell, G, min_radius);
  // FindSphvds_fastV2(limit, rawfname, nshell, G, min_radius);

  return (0);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int numprocs = 1;
  int proc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  cout << "MPI: proc " << proc << " of " << numprocs;

  // initialize instance
  VariablesManager *vars = VariablesManager::getInstance();
  // load variables
  loadConfigurationVars(argc, argv);

  NCORES = omp_get_max_threads();
  fprintf(stdout, "\n ====>>>> Identificador corriendo en %d core(s) \n", (int)NCORES);
  fflush(stdout);

  string fmt = vars->valueFromKey(K_FILEFMT);
  string fname_tracer = vars->valueFromKey(K_TRSFILE);
  string sph_file = vars->valueFromKey(K_SPHFILE);
  int num_files = Converter::stringToInt(vars->valueFromKey(K_NUM_FILE));
  radmax = Converter::stringToFloat(vars->valueFromKey(K_MAXRADIUS));
  float densth = Converter::stringToFloat(vars->valueFromKey(K_DENSTH));
  // underdensity parameter (typical values -0.9 o -0.8)
  float min_radius = Converter::stringToFloat(vars->valueFromKey(K_MINRADIUS));

  if (!fmt.compare("ASCII") || !fmt.compare("Minerva")) {
    float boxlength = Converter::stringToFloat(vars->valueFromKey(K_BOXSIZE));
    masslim = Converter::stringToFloat(vars->valueFromKey(K_MASSMIN));
    cout << "minimum mass= " << masslim << endl;
    boxsize.x = boxlength;
    boxsize.y = boxlength;
    boxsize.z = boxlength;
  }
  if (!fmt.compare("HDF5_SUBFIND_GROUPS")) {
    masslim = Converter::stringToFloat(vars->valueFromKey(K_MASSMIN));
    cout << "minimum mass= " << masslim << endl;
  }

  if (!fmt.compare("HDF5_SUBFIND_SUBHALOS")) {
    masslim = Converter::stringToFloat(vars->valueFromKey(K_MASSMIN));
    cout << "minimum mass= " << masslim << endl;
  }
  cout << "Reading Header " << endl;
  ReadHeaderTracers(fmt, fname_tracer, num_files);
  cout << "Reading Tracers" << endl;
  ReadTracers(fmt, fname_tracer, num_files);
  cout << "Running spherical void finder" << endl;
  cout.flush();

  svf_mpi_start(densth, min_radius, radmax, proc);
  // gather results on rank 0
  auto all_voids = gather_sphvds(proc, numprocs, Sphvd);
  Sphvd = all_voids;

  if (proc == 0) {
    Grid *Gv;
    long long int Ngrid = find_ngrid((long long int)Tracer.size(), min_radius, boxsize.x);
    Gv = new Grid(Ngrid, Sphvd.size(), boxsize, 2 * radmax);
    (*Gv).build(&Sphvd);
    CleanSphvds_Tol0(Gv, radmax);
    Write_sv(sph_file);
  }
  MPI_Finalize();

  return 0;
}
