#include <omp.h>
#include <unistd.h>

#include <cstdlib> /* system*/
#include <iostream>
#include <string>

#include "allvars.h"
#include "clean_duplicates.h"
#include "colors.h"
#include "compute_intersecs.h"
#include "constants.h" //  K_OUTFILE_INTERSECS K_THRSLD K_CRTP K_EPS
#include "grid.h"
#include "io_lib/io.h"
#include "popcorn.h"
#include "svf.h"
// POPFILE K_POPFILE_COMPUTE_INTERSECS
#include "configuration_lib/converter.h"
#include "configuration_lib/variable_control.h"
#include "configuration_lib/variables_manager.h"
#include "configuration_lib/variables_tools.h"

using namespace std;
using namespace grid;
using namespace iate_variables;
void cartel_inicio();

int main(int argc, char *argv[]) {

  cartel_inicio();
  NCORES = omp_get_max_threads(); // omp_get_num_threads();

  fprintf(stdout, "\n ====>>>> Identificador corriendo en %d core(s) \n", (int)NCORES);
  // initialize instance
  VariablesManager *vars = VariablesManager::getInstance();
  // load variables
  loadConfigurationVars(argc, argv);

  conf_var_file confile = VariableControl::loader_file_var(vars);
  conf_var confglobal = VariableControl::loader_global_var(vars, confile.fmt);

  radmax = confglobal.radmax;
  boxsize.x = confglobal.boxsize.x;
  boxsize.y = confglobal.boxsize.y;
  boxsize.z = confglobal.boxsize.z;

  VariableControl::show_file_var(confile);
  VariableControl::show_conf_var(confglobal);

  ReadTracers(confile.fmt, confile.tracer_file, confile.num_files);

  if (access(confile.sph_file.c_str(), F_OK) != -1) {
    // read into Sphvd
    Read_sv(confile.sph_file, confglobal.boxsize);
  } else {
    svf_start(confile.densth, confile.sph_rawfile, confile.min_radius);

    {
      vector<Sphere> copy;
      for (long unsigned int i = 0; i < Sphvd.size(); i++) {
        if (Sphvd[i].ToF) {
          copy.push_back(Sphvd[i]);
        }
      }
      Sphvd.clear();
      for (long unsigned int i = 0; i < copy.size(); i++)
        Sphvd.push_back(copy[i]);
      copy.clear();

      NVOID = Sphvd.size();
    }
    // Escritura void sphericos
    Write_sv(confile.sph_file);
  }

  popcorn_start(confile.min_radius, confile.eps, confile.densth, 0, 1);
  if (confile.aux_files)
    Write_raw_popcorn(confile.rawpop_file);

  compute_intersecs_start();
  if (confile.aux_files)
    Write_pairs(confile.pairs_file);

  vector<bool> erase;
  clean_duplicates_start(confile.overlap, &erase);
  Write_clean_popcorn(confile.pop_file, &erase);

  Sphvd.clear();
  return 0;
}

void cartel_inicio() {
  cout << "\n\n";
  cout << RED_FONT;
  cout << "   ██████╗  ██████╗ ██████╗  ██████╗ ██████╗ ██████╗ ███╗   ██╗" << endl;
  cout << "   ██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔═══██╗██╔══██╗████╗  ██║" << endl;
  cout << "   ██████╔╝██║   ██║██████╔╝██║     ██║   ██║██████╔╝██╔██╗ ██║" << endl;
  cout << "   ██╔═══╝ ██║   ██║██╔═══╝ ██║     ██║   ██║██╔══██╗██║╚██╗██║" << endl;
  cout << "   ██║     ╚██████╔╝██║     ╚██████╗╚██████╔╝██║  ██║██║ ╚████║" << endl;
  cout << "   ╚═╝      ╚═════╝ ╚═╝      ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝" << endl;
  cout << "\n";
  cout << "   ███████████████████████████████████████████████████████████████╗" << endl;
  cout << "   ╚══════════════════════════════════════════════════════════════╝" << endl;
  cout << RESTORE << WHITE_FONT << "\n";
}
