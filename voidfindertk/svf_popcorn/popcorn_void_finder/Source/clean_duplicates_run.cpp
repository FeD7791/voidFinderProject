#include <string>

#include "allvars.h"
#include "clean_duplicates.h"
#include "constants.h" // K_OUTFILE_INTERSECS K_THRSLD K_CRTP K_EPS
#include "io_lib/io.h"

// POPFILE K_POPFILE_COMPUTE_INTERSECS
#include "configuration_lib/converter.h"
#include "configuration_lib/variables_manager.h"
#include "configuration_lib/variables_tools.h"

using namespace std;
using namespace iate_variables;

int main(int argc, char *argv[]) {
  // initialize instance
  double pirad = 3.141592653589793238462643383279;
  VariablesManager *vars = VariablesManager::getInstance();
  // load variables
  loadConfigurationVars(argc, argv);
  string rawpop_file = vars->valueFromKey(K_RAWPOPFILE);
  string pairs_file = vars->valueFromKey(K_PAIRSFILE);
  string pop_file = vars->valueFromKey(K_POPFILE);
  string sph_file = vars->valueFromKey(K_SPHFILE);

  string fname_tracer = vars->valueFromKey(K_TRSFILE);
  string fmt = vars->valueFromKey(K_FILEFMT);
  float overlap = Converter::stringToFloat(vars->valueFromKey(K_MINRADIUS));
  overlap = 4.0 / 3.0 * pirad * overlap * overlap * overlap;

  int num_files = Converter::stringToInt(vars->valueFromKey(K_NUM_FILE));

  if (!fmt.compare("ASCII")) {
    float boxlength = Converter::stringToFloat(vars->valueFromKey(K_BOXSIZE));
    boxsize.x = boxlength;
    boxsize.y = boxlength;
    boxsize.z = boxlength;
  }

  ReadHeaderTracers(fmt, fname_tracer, num_files);
  Read_sv(sph_file, boxsize);
  Read_raw_popcorn(rawpop_file);
  Read_pairs(pairs_file);
  vector<bool> erase;
  clean_duplicates_start(overlap, &erase);
  Write_clean_popcorn(pop_file, &erase);
  return 0;
}
