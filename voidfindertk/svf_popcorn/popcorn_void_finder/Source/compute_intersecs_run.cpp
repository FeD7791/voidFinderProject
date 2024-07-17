#include "allvars.h"
#include "compute_intersecs.h"
#include "constants.h" // K_OUTFILE_INTERSECS K_THRSLD K_CRTP K_EPS
#include "io_lib/io.h"
#include <string>
// POPFILE K_POPFILE_COMPUTE_INTERSECS
#include "configuration_lib/converter.h"
#include "configuration_lib/variables_manager.h"
#include "configuration_lib/variables_tools.h"
using namespace std;
using namespace iate_variables;
int main(int argc, char *argv[]) {
  // initialize instance
  VariablesManager *vars = VariablesManager::getInstance();
  // load variables
  loadConfigurationVars(argc, argv);
  string rawpop_file = vars->valueFromKey(K_RAWPOPFILE);
  string pairs_file = vars->valueFromKey(K_PAIRSFILE);

  string fname_tracer = vars->valueFromKey(K_TRSFILE);
  string fmt = vars->valueFromKey(K_FILEFMT);
  int num_files = Converter::stringToInt(vars->valueFromKey(K_NUM_FILE));

  if (!fmt.compare("ASCII")) {
    float boxlength = Converter::stringToFloat(vars->valueFromKey(K_BOXSIZE));
    boxsize.x = boxlength;
    boxsize.y = boxlength;
    boxsize.z = boxlength;
  }

  ReadHeaderTracers(fmt, fname_tracer, num_files);

  Read_raw_popcorn(rawpop_file);
  compute_intersecs_start();
  Write_pairs(pairs_file);
  return 0;
}
