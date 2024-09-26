/*
 * Copyright 2021 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

#include "variable_control.h"

#include <iostream>
#include <string>

#include "../allvars.h"
#include "../constants.h"
#include "../mensaje.h"
#include "converter.h"
#include "variables_manager.h"

using std::cout;
using std::string;

namespace iate_variables {

/* static */
conf_var_file VariableControl::loader_file_var(VariablesManager *vars) {
  conf_var_file conf;

  conf.pop_file = vars->valueFromKey(K_POPFILE, false);
  conf.rawpop_file = vars->valueFromKey(K_RAWPOPFILE, false);
  conf.sph_file = vars->valueFromKey(K_SPHFILE, false);
  conf.tracer_file = vars->valueFromKey(K_TRSFILE, false);
  conf.pairs_file = vars->valueFromKey(K_PAIRSFILE, false);
  conf.fmt = vars->valueFromKey(K_FILEFMT, false);
  conf.sph_rawfile = vars->valueFromKey(K_RAWSPHFILE, false);

  string masslim_str = vars->valueFromKey(K_MASSMIN, false);
  if (masslim_str == "") {
    conf.masslim = 0.0f;
  } else {
    conf.masslim = Converter::stringToFloat(masslim_str);
  }

  string num_file_str = vars->valueFromKey(K_NUM_FILE, false);

  if (num_file_str == "") {
    conf.num_files = 0;
  } else {
    conf.num_files = Converter::stringToInt(num_file_str);
  }

  string aux_files_str = vars->valueFromKey(K_AUXFILES, false);

  if (aux_files_str == "") {
    conf.aux_files = false;
  } else {
    conf.aux_files = Converter::stringToBool(aux_files_str);
  }

  string densth_str = vars->valueFromKey(K_DENSTH, false);
  if (densth_str == "") {
    conf.densth = 0.0f;
  } else {
    conf.densth = Converter::stringToFloat(densth_str);
    // underdensity parameter (typical values -0.9 o -0.8)
  }

  string eps_str = vars->valueFromKey(K_EPS, false);

  if (eps_str == "") {
    conf.eps = 0.0f;
  } else {
    conf.eps = Converter::stringToDouble(eps_str);
  }

  string min_radius_str = vars->valueFromKey(K_MINRADIUS, false);

  if (min_radius_str == "") {
    conf.min_radius = 0.0f;
  } else {
    conf.min_radius = Converter::stringToFloat(min_radius_str);
  }

  string max_radius_str = vars->valueFromKey(K_MAXRADIUS, false);

  if (max_radius_str == "") {
    conf.max_radius = 0.0f;
  } else {
    conf.max_radius = Converter::stringToFloat(max_radius_str);
  }

  return conf;
}

/* static */
conf_var VariableControl::loader_global_var(VariablesManager *vars, string fmt) {
  conf_var conf;
  string boxsize_str = vars->valueFromKey(K_BOXSIZE, false);

  float boxsize_float = Converter::stringToFloat(boxsize_str);
  if (!fmt.compare("ASCII")) {
    conf.boxsize.x = boxsize_float;
    conf.boxsize.y = boxsize_float;
    conf.boxsize.z = boxsize_float;
  } else {
    if (boxsize_float != 0)
      cout << "\n" << MSJ_WARNING_ASCII << boxsize_str << "\n";
  }

  string radmax_str = vars->valueFromKey(K_RADMAX, false);
  if (radmax_str == "") {
    conf.radmax = 0.0f;
  } else {
    conf.radmax = Converter::stringToFloat(radmax_str);
  }

  return conf;
}

void VariableControl::show_file_var(conf_var_file conf) {
  cout << MSJ_VAR_FILEVAR << "\n";
  cout << MSJ_VAR_POPFILE << conf.pop_file << "\n";
  cout << MSJ_VAR_RAWPOPFILE << conf.rawpop_file << "\n";
  cout << MSJ_VAR_SPHFILE << conf.sph_file << "\n";
  cout << MSJ_VAR_TRACERFILE << conf.tracer_file << "\n";
  cout << MSJ_VAR_PAIRFILE << conf.pairs_file << "\n";
  cout << MSJ_VAR_FMT << conf.fmt << "\n";
  cout << MSJ_VAR_SPHRAWFILE << conf.sph_rawfile << "\n";
  cout << MSJ_VAR_NUMFILE << conf.num_files << "\n";
  cout << MSJ_VAR_AUXFILE << conf.aux_files << "\n";
  cout << MSJ_VAR_DENSTH << conf.densth << "\n";
  cout << MSJ_VAR_EPS << conf.eps << "\n";
  cout << MSJ_VAR_MINRADIUS << conf.min_radius << "\n";
  cout << MSJ_VAR_MAXRADIUS << conf.max_radius << "\n";
}

void VariableControl::show_conf_var(conf_var conf) {
  cout << MSJ_VAR_CONFVAR << "\n";
  cout << MSJ_VAR_BOXSIZE << conf.boxsize.x << "\n";
  cout << MSJ_VAR_BOXSIZE << conf.boxsize.y << "\n";
  cout << MSJ_VAR_BOXSIZE << conf.boxsize.z << "\n";
  cout << MSJ_VAR_RADMAX << conf.radmax << "\n";
}

} // namespace iate_variables
