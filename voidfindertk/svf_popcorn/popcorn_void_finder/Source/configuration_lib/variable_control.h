/*
 * Copyright 2021 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */
#ifndef IATE_VARIABLES_VARIABLESLOADER_H_
#define IATE_VARIABLES_VARIABLESLOADER_H_

#include <string>

#include "../allvars.h"
#include "variables_manager.h"

namespace iate_variables {

class VariableControl {

public:
  /* Read
   *
   *  param vars =
   *
   *  returns
   */

  static conf_var_file loader_file_var(VariablesManager *vars);
  static conf_var loader_global_var(VariablesManager *vars, std::string fmt);
  static void show_file_var(conf_var_file conf);
  static void show_conf_var(conf_var conf);

}; // class

} // namespace iate_variables

#endif // IATE_VARIABLES_VARIABLESLOADER_H_