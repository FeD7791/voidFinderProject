/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

// Collection of several small tools

#include "variables_tools.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include "string_tools.h"
#include "variables_loader.h"
#include "variables_manager.h"

// #define VARPATH "configuration/"
#define KEYVARFILE "config"
#define KEYNOCONFIGFILE "no"
#define VARFILE "configuration/vars.conf"

using namespace std;
using namespace iate_variables;

void loadConfigurationVars(int argc, char *argv[]) {
  // string configpath(VARPATH);
  string configfile;
  stringstream config;

  VariablesManager *vars = VariablesManager::getInstance();
  vars->addMap(VariablesLoader::loadFromArgv(argc, argv));
  /*the value is requested in unsafe mode.
    It does  not throw exception in case of not found the key*/
  configfile = vars->valueFromKey(KEYVARFILE, false);

  if (configfile.empty()) {
    // config << configpath << string(VARFILE);
    config << string(VARFILE);
    // load filename from argv
    vars->addMap(VariablesLoader::loadFromFile(config.str()));
  } else {
    configfile = clean_left(clean_right(configfile));
    if (configfile.compare(KEYNOCONFIGFILE) != 0) {
      // config << configpath << configfile;
      config << configfile;
      vars->addMap(VariablesLoader::loadFromFile(config.str()));
    } // argv -->config=no : do not load var
  }
}
