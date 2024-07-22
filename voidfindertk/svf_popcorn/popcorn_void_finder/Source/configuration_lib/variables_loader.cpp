/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

#include "variables_loader.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "string_tools.h"

using namespace std;

namespace iate_variables {

#define DEFAULT_SEPARATOR '='

/* static */
map<string, string> VariablesLoader::loadFromFile(string nameFile, char separator, bool loadDataRaw) {
  string line;
  string lineAux;
  string key;
  string value;
  map<string, string> mapVar;

  try {
    ifstream file_in(nameFile.c_str());
    while (getline(file_in, line)) {
      lineAux = clean_left(clean_right(line));
      if (line.empty()) {
        lineAux = "#";
      }
      if (lineAux.at(0) != '#') {
        stringstream ss(line);
        getline(getline(ss, key, separator) >> std::ws, value);
        if (!loadDataRaw) {
          key = clean_left(clean_right(key));
          value = clean_left(clean_right(value));
        }
        mapVar[key] = value;
      }
    }
  } catch (const ifstream::failure &e) {
    cout << "Exception opening/reading " << nameFile << " file. \n";
  }
  return mapVar;
}

/*static */
map<string, string> VariablesLoader::loadFromFile(string nameFile) {
  return VariablesLoader::loadFromFile(nameFile, DEFAULT_SEPARATOR, false);
}

/*static*/
map<string, string> VariablesLoader::loadFromArgv(int argc, char *argv[]) {
  string key;
  string value;
  map<string, string> mapVar;

  for (int i = 1; i < argc; ++i) {
    istringstream param(argv[i]);
    if (getline(getline(param, key, DEFAULT_SEPARATOR) >> std::ws, value)) {
      key = clean_left(clean_right(key));
      value = clean_left(clean_right(value));
      mapVar[key] = value;
    }
  }
  return mapVar;
}

} // namespace iate_variables
