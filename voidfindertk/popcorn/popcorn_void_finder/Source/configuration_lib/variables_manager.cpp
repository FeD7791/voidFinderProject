/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

#include "variables_manager.h"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

using namespace std;

namespace iate_variables {

/* static */
VariablesManager *VariablesManager::getInstance(void) {
  if (instance == 0) {
    instance = new VariablesManager();
  }
  return instance;
}

string VariablesManager::valueFromKey(string key) { return this->valueFromKey(key, true); }

string VariablesManager::valueFromKey(string key, bool safeMode) {
  string value;
  map<string, string>::iterator it = _mapVars.find(key);

  if (it == _mapVars.end()) {
    if (safeMode) {
      throw invalid_argument("VariableManager invalid key: " + key);
    } else {
      value = string("");
    }
  } else {
    value = string(it->second);
    if (value.size() < 1) {
      if (safeMode) {
        throw invalid_argument("VariableManager invalid key: " + key);
      } else {
        value = string("");
      }
    }
  }
  return value;
}

void VariablesManager::addVariable(string key, string value) {
  _mapVars.insert(_mapVars.end(), pair<string, string>(key, value));
}

void VariablesManager::updateVariable(string key, string value) { _mapVars[key] = value; }

void VariablesManager::clear() { _mapVars.clear(); }

void VariablesManager::addMap(map<string, string> externMap) { _mapVars.insert(externMap.begin(), externMap.end()); }

size_t VariablesManager::size() { return _mapVars.size(); }

map<string, string> VariablesManager::getMap() {
  map<string, string> mapAux;
  mapAux = _mapVars;
  return mapAux;
}

void VariablesManager::showVariables() {
  map<string, string>::iterator iterator;
  cout << "Vars: \n";
  for (iterator = _mapVars.begin(); iterator != _mapVars.end(); iterator++) {
    cout << "key: " << iterator->first << " value: " << iterator->second << " \n";
  }
}

/* instance will be initialized on demand.*/
VariablesManager *VariablesManager::instance = 0;

} // namespace iate_variables
