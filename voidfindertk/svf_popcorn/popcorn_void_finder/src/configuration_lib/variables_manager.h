/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

// A singleton with the variables.
//
// Usage:
//
// // initialize instance
// VariablesManager *vars = VariablesManager::getInstance();
//
// // add variable
// vars->addVariable("key","value");
//
// // get variable
// std::cout << vars2->valueFromKey("key")
//
// // load map from a configuration file
// vars->addMap( VariablesLoader::loadFromFile("popcorn.conf"));
//

#ifndef IATE_VARIABLES_VARIABLESMANAGER_H_
#define IATE_VARIABLES_VARIABLESMANAGER_H_

#include <map>
#include <string>

namespace iate_variables {

class VariablesManager {

public:
  /* method to create instance of singleton class */
  static VariablesManager *getInstance(void);

  /* get a value from the key in safe mode
     Safe mode: throw exception by key not found or if the value is empty*/
  std::string valueFromKey(std::string key, bool safeMode);

  /* get a value from the keyin safe mode.*/
  std::string valueFromKey(std::string key);

  /* add variable to the map */
  void addVariable(std::string key, std::string value);

  /* modify the value of a variable */
  void updateVariable(std::string key, std::string value);

  /* Removes all variables */
  void clear();

  /* Add variables to the manager from map */
  void addMap(std::map<std::string, std::string> externMap);

  /* Returns the number of variables */
  size_t size();

  /* Return a copy of the map (all )*/
  std::map<std::string, std::string> getMap();

  /* Shows the value of the variables by the standard output */
  void showVariables();

private:
  /* Here will be the instance stored */
  static VariablesManager *instance;

  /* Variables container */
  std::map<std::string, std::string> _mapVars;

  /* Constructor protected */
  VariablesManager() {}
};

} // namespace iate_variables

#endif // IATE_VARIABLES_VARIABLESMANAGER_H_
