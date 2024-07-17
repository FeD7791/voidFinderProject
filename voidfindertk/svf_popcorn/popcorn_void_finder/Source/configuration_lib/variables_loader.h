/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

// Load the variables from a file and return a map<key,value>
// The format of the reading file should be:
// -key"separator character"value
// -Each pair must be in a new line
//
//  Example:
//   key1=value1
//   key2=value2
//   ..........
//
// Usage:
//
// map<string,string> varMap = VariablesLoader::loadFromFile("file_name.conf")
//

#ifndef IATE_VARIABLES_VARIABLESLOADER_H_
#define IATE_VARIABLES_VARIABLESLOADER_H_

#include <map>
#include <string>

namespace iate_variables {

class VariablesLoader {

public:
  /* Reads variables from file with key-separator-value format
   *
   *  param nameFile = full path of the name file to load
   *  param separator = separator character between key and value
   *  param loadDataRaw = flag to enable the cleaning of the string
   *
   *  returns one map with format map(string = key, string = value)
   */
  static std::map<std::string, std::string> loadFromFile(std::string nameFile, char separator, bool loadDataRaw);
  /* Separator used '=' */
  static std::map<std::string, std::string> loadFromFile(std::string nameFile);

  /* Load variables from command line arguments
   *  param argv is an array of c-string pointers
   *  param number of c-string pointers pointed to by argv
   *
   *  returns one map with format map(string=key)
   *
   *  Format of  c-string pointers "key=value"
   *  another format is ignored
   */
  static std::map<std::string, std::string> loadFromArgv(int argc, char *argv[]);

}; // class VariablesLoader

} // namespace iate_variables

#endif // IATE_VARIABLES_VARIABLESLOADER_H_
