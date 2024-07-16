/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

// A singleton with the variables.
//
// Usage:
//

#ifndef IATE_VARIABLES_CONVERTER_H_
#define IATE_VARIABLES_CONVERTER_H_

#include <iostream>
#include <map>
#include <string>

namespace iate_variables {

class Converter {
public:
  /* Parses the string 'value' interpreting its contents as an integer number,
   * which is returned as a valueof type int*/
  static int stringToInt(std::string value);

  /* Parses the string 'value' interpreting its contents as an single precision
   * floating point number, which is returned as a valueof type float*/
  static float stringToFloat(std::string value);

  /* Parses the string 'value' interpreting its contents as an double precision
   * floating point number, which is returned as a valueof type double*/
  static double stringToDouble(std::string value);

  /* Parses the string 'value' interpreting its contents as an boolean value
   * , which is returned as a valueof type bool*/
  static bool stringToBool(std::string value);
};

} // namespace iate_variables

#endif // IATE_VARIABLES_CONVERTER_H_
