/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

#include "converter.h"

#include <algorithm> // std::transform
#include <iostream>
#include <sstream>
#include <string>

#include "string_tools.h"

using std::string;

namespace iate_variables {

/* static */
int Converter::stringToInt(string value) {
  float result;
  value = clean_left(clean_right(value));
  std::istringstream iss_value(value);
  iss_value >> result;
  return result;
}

/* static */
float Converter::stringToFloat(string value) {
  float result;
  value = clean_left(clean_right(value));
  std::istringstream iss_value(value);
  iss_value >> result;
  return result;
}

/* static */
double Converter::stringToDouble(string value) {
  double result;
  value = clean_left(clean_right(value));
  std::istringstream iss_value(value);
  iss_value >> result;
  return result;
}

/* static */
bool Converter::stringToBool(string value) {
  bool result;
  value = clean_left(clean_right(value));

  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  std::istringstream aux(value);

  if (aux.str() == string("false")) {
    result = false;
  } else if (aux.str() == string("true")) {
    result = true;
  } else {
    string site("Convert::stringToBool ");
    string message("value to convert boolean is not correct: ");
    throw std::invalid_argument(site + message + value);
  }
  return result;
}

} // namespace iate_variables
