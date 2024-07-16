/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */
#include "string_tools.h"

#include <algorithm>
#include <iostream>
#include <string>
using std::string;

const string K_WHITE_SPACE = " \n\r\t\f\v";

string clean_left(const string &text) {
  string result = "";
  if (!text.empty()) {
    size_t startString = text.find_first_not_of(K_WHITE_SPACE);
    if (startString != std::string::npos) {
      result = text.substr(startString);
    }
  }
  return result;
}

string clean_right(const string &text) {
  string result = "";
  if (!text.empty()) {
    size_t endString = text.find_last_not_of(K_WHITE_SPACE);
    if (endString != std::string::npos) {
      result = text.substr(0, endString + 1);
    }
  }
  return result;
}
