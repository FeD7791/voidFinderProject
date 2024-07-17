/*
 * Copyright 2019 Popcorn Authors. All Rights Reserved.
 *
 * License that can be found in the LICENSE file.
 */

// Two functions for remove characters ( \n\r\t\f\v) to
// right and left of the allowed value.

#ifndef STRING_TOOLS_C
#define STRING_TOOLS_C

#include <string>

std::string clean_left(const std::string &text);

std::string clean_right(const std::string &text);

#endif // STRING_TOOLS_C
