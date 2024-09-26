#include "converter.h"
#include "string_tools.h"
#include "variables_loader.h"
#include "variables_manager.h"
#include <iostream>
#include <math.h> /* fabs */
#include <string>

using namespace std;
using namespace iate_variables;

void pruebaConfigVar();

// TEST STRING_TOOLS.H

bool exit_success__string_tools__clean_left(string strRaw, string strFormat) {
  bool result = false;
  string stringClean = clean_left(string(" /t/n/r/f/v value "));
  if (stringClean.compare("value ")) {
    result = true;
  }
  return result;
}

bool exit_success__string_tools__clean_right(string strRaw, string strFormat) {
  bool result = false;
  string stringClean = clean_right(string(" value /n/t/r/f/v"));
  if (stringClean.compare(" value")) {
    result = true;
  }
  return result;
}

// TEST VARIABLES_LOADER

bool map_comp(map<string, string> map1, map<string, string> map2) {
  return map1.size() == map2.size() && std::equal(map1.begin(), map1.end(), map2.begin());
}

bool exit_success__variables_loader(string nameFile, map<string, string> correctMap) {

  bool result = false;
  try {
    map<string, string> mapVars = VariablesLoader::loadFromFile(nameFile);

    result = map_comp(mapVars, correctMap);
  } catch (exception &e) {
    cout << "A exception was caught, with message '" << e.what() << "'\n";
    throw;
  }

  return result;
}

// TEST VARIABLE_MANAGER

bool exit_success__variables_manager__load(string nameFile, map<string, string> correctMap) {

  bool result = false;
  try {
    map<string, string> mapVars = VariablesLoader::loadFromFile(nameFile);
    VariablesManager *vars = VariablesManager::getInstance();
    vars->addMap(VariablesLoader::loadFromFile("vars.conf"));

    VariablesManager *vars2 = VariablesManager::getInstance();
    map<string, string> mapAux = vars2->getMap();

    result = map_comp(mapAux, correctMap);
  } catch (exception &e) {
    cout << "A exception was caught, with message '" << e.what() << "'\n";
    throw;
  }

  return result;
}

// TEST CONVERTER
bool exit_success__coverter__stringToInt(string stringValue, int intValue) {
  bool result = false;
  int resultValue = 0;
  try {
    resultValue = Converter::stringToInt(stringValue);
    if (resultValue == intValue) {
      result = true;
    }
  } catch (exception &e) {
    cout << "A exception was caught, with message '" << e.what() << "'\n";
    throw;
  }

  return result;
}

bool exit_success__coverter__stringToFloat(string stringValue, float floatValue) {
  bool result = false;
  float resultValue = 0;
  try {
    resultValue = Converter::stringToFloat(stringValue);
    if (resultValue == floatValue) {
      result = true;
    }
  } catch (exception &e) {
    cout << "A exception was caught, with message '" << e.what() << "'\n";
    throw;
  }

  return result;
}

bool exit_success__coverter__stringToDouble(string stringValue, double long doubleValue) {
  bool result = false;
  double long resultValue = 0;
  try {
    resultValue = Converter::stringToDouble(stringValue);
    // conversion loses presicion
    if (fabs(resultValue - doubleValue) < 0.01) {
      result = true;
    }
  } catch (exception &e) {
    cout << "A exception was caught, with message '" << e.what() << "'\n";
    throw;
  }
  return result;
}

bool exit_success__coverter__stringToBool(string stringValue, bool boolValue) {
  bool result = false;
  bool resultValue = 0;
  try {
    resultValue = Converter::stringToBool(stringValue);
    if (resultValue == boolValue) {
      result = true;
    }
  } catch (exception &e) {
    cout << "A exception was caught, with message '" << e.what() << "'\n";
    throw std::invalid_argument("");
  }
  return result;
}

// TEST
void showResultTest(int countFunction, string nameFunction, bool result) {
  if (result) {
    cout << "  " << countFunction << " " << nameFunction << " \033[1;32m OK\033[0m\n";
  } else {
    cout << " " << countFunction << " " << nameFunction << " \033[1;31m FAILURE\033[0m\n";
  }
}

int main() {

  /*
  VariablesManager *vars = VariablesManager::getInstance();
  vars->addMap( VariablesLoader::loadFromFile("vars.conf"));

  VariablesManager *vars2 = VariablesManager::getInstance();
  vars2->showVariables();

  std::cout << "\n \n";
  delete vars2;*/

  cout << "\n\nFunction:\n\n";

  bool result = false;
  int count = 0;

  // string_tools
  cout << "  string_tools:\n";
  result = exit_success__string_tools__clean_left(string(" /t/n/r/f/v value "), string("value "));
  count = count + 1;
  showResultTest(count, string("exit_success__string_tools__clean_left"), result);

  result = exit_success__string_tools__clean_right(string("  value /t/n/r/f/v "), string(" value "));
  count = count + 1;
  showResultTest(count, string("exit_success__string_tools__clean_right"), result);

  //  variables_loader
  cout << "\n  variables_loader:\n";
  map<string, string> dic;
  dic["key1"] = "value1";
  dic["key2"] = "value2";
  dic["key3"] = "value3";
  dic["key4"] = "value4";
  dic["key5"] = "value5";
  dic["key6"] = "value6";
  dic["keyEnd"] = "valueEnd";

  result = exit_success__variables_loader("vars.conf", dic);
  count = count + 1;
  showResultTest(count, string("exit_success__variables_loader"), result);

  dic["key3"] = "valorErroneo";
  result = exit_success__variables_loader("vars.conf", dic);
  count = count + 1;
  showResultTest(count, string("exit_failure__variables_loader"), result);

  // variables_manager
  cout << "\n  variables_manager:\n";
  dic["key3"] = "value3";
  result = exit_success__variables_manager__load("vars.conf", dic);
  count = count + 1;
  showResultTest(count, string("exit_success__variables_manager__load"), result);

  // converter
  cout << "\n converter: \n";
  // int
  result = exit_success__coverter__stringToInt(string("2"), 2);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToInt"), result);

  result = exit_success__coverter__stringToInt(string("2"), 4);
  count = count + 1;
  showResultTest(count, string("exit_failure__coverter__stringToInt"), result);

  result = exit_success__coverter__stringToInt(string("  2332  "), 2332);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToInt"), result);

  result = exit_success__coverter__stringToInt(string("5  2332 5 "), 2332);
  count = count + 1;
  showResultTest(count, string("exit_failure__coverter__stringToInt"), result);

  // float
  result = exit_success__coverter__stringToFloat(string("2.654"), 2.654);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToFloat"), result);

  result = exit_success__coverter__stringToFloat(string(" 5.66781   "), 5.66782);
  count = count + 1;
  showResultTest(count, string("exit_failure__coverter__stringToFloat"), result);

  result = exit_success__coverter__stringToFloat(string(" 56.76666788  "), 56.76666788);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToFloat"), result);

  result = exit_success__coverter__stringToFloat(string("5 45.65   5 "), 45.65);
  count = count + 1;
  showResultTest(count, string("exit_failure__coverter__stringToFloat"), result);

  // double
  double long double_aux = 324.656;
  result = exit_success__coverter__stringToDouble(string("324.654"), double_aux);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToDouble"), result);

  result = exit_success__coverter__stringToDouble(string(" 100005.6682    "), 100006.6682);
  count = count + 1;
  showResultTest(count, string("exit_failure__coverter__stringToDouble"), result);

  result = exit_success__coverter__stringToDouble(string(" 56.77811919  "), 56.77811919);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToDouble"), result);

  result = exit_success__coverter__stringToDouble(string(" 45.65   5 "), 45.665);
  count = count + 1;
  showResultTest(count, string("exit_failure__coverter__stringToDouble"), result);

  // bool
  result = exit_success__coverter__stringToBool(string("true"), true);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToBool"), result);

  result = exit_success__coverter__stringToBool(string("True"), true);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToBool"), result);

  result = exit_success__coverter__stringToBool(string("false"), false);
  count = count + 1;
  showResultTest(count, string("exit_success__coverter__stringToBool"), result);

  try {
    count = count + 1;
    result = exit_success__coverter__stringToBool(string("FALSE"), false);
    showResultTest(count, string("exit_success__coverter__stringToBool"), result);
  } catch (exception &e) {
    cout << " \033[1;33m ";
    cout << "\nexception in:  exit_failure__coverter__stringToBool "
         << "  test " << count << "param '" << string(" FALSE ") << "' param "
         << "false";
    cout << "\033[0m\n";
  }

  count = count + 1;
  try {
    result = exit_success__coverter__stringToBool(string(" FALSE "), true);
    showResultTest(count, string("exit_failure__coverter__stringToBool"), result);
  } catch (exception &e) {
    cout << " \033[1;33m ";
    cout << "\nexception in:  exit_failure__coverter__stringToBool "
         << "  test " << count << "param '" << string(" FALSE ") << "' param "
         << "true";
    cout << "\033[0m\n";
  }

  count = count + 1;
  try {
    result = exit_success__coverter__stringToBool(string(" falsE "), false);
    showResultTest(count, string("exit_success__coverter__stringToBool"), result);
  } catch (exception &e) {
    cout << " \033[1;33m ";
    cout << "\nexception in:  exit_failure__coverter__stringToBool "
         << "  test " << count << "param '" << string(" falsE ") << "' param "
         << "false";
    cout << "\033[0m\n";
  }

  count = count + 1;
  try {
    result = exit_success__coverter__stringToBool(string("5 false  5 "), false);
    showResultTest(count, string("exit_failure__coverter__stringToBool"), result);
  } catch (exception &e) {
    cout << " \033[1;33m ";
    cout << "  " << count << " exception: exit_failure__coverter__stringToBool"
         << "  test "
         << "param '" << string("5 false  5 ") << "' param "
         << "false";
    cout << "\033[0m\n";
  }

  // end
  cout << "\n\n";
}
