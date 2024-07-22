#include "arvo_variables_tools.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>

#include "arvo_functions.h"
#include "constants.h"

using namespace std;
namespace arvo {

void SetInd(const int i, const int val, int **ind_temp, int *sizeInd_temp) {
  static string message = "Error (re)allocating neighborsIndices. Exiting.\n";
  reallocVarPointer(i, 0, 1, val, ind_temp, sizeInd_temp, message);
}

void SetNI(const int i, const int val, int **neighborsIndices_temp, int *sizeNeighborsIndices_temp) {
  static string message = "Error (re)allocating neighborsIndices. Exiting.\n";
  reallocVarPointer(i, 0, 1, val, neighborsIndices_temp, sizeNeighborsIndices_temp, message);
}
} // namespace arvo
