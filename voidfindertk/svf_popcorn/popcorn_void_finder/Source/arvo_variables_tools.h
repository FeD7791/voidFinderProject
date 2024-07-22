#ifndef ARVO_VARIABLES_TOOL_H_
#define ARVO_VARIABLES_TOOL_H_

#include <string>

namespace arvo {

#define ARRAY_INCREMENT 10 // by this size are arrays resized when necessary

void SetInd(const int i, const int val, int **ind_temp, int *sizeInd_temp);

void SetNI(const int i, const int val, int **neighborsIndices, int *sizeNeighborsIndices);

template <typename T>
void reallocVarPointer(const int i, const int j, int multiplier, const T val, T **id_temp, int *sizeId_temp,
                       const std::string &message) {
  T **id = id_temp;
  int *size = sizeId_temp;

  if (i < *size) {
    *((*id) + i * multiplier + j) = val;
    return;
  }
  while (*size <= i)
    *size += ARRAY_INCREMENT;

  *id = (T *)realloc(*id, (*size) * multiplier * sizeof(T));
  if (*id == NULL) {
    printf("%s\n", message.c_str());
    std::bad_alloc exception;
    throw exception;
  }
  *((*id) + i * multiplier + j) = val;
}

} // namespace arvo
#endif
