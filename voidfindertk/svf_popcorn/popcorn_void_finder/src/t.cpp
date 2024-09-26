
#include "t.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

int CountLines(char *filename) {
  FILE *fp;
  int count = 0;
  char string[256];

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stdout, "ERROR!! %s no existe \n", filename);
    exit(EXIT_FAILURE);
  }

  while (fgets(string, 256, fp))
    count++;
  fclose(fp);

  return count;
}

void Time(clock_t ti, int N) {
  clock_t tf;
  float tseg, tmin;

  tf = clock();

  tseg = (float)(tf - ti) / CLOCKS_PER_SEC / (float)N;
  tmin = tseg / 60.0;

  fprintf(stdout, " | -> tiempo = %f min. \n", tmin);
}

void Progress(int par, int tot) {
  float prog;

  prog = (float)par / (float)tot * 100.0;
  fprintf(stdout, " | %4.1f %s \r", prog, "%");
  fflush(stdout);
}

long long int find_ngrid(long long int ntracers, float radmin, float boxsize) {
  float v;

  // in order to have approx 10 particles by gridi for low density samples
  if (ntracers < 16777216) {
    v = ((float)ntracers) / 10.0;
    v = pow(v, 0.33333333);
  } else {
    v = ((float)boxsize) / radmin;
  }

  // closest power of 2 integer
  long long int ngrid = pow(2, ceil(log(v) / log(2)));
  while (radmin < boxsize / (float)ngrid) {
    ngrid = ngrid * 2;
  }

  return (ngrid);
}
