#include "io_tracers_stream.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "../allvars.h"
#include "../grid.h"
#include "../t.h"

using namespace grid;

void ReadTracers_stream(std::string fname) {
  clock_t t;
  halo pp;
  FILE *f1;

  fprintf(stdout, "\n Lectura de trazadores\n");
  fflush(stdout);
  fprintf(stdout, " | File = %s\n", fname.c_str());
  fflush(stdout);
  t = clock();

  f1 = fopen(fname.c_str(), "r");

  fread(&NTRAC, sizeof(int), 1, f1);

  float boxsizex, boxsizey, boxsizez;

  fread(&boxsizex, sizeof(float), 1, f1);
  fread(&boxsizey, sizeof(float), 1, f1);
  fread(&boxsizez, sizeof(float), 1, f1);

  // std::vector<float> vox{boxsizex,boxsizey,boxsizez};
  // std::vector<float>::iterator result;
  // result = std::max_element(vox.begin(), vox.end());
  // boxsize= *result;

  boxsize.x = boxsizex;
  boxsize.y = boxsizey;
  boxsize.z = boxsizez;

  while (1) {
    if (fread(&(pp.x), sizeof(float), 1, f1) == 0)
      break;
    fread(&(pp.y), sizeof(float), 1, f1);
    fread(&(pp.z), sizeof(float), 1, f1);
    fread(&(pp.vx), sizeof(float), 1, f1);
    fread(&(pp.vy), sizeof(float), 1, f1);
    fread(&(pp.vz), sizeof(float), 1, f1);
    pp.x = periodicity(pp.x, boxsize.x);
    pp.y = periodicity(pp.y, boxsize.y);
    pp.z = periodicity(pp.z, boxsize.z);
    Tracer.push_back(pp);
  }

  if (NTRAC != (long long int)Tracer.size()) {
    std::cout << "ERROR at reading STREAM File " << NTRAC << " " << Tracer.size() << std::endl;
    exit(EXIT_FAILURE);
  };
  NMEAN = (float)NTRAC / (boxsizex * boxsizey * boxsizez);

  fprintf(stdout, " | Numero de trazadores    = %lld \n", NTRAC);
  fprintf(stdout, " | Densidad media [h/Mpc]³ = %f \n", NMEAN);

  fclose(f1);
  Time(t, 1);
}

void ReadHeaderTracers_stream(std::string fname) {
  FILE *f1 = fopen(fname.c_str(), "r");

  fread(&NTRAC, sizeof(int), 1, f1);

  float boxsizex, boxsizey, boxsizez;

  fread(&boxsizex, sizeof(float), 1, f1);
  fread(&boxsizey, sizeof(float), 1, f1);
  fread(&boxsizez, sizeof(float), 1, f1);

  // std::vector<float> vox{boxsizex,boxsizey,boxsizez};
  // std::vector<float>::iterator result;
  // result = std::max_element(vox.begin(), vox.end());
  // boxsize= *result;

  boxsize.x = boxsizex;
  boxsize.y = boxsizey;
  boxsize.z = boxsizez;

  fclose(f1);
}
