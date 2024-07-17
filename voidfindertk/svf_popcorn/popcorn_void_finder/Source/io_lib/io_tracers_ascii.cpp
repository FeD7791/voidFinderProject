#include "io_tracers_ascii.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "../allvars.h"
#include "../grid.h"
#include "../t.h"

using namespace grid;

void ReadTracers_ascii(std::string fname) {
  std::ifstream file_tracer;
  clock_t t;
  halo pp;

  fprintf(stdout, "\n Lectura de trazadores\n");
  fflush(stdout);
  fprintf(stdout, " | File = %s\n", fname.c_str());
  fflush(stdout);
  t = clock();

  file_tracer.open(fname.c_str());

  while (file_tracer >> pp) {
    pp.x = periodicity(pp.x, boxsize.x);
    pp.y = periodicity(pp.y, boxsize.y);
    pp.z = periodicity(pp.z, boxsize.z);
    Tracer.push_back(pp);
  }
  NTRAC = Tracer.size();
  NMEAN = (float)NTRAC / (boxsize.x * boxsize.y * boxsize.z);

  fprintf(stdout, " | Numero de trazadores    = %lld \n", NTRAC);
  fprintf(stdout, " | Densidad media [h/Mpc]³ = %f \n", NMEAN);

  Time(t, 1);
}

void ReadHeaderTracers_ascii(std::string fname) {
  std::ifstream file_tracer;

  file_tracer.open(fname.c_str());

  std::string line;
  NTRAC = 0;
  while (getline(file_tracer, line))
    NTRAC++;
  file_tracer.close();
}
