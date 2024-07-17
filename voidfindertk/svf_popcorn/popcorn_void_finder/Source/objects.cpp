#include "objects.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "constants.h" // K_SPHFILE K_OUTFILE K_TRSFILE K_TRSFILE K_THRSLD
                       // K_CRTP K_MINRADIUS K_EPS K_BOXSIZE

using namespace std;
namespace grid {

//  class Sphere

// ordenamiento por radios
bool operator<(const Sphere &str1, const Sphere &str2) { return (str1.radius < str2.radius); }

bool operator>(const Sphere &str1, const Sphere &str2) { return (str1.radius > str2.radius); }

std::ostream &operator<<(std::ostream &output, const Sphere &D) {
  output << std::scientific << std::setprecision(6) << D.x << " " << std::scientific << std::setprecision(6) << D.y
         << " " << std::scientific << std::setprecision(6) << D.z << " " << std::scientific << std::setprecision(6)
         << D.radius << " " << std::scientific << std::setprecision(6) << D.vol << " " << D.lvl << std::endl;
  return output;
}

std::ifstream &operator>>(std::ifstream &input, Sphere &D) {
  input >> D.x >> D.y >> D.z >> D.radius >> D.vol >> D.lvl;
  return input;
}

Sphere::Sphere(float xc, float yc, float zc, float deltac, float volc) {
  x = xc;
  y = yc;
  z = zc;
  delta = deltac;
  vol = volc;
  id = -1;
  erase = false;
};

// class halo

ifstream &operator>>(ifstream &input, halo &D) {
  float npf;
  input >> npf >> D.x >> D.y >> D.z >> D.vx >> D.vy >> D.vz;
  D.np = (int)npf;
  D.mass = npf;
  D.visit = false;
  return input;
}

ostream &operator<<(ostream &output, const halo &D) {
  output << D.np << " " << std::scientific << std::setprecision(6) << D.x << " " << std::scientific
         << std::setprecision(6) << D.y << " " << std::scientific << std::setprecision(6) << D.z << " "
         << std::scientific << std::setprecision(6) << D.vx << " " << std::scientific << std::setprecision(6) << D.vy
         << " " << std::scientific << std::setprecision(6) << D.vz << endl;
  return output;
}

//  class cell

cell::cell() { ix = iy = iz = 0; };

cell::cell(int iix, int iiy, int iiz) {
  int i;
  ix = iix;
  iy = iiy;
  iz = iiz;
  for (i = 0; i < 27; i++)
    next[i] = -1;
};

// class popcorn
popcorn::popcorn(const popcorn &D) {
  id = D.id;
  nmem = D.nmem;
  vol = D.vol;
  npart = D.npart;
  isCave = D.isCave;

  for (long unsigned int i = 0; i < D.membs.size(); i++)
    membs.push_back(D.membs[i]);

  for (long unsigned int i = 0; i < D.in_halos.size(); i++)
    in_halos.push_back(D.in_halos[i]);

}; // copy constructor

ostream &operator<<(ostream &output, const popcorn &D) {

  output << D.id << " " << D.nmem << " " << std::scientific << std::setprecision(6) << D.vol << " " << D.npart << " "
         << D.isCave << endl;
  for (long unsigned int i = 0; i < D.membs.size(); i++)
    output << (D.membs[i]);

  for (long unsigned int i = 0; i < D.in_halos.size(); i++)
    output << (D.in_halos[i]) << endl;
  return output;
}

ifstream &operator>>(ifstream &input, popcorn &D) {
  Sphere pp;
  int i;
  int iiid;

  input >> D.id >> D.nmem >> D.vol >> D.npart >> D.isCave;
  for (i = 0; i < D.nmem; i++) {
    input >> pp;
    D.membs.push_back(pp);
  }

  for (i = 0; i < D.npart; i++) {
    input >> iiid;
    D.in_halos.push_back(iiid);
  }
  return input;
}

// ordenamiento por volumen
bool operator<(const popcorn &str1, const popcorn &str2) { return (str1.vol < str2.vol); }

bool operator>(const popcorn &str1, const popcorn &str2) { return (str1.vol > str2.vol); }

//  class vecv

vecv::vecv() { v.reserve(1); }

// class float_in

float_in::float_in() {
  d = 0.0;
  id = 0;
  cnt = 0;
}

// ordenamiento por distancias
bool operator<(const float_in &str1, const float_in &str2) { return (str1.d < str2.d); }

bool operator>(const float_in &str1, const float_in &str2) { return (str1.d > str2.d); }

struct greater {
  template <class T> bool operator()(T const &a, T const &b) const { return a > b; }
};
} // namespace grid
