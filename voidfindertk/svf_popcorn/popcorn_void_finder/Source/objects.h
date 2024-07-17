
#ifndef OBJECT_H_
#define OBJECT_H_

#include <fstream>
#include <iostream>
#include <vector>

struct float3 {
  float x;
  float y;
  float z;
};

namespace grid {

class Sphere {
public:
  float x = 0.0;
  float y = 0.0;
  float z = 0.0;
  float radius = 0.0;
  float delta = 0.0;
  float vol = 0.0;
  int lvl = -1;
  int id = -1;
  bool erase = false;
  bool ToF = false;
  int np = -1;
  int irad = -1;

  // ordenamiento por radios
  friend bool operator<(const Sphere &str1, const Sphere &str2);

  friend bool operator>(const Sphere &str1, const Sphere &str2);

  friend std::ostream &operator<<(std::ostream &output, const Sphere &D);

  friend std::ifstream &operator>>(std::ifstream &input, Sphere &D);

  bool operator==(const Sphere &other) const {
    return (x == other.x) && (y == other.y) && (z == other.z) && (radius == other.radius) && (delta == other.delta) &&
           (vol == other.vol) && (lvl == other.lvl) && (id == other.id) && (erase == other.erase) &&
           (ToF == other.ToF) && (np == other.np) && (irad == other.irad);
  }

  Sphere() = default;
  Sphere(float xc, float yc, float zc, float deltac, float volc);
};

//// class rootvds
//// {
////  public:
////   float x,y,z,delta,vol;
////   int id;
////   bool erase;
////   float radius;
////   //constructor
////   rootvds(float xc=0.0, float yc=0.0, float zc=0.0, float deltac=0.0, float
/// volc=0.0);
////
////   friend ifstream &operator>>( ifstream  &input, rootvds &D );
////
////   friend ostream &operator<<( ostream &output, const rootvds &D );
////
////   //ordenamiento por delta
////  // friend bool operator < (const rootvds& str1, const rootvds& str2);
////
////  // friend bool operator > (const rootvds& str1, const rootvds& str2);
//// };
////
struct par {
  int idmin, idmax;
  int id_rad_min;
  int id_rad_max;
  float root_rad_min, root_rad_max;
  float volint, minvol, pop_Delta;
  float joinvol;
  int npart1, npart2, common;
  float vol1, vol2;
  bool isCave;
};
class halo {
public:
  int np;
  float x, y, z;
  float vx, vy, vz;
  bool visit;
  long long int gbin;
  float mass;

  // int id;

  friend std::ifstream &operator>>(std::ifstream &input, halo &D);

  friend std::ostream &operator<<(std::ostream &output, const halo &D);
};

class cell {
public:
  int ix, iy, iz;
  // vector con las particulas
  std::vector<int> in;
  int next[27];
  cell();

  cell(int iix, int iiy, int iiz);
};

class popcorn {
public:
  int id = -1;
  int nmem = 0;
  std::vector<Sphere> membs;
  std::vector<int> in_halos;
  double vol = 0.0;
  int npart = 0;
  int isCave = -1;

  popcorn() = default;

  popcorn(const popcorn &D); // copy constructor

  friend std::ostream &operator<<(std::ostream &output, const popcorn &D);

  friend std::ifstream &operator>>(std::ifstream &input, popcorn &D);

  // ordenamiento por volumen
  friend bool operator<(const popcorn &str1, const popcorn &str2);

  friend bool operator>(const popcorn &str1, const popcorn &str2);

  bool operator==(const popcorn &other) const {
    return (id == other.id) && (nmem == other.nmem) && (membs == other.membs) && (in_halos == other.in_halos) &&
           (vol == other.vol) && (npart == other.npart) && (isCave == other.isCave);
  }
};

class vecv {
public:
  std::vector<int> v;
  vecv();
};

class float_in {
public:
  float d;
  int id;
  unsigned long long int cnt;
  std::vector<int> vec;

  float_in();

  // ordenamiento por distancias
  friend bool operator<(const float_in &str1, const float_in &str2);

  friend bool operator>(const float_in &str1, const float_in &str2);
};

struct small_grid_index_first {
  bool operator()(halo const &a, halo const &b) const { return a.gbin < b.gbin; }
};
} // namespace grid

#endif
