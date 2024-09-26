#ifndef ALLVARS_H_
#define ALLVARS_H_

#include <string>
#include <vector>

#include "objects.h"

// Constantes
#ifndef PI
#define PI 3.14159265358979323846264
#endif
#define TWOPI (2.0 * PI)
#define HALFPI (0.5 * PI)

// Parametros identificacion
#define NRAND 100 // Cantidad de pasos en recentrado

extern long long int NTRAC; // Numero de trazadores
extern float NMEAN;         // Densidad media de trazadores
extern int NVOID;           // Numero de candidatos a void
extern float radmax;

extern float NCORES;
extern float masslim;

extern double Redshift;
extern double HubbleParam;
// Trazadores
struct vorocell {
  float Cen[3];
  float Delta;
  float Volume;
};

extern struct float3 boxsize; // estructura del paralelepipedo del box
extern struct vorocell *Tessell;

extern std::vector<grid::halo> Tracer;

extern std::vector<grid::Sphere> Sphvd;
extern std::vector<grid::Sphere> seed;
extern std::vector<grid::popcorn> pop_voids; // the output: popcorn voids

extern std::vector<grid::par> pairs;

struct conf_var {
  float3 boxsize;
  float radmax;
};
extern struct conf_var global_conf;

struct conf_var_file {
  std::string pop_file;
  std::string rawpop_file;
  std::string sph_file;
  std::string tracer_file;
  std::string voro_file;
  std::string pairs_file;
  std::string fmt;
  std::string sph_rawfile;
  int num_files;
  bool aux_files;
  float densth;
  float overlap;
  double eps;
  float masslim;
  float min_radius;
  float max_radius;
};
extern struct conf_var_file global_conf_file;

#endif // ALLVARS_H_
