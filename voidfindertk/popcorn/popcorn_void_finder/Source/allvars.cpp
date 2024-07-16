#include "allvars.h"

// struct tracers   *Tracer;
struct vorocell *Tessell;
std::vector<grid::halo> Tracer;
std::vector<grid::Sphere> Sphvd;
std::vector<grid::Sphere> seed;
std::vector<grid::popcorn> pop_voids; // the output: popcorn voids
std::vector<grid::par> pairs;
int NVOID;
long long int NTRAC;
float NMEAN;
float NCORES;
float masslim;

struct float3 boxsize; // estructura del paralelepipedo del box

float radmax;
struct conf_var conf;         //
struct conf_var_file confile; // TODO: Primer paso reemplazar el struct global
                              // por algunas variables globales
                              //       segundo paso eliminar la struct global
                              //       por envios de variables.
double Redshift;
double HubbleParam;
