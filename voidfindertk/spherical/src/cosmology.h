
#define WORKSIZE 100000

struct cosmoparam {
  double OmM; 
  double OmL; 
  double OmK;
  double Hub;
};

double evolution_param(double, struct cosmoparam *);
double comoving_distance_integ(double, void *); 
double comoving_distance(double, struct cosmoparam *);
double angular_distance(double, struct cosmoparam *);

