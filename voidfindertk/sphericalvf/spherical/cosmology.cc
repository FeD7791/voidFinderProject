
#include "allvars.h"
#include "cosmology.h"
#include "gsl_integration.h"
#include "gsl_math.h"

double evolution_param(double z, struct cosmoparam *C)
{
  double q = C->OmM*(1.0 + z)*(1.0 + z)*(1.0 + z)
           + C->OmK*(1.0 + z)*(1.0 + z) 
           + C->OmL;

  return(sqrt(q));
}

double comoving_distance_integ(double z, void *p) 
{
  struct cosmoparam *C = (struct cosmoparam *)p;

  double q = 1.0/evolution_param(z,C);
  
  return(q);
}

double comoving_distance(double z, struct cosmoparam *C)
{
  double result,abserr;
  gsl_function F;
  gsl_integration_workspace *workspace;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &comoving_distance_integ;
  F.params = C;
  
  gsl_integration_qag(&F,0.0,z,0,1.0e-8,WORKSIZE,GSL_INTEG_GAUSS41,workspace,&result,&abserr);

  double distance = (CVEL/C->Hub)*result;

  gsl_integration_workspace_free(workspace);

  return(distance);
}

double angular_distance(double z, struct cosmoparam *C)
{

  double dcom = comoving_distance(z,C);
  double dh = CVEL/C->Hub;
  double dang;

  if (C->OmK == 0.0) {
     dang = dcom;
  } else if (C->OmK > 0.0) {
     dang = (dh/sqrt(C->OmK))*sinh(sqrt(C->OmK)*dcom/dh);	  
  } else if (C->OmK < 0.0) {
     dang = (dh/sqrt(fabs(C->OmK)))*sin(sqrt(fabs(C->OmK))*dcom/dh);	  
  }

  return(dang/(1.0 + z));  
}	

