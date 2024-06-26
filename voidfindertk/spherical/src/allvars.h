
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <vector>
#include "voro++.hh"

using namespace std;
using namespace voro;

#define MAXCHAR 500

// Some constants

#define PI   3.141592653589793238
#define CVEL  299792.469 

extern double RadIncrement;
extern int    NumRanWalk;
extern double BoxSize;          
extern double MaxRadiusSearch;  
extern double ProxyGridSize;  
extern double FracRadius;       
extern double DeltaThreshold;   
extern double DeltaSeed;        
extern double OverlapTol;   
extern int    FormatTracers;
extern int    NumFiles;
extern char   FileTracers[MAXCHAR];      
extern char   FileVoids[MAXCHAR];        
extern int    OMPcores;         
extern int    RSDist;           
extern double Redshift;         
extern double OmegaMatter;
extern double OmegaLambda;
extern double Hubble;           
extern int    GDist;           
extern double FidOmegaMatter;        
extern double FidOmegaLambda;        
extern double FidHubble;        
extern int    WriteProfiles;    
extern double MinProfileDist;   
extern double MaxProfileDist;   
extern int    NumProfileBins;   
extern char   PathProfiles[MAXCHAR];     
extern double InnerShellVel;  
extern double OuterShellVel;  
extern double ScalePos;
extern double ScaleVel;
extern int    RunFlag;

extern vector <double> StepTime;
extern vector <string> StepName;   

extern int    NumTrac;
extern double MeanNumTrac;
extern double MeanSeparation;
extern int    NumVoid;
extern double LBox[3];
extern FILE   *logfile;

// Tracers

struct tracers {
  float Pos[3];
  float Vel[3];
  float Cen[3];
  float Delta;
  float Volume;
  float size;
}; 
extern vector <tracers> Tracer;

// Voids

struct voids {
  float Rad;
  float Rini;
  float Ini[3];
  float Pos[3];
  float Vel[3];
  float Dtype;
  float Delta;
  float Poisson;
  float Dist4;
  bool  ToF;  
  int   Nran;
};
//This struct is to put filtered values from struct voids. If all the data is needed just use struct voids instead
struct p_voids{
  int n_voids;
  float Rad;
  float Pos[3];
  float Vel[3];
  float Dtype;
  float Delta;
  float Poisson;
  int Nran;
  
};
extern vector <voids> Void;
//vector void for python
///ADDING RETURN
   struct voidArray{	
      struct p_voids voids_arr[1];
	
   };
///

struct InputParams
{
    double RadIncrement;
    double BoxSize;          
    double MaxRadiusSearch;  
    double ProxyGridSize;  
    double FracRadius;       
    double DeltaThreshold;   
    double DeltaSeed;        
    double OverlapTol;   
    double Redshift;         
    double OmegaMatter;
    double OmegaLambda;
    double Hubble;  
    double FidOmegaMatter;        
    double FidOmegaLambda;        
    double FidHubble;  
    double MinProfileDist;   
    double MaxProfileDist;   
    double ScalePos;
    double ScaleVel;
    double InnerShellVel;  
    double OuterShellVel; 

    int    FormatTracers;
    int    NumFiles;
    int    NumRanWalk;
    int    OMPcores;         
    int    RSDist;           
    int    GDist;                     
    int    WriteProfiles;        
    int    NumProfileBins;       
    //int    RunFlag;
};

