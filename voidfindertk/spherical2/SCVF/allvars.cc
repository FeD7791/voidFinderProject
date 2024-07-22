
#include "allvars.h"

double RadIncrement;
int    NumRanWalk;
int    RandomSeed;
double BoxSize;          
double MaxRadiusSearch;  
double ProxyGridSize;  
double FracRadius;       
double DeltaThreshold;   
double DeltaSeed;        
double OverlapTol;      
int    FormatTracers;
int    NumFiles;
char   FileTracers[MAXCHAR];      
char   FileVoids[MAXCHAR];        
int    OMPcores;         
int    RSDist;           
double Redshift;         
double OmegaMatter;
double OmegaLambda;
double Hubble;           
int    GDist;           
double FidOmegaMatter;        
double FidOmegaLambda;        
double FidHubble;        
int    WriteProfiles;    
double MinProfileDist;   
double MaxProfileDist;   
int    NumProfileBins;   
char   PathProfiles[MAXCHAR];     
double InnerShell;  
double OuterShell;  
double ScalePos;
double ScaleVel;
int    RunFlag;

vector <double> StepTime;
vector <string> StepName;

int    NumTrac;
double MeanNumTrac;
double MeanSeparation;
int    NumVoid;
double LBox[3];
FILE   *logfile;

vector <tracers> Tracer;
vector <voids> Void;

