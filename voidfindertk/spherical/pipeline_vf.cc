
#include "allvars.h"
#include "io.h"
#include "voronoi.h" 
#include "finder.h"
#include "velocity.h"
#include "profiles.h"
#include "tools.h"




void print_void_array(struct voidArray* va) {
  if (va == NULL) {
    printf("voidArray is NULL.\n");
    return;
  }

  printf("Rad: %3.3f\n", va->voids_arr[0].Rad);
  printf("Pos: (%3.3f, %f, %f)\n", va->voids_arr[0].Pos[0], va->voids_arr[0].Pos[1], va->voids_arr[0].Pos[2]);
  //printf("Vel: (%3.3f, %f, %f)\n", va->voids_arr[0].Vel[0], va->voids_arr[0].Vel[1], va->voids_arr[0].Vel[2]);
  //printf("Delta: %f\n", va->voids_arr[0].Delta);
  //printf("Dtype: %d\n", va->voids_arr[0].Dtype);
  //printf("Poisson: %d\n", va->voids_arr[0].Poisson);
  //printf("Dist4: %f\n", va->voids_arr[0].Dist4);
  //printf("Nran: %d\n", va->voids_arr[0].Nran);

  // Add a newline for visual clarity:
  printf("\n");
}



int main(){
    return 0;
}



std::vector <tracers>  load_tracer(
    double *arr_x,
    double *arr_y,
    double *arr_z,
    double *arr_vx,
    double *arr_vy,
    double *arr_vz,
    double *arr_m,
    int size){
    
    

    NumTrac = 0;  
    
    std::vector <tracers> vector_tracer;
    for (int i=0; i<size; i++) {

       

       vector_tracer.push_back(tracers());

       
       vector_tracer.at(i).Pos[0] = arr_x[i];
       vector_tracer.at(i).Pos[1] = arr_y[i];
       vector_tracer.at(i).Pos[2] = arr_z[i];
       vector_tracer.at(i).Vel[0] = arr_vx[i];
       vector_tracer.at(i).Vel[1] = arr_vy[i];
       vector_tracer.at(i).Vel[2] = arr_vz[i];
       
       
       NumTrac++;
   }
   
   

   


   return vector_tracer;
}



void tracer_verification(std::vector<tracers>&input_Tracer){

   

   
   
   clock_t t = clock();
   float xmax = 0.0;
   float ymax = 0.0;
   float zmax = 0.0;
   
   float diff = 0.999;
   int NumTrac = input_Tracer.size();


   
   for (int i=0; i<NumTrac; i++) {
       if (input_Tracer[i].Pos[0] > xmax) xmax = input_Tracer[i].Pos[0];	   
       if (input_Tracer[i].Pos[1] > ymax) ymax = input_Tracer[i].Pos[1];	   
       if (input_Tracer[i].Pos[2] > zmax) zmax = input_Tracer[i].Pos[2];
       if (input_Tracer.at(i).Pos[0] > xmax) xmax = input_Tracer.at(i).Pos[0];	   
       if (input_Tracer.at(i).Pos[1] > ymax) ymax = input_Tracer.at(i).Pos[1];	   
       if (input_Tracer.at(i).Pos[2] > zmax) zmax = input_Tracer.at(i).Pos[2];	   
   }
   
   
   

   if (xmax/BoxSize < diff || ymax/BoxSize < diff || zmax/BoxSize < diff) {
      fprintf(stdout,"\n Error. Wrong BoxSize? - MAX = (%f,%f,%f) \n",xmax,ymax,zmax);
      fflush(stdout);
      exit(EXIT_FAILURE);	
   } 
   if (xmax/BoxSize < diff || ymax/BoxSize < diff || zmax/BoxSize < diff) {
      printf("\n Error. Wrong BoxSize? - MAX = (%f,%f,%f) \n",xmax,ymax,zmax);
      
      exit(EXIT_FAILURE);	
   } 
   LBox[0] = BoxSize;
   LBox[1] = BoxSize;
   LBox[2] = BoxSize;

   

   if (RSDist == 1) redshift_space_distortions();
   if (GDist == 1) geometrical_distortions();

   double Volume = LBox[0]*LBox[1]*LBox[2];
   MeanNumTrac = (double)NumTrac/Volume;
   MeanSeparation = cbrt(Volume/(double)NumTrac);

   fprintf(logfile,"\n READING TRACERS \n");
 
   fprintf(logfile," | Number of tracers = %d \n",NumTrac);
   fprintf(logfile," | Size of the box: x-axis = %f \n",LBox[0]);
   fprintf(logfile," |                  y-axis = %f \n",LBox[1]);
   fprintf(logfile," |                  z-axis = %f \n",LBox[2]);
   fprintf(logfile," | Mean number density [h³/Mpc³] = %e \n",MeanNumTrac);
   fprintf(logfile," | Mean separation [Mpc/h] = %e \n",MeanSeparation);

   


	


}



void load_input_file()
{



  OMPcores= 8;  
  BoxSize = 500;
  MaxRadiusSearch = 40.0;
  ProxyGridSize =     5.0 ;
  DeltaThreshold =   -0.9 ;
  DeltaSeed =        -0.7 ;
  OverlapTol =        0 ;

  FormatTracers=     0 ;
  NumFiles=          32 ;
  
  //strcpy(FileTracers, "data/galaxies.dat" );
  strcpy(FileVoids, "voids_file.dat");

  ScalePos = 1 ;
  ScaleVel = 1 ;

  FracRadius = 0.5 ;
  NumRanWalk = 75  ;
  RadIncrement = 0 ;

  RSDist = 0;
  Redshift = 0.99   ;
  OmegaMatter = 0.25;
  OmegaLambda = 0.75;
  Hubble = 0.73;

  GDist = 0;
  FidOmegaMatter = 0.2;  
  FidOmegaLambda = 0.8 ;
  FidHubble = 0.7;

  WriteProfiles = 0 ;
  MinProfileDist = 0.5;
  MaxProfileDist = 3.0;
  NumProfileBins = 100;
  //PathProfiles= "data/profiles/";
  strcpy(PathProfiles, "data/profiles/");

  InnerShellVel = 0.8;
  OuterShellVel = 1.2;
  
}

int ToF_void_counter(int n_void){
//Not all voids finded are good candidates, only the one with Void[j]=True are, this function counts how many of this true voids are
	int counter = 0;
	for(int j=0; j<n_void; j++){
	if (Void[j].ToF) {
	counter = counter + 1;
}
	}
return counter;
}



extern "C"{

struct voidArray* arr_void_creator(int n_eff_voids){
   

   struct voidArray* va = (struct voidArray*)malloc((n_eff_voids)*sizeof(struct voidArray));////

   if (va == NULL) {
        // Handle memory allocation failure
        return NULL;
    }
	

int counter = 0;
for (int i=0; i<NumVoid; i++) {
	if (Void[i].ToF) {
	    va->voids_arr[counter].n_voids = n_eff_voids;
            va->voids_arr[counter].Rad = Void[i].Rad;
            va->voids_arr[counter].Pos[0] = Void[i].Pos[0];
            va->voids_arr[counter].Pos[1] = Void[i].Pos[1];
            va->voids_arr[counter].Pos[2] = Void[i].Pos[2];
	    va->voids_arr[counter].Vel[0] = Void[i].Vel[0];
	    va->voids_arr[counter].Vel[1] = Void[i].Vel[1];
	    va->voids_arr[counter].Vel[2] = Void[i].Vel[2];
	    va->voids_arr[counter].Dtype = Void[i].Dtype;
	    va->voids_arr[counter].Delta = Void[i].Delta;
	    va->voids_arr[counter].Poisson = Void[i].Poisson;
	    va->voids_arr[counter].Nran = Void[i].Nran;
	    

	counter = counter + 1;
        }


}

//va->n_voids = counter;//////
return va;
}


struct voidArray* execute_void_finder(
    double *arr_x,
    double *arr_y,
    double *arr_z,
    double *arr_vx,
    double *arr_vy,
    double *arr_vz,
    double *arr_m , int size) 

{

clock_t t = clock();
logfile = fopen("log_file.txt", "a");

   load_input_file();



   if (Redshift == 0.0 && GDist == 1) {
      fprintf(stdout,"\nError. Geometrical distortions not available for z = 0\n");
      exit(EXIT_FAILURE);
   }

   omp_set_num_threads(OMPcores);
   fprintf(stdout,"\n ====>>>> Void finder runnning in %d core(s) <<<<==== \n",OMPcores);



   fprintf(stdout,"\nReading tracers... ");fflush(stdout);
   
   Tracer = load_tracer(
     arr_x,
     arr_y,
     arr_z,
     arr_vx,
     arr_vy,
     arr_vz,
     arr_m,
     size
   );
	
   
   tracer_verification(Tracer);
	StepName.push_back("Reading tracers");

	StepTime.push_back(get_time(t,1));
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing Voronoi tessellation... ");fflush(stdout); 
   compute_voronoi();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nSearching candidates... ");fflush(stdout);
   find_void_candidates(); 
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nPerforming void identification... ");fflush(stdout);
   find_voids();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nCleaning void catalogue... ");fflush(stdout);
   clean_voids();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void velocities... ");fflush(stdout);
   compute_velocity();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void profiles... ");fflush(stdout);
   compute_profiles();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nWrinting void catalogue... ");fflush(stdout);
   write_voids();
   fprintf(stdout,"Done.\n\n");fflush(stdout);
	
	int n_eff_voids = ToF_void_counter(NumVoid); //In order to not have garbage in memory we use n_eff_voids tho allocate just the right space in memory
	struct voidArray* va = (struct voidArray*)malloc((n_eff_voids)*sizeof(struct voidArray));
	
	
  	
	//va = arr_void_creator(NumVoid);
	va = arr_void_creator(n_eff_voids);
	
	
   
	time_resume();
	Tracer.clear();
   	Void.clear();

	fclose(logfile);

	return va;

   
 

}
}








