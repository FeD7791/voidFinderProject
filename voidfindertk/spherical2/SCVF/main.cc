
#include "allvars.h"
#include "io.h"
#include "voronoi.h" 
#include "finder.h"
#include "velocity.h"
#include "profiles.h"
#include "tools.h"

int main(int argc, char **argv) 
{
   
   if (argc < 2) {
       fprintf(stdout, "\n Error. Missing input file and flags.\n");
       fprintf(stdout, "./main.x <input_param> [<run flag>] [...] \n\n");
       exit(EXIT_FAILURE);
   } else if (argc == 2) {
       RunFlag = 0;	   
   } else {
       sscanf(argv[2], "%d", &RunFlag);   
   }

   read_input_file(argv[1]);

   if (RunFlag == 1) {
      if (argc < 4) {
	 fprintf(stdout, "\n Error. Missing void ID.\n");
         fprintf(stdout, "./main.x <input_param> 1 <void_ID> \n\n");
         exit(EXIT_FAILURE);
      }	      
      int voidID; 
      sscanf(argv[3], "%d", &voidID);
      bin2ascii_profile(voidID); 
      exit(EXIT_SUCCESS);
   } 

   if (Redshift == 0.0 && GDist == 1) {
      fprintf(stdout,"\nError. Geometrical distortions not available for z = 0\n");
      exit(EXIT_FAILURE);
   }

   omp_set_num_threads(OMPcores);
   fprintf(stdout,"\n ====>>>> Void finder runnning in %d core(s) <<<<==== \n",OMPcores);

   fprintf(stdout,"\nReading tracers... ");fflush(stdout);
   read_tracers();
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

   fprintf(stdout,"\nComputing void velocities and mass centers... ");fflush(stdout);
   compute_velocity_and_masscenter();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nComputing void profiles... ");fflush(stdout);
   compute_profiles();
   fprintf(stdout,"Done.\n");fflush(stdout);

   fprintf(stdout,"\nWrinting void catalogue... ");fflush(stdout);
   write_voids();
   fprintf(stdout,"Done.\n\n");fflush(stdout);

   time_resume();
   
   Tracer.clear();
   Void.clear();
   fclose(logfile);

   return(0);
}
