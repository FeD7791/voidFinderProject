
#include "allvars.h"
#include "voronoi.h"
#include "tools.h"
#include "grid.h"

void compute_voronoi()
{
   int            G,i,j,k,p,N,l,id,count,NumGrid;
   struct grid    *GridList;
   vector <int>   IDs,indx;
   double         min[3],max[3],ref[3],GridSize[3];   
   double         xp[3],xc[3],rr,Vol;
   bool           check;
   voronoicell    cell;
   container      *con;
   c_loop_order   *clo;
   particle_order *po;
   FILE           *fd;
   clock_t        t;

   fprintf(logfile,"\n COMPUTING VORONOI TESSELLATION \n");
  
   t = clock();

   NumGrid = (int)round(cbrt((double)NumTrac/MeanPartPerGrid));
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   build_grid_list(Tracer,NumTrac,GridList,NumGrid,GridSize,true);

   Vol = LBox[0]*LBox[1]*LBox[2];

   #pragma omp parallel for default(none) schedule(static)           \
    shared(Vol,stdout,GridSize,NumGrid,LBox,Tracer,GridList,NumTrac) \
    private(l,N,k,j,i,ref,min,max,xp,xc,G,count,p,id,IDs,rr,check,   \
            cell,con,clo,po,indx)                                    \

   for (l=0; l<NumGrid*NumGrid*NumGrid; l++) {

       // Count number of particles
        
       N = 0;
       for (j=0; j<27; j++) {
           i = GridList[l].Neighbour[j];	  
           N += GridList[i].NumMem;
       }
   
       // Global parameters

       G = (int)cbrt((double)N/5.0);

       indx = index_3d(l,NumGrid);
   
       for (k=0; k<3; k++) {
           ref[k] = (double)(indx[k]  )*GridSize[k];
           min[k] = (double)(indx[k]-1)*GridSize[k];
           max[k] = (double)(indx[k]+2)*GridSize[k];	      
       }

       con = new container(min[0],max[0],min[1],max[1],min[2],max[2],G,G,G,false,false,false,8);
       po = new particle_order();
       count = 0;

       // Load particles
       
       for (j=0; j<27; j++) {

           i = GridList[l].Neighbour[j];

           for (p=0; p<GridList[i].NumMem; p++) {
               
               id = GridList[i].Member[p];

     	       for (k=0; k<3; k++) {
     	           xp[k] = (double)Tracer[id].Pos[k];
                   if (xp[k] - ref[k] >  0.5*LBox[k]) xp[k] -= LBox[k];	      
                   if (xp[k] - ref[k] < -0.5*LBox[k]) xp[k] += LBox[k];	      
     	       }

     	       if (i == l) 
                  (*con).put(*po,count,xp[0],xp[1],xp[2]);
     	       else
                  (*con).put(count,xp[0],xp[1],xp[2]);
            
               IDs.push_back(id);
               count++;
           }
       }

       // Create loop class
       
       clo = new c_loop_order(*con,*po);

       // Compute Voronoi cells (only those in *po)

       if ((*clo).start()) do {

          // Id and position of the tracer
          (*clo).pos(i,xp[0],xp[1],xp[2],rr);	 

          check = (*con).compute_cell(cell,*clo);
          if (!check) {
             fprintf(stdout,"\n Error. Voronoi cell can not be computed: i = %d \n",i);
     	     fflush(stdout);
             exit(EXIT_FAILURE);	 
          } 
          
          // Centroid position
          cell.centroid(xc[0],xc[1],xc[2]);

          for (k=0; k<3; k++)  
     	      Tracer[IDs[i]].Cen[k] = (float)periodic_position(xp[k] + xc[k],LBox[k]);
         
          // Volume of the cell
          Tracer[IDs[i]].Volume = (float)cell.volume(); 
       
          // Delta
          Tracer[IDs[i]].Delta = (float)((Vol/Tracer[IDs[i]].Volume)/(double)(NumTrac) - 1.0);

       } while ((*clo).inc());

       (*con).clear();
       delete con;
       delete po;
       delete clo;
       IDs.clear();
   }

   free_grid_list(GridList,NumGrid);

   StepName.push_back("Computing Voronoi tessellation");
   StepTime.push_back(get_time(t,OMPcores));

}
