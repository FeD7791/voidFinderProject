
#include "allvars.h"
#include "velocity.h"
#include "grid.h"
#include "tools.h"

void compute_velocity()
{
   int          i,k,ic,jc,kc,in,m,NumGrid;
   int          ii,jj,kk,l,next,Counter;
   double       xc[3],xt[3],dx[3],vt[3],GridSize[3];
   double       dist,Radius,PLUS,MaxDist,MinDist,GAP;
   struct grid  *GridList;
   int          NumQuery;
   struct query Query;
   clock_t      t;

   fprintf(logfile,"\n COMPUTING VOID BULK VELOCITIES \n");
   t = clock();

   NumGrid = (int)round(cbrt((double)NumTrac/100.0));
   if (NumGrid < 50) NumGrid = 50;
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   build_grid_list(Tracer,NumTrac,GridList,NumGrid,GridSize,false);

   // Selecciono grides

   MinDist = 0.0;
   MaxDist = 0.0;
   for (i=0; i<NumVoid; i++) {
       if (!Void[i].ToF) continue;
       if (Void[i].Rad > MaxDist) MaxDist = Void[i].Rad;
   }
   MaxDist *= 1.5;
   query_grid(&Query,GridSize,MinDist,MaxDist);
   NumQuery = Query.i.size();
   GAP = 0.5*sqrt(3.0)*max_grid_size(GridSize);
  
   fprintf(logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
   fflush(logfile);

   #pragma omp parallel for default(none) schedule(dynamic)      \
    shared(NumVoid,Void,Tracer,NumQuery,Query,LBox,InnerShellVel,\
           OuterShellVel,NumGrid,GridSize,GridList,GAP)          \
   private(i,l,k,m,Radius,xc,ic,jc,kc,ii,jj,kk,next,dx,xt,dist,  \
           Counter,vt,PLUS,in)

   for (i=0; i<NumVoid; i++) {
       
       if (!Void[i].ToF) continue;

       Counter = 0;
       PLUS = 0.0;
       Radius = Void[i].Rad;
       for (k=0; k<3; k++) {
           xc[k] = (double)Void[i].Pos[k];
           Void[i].Vel[k] = 0.0;
       }

       ic = (int)(xc[0]/GridSize[0]);
       jc = (int)(xc[1]/GridSize[1]);
       kc = (int)(xc[2]/GridSize[2]);

       do {

          for (in=0; in<NumQuery; in++) {
	    
	      ii = Query.i[in]; 
	      jj = Query.j[in]; 
	      kk = Query.k[in]; 

	      dist = (double)(ii*ii)*(GridSize[0]*GridSize[0])
	           + (double)(jj*jj)*(GridSize[1]*GridSize[1])
	           + (double)(kk*kk)*(GridSize[2]*GridSize[2]);
	      dist = sqrt(dist);

	      if (dist > 1.5*Radius+GAP) continue;

              ii = periodic_grid(ii + ic,NumGrid); 
	      jj = periodic_grid(jj + jc,NumGrid); 
	      kk = periodic_grid(kk + kc,NumGrid); 	  

              l = index_1d(ii,jj,kk,NumGrid);

              if (GridList[l].NumMem == 0) continue;

              for (m=0; m<GridList[l].NumMem; m++) {
		
	          next = GridList[l].Member[m];	

                  for (k=0; k<3; k++) {
                      xt[k] = (double)Tracer[next].Pos[k];	 
                      vt[k] = (double)Tracer[next].Vel[k];
                      dx[k] = periodic_delta(xc[k] - xt[k],LBox[k]);
                  }

                  dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
                  dist /= Radius;

                  if (dist > InnerShellVel-PLUS && dist < OuterShellVel+PLUS) {
                     Void[i].Vel[0] += vt[0];
                     Void[i].Vel[1] += vt[1];
                     Void[i].Vel[2] += vt[2];
                     Counter++;	 
                  }
              } 
          }

          PLUS += 0.05;     

       } while (Counter == 0);

       Void[i].Vel[0] /= (double)Counter;
       Void[i].Vel[1] /= (double)Counter;
       Void[i].Vel[2] /= (double)Counter; 
   }
  
   free_grid_list(GridList,NumGrid);
   free_query_grid(&Query);

   StepName.push_back("Computing velocities");
   StepTime.push_back(get_time(t,OMPcores));
}
