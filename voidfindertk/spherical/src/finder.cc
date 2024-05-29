#include "allvars.h"
#include "finder.h"
#include "grid.h"
#include "qsort.h"
#include "tools.h"

void find_void_candidates()
{
  int     i;
  clock_t t;
   
  fprintf(logfile,"\n SELECTING UNDERDENSE REGIONS \n");
  t = clock();

  NumVoid = 0;
  for (i=0; i<NumTrac; i++) {
      if (Tracer[i].Delta <= DeltaSeed) {
         Void.push_back(voids());
         Void[NumVoid].Pos[0] = Tracer[i].Cen[0];	 
         Void[NumVoid].Pos[1] = Tracer[i].Cen[1];	 
         Void[NumVoid].Pos[2] = Tracer[i].Cen[2];
	 Void[NumVoid].Ini[0] = Tracer[i].Cen[0];	 
         Void[NumVoid].Ini[1] = Tracer[i].Cen[1];	 
         Void[NumVoid].Ini[2] = Tracer[i].Cen[2];
         Void[NumVoid].Rini = 1.5*cbrt(0.75*(double)Tracer[i].Volume/PI);
         Void[NumVoid].Rad = 0.0;
         Void[NumVoid].ToF = false;
	 Void[NumVoid].Nran = 0;
         NumVoid++;	 
      }
  }

  fprintf(logfile," | Void candidates = %d \n",NumVoid);
  
  StepName.push_back("Finding centers");
  StepTime.push_back(get_time(t,1));

}

void find_voids() 
{
  struct grid    *GridList;
  int            NumCores,NumGrid;
  int            iv,CheckRan,TotRan,ir,next,l,kappa,p,in;
  int            ic,jc,kc,ii,jj,kk,k,Nsort,m;
  double         Radius,BiggestRadius,lambda,MinDist,MaxDist;
  double         dx[3],xr[3],xc[3],dist,GridSize[3];
  double         the,phi,rad,Volume,Delta;
  vector <float> val;
  struct sort    *SortArr;
  bool           done;
  clock_t        t;

  fprintf(logfile,"\n VOID IDENTIFICATION \n");
  t = clock();

  srand(time(NULL));
  
  NumGrid = (int)(BoxSize/ProxyGridSize);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
  build_grid_list(Tracer,NumTrac,GridList,NumGrid,GridSize,false);

// int          NumShell = (int)round(0.5*MaxRadiusSearch/max_grid_size(GridSize));
  int          NumShell = (int)round(MaxRadiusSearch/max_grid_size(GridSize));
  int          NumQuery[NumShell];
  struct query Query[NumShell];

  // Selecciono vecinos

  if (OMPcores > NumShell) 
     NumCores = NumShell;
  else
     NumCores = OMPcores;	

  #pragma omp parallel for default(none) num_threads(NumCores)    \
   shared(NumShell,NumQuery,Query,GridSize,stdout,MaxRadiusSearch)\
   private(p,MinDist,MaxDist)

  for (p=0; p<NumShell; p++) {
      MinDist = MaxRadiusSearch/(double)NumShell*(double)p;  
      MaxDist = MaxRadiusSearch/(double)NumShell*(double)(p+1);  
      query_grid(&Query[p],GridSize,MinDist,MaxDist);
      NumQuery[p] = Query[p].i.size();
  } 

  for (p=0; p<NumShell; p++) {
      MinDist = MaxRadiusSearch/(double)NumShell*(double)p;  
      MaxDist = MaxRadiusSearch/(double)NumShell*(double)(p+1); 
      fprintf(logfile," | Shell N° %2d: MinDist - MaxDist = %5.2f - %5.2f [Mpc/h], %5d grids (Overlap = %f) \n",
		      p,MinDist,MaxDist,NumQuery[p],0.5*sqrt(3.0)*max_grid_size(GridSize));
  }
  fflush(logfile);

  #pragma omp parallel for default(none) schedule(static)                     \
   shared(NumVoid,MeanNumTrac,Void,Tracer,Query,NumQuery,NumShell,LBox,stdout,\
          MaxRadiusSearch,FracRadius,DeltaThreshold,NumGrid,GridSize,GridList,\
	  OMPcores,RadIncrement,NumRanWalk)                                   \
   private(iv,ir,ic,jc,kc,ii,jj,kk,xc,xr,dx,l,Radius,BiggestRadius,next,k,    \
           dist,val,kappa,SortArr,Nsort,the,phi,rad,Volume,Delta,lambda,p,m,  \
	   MinDist,MaxDist,done,in,CheckRan,TotRan)

  for (iv=0; iv<NumVoid; iv++) {

      BiggestRadius = 0.1; 
      TotRan = 0;
      CheckRan = 0;

      do {

	  TotRan++;
	  CheckRan++;

	  if (TotRan == 1) {

             for (k=0; k<3; k++) 
		 xr[k] = 0.0; 

	  } else {

	     the = acos(2.0*random_number() - 1.0);
	     phi = 2.0*PI*random_number();
	     rad = (double)Void[iv].Rad;

	     if (rad == 0.0) 
	        rad = (double)Void[iv].Rini*random_number();
	     else  
	        rad *= FracRadius*random_number();

	     xr[0] = rad*sin(the)*cos(phi);
	     xr[1] = rad*sin(the)*sin(phi);
	     xr[2] = rad*cos(the);
	  }

	  if (Void[iv].Rad == 0.0) {

             for (k=0; k<3; k++) 
		 xc[k] = periodic_position((double)Void[iv].Ini[k] + xr[k],LBox[k]);

	  } else {

             for (k=0; k<3; k++) {
	         xc[k] = periodic_position((double)Void[iv].Pos[k] + xr[k],LBox[k]);
	         dx[k] = periodic_delta(xc[k] - (double)Void[iv].Ini[k],LBox[k]);
	     }
	     dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	     if (dist > Void[iv].Rad) // avoid big migration 
	        for (k=0; k<3; k++) 
	            xc[k] = periodic_position((double)Void[iv].Ini[k] + xr[k],LBox[k]); 	  
	  }

	  ic = (int)(xc[0]/GridSize[0]);
	  jc = (int)(xc[1]/GridSize[1]);
	  kc = (int)(xc[2]/GridSize[2]);

	  done = false;
	  p = 0;

	  do {

             MinDist = MaxRadiusSearch/(double)NumShell*(double)p;		  
             MaxDist = MaxRadiusSearch/(double)NumShell*(double)(p+1);		  

	     for (in=0; in<NumQuery[p]; in++) {

	         ii = periodic_grid(Query[p].i[in] + ic,NumGrid); 
	         jj = periodic_grid(Query[p].j[in] + jc,NumGrid); 
	         kk = periodic_grid(Query[p].k[in] + kc,NumGrid); 

		 l = index_1d(ii,jj,kk,NumGrid);

		 if (GridList[l].NumMem == 0) continue;

		 for (m=0; m<GridList[l].NumMem; m++) {
		
	             next = GridList[l].Member[m];		 

		     for (k=0; k<3; k++) 
	                 dx[k] = periodic_delta(xc[k] - (double)Tracer[next].Pos[k],LBox[k]);

                     dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	             if (dist >= MinDist && dist < MaxDist) 
	                val.push_back(dist);	 

		 }

	     }

	     Nsort = val.size();
	     SortArr = (struct sort *) malloc(Nsort*sizeof(struct sort));
	     for (k=0; k<Nsort; k++) {
	         SortArr[k].val = val[k];
                 SortArr[k].ord = k;	      
	     }
             qsort(SortArr,0,Nsort-1);

	     Radius = 0.5*(SortArr[Nsort-2].val + SortArr[Nsort-1].val);    
   	     Volume = (4.0/3.0)*PI*Radius*Radius*Radius;
	     Delta = (double)(Nsort-1)/Volume/MeanNumTrac - 1.0;

	     if (Delta < DeltaThreshold) {
	        p++;    
		free(SortArr); 
	     } else {
		done = true;
		val.clear();
	     }            

	  } while(!done && p != NumShell);

	  if (!done && p == NumShell) {
             fprintf(stdout,"\n Error. Increase MaxRadiusSearch.\n");
	     fflush(stdout);
             exit(EXIT_FAILURE);	     
	  }

	  for (ir=0; ir<Nsort-3; ir++) {
	
	      Radius = 0.5*(SortArr[ir].val + SortArr[ir+1].val);    
   	      Volume = (4.0/3.0)*PI*Radius*Radius*Radius;
	      Delta = (double)(ir+1)/Volume/MeanNumTrac - 1.0;

              if (Delta < DeltaThreshold && Radius > BiggestRadius) {
		 
		 if (Radius/BiggestRadius - 1.0 >= RadIncrement) CheckRan = 0;
		
	         Void[iv].Rad = (float)Radius;
		 Void[iv].Delta = (float)Delta;
	         Void[iv].Pos[0] = (float)xc[0];	       
	         Void[iv].Pos[1] = (float)xc[1];	       
	         Void[iv].Pos[2] = (float)xc[2];
	         Void[iv].ToF = true;
		 kappa = ir + 1; 
	         BiggestRadius = Radius;
		 Void[iv].Dist4 = (SortArr[ir+3].val - SortArr[ir].val) / Void[iv].Rad;
	      
	      } /* Fin lazo Dcum < DeltaThreshold */

	  } /* Fin lazo bines */

	  free(SortArr);
 
      } while (CheckRan < NumRanWalk); /* Fin lazo random */      

      if (Void[iv].ToF) {
	 lambda = (4.0/3.0)*PI*pow((double)Void[iv].Rad,3)*MeanNumTrac;
         Void[iv].Poisson = (double)kappa*log(lambda) - lambda - ln_factorial(kappa); 	 
	 Void[iv].Nran = TotRan;
      }    

  } /* Fin lazo voids */

  for (p=0; p<NumShell; p++) 
      free_query_grid(&Query[p]);	  
  free_grid_list(GridList,NumGrid);

  StepName.push_back("Finding voids");
  StepTime.push_back(get_time(t,OMPcores));

}

void clean_voids()
{
  int          i,k,l,m,indx,next,in,it;
  int          ii,jj,kk,ic,jc,kc,NumTrueVoid;
  int          NumQuery,NumGrid;
  double       dist,Ri,Rj,Vi,Vj,Vij,MinDist,MaxDist;
  double       xi[3],xj[3],dx[3],GridSize[3];
  struct query Query;
  struct sort  *SortArr; 
  struct grid  *GridList;
  clock_t      t;

  fprintf(logfile,"\n CLEANING VOID CATALOGUE BY OVERLAP (TOL = %4.2f) \n",OverlapTol);
  t = clock();

  NumGrid = (int)cbrt((double)NumVoid/10.0);
  GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid)); 
  build_grid_list(Void,NumVoid,GridList,NumGrid,GridSize,false);
  
  NumTrueVoid = 0;
  for (i=0; i<NumVoid; i++)
      if (Void[i].ToF) 
	 NumTrueVoid++;	  
  
  fprintf(logfile," | Number of true voids before cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)NumVoid*100.0);
  fflush(logfile);

  SortArr = (struct sort *) malloc(NumTrueVoid*sizeof(struct sort));

  it = 0;
  for (i=0; i<NumVoid; i++) {
      if (Void[i].ToF) {	      
         SortArr[it].val = Void[i].Rad;
         SortArr[it].ord = i;
	 it++;
      }     
  }

  qsort(SortArr,0,NumTrueVoid-1);

  // Selecciono vecinos
  
  MinDist = 0.0;
  MaxDist = 2.0*MaxRadiusSearch;  

  query_grid(&Query,GridSize,MinDist,MaxDist);
  NumQuery = Query.i.size();
  
  fprintf(logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
  fflush(logfile);

  for (i=NumTrueVoid-1; i>=0; i--) {

      indx = SortArr[i].ord;
      if (!Void[indx].ToF) continue;

      xi[0] = (double)Void[indx].Pos[0];      
      xi[1] = (double)Void[indx].Pos[1];      
      xi[2] = (double)Void[indx].Pos[2];
      Ri = (double)Void[indx].Rad;
      Vi = (4.0/3.0)*PI*Ri*Ri*Ri;

      ic = (int)(xi[0]/GridSize[0]);
      jc = (int)(xi[1]/GridSize[1]);
      kc = (int)(xi[2]/GridSize[2]);
    
      for (in=0; in<NumQuery; in++) {

   	  ii = periodic_grid(Query.i[in] + ic,NumGrid); 
	  jj = periodic_grid(Query.j[in] + jc,NumGrid); 
	  kk = periodic_grid(Query.k[in] + kc,NumGrid); 

	  l = index_1d(ii,jj,kk,NumGrid);

          if (GridList[l].NumMem == 0) continue;

	  for (m=0; m<GridList[l].NumMem; m++) {
		
	      next = GridList[l].Member[m];		

	      if (Void[next].ToF && next != indx) {

                 xj[0] = (double)Void[next].Pos[0];      
                 xj[1] = (double)Void[next].Pos[1];      
                 xj[2] = (double)Void[next].Pos[2];
                 Rj = (double)Void[next].Rad;
	         Vj = (4.0/3.0)*PI*Rj*Rj*Rj;
    
	         for (k=0; k<3; k++) 
	             dx[k] = periodic_delta(xi[k] - xj[k],LBox[k]);

                 dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	         // Volumen de la interseccion, normalizado a Vj
                 // Solo vale para d <= Ri+Rj

	         if (dist > Ri + Rj) continue;

	         Vij  = (PI/12.0/dist);
                 Vij *= pow(Ri + Rj - dist,2);
	         Vij *= (pow(dist,2) + 2.0*dist*(Ri + Rj) - 3.0*pow(Ri - Rj,2));
	         Vij /= Vj;

	         if (Vij > OverlapTol) Void[next].ToF = false;
    
	      } 
    	  }
      }
  }

  NumTrueVoid = 0;
  for (i=0; i<NumVoid; i++)
      if (Void[i].ToF) 
	 NumTrueVoid++;	  
  
  fprintf(logfile," | Number of true void after cleaning = %d (%4.2f %)\n",
		  NumTrueVoid,(double)NumTrueVoid/(double)NumVoid*100.0);

  free(SortArr);
  free_grid_list(GridList,NumGrid);
  free_query_grid(&Query);

  StepName.push_back("Cleaning voids"); 
  StepTime.push_back(get_time(t,1));

}

