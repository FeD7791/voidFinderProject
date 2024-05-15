
#include "allvars.h"
#include "grid.h"
#include "tools.h"
#include "profiles.h"

void compute_profiles()
{
   int            i,k,ic,jc,kc,l,ii,jj,kk,next,ibin,in,m,NumGrid;
   double         xc[3],xt[3],dx[3],vt[3],dist,GridSize[3];
   double         dR,MinDist,MaxDist,GAP;
   double         VRad,Vol,DeltaMax;
   float          Radius;
   vector <int>   Indx;
   struct profile Prof[NumProfileBins];  
   char           TxtFile[MAXCHAR],BinFile[MAXCHAR];
   int            NumQuery;
   struct query   Query;
   struct grid    *GridList;
   FILE           *ftxt,*fbin;
   clock_t        t;

   fprintf(logfile,"\n COMPUTING VOID PROFILES \n");
   t = clock();

   NumGrid = (int)(BoxSize/ProxyGridSize);
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   build_grid_list(Tracer,NumTrac,GridList,NumGrid,GridSize,false);

   // Only for true voids

   for (i=0; i<NumVoid; i++) 
       if (Void[i].ToF) 
	  Indx.push_back(i);       

   dR = (MaxProfileDist - MinProfileDist)/(double)NumProfileBins;
   
   // Selecciono grides

   MinDist = 0.0;
   MaxDist = 0.0;
   for (i=0; i<NumVoid; i++) {
       if (!Void[i].ToF) continue;
       if (Void[i].Rad > MaxDist) MaxDist = Void[i].Rad;
   }
   MaxDist *= MaxProfileDist;
   query_grid(&Query,GridSize,MinDist,MaxDist);
   NumQuery = Query.i.size();
   GAP = 0.5*sqrt(3.0)*max_grid_size(GridSize);
   
   fprintf(logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
   fflush(logfile);

   if (WriteProfiles == 2) {
      sprintf(BinFile,"%s/profiles.bin",PathProfiles);
      fbin = safe_open(BinFile,"w");
   }

   #pragma omp parallel for default(none) schedule(dynamic)                   \
    shared(NumVoid,Void,Tracer,NumQuery,Query,dR,NumGrid,GridSize,GridList,   \
           Indx,MeanNumTrac,LBox,NumProfileBins,MinProfileDist,MaxProfileDist,\
	   WriteProfiles,PathProfiles,GAP,fbin,BinFile)                       \
    private(i,m,k,ii,jj,kk,l,Radius,ic,jc,kc,xc,xt,dx,vt,next,Prof,dist,VRad, \
	    ibin,DeltaMax,Vol,ftxt,in,TxtFile)

   for (i=0; i<(int)Indx.size(); i++) {
       
       for (k=0; k<NumProfileBins; k++) {
           Prof[k].DeltaDiff = 0.0;
           Prof[k].DeltaCum = 0.0;
           Prof[k].Velocity = 0.0;
       }
       
       Radius = Void[Indx[i]].Rad;
       for (k=0; k<3; k++) 
	   xc[k] = Void[Indx[i]].Pos[k];
       Void[Indx[i]].Dtype = 0.0;
      
       ic = (int)(xc[0]/GridSize[0]);
       jc = (int)(xc[1]/GridSize[1]);
       kc = (int)(xc[2]/GridSize[2]);

       for (in=0; in<NumQuery; in++) {

       	   ii = Query.i[in]; 
	   jj = Query.j[in]; 
	   kk = Query.k[in]; 

	   dist = (double)(ii*ii)*(GridSize[0]*GridSize[0])
	        + (double)(jj*jj)*(GridSize[1]*GridSize[1])
	        + (double)(kk*kk)*(GridSize[2]*GridSize[2]);
	   dist = sqrt(dist);

	   if (dist > MaxProfileDist*Radius + 2.0*GAP) continue;

       	   ii = periodic_grid(ii + ic, NumGrid); 
	   jj = periodic_grid(jj + jc, NumGrid); 
	   kk = periodic_grid(kk + kc, NumGrid); 	       
	   
	   l = index_1d(ii,jj,kk,NumGrid);

           if (GridList[l].NumMem == 0) continue;

           for (m=0; m<GridList[l].NumMem; m++) {
		
	       next = GridList[l].Member[m];	

	       for (k=0; k<3; k++) {
	           xt[k] = (double)Tracer[next].Pos[k];	 
	           vt[k] = (double)(Tracer[next].Vel[k] - Void[Indx[i]].Vel[k]);
	           dx[k] = periodic_delta(xt[k] - xc[k],LBox[k])/Radius;
	       }

               dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	       if (dist > MinProfileDist && dist < MaxProfileDist) {
           
	          ibin = (int)((dist - MinProfileDist)/dR);

	          VRad = vt[0]*dx[0] + vt[1]*dx[1] + vt[2]*dx[2];
	          VRad /= dist;

	          Prof[ibin].Velocity += VRad;
	          Prof[ibin].DeltaDiff += 1.0; 

	       }
	   }
       }
       
       for (k=0; k<NumProfileBins; k++) {
	   if (Prof[k].DeltaDiff < 3.0) {
	       Prof[k].DeltaDiff = 0.0;
               Prof[k].Velocity = 0.0;
           } else {
               Prof[k].Velocity /= Prof[k].DeltaDiff;	   
	   }
           for (kk=0; kk<=k; kk++) 
	       Prof[k].DeltaCum += Prof[kk].DeltaDiff;
       }

       DeltaMax = -1.0;
       for (k=0; k<NumProfileBins; k++) {

	   Prof[k].Ri = (float)(k    )*dR + MinProfileDist;
	   Prof[k].Rm = (float)(k+0.5)*dR + MinProfileDist;
	   Prof[k].Rs = (float)(k+1.0)*dR + MinProfileDist;

	   Prof[k].Ri *= Radius;
	   Prof[k].Rm *= Radius;
	   Prof[k].Rs *= Radius;

	   Vol = (4.0/3.0)*PI*(pow(Prof[k].Rs,3) - pow(Prof[k].Ri,3));
	   Prof[k].DeltaDiff = Prof[k].DeltaDiff/Vol/MeanNumTrac - 1.0;

	   Vol = (4.0/3.0)*PI*pow(Prof[k].Rs,3);
	   Prof[k].DeltaCum = Prof[k].DeltaCum/Vol/MeanNumTrac - 1.0;

	   if (Prof[k].Rs < 2.0*Radius || Prof[k].Rs > 3.0*Radius) continue;

	   if (Prof[k].DeltaCum > DeltaMax) DeltaMax = Prof[k].DeltaCum;

       }

       Void[Indx[i]].Dtype = DeltaMax;
       
       if (WriteProfiles == 1) {
          sprintf(TxtFile,"%s/profile_void_%d.dat",PathProfiles,i);
          ftxt = safe_open(TxtFile,"w");
          for (k=0; k<NumProfileBins; k++)
	      fprintf(ftxt,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",
			   Prof[k].Ri,Prof[k].Rm,Prof[k].Rs,Prof[k].DeltaDiff,
			   Prof[k].DeltaCum,Prof[k].Velocity,Radius);  
          fclose(ftxt);
       }
       
       if (WriteProfiles == 2) {
          #pragma omp critical
	  {
	     fwrite(&i,sizeof(int),1,fbin);
             fwrite(&Radius,sizeof(float),1,fbin);
             fwrite(&NumProfileBins,sizeof(int),1,fbin);
	     for (k=0; k<NumProfileBins; k++) {
	         fwrite(&Prof[k].Ri,sizeof(float),1,fbin);
	         fwrite(&Prof[k].Rm,sizeof(float),1,fbin);
	         fwrite(&Prof[k].Rs,sizeof(float),1,fbin);
	         fwrite(&Prof[k].DeltaDiff,sizeof(float),1,fbin);
	         fwrite(&Prof[k].DeltaCum,sizeof(float),1,fbin);
	         fwrite(&Prof[k].Velocity,sizeof(float),1,fbin);
	     }
	  }
       }

       Void[Indx[i]].Dtype = DeltaMax;
   } 
   
   if (WriteProfiles == 2) fclose(fbin);

   Indx.clear();
   free_query_grid(&Query);
   free_grid_list(GridList,NumGrid);
   
   StepName.push_back("Computing profiles");
   StepTime.push_back(get_time(t,OMPcores));

}

void bin2ascii_profile(int voidID)
{
   int            SkipBlock,iv,k;
   float          Radius;  
   struct profile Prof[NumProfileBins];
   FILE           *fbin,*fout;
   char           filename[MAXCHAR];
   
   SkipBlock = sizeof(float) + sizeof(int) + 6*NumProfileBins*sizeof(float);
   sprintf(filename,"%s/profiles.bin",PathProfiles);  
   fbin = safe_open(filename,"r");

   do {
   
      fread(&iv,sizeof(int),1,fbin);

      if (iv+1 != voidID)  
         fseek(fbin,SkipBlock,SEEK_CUR);
      else {
         fread(&Radius,sizeof(float),1,fbin);
         fread(&NumProfileBins,sizeof(int),1,fbin);
	 for (k=0; k<NumProfileBins; k++) {
	     fread(&Prof[k].Ri,sizeof(float),1,fbin);
	     fread(&Prof[k].Rm,sizeof(float),1,fbin);
	     fread(&Prof[k].Rs,sizeof(float),1,fbin);
	     fread(&Prof[k].DeltaDiff,sizeof(float),1,fbin);
	     fread(&Prof[k].DeltaCum,sizeof(float),1,fbin);
	     fread(&Prof[k].Velocity,sizeof(float),1,fbin);
	 }
      }
   } while (iv+1 != voidID);
   fclose(fbin);

   fprintf(logfile,"\n Writting profile for void %d, with radius %5.3f [Mpc/h] \n",voidID,Radius);

   sprintf(filename,"%s/profile_void_%d.dat",PathProfiles,voidID);
   fout = safe_open(filename,"w");

   for (k=0; k<NumProfileBins; k++)
       fprintf(fout,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",
		     Prof[k].Ri,Prof[k].Rm,Prof[k].Rs,Prof[k].DeltaDiff,
		     Prof[k].DeltaCum,Prof[k].Velocity,Radius);  
   fclose(fout);

}
