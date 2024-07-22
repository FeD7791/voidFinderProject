
#include "allvars.h"
#include "tools.h"
#include "grid.h"

template <class T>
void build_grid_list(T Points, int NumPoints, struct grid *GridList, int NumGrid, double *GridSize, bool compute_neigh)
{

  fprintf(logfile," | Creating grid-list \n");

  for (int k=0; k<3; k++) 
      GridSize[k] = LBox[k]/(double)NumGrid;    

  fprintf(logfile," | Number of grids = %d \n",NumGrid);
  fprintf(logfile," | Grid sizes: x = %f [Mpc/h] \n",GridSize[0]);
  fprintf(logfile," |             y = %f [Mpc/h] \n",GridSize[1]);
  fprintf(logfile," |             z = %f [Mpc/h] \n",GridSize[2]);
  fflush(logfile);

  for (int l=0; l<NumGrid*NumGrid*NumGrid; l++)
      GridList[l].NumMem = 0;
   
  for (int p=0; p<NumPoints; p++) {

      int i = (int)(Points[p].Pos[0]/GridSize[0]);
      int j = (int)(Points[p].Pos[1]/GridSize[1]);
      int k = (int)(Points[p].Pos[2]/GridSize[2]);

      if (i == NumGrid) i--;
      if (j == NumGrid) j--;
      if (k == NumGrid) k--;

      int l = index_1d(i,j,k,NumGrid);

      GridList[l].NumMem++;
  }

  for (int l=0; l<NumGrid*NumGrid*NumGrid; l++) {
      GridList[l].Member = (int *) malloc(GridList[l].NumMem*sizeof(int));	  
      GridList[l].NumMem = 0;
  }

  for (int p=0; p<NumPoints; p++) {

      int i = (int)(Points[p].Pos[0]/GridSize[0]);
      int j = (int)(Points[p].Pos[1]/GridSize[1]);
      int k = (int)(Points[p].Pos[2]/GridSize[2]);

      if (i == NumGrid) i--;
      if (j == NumGrid) j--;
      if (k == NumGrid) k--;

      int l = index_1d(i,j,k,NumGrid);

      GridList[l].Member[GridList[l].NumMem] = p;
      GridList[l].NumMem++;

  }

  if (!compute_neigh) return;

  int nv;

  for (int l=0; l<NumGrid*NumGrid*NumGrid; l++) 
      GridList[l].Neighbour = (int *) malloc(27*sizeof(int));	  

  for (int i=0; i<NumGrid; i++) {
      for (int j=0; j<NumGrid; j++) {
          for (int k=0; k<NumGrid; k++) {

	      nv = 0;
              int l = index_1d(i,j,k,NumGrid);

	      for (int it=i-1; it<=i+1; it++) {
		  int ii = periodic_grid(it,NumGrid);

	          for (int jt=j-1; jt<=j+1; jt++) {
		      int jj = periodic_grid(jt,NumGrid);    

	              for (int kt=k-1; kt<=k+1; kt++) {
		          int kk = periodic_grid(kt,NumGrid);    
              
			  int ll = index_1d(ii,jj,kk,NumGrid);
			  GridList[l].Neighbour[nv] = ll;
		          nv++;

		      }
		  }
	      }	  
	  }
      }
  }

  return;
}

void query_grid(struct query *Query, double *GridSize, double R1, double R2)
{
  int    i,j,k,NG[3];
  int    i1,i2,j1,j2,k1,k2;
  double x1,x2,y1,y2,z1,z2;
  double MinDist,MaxDist,GAP,d1,d2;

  GAP = 0.5*sqrt(3.0)*max_grid_size(GridSize);  
  MinDist = R1 - GAP;
  MaxDist = R2 + GAP;
  if (MinDist < 0.0) MinDist = 0.0;

  for (k=0; k<3; k++)
      NG[k] = (int)ceil(MaxDist/GridSize[k]);  

  for (i=-NG[0]; i<=NG[0]; i++) {

      i1 = i2 = i; 
      i >= 0 ? i1++ : i2++;

      for (j=-NG[1]; j<=NG[1]; j++) {
     
          j1 = j2 = j;
          j >= 0 ? j1++ : j2++;

	  for (k=-NG[2]; k<=NG[2]; k++) {

	      k1 = k2 = k;
              k >= 0 ? k1++ : k2++;	      
     
	      x1 = ((double)i1 - 0.5)*GridSize[0];
	      x2 = ((double)i2 - 0.5)*GridSize[0];
	      y1 = ((double)j1 - 0.5)*GridSize[1];
	      y2 = ((double)j2 - 0.5)*GridSize[1];
	      z1 = ((double)k1 - 0.5)*GridSize[2];
	      z2 = ((double)k2 - 0.5)*GridSize[2];

	      d1 = sqrt(x1*x1 + y1*y1 + z1*z1);
	      d2 = sqrt(x2*x2 + y2*y2 + z2*z2);
	      
              if ((d1 > MinDist && d1 < MaxDist) || (d2 > MinDist && d2 < MaxDist)) {
                 (*Query).i.push_back(i);     
                 (*Query).j.push_back(j);     
                 (*Query).k.push_back(k);     
              }
          }
      }
  }  

}

void free_grid_list(struct grid *GridList, int NG)
{
   int i;

   for (i=0; i<NG*NG*NG; i++)  
       free(GridList[i].Member);
   free(GridList);
}

void free_query_grid(struct query *Query)
{
   (*Query).i.clear();     
   (*Query).j.clear();     
   (*Query).k.clear();
}

double max_grid_size(double *GridSize)
{
  double max = 0.0;
  int    k;

  for (k=0; k<3; k++) 
      if (GridSize[k] > max) 
	 max = GridSize[k];	
  
  return max;  
}

template void build_grid_list<vector <tracers>>(vector <tracers> Points, int NumPoints, struct grid *GridList, 
		                             int NumGrid, double *GridSize, bool compute_neigh);
template void build_grid_list<vector <voids>>(vector <voids> Points, int NumPoints, struct grid *GridList, 
		                            int NumGrid, double *GridSize, bool compute_neigh);
