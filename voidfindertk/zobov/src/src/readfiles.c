#include <stdio.h>
#include <stdlib.h>
#include "voz.h"

/* Positions */
/* Returns number of particles read */
int posread(char *posfile, realT ***p, realT fact) {

  FILE *pos;
  int np,dum,d,i;
  realT xmin,xmax,ymin,ymax,zmin,zmax;
  realT *ptemp;

  pos = fopen(posfile, "rb");
  if (pos == NULL) {
    printf("Unable to open position file %s\n\n",posfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
   fread(&np,1, sizeof(int),pos); 

  /* Allocate the arrays */
  (*p) = (realT **)malloc(np*sizeof(realT *));
  ptemp = (realT *)malloc(np*sizeof(realT));

  printf("np = %d\n",np);

  /* Fill the arrays */

  for (i=0; i<np; i++) {
    (*p)[i] = (realT *)malloc(3*sizeof(realT));
    if ((*p)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
  }

  fread(ptemp,np,sizeof(realT),pos);   
  for (i=0; i<np; i++) (*p)[i][0] = ptemp[i];   

  fread(ptemp,np,sizeof(realT),pos);
  for (i=0; i<np; i++) (*p)[i][1] = ptemp[i];

  fread(ptemp,np,sizeof(realT),pos);
  for (i=0; i<np; i++) (*p)[i][2] = ptemp[i];
   

  fclose(pos);
  free(ptemp);

  /* Get into physical units (Mpc/h) */
  
  for (i=0; i<np; i++) DL (*p)[i][d] *= fact;


  /* Test range -- can comment out */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
    if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
    if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
  }
  printf("np: %d, x: %g,%g; y: %g,%g; z: %g,%g\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  return(np);
}

/* Velocities */
/* Returns number of particles read */
int velread(char *velfile, realT ***v, realT fact) {

  FILE *vel;
  int np,dum,d,i;
  realT xmin,xmax,ymin,ymax,zmin,zmax;

  vel = fopen(velfile, "rb");
  if (vel == NULL) {
    printf("Unable to open velocity file %s\n\n",velfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
   fread(&np,1, sizeof(int),vel); 

  /* Allocate the arrays */
  (*v) = (realT **)malloc(3*sizeof(realT*));
  for (i=0;i<3;i++) (*v)[i] = (realT *)malloc(np*sizeof(realT));

  /* Fill the arrays */
  fread((*v)[0],np,sizeof(realT),vel); 
  fread((*v)[1],np,sizeof(realT),vel); 
  fread((*v)[2],np,sizeof(realT),vel); 

  fclose(vel);

  /* Convert from code units into physical units (km/sec) */
  
  for (i=0; i<np; i++) DL (*v)[d][i] *= fact;

  /* Test range -- can comment out */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*v)[0][i] < xmin) xmin = (*v)[0][i]; if ((*v)[0][i] > xmax) xmax = (*v)[0][i];
    if ((*v)[1][i] < ymin) ymin = (*v)[1][i]; if ((*v)[1][i] > ymax) ymax = (*v)[1][i];
    if ((*v)[2][i] < zmin) zmin = (*v)[2][i]; if ((*v)[2][i] > zmax) zmax = (*v)[2][i];
  }
  printf("vx: %g,%g; vy: %g,%g; vz: %g,%g\n",xmin,xmax, ymin,ymax, zmin,zmax);
  
  return(np);
}
