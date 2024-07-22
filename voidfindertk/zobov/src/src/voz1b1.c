#include "voz.h"

int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs);
int vorvol (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints, realT *vol);
int posread(char *posfile, realT ***p, realT fact);

int main(int argc, char *argv[]) {
  int exitcode;
  int i, j, np;
  realT **r;
  coordT rtemp[3], *parts;
  coordT deladjs[3*MAXVERVER], points[3*MAXVERVER];
  pointT intpoints[3*MAXVERVER];
  FILE *pos, *out;
  char *posfile, outfile[80], *suffix;
  PARTADJ *adjs;
  realT *vols;
  realT predict, xmin,xmax,ymin,ymax,zmin,zmax;
  int *orig;
  
  int isitinbuf;
  char isitinmain, d;
  int numdiv, nvp, nvpall, nvpbuf;
  realT width, width2, totwidth, totwidth2, bf, s, g;
  realT border, boxsize;
  realT c[3];
  int b[3];
  realT totalvol;

  if (argc != 9) {
    printf("Wrong number of arguments.\n");
    printf("arg1: position file\n");
    printf("arg2: border size\n");
    printf("arg3: box size\n");
    printf("arg4: suffix\n");
    printf("arg5: number of divisions\n");
    printf("arg6-8: b[0-2]\n\n");
    exit(0);
  }
  posfile = argv[1];
  if (sscanf(argv[2],"%"vozRealSym,&border) != 1) {
    printf("That's no border size; try again.\n");
    exit(0);
  }
  if (sscanf(argv[3],"%"vozRealSym,&boxsize) != 1) {
    printf("That's no boxsize; try again.\n");
    exit(0);
  }
  suffix = argv[4];
  if (sscanf(argv[5],"%d",&numdiv) != 1) {
    printf("%s is no number of divisions; try again.\n",argv[5]);
    exit(0);
  }
  if (numdiv == 1) {
    printf("Only using one division; should only use for an isolated segment.\n");
  }
  if (numdiv < 1) {
    printf("Cannot have a number of divisions less than 1.  Resetting to 1.\n");
    numdiv = 1;
  }
  if (sscanf(argv[6],"%d",&b[0]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  if (sscanf(argv[7],"%d",&b[1]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  if (sscanf(argv[8],"%d",&b[2]) != 1) {
    printf("That's no b index; try again.\n");
    exit(0);
  }
  
  /* Boxsize should be the range in r, yielding a range 0-1 */
  np = posread(posfile,&r,1./boxsize);
  printf("%d particles\n",np);fflush(stdout);
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
    if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
    if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];
  }
  printf("np: %d, x: %g,%g; y: %g,%g; z: %g,%g\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  width = 1./(realT)numdiv;
  width2 = 0.5*width;
  if (border > 0.) bf = border;
  else bf = 0.1;
      /* In units of 0-1, the thickness of each subregion's buffer*/
  totwidth = width+2.*bf;
  totwidth2 = width2 + bf;
  
  s = width/(realT)NGUARD;
  if ((bf*bf - 2.*s*s) < 0.) {
    printf("bf = %g, s = %g.\n",bf,s);
    printf("Not enough guard points for given border.\nIncrease guards to >= %g\n.",
	   sqrt(2.)*width/bf);
    exit(0);
  }
  g = (bf / 2.)*(1. + sqrt(1 - 2.*s*s/(bf*bf)));
  printf("s = %g, bf = %g, g = %g.\n",s,bf,g);
  
  fflush(stdout);

  adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
  if (adjs == NULL) {
    printf("Unable to allocate adjs\n");
    exit(0);
  }
  
  DL c[d] = ((realT)b[d]+0.5)*width;
  printf("c: %g,%g,%g\n",c[0],c[1],c[2]);
  /* Assign temporary array*/
  nvpbuf = 0; /* Number of particles to tesselate, including
		 buffer */
  nvp = 0; /* Without the buffer */
  for (i=0; i<np; i++) {
    isitinbuf = 1;
    isitinmain = 1;
    DL {
      rtemp[d] = (realT)r[i][d] - (realT)c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinbuf = isitinbuf && (fabs(rtemp[d]) < totwidth2);
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
  
    if (isitinbuf) nvpbuf++;
    if (isitinmain) nvp++;
  }
  
  nvpbuf += 6*(NGUARD+1)*(NGUARD+1); /* number of guard
					points */

  parts = (coordT *)malloc(3*nvpbuf*sizeof(coordT));
  orig = (int *)malloc(nvpbuf*sizeof(int));

  if (parts == NULL) {
    printf("Unable to allocate parts\n");
    fflush(stdout);
  }
  if (orig == NULL) {
    printf("Unable to allocate orig\n");
    fflush(stdout);
  }

  nvp = 0; nvpall = 0; /* nvp = number of particles without buffer */
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np; i++) {
    isitinmain = 1;
    DL {
      rtemp[d] = r[i][d] - c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinmain = isitinmain && (fabs(rtemp[d]) <= width2);
    }
    if (isitinmain) {
      parts[3*nvp] = rtemp[0];
      parts[3*nvp+1] = rtemp[1];
      parts[3*nvp+2] = rtemp[2];
      orig[nvp] = i;
      nvp++;
      if (rtemp[0] < xmin) xmin = rtemp[0];
      if (rtemp[0] > xmax) xmax = rtemp[0];
      if (rtemp[1] < ymin) ymin = rtemp[1];
      if (rtemp[1] > ymax) ymax = rtemp[1];
      if (rtemp[2] < zmin) zmin = rtemp[2];
      if (rtemp[2] > zmax) zmax = rtemp[2];
    }
  }
  printf("nvp = %d\n",nvp);
  printf("x: %g,%g; y: %g,%g; z:%g,%g\n",xmin,xmax,ymin,ymax,zmin,zmax);
  nvpbuf = nvp;
  for (i=0; i<np; i++) {
    isitinbuf = 1;
    DL {
      rtemp[d] = r[i][d] - c[d];
      if (rtemp[d] > 0.5) rtemp[d] --;
      if (rtemp[d] < -0.5) rtemp[d] ++;
      isitinbuf = isitinbuf && (fabs(rtemp[d])<totwidth2);
    }
    if ((isitinbuf > 0) &&
	((fabs(rtemp[0])>width2)||(fabs(rtemp[1])>width2)||(fabs(rtemp[2])>width2))) {
      
      /*printf("%3.3f ",sqrt(rtemp[0]*rtemp[0] + rtemp[1]*rtemp[1] +
	rtemp[2]*rtemp[2]));
	printf("|%2.2f,%2.2f,%2.2f,%g,%g",r[i][0],r[i][1],r[i][2],width2,totwidth2);*/
      parts[3*nvpbuf] = rtemp[0];
      parts[3*nvpbuf+1] = rtemp[1];
      parts[3*nvpbuf+2] = rtemp[2];
      orig[nvpbuf] = i;

      nvpbuf++;
      if (rtemp[0] < xmin) xmin = rtemp[0];
      if (rtemp[0] > xmax) xmax = rtemp[0];
      if (rtemp[1] < ymin) ymin = rtemp[1];
      if (rtemp[1] > ymax) ymax = rtemp[1];
      if (rtemp[2] < zmin) zmin = rtemp[2];
      if (rtemp[2] > zmax) zmax = rtemp[2];
    }
  }
  printf("nvpbuf = %d\n",nvpbuf);
  printf("x: %g,%g; y: %g,%g; z:%g,%g\n",xmin,xmax,ymin,ymax,zmin,zmax);
  nvpall = nvpbuf;
  predict = pow(width+2.*bf,3)*(realT)np;
  printf("There should be ~ %g points; there are %d\n",predict,nvpbuf);

  for (i=0;i<np;i++) free(r[i]);
  free(r);
  
  /* Add guard points */
  for (i=0; i<NGUARD+1; i++) {
    for (j=0; j<NGUARD+1; j++) {
      /* Bottom */
      parts[3*nvpall]   = -width2 + (realT)i * s;
      parts[3*nvpall+1] = -width2 + (realT)j * s;
      parts[3*nvpall+2] = -width2 - g;
      nvpall++;
      /* Top */
      parts[3*nvpall]   = -width2 + (realT)i * s;
      parts[3*nvpall+1] = -width2 + (realT)j * s;
      parts[3*nvpall+2] = width2 + g;
      nvpall++;
    }
  }
  for (i=0; i<NGUARD+1; i++) { /* Don't want to overdo the corners*/
    for (j=0; j<NGUARD+1; j++) {
      parts[3*nvpall]   = -width2 + (realT)i * s;
      parts[3*nvpall+1] = -width2 - g;
      parts[3*nvpall+2] = -width2 + (realT)j * s;
      nvpall++;
      
      parts[3*nvpall]   = -width2 + (realT)i * s;
      parts[3*nvpall+1] = width2 + g;
      parts[3*nvpall+2] = -width2 + (realT)j * s;
      nvpall++;
    }
  }
  for (i=0; i<NGUARD+1; i++) {
    for (j=0; j<NGUARD+1; j++) {
      parts[3*nvpall]   = -width2 - g;
      parts[3*nvpall+1] = -width2 + (realT)i * s;
      parts[3*nvpall+2] = -width2 + (realT)j * s;
      nvpall++;
      
      parts[3*nvpall]   = width2 + g;
      parts[3*nvpall+1] = -width2 + (realT)i * s;
      parts[3*nvpall+2] = -width2 + (realT)j * s;
      nvpall++;
    }
  }
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=nvpbuf;i<nvpall;i++) {
    if (parts[3*i] < xmin) xmin = parts[3*i];
    if (parts[3*i] > xmax) xmax = parts[3*i];
    if (parts[3*i+1] < ymin) ymin = parts[3*i+1];
    if (parts[3*i+1] > ymax) ymax = parts[3*i+1];
    if (parts[3*i+2] < zmin) zmin = parts[3*i+2];
    if (parts[3*i+2] > zmax) zmax = parts[3*i+2];
  }
  
  printf("Added guard points to total %d points (should be %d)\n",nvpall,
	 nvpbuf + 6*(NGUARD+1)*(NGUARD+1));
  printf("x: %g,%g; y: %g,%g; z:%g,%g\n",xmin,xmax,ymin,ymax,zmin,zmax);
  
  /* Do tesselation*/
  printf("File read.  Tessellating ...\n"); fflush(stdout);
  exitcode = delaunadj(parts, nvp, nvpbuf, nvpall, &adjs);
  
  /* Now calculate volumes*/
  printf("Now finding volumes ...\n"); fflush(stdout);
  vols = (realT *)malloc(nvp*sizeof(realT));
  
  for (i=0; i<nvp; i++) { /* Just the original particles
			     Assign adjacency coordinate array*/
    /* Volumes */
    for (j = 0; j < adjs[i].nadj; j++)
      DL {
	deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d];
	if (deladjs[3*j+d] < -0.5) deladjs[3*j+d]++;
	if (deladjs[3*j+d] > 0.5) deladjs[3*j+d]--;
      }
    
    exitcode = vorvol(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
    vols[i] *= (realT)np;
    if (i % 1000 == 0)
      printf("%d: %d, %g\n",i,adjs[i].nadj,vols[i]);
  }

  /* Get the adjacencies back to their original values */

  for (i=0; i<nvp; i++)
    for (j = 0; j < adjs[i].nadj; j++)
      adjs[i].adj[j] = orig[adjs[i].adj[j]];
  
  totalvol = 0.;
  for (i=0;i<nvp; i++) {
    totalvol += (realT)vols[i];
  }
  printf("Average volume = %g\n",totalvol/(realT)nvp);
  
  /* Now the output!
     First number of particles */
  sprintf(outfile,"part.%s.%02d.%02d.%02d",suffix,b[0],b[1],b[2]);

  printf("Output to %s\n\n",outfile);
  out = fopen(outfile,"w");
  if (out == NULL) {
    printf("Unable to open %s\n",outfile);
    exit(0);
  }
  fwrite(&np,1, sizeof(int),out);
  fwrite(&nvp,1, sizeof(int),out);
  printf("nvp = %d\n",nvp);

  /* Tell us where the original particles were */
  fwrite(orig,sizeof(int),nvp,out);
  /* Volumes*/
  fwrite(vols,sizeof(realT),nvp,out);
  /* Adjacencies */
  for (i=0;i<nvp;i++) {
    fwrite(&(adjs[i].nadj),1,sizeof(int),out);
    if (adjs[i].nadj > 0)
      fwrite(adjs[i].adj,adjs[i].nadj,sizeof(int),out);
    else printf("0");
  }
  fclose(out);
  
  return(0);
}
