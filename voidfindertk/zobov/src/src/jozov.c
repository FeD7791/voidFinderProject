/* jovoz.c by Mark Neyrinck */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "voz.h"

#define NLINKS 300 /* Number of possible links with the same rho_sl */
#define FF fflush(stdout)

typedef struct Particle {
  realT dens;
  int nadj;
  int *adj;
} PARTICLE;

typedef struct Zone {
  int core; /* Identity of peak particle */
  int np; /* Number of particles in zone */
  int npjoin; /* Number of particles in the joined void */
  int nadj; /* Number of adjacent zones */
  int nhl; /* Number of zones in final joined void */
  realT leak; /* Volume of last leak zone*/
  int *adj; /* Each adjacent zone, with ... */
  realT *slv; /* Smallest Linking Volume */
  realT denscontrast; /* density contrast */
  double vol; /* Total volume of all particles in the zone */
  double voljoin; /* Total volume of all particles in the joined void */
} ZONE;

typedef struct ZoneT {
  int nadj; /* Number of zones on border */
  int *adj; /* Each adjacent zone, with ... */
  realT *slv; /* Smallest Linking Volume */
} ZONET;

void findrtop(double *a, int na, int *iord, int nb);

int main(int argc,char **argv) {

  FILE *adj, *vol, *zon, *zon2, *txt;
  PARTICLE *p;
  ZONE *z;
  ZONET *zt;
  char *adjfile, *volfile, *zonfile, *zonfile2, *txtfile;
  int i, j,k,l, h, h2,hl,n,np, np2, nzones, nhl, nhlcount, nhl2;
  int *jumper, *jumped, *numinh;
  int *zonenum, *zonelist, *zonelist2;
  int link[NLINKS], link2, nl;
  realT lowvol, voltol, prob;

  int q,q2;

  int za, nin;
  int testpart;
  char already, interior, *inyet, *inyet2, added, beaten;
  int *nm, **m;

  realT maxvol, minvol;
  double *sorter, e1,maxdenscontrast;
  int *iord;

  e1 = exp(1.)-1.;

  if (argc != 7) {
    printf("Wrong number of arguments.\n");
    printf("arg1: adjacency file\n");
    printf("arg2: volume file\n");
    printf("arg3: output file containing particles in each zone\n");
    printf("arg4: output file containing zones in each void\n");
    printf("arg5: output text file\n");
    printf("arg6: Density threshold (0 for no threshold)\n\n");
    exit(0);
  }
  adjfile = argv[1];
  volfile = argv[2];
  zonfile = argv[3];
  zonfile2 = argv[4];
  txtfile = argv[5];
  if (sscanf(argv[6],"%"vozRealSym,&voltol) == 0) {
    printf("Bad density threshold.\n");
    exit(0);
  }
  if (voltol <= 0.) {
    printf("Proceeding without a density threshold.\n");
    voltol = 1e30;
  }

  adj = fopen(adjfile, "r");
  if (adj == NULL) {
    printf("Unable to open %s\n",adjfile);
    exit(0);
  }
  fread(&np,1, sizeof(int),adj);
  
  printf("adj: %d particles\n", np);
  FF;

  p = (PARTICLE *)malloc(np*sizeof(PARTICLE));
  /* Adjacencies*/
  for (i=0;i<np;i++) {
    fread(&p[i].nadj,1,sizeof(int),adj); 
    /* The number of adjacencies per particle */
    if (p[i].nadj > 0)
      p[i].adj = (int *)malloc(p[i].nadj*sizeof(int));
    p[i].nadj = 0; /* Temporarily, it's an adj counter */
  }
  for (i=0;i<np;i++) {
    fread(&nin,1,sizeof(int),adj);
    if (nin > 0)
      for (k=0;k<nin;k++) {
	fread(&j,1,sizeof(int),adj);
	if (j < np) {
	  /* Set both halves of the pair */
	  p[i].adj[p[i].nadj] = j;
	  p[j].adj[p[j].nadj] = i;
	  p[i].nadj++; p[j].nadj++;
	} else {
	  printf("%d: adj = %d\n",i,j);
	}
      }
  }
  fclose(adj);

  /* Check that we got all the pairs */
  adj = fopen(adjfile, "r");
  fread(&np,1, sizeof(int),adj);
  for (i=0;i<np;i++) {
    fread(&nin,1,sizeof(int),adj); /* actually nadj */
    if (nin != p[i].nadj) {
      printf("We didn't get all of %d's adj's; %d != %d.\n",i,nin,p[i].nadj);
      /*exit(0);*/
    }
  }
  fclose(adj);

  /* Volumes */
  vol = fopen(volfile, "r");
  if (vol == NULL) {
    printf("Unable to open volume file %s.\n\n",volfile);
    exit(0);
  }
  fread(&np2,1, sizeof(int),adj);
  if (np != np2) {
    printf("Number of particles doesn't match! %d != %d\n",np,np2);
    exit(0);
  }
  printf("%d particles\n", np);
  FF;
  for (i=0;i<np;i++) {
    fread(&p[i].dens,1,sizeof(realT),vol);
    if ((p[i].dens < 1e-30) || (p[i].dens > 1e30)) {
      printf("Whacked-out volume found, of particle %d: %g\n",i,p[i].dens);
      p[i].dens = 1.;
    }
    p[i].dens = 1./p[i].dens; /* Get density from volume */
  }
  fclose(vol);

  jumped = (int *)malloc(np*sizeof(int));  
  jumper = (int *)malloc(np*sizeof(int));
  numinh = (int *)malloc(np*sizeof(int));

  /* find jumper */
  for (i = 0; i < np; i++) {
    minvol = p[i].dens; jumper[i] = -1;
    for (j=0; j< p[i].nadj; j++) {
      if (p[p[i].adj[j]].dens < minvol) {
	jumper[i] = p[i].adj[j];
	minvol = p[jumper[i]].dens;
      }
    }
    numinh[i] = 0;
  }

  printf("About to jump ...\n"); FF;
  
  /* Jump */
  for (i = 0; i < np; i++) {
    jumped[i] = i;
    while (jumper[jumped[i]] > -1)
      jumped[i] = jumper[jumped[i]];
    numinh[jumped[i]]++;
  }
  printf("Post-jump ...\n"); FF;
  
  nzones = 0;
  for (i = 0; i < np; i++)
    if (numinh[i] > 0) nzones++;
  printf("%d initial zones found\n", nzones);

  z = (ZONE *)malloc(nzones*sizeof(ZONE));
  if (z == NULL) {
    printf("Unable to allocate z\n");
    exit(0);
  }
  zt = (ZONET *)malloc(nzones*sizeof(ZONET));
  if (zt == NULL) {
    printf("Unable to allocate zt\n");
    exit(0);
  }
  
  for (h=0;h<nzones;h++) {
    zt[h].nadj = 0;
  }
  zonenum = (int *)malloc(np*sizeof(int));
  if (zonenum == NULL) {
    printf("Unable to allocate zonenum\n");
    exit(0);
  }

  h = 0;
  for (i = 0; i < np; i++)
    if (numinh[i] > 0) {
      z[h].core = i;
      zonenum[i] = h;
      h++;
    } else {
      zonenum[i] = -1;
    }
 
  /* Count border particles */
  for (i = 0; i < np; i++)
    for (j = 0; j < p[i].nadj; j++) {
      testpart = p[i].adj[j];
      if (jumped[i] != jumped[testpart])
	zt[zonenum[jumped[i]]].nadj++;
    }
  
  for (h=0;h<nzones;h++) {
    zt[h].adj = (int *)malloc(zt[h].nadj*sizeof(int));
    if (zt[h].adj == NULL) {
      printf("Unable to allocate %d adj's of zone %d\n",zt[h].nadj,h);
      exit(0);
    }
    zt[h].slv = (realT *)malloc(zt[h].nadj*sizeof(realT));
    if (zt[h].slv == NULL) {
      printf("Unable to allocate %d slv's of zone %d\n",zt[h].nadj,h);
      exit(0);
    }
    zt[h].nadj = 0;
  }

  /* Find "weakest links" */
  for (i = 0; i < np; i++) {
    h = zonenum[jumped[i]];
    for (j = 0; j < p[i].nadj; j++) {
      testpart = p[i].adj[j];
      if (h != zonenum[jumped[testpart]]) {
	if (p[testpart].dens > p[i].dens) {
	  /* there could be a weakest link through testpart */
	  already = 0;
	  for (za = 0; za < zt[h].nadj; za++)
	    if (zt[h].adj[za] == zonenum[jumped[testpart]]) {
	      already = 1;
	      if (p[testpart].dens < zt[h].slv[za]) {
		zt[h].slv[za] = p[testpart].dens;
	      }
	    }
	  if (already == 0) {
	    zt[h].adj[zt[h].nadj] = zonenum[jumped[testpart]];
	    zt[h].slv[zt[h].nadj] = p[testpart].dens;
	    zt[h].nadj++;
	  }
	} else { /* There could be a weakest link through i */
	  already = 0;
	  for (za = 0; za < zt[h].nadj; za++)
	    if (zt[h].adj[za] == zonenum[jumped[testpart]]) {
	      already = 1;
	      if (p[i].dens < zt[h].slv[za]) {
		zt[h].slv[za] = p[i].dens;
	      }
	    }
	  if (already == 0) {
	    zt[h].adj[zt[h].nadj] = zonenum[jumped[testpart]];
	    zt[h].slv[zt[h].nadj] = p[i].dens;
	    zt[h].nadj++;
	  }
	}
      }
    }
  }
  printf("Found zone adjacencies\n"); FF;

  /* Free particle adjacencies */
  for (i=0;i<np; i++) free(p[i].adj);

  /* Use z instead of zt */
  for (h=0;h<nzones;h++) {
    /*printf("%d ",zt[h].nadj);*/
    z[h].nadj = zt[h].nadj;
    z[h].adj = (int *)malloc(zt[h].nadj*sizeof(int));
    z[h].slv = (realT *)malloc(zt[h].nadj*sizeof(realT));
    for (za = 0; za<zt[h].nadj; za++) {
      z[h].adj[za] = zt[h].adj[za];
      z[h].slv[za] = zt[h].slv[za];
    }
    free(zt[h].adj);
    free(zt[h].slv);
    z[h].np = numinh[z[h].core];
  }
  free(zt);
  free(numinh);

  m = (int **)malloc(nzones*sizeof(int *));
  /* Not in the zone struct since it'll be freed up (contiguously, we hope)
     soon */
  nm = (int *)malloc(nzones*sizeof(int));
  for (h=0; h<nzones; h++) {
    m[h] = (int *)malloc(z[h].np*sizeof(int));
    nm[h] = 0;
    z[h].vol = 0.;
  }
  for (i=0; i<np; i++) {
    h = zonenum[jumped[i]];
    if (i == z[h].core) {
      m[h][nm[h]] = m[h][0];
      m[h][0] = i; /* Put the core particle at the top of the list */
    } else {
      m[h][nm[h]] = i;
    }
    z[h].vol += 1.0/(double)p[i].dens;
    nm[h] ++;
  }
  free(nm);


  zon = fopen(zonfile,"w");
  if (zon == NULL) {
    printf("Problem opening zonefile %s.\n\n",zonfile);
    exit(0);
  }
  fwrite(&np,1,sizeof(int),zon);
  fwrite(&nzones,1,sizeof(int),zon);
  for (h=0; h<nzones; h++) {
    fwrite(&(z[h].np),1,sizeof(int),zon);
    fwrite(m[h],z[h].np,sizeof(int),zon);
    free(m[h]);
  }
  free(m);
  close(zon);

  inyet = (char *)malloc(nzones*sizeof(char));
  inyet2 = (char *)malloc(nzones*sizeof(char));
  zonelist = (int *)malloc(nzones*sizeof(int));
  zonelist2 = (int *)malloc(nzones*sizeof(int));
  sorter = (double *)malloc((nzones+1)*sizeof(double));

  for (h = 0; h< nzones; h++) {
    inyet[h] = 0;
    inyet2[h] = 0;
  }

  nhl = 0; 

  maxvol = 0.;
  minvol = BF;
  maxdenscontrast = 0.;
  for(i=0;i<np; i++){
    if (p[i].dens > maxvol) maxvol = p[i].dens;
    if (p[i].dens < minvol) minvol = p[i].dens;
  }
  printf("Densities range from %e to %e.\n",minvol,maxvol);FF;

  zon2 = fopen(zonfile2,"w");
  if (zon2 == NULL) {
    printf("Problem opening zonefile %s.\n\n",zonfile2);
    exit(0);
  }
  fwrite(&nzones,1,sizeof(int),zon2);

  for (h = 0; h<nzones; h++) {
    nhlcount = 0;
    for (hl = 0; hl < nhl; hl++)
      inyet[zonelist[hl]] = 0;

    zonelist[0] = h;
    inyet[h] = 1;
    nhl = 1;
    z[h].npjoin = z[h].np;
    do {
      /* Find the lowest-volume (highest-density) adjacency */
      lowvol = BF; nl = 0; beaten = 0;
      for (hl = 0; hl < nhl; hl++) {
	h2 = zonelist[hl];
	if (inyet[h2] == 1) { /* If it's not already identified as 
				 an interior zone, with inyet=2 */
	  interior = 1;
	  for (za = 0; za < z[h2].nadj; za ++) {
	    if (inyet[z[h2].adj[za]] == 0) {
	      interior = 0;
	      if (z[h2].slv[za] == lowvol) {
		link[nl] = z[h2].adj[za];
		nl ++;
		if (nl == NLINKS) {
		  printf("Too many links with the same rho_sl!  Increase NLINKS from %d\n",nl);
		  exit(0);
		}
	      }
	      if (z[h2].slv[za] < lowvol) {
		lowvol = z[h2].slv[za];
		link[0] = z[h2].adj[za];
		nl = 1;
	      }
	    }
	  }
	  if (interior == 1) inyet[h2] = 2; /* No bordering exter. zones */
	}
      }

      if (nl == 0) {
	beaten = 1;
	z[h].leak = maxvol;
	continue;
      }
	
      if (lowvol > voltol) {
	beaten = 1;
	z[h].leak = lowvol;
	continue;
      }

      for (l=0; l < nl; l++)
	if (p[z[link[l]].core].dens < p[z[h].core].dens)
	  beaten = 1;
      if (beaten == 1) {
	z[h].leak = lowvol;
	continue;
      }
      /* Add everything linked to the link(s) */
      nhl2 = 0;
      for (l=0; l < nl; l++) {
	if (inyet2[link[l]] == 0) {
	  zonelist2[nhl2] = link[l];
	  inyet2[link[l]] = 1;
	  nhl2 ++;
	  added = 1;
	  while ((added == 1) && (beaten == 0)) {
	    added = 0;
	    for (hl = 0; (hl < nhl2) && (beaten == 0); hl++) {
	      h2 = zonelist2[hl];
	      if (inyet2[h2] == 1) {
		interior = 1; /* Guilty until proven innocent */
		for (za = 0; za < z[h2].nadj; za ++) {
		  link2 = z[h2].adj[za];
		  if ((inyet[link2]+inyet2[link2]) == 0) {
		    interior = 0;
		    if (z[h2].slv[za] <= lowvol) {
		      if (p[z[link2].core].dens < p[z[h].core].dens) {
			beaten = 1;
			break;
		      }
		      zonelist2[nhl2] = link2;
		      inyet2[link2] = 1;
		      nhl2++;
		      added = 1;
		    }
		  }
		}
		if (interior == 1) inyet2[h2] = 2;
	      }
	    }
	  }
	}
      }
      for (hl = 0; hl < nhl2; hl++)
	inyet2[zonelist2[hl]] = 0;
      
      /* See if there's a beater */
      if (beaten == 1) {
	z[h].leak = lowvol;
      } else {
	for (h2 = 0; h2 < nhl2; h2++) {
	  zonelist[nhl] = zonelist2[h2];
	  inyet[zonelist2[h2]] = 1;
	  nhl++;
	  z[h].npjoin += z[zonelist2[h2]].np;
	}
      }
      if (nhl/10000 > nhlcount) {
	nhlcount = nhl/10000;
	printf(" %d",nhl); FF;
      }
    } while((lowvol < BF) && (beaten == 0));
    
    z[h].denscontrast = z[h].leak/p[z[h].core].dens;
    if (z[h].denscontrast < 1.) z[h].denscontrast = 1.;
    
    /* find biggest denscontrast */
    if (z[h].denscontrast > maxdenscontrast) {
      maxdenscontrast = (double)z[h].denscontrast;
    }

    /* Don't sort; want the core zone to be first */

    if (nhlcount > 0) { /* Outputs the number of zones in large voids */
      printf(" h%d:%d\n",h,nhl);
      FF;
    }
    /* Calculate volume */
    z[h].voljoin = 0.;
    for (q = 0; q<nhl; q++) {
      z[h].voljoin += z[zonelist[q]].vol;
    }

    z[h].nhl = nhl;

    fwrite(&nhl,1,sizeof(int),zon2);
    fwrite(zonelist,nhl,sizeof(int),zon2);
  }
  fclose(zon2);

  printf("Maxdenscontrast = %g.\n",maxdenscontrast);

  /* Assign sorter by probability (could use volume instead) */
  for (h=0; h< nzones; h++)
    sorter[h] = (double)z[h].denscontrast;
    
  /* Text output file */

  printf("about to sort ...\n");FF;

  iord = (int *)malloc(nzones*sizeof(int));

  findrtop(sorter, nzones, iord, nzones);

  txt = fopen(txtfile,"w");
  fprintf(txt,"%d particles, %d voids.\n", np, nzones);
  fprintf(txt,"Void# FileVoid# CoreParticle CoreDens ZoneVol Zone#Part Void#Zones VoidVol Void#Part VoidDensContrast VoidProb\n");
  for (h=0; h<nzones; h++) {
    i = iord[h];
    prob = exp(-5.12*(z[i].denscontrast-1.) - 0.8*pow(z[i].denscontrast-1.,2.8));
    fprintf(txt,"%d %d %d %e %e %d %d %e %d %f %6.2e\n",
	    h+1, i, z[i].core, p[z[i].core].dens, z[i].vol, z[i].np, z[i].nhl, z[i].voljoin, z[i].npjoin, z[i].denscontrast, prob);

  } /* h+1 to start from 1, not zero */
  fclose(txt);

  return(0);
}
