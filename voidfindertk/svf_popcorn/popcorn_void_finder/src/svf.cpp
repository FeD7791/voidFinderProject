#include "svf.h"

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "allvars.h"
#include "finder.h"
#include "grid.h"
#include "io_lib/io.h"
#include "t.h"

using namespace std;
using namespace grid;

#define likely(x) __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)

int svf_start(float limit, std::string rawfname, float min_radius) {
  Grid *G, *Gv;
  long long int Ngrid = find_ngrid((long long int)Tracer.size(), min_radius, boxsize.x);
  // int Ngrid=12;
  float maxscale;
  int nshell = 30;
  float rhobar;

  cout << "using a grid of " << Ngrid << "**3 boxcells" << endl;
  cout << "number of tracers: " << NTRAC << " in  a box of " << boxsize.x << endl;
  cout.flush();

  //   ofstream ftest;
  //   ftest.open ("test.dat");
  //   for(int ii=0; ii < Tracer.size();ii++)
  //  {
  //       if(Tracer[ii].z<boxsize.z/10.0)
  //       {
  //           ftest << Tracer[ii].x << " ";
  //           ftest << Tracer[ii].y << " ";
  //           ftest << Tracer[ii].z << endl;
  //       }
  //  }
  //   ftest.close();
  // exit(-1);

  // cargo el grid en el que está cada tracer
  for (unsigned long int ii = 0; ii < Tracer.size(); ii++) {
    long long int i = (long long int)(Tracer[ii].x / boxsize.x * Ngrid);
    if (unlikely(i > (Ngrid - 1)))
      i = Ngrid - 1;
    long long int j = (long long int)(Tracer[ii].y / boxsize.y * Ngrid);
    if (unlikely(j > (Ngrid - 1)))
      j = Ngrid - 1;
    long long int k = (long long int)(Tracer[ii].z / boxsize.z * Ngrid);
    if (unlikely(k > (Ngrid - 1)))
      k = Ngrid - 1;
    long long int ind = i + Ngrid * (j + Ngrid * k);
    Tracer[ii].gbin = ind;
    // Tracer[ii].id=ii;
  }

  if (Tracer.size() < 16777216) {
    cout << "sorting particles by grid..." << endl;
    cout.flush();
    sort(Tracer.begin(), Tracer.end(), small_grid_index_first());
    cout << "done." << endl;
  }

  maxscale = boxsize.x / 10.0;
  NMEAN = ((double)Tracer.size() / (boxsize.x * boxsize.y * boxsize.z));

  // Construyo grid para identificar voids y calculo de velocidades
  G = new Grid(Ngrid, Tracer.size(), boxsize, maxscale);

  cout << "placing particles into the grid... " << endl;
  cout.flush();
  (*G).build(&Tracer);
  cout << "done." << endl;
  cout.flush();

  // cout << "array order: ";
  // int init=(*G).next[(*G).first[2]];
  // int iend=(*G).next[(*G).first[2]]+(*G).cnt[2];
  // for(int i=init ; i < iend; i++)
  //{
  //	cout << Tracer[i].id << " ";
  //}
  // cout << endl;
  // cout.flush();

  // cout << "linked list order: ";
  // int fst=(*G).next[(*G).first[2]];
  // int nxt=fst;
  // do
  //{
  //	cout << Tracer[nxt].id << " ";
  //	nxt=(*G).next[nxt];
  //	if(nxt==fst)break;
  //}while(1);
  // cout << endl;
  // cout.flush();

  // exit(0);

  // cell_state **status_rad;
  // status_rad =compute_steps(G,nshell);
  // cout << "done." << endl; cout.flush();

  rhobar = (float)Tracer.size() / (boxsize.x * boxsize.y * boxsize.z);
  if (access(rawfname.c_str(), F_OK) != -1) {
    cout << "leyendo raw" << endl;
    FILE *fd;
    Sphere dumm;
    fd = fopen(rawfname.c_str(), "r");
    int coun = 0;
    while (fread(&(dumm.id), sizeof(int), 1, fd) == 1) {
      fread(&(dumm.radius), sizeof(float), 1, fd);
      // fread(&(dumm.Ref    ),sizeof(float),1,fd);
      fread(&(dumm.x), sizeof(float), 1, fd);
      fread(&(dumm.y), sizeof(float), 1, fd);
      fread(&(dumm.z), sizeof(float), 1, fd);
      fread(&(dumm.delta), sizeof(float), 1, fd);
      fread(&(dumm.ToF), sizeof(bool), 1, fd);
      fread(&(dumm.np), sizeof(int), 1, fd);
      fread(&(dumm.irad), sizeof(int), 1, fd);
      if (dumm.ToF)
        coun++;
      Sphvd.push_back(dumm);
    }
    fclose(fd);
    cout << "A encontrar: " << coun << endl;
    NVOID = Sphvd.size();
  } else {

    // esta funcion deberia llamarse diferente. ya que sirve para
    // darle una escala de densidad a cada grid del box, no es exactamente
    // una busqueda de centros (de hecho la información del centro estaría
    // contenida en el índice si no se elimina niguno (como pasa el 100% de
    // los casos, pero deberíamos reescribir esa parte para asegurarlo si
    // removemos los campos .x .y .z y los .Ini)
    // FindCenters_fast(limit,rhobar,nshell,G, &Tracer,status_rad,min_radius);
    FindCenters_fast(limit, rhobar, nshell, G, min_radius);

    // FindSphvds_fast(limit,rawfname,nshell, G, status_rad);
    // FindSphvds_fastV2(limit,rawfname,nshell, G, status_rad,min_radius);
    FindSphvds_fastV2(limit, rawfname, nshell, G, min_radius);
  }
  // GridV();
  Gv = new Grid(Ngrid, Sphvd.size(), boxsize, 2 * maxscale);
  (*Gv).build(&Sphvd);

  // cell_state **status_rad2;
  // int nshell2 = 2 * nshell;
  // status_rad2 =compute_steps(Gv,nshell2);
  CleanSphvds_Tol0(Gv, maxscale);
  // CleanSphvds(Gv, maxscale, nshell2);

  // FreeGridV();
  // FreeGridT();

  // Calculo de velocidad del void
  // ComputeVelocity();

  // Calculo de perfiles
  // ComputeProfiles(fname_tracer);

  // Libero memoria

  return (0);
}
