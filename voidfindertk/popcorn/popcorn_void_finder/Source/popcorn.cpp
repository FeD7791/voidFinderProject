#include "popcorn.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// global data
#include "allvars.h"
#include "colors.h"
#include "t.h"

#include "arvo_functions.h"
#include "grid.h"
#include "io_lib/io.h"
#include "objects.h"

using namespace std;
using namespace grid;

struct maxball {
  double volret = -1.0;
  int imax = -1;
  vector<int> in_halos;
  int nhret = -1;
  bool stat = false;
};

void popcorn_work(int cnt, int nshell, float min_radius, double eps, float limit, Grid *G) {

  //// mean numerical density of matter tracers
  float rhobar = Tracer.size() / boxsize.x / boxsize.y / boxsize.z;

  // nuevo pop void sin miembros
  popcorn *seed = new popcorn();
  // copiando las coordendas del candidato esférico a void tipo pop
  Sphere sph;
  sph.x = Sphvd[cnt].x;
  sph.y = Sphvd[cnt].y;
  sph.z = Sphvd[cnt].z;
  sph.radius = Sphvd[cnt].radius;
  sph.lvl = 0;              // esfera base del pop
  sph.id = 0;               // esfera base del pop
  seed->id = Sphvd[cnt].id; // mismo id que el del void esférico

  // ids of tracers inside the popcorn
  vector<int> in_halos;
  int nhret = init_flags(&sph, G, &Tracer, &in_halos);
  // ACA marco las particulas internas como visit y lo meto de una al void
  // esferico el seed, luego la funcion inflate pasa a solo tratar con
  // seed pops que tienen al menos una esfera

  double volret = 4.0 / 3.0 * PI * sph.radius * sph.radius * sph.radius;
  double convergence =
      1.0; // incremento relativo del volumen (debe ser siempre decreciente para asegurar convergencia uniforme)
  seed->vol = volret;         // actualizo  el volumen
  sph.vol = seed->vol;        // volumen inicial
  seed->npart = nhret;        // actualizando numero de particulas
  seed->membs.push_back(sph); // cargando esfera aceptada
  seed->nmem = 1;             // cargando esfera aceptada

  seed->in_halos.insert(seed->in_halos.end(), in_halos.begin(), in_halos.end());

  // boolean flags for matter tracers
  vector<bool> visit(Tracer.size(), false);
  for (long unsigned int i = 0; i < in_halos.size(); i++) {
    visit[in_halos[i]] = true;
  }

  vector<Sphere> chunck;
  vector<Sphere> level;

  level.push_back(sph);
  while (!level.empty()) {
    if (sph.radius < min_radius) {
      level.clear();
      break; // radios por debajo del radio de shotnoise no se pochoclean
    }

    chunck.clear(); // donde voy a guardar la lista de pochocleados
    // pochocleando el popcorn
    for (auto its = level.begin(); its != level.end(); ++its) { // iterando sobre el level
      // subarray de pochoclos de cada bola del level
      vector<Sphere> chucky = pochoclea(*its, &(seed->membs), boxsize, eps);

      for (long unsigned int e = 0; e < chucky.size(); e++) {
        chucky[e].lvl = its->lvl + 1; // subo un lvl
      }
      chunck.insert(chunck.end(), chucky.begin(), chucky.end());
    }

    level.clear();
    // buscando en orden de volumen máximo de contribucion las esferas del lvl que  amplien el volumen del pop
    // vamos lvl por lvl
    while (!chunck.empty()) {

      maxball maxi;
      vector<maxball> arr(chunck.size());
#pragma omp parallel default(none) shared(chunck, seed, G, nshell, visit, limit, rhobar, Tracer, min_radius, arr)
      {
#pragma omp single
        {
          for (long unsigned int ii = 0; ii < chunck.size(); ii++) {
#pragma omp task default(none) firstprivate(ii)                                                                        \
    shared(seed, G, nshell, visit, limit, rhobar, Tracer, min_radius, arr, chunck)
            {
              arr[ii].stat = false;
              arr[ii].in_halos.clear();
              arr[ii].stat = inflate_fast(&(chunck[ii]), seed, G, nshell, &visit, &(arr[ii].in_halos),
                                          &(arr[ii].volret), &(arr[ii].nhret), limit, rhobar, Tracer, min_radius);
            }
          }
        }
      }

      maxi.imax = -1;
      maxi.volret = -1.E12;
      for (long unsigned int ii = 0; ii < chunck.size(); ii++) {
        if ((arr[ii].stat) && (arr[ii].volret > maxi.volret)) {
          maxi.imax = ii;
          maxi.volret = arr[ii].volret;
          maxi.in_halos = arr[ii].in_halos;
          maxi.nhret = arr[ii].nhret;
        }
      }

      if (maxi.imax >= 0) {
        if (maxi.volret < 0.0) {
          cout << "volumen negativo " << endl;
          exit(-1);
        } // exception!

        /// NUEVO CRITERIO DE CONVERGENCIA
        // la serie de correciones debe ser decreciente
        double new_increment = (maxi.volret - (seed->vol)) / (maxi.volret);
        if (new_increment > convergence) {
          maxi.volret = -1.0;
          maxi.imax = -1;
          maxi.in_halos.clear();
          maxi.nhret = -1;
          maxi.stat = false;
          break;
        } else {
          convergence = new_increment;
          seed->vol = maxi.volret;                  // actualizo  el volumen
          seed->npart = maxi.nhret;                 // actualizando numero de particulas
          chunck[maxi.imax].id = seed->nmem;        // id de la esfera aceptada dentro del pop
          chunck[maxi.imax].vol = seed->vol;        // volumen al momento de cargar la esfera
          seed->membs.push_back(chunck[maxi.imax]); // cargando esfera aceptada
          seed->nmem++;                             // numero de esferas aceptadas
          seed->in_halos.insert(seed->in_halos.end(), maxi.in_halos.begin(), maxi.in_halos.end());
          for (long unsigned int ii = 0; ii < maxi.in_halos.size(); ii++) {
            visit[maxi.in_halos[ii]] = true;
          }

          level.push_back(chunck[maxi.imax]);       // anotandola para que se pochoclee en el siguiente lvl
          chunck.erase(chunck.begin() + maxi.imax); // sacando de la lista de candidatos al que agregamos
          // ahora se vuelve a revisar la lista de candidatos de este lvl buscando el siguiente máximo
          // antes de volver a pochoclear el lvl que sigue
          maxi.volret = -1.0;
          maxi.imax = -1;
          maxi.in_halos.clear();
          maxi.nhret = -1;
          maxi.stat = false;
        }

      } else {
        maxi.volret = -1.0;
        maxi.imax = -1;
        maxi.in_halos.clear();
        maxi.nhret = -1;
        maxi.stat = false;
        // ninguna de las esperas candidatos dio un stat positivo
        break;
      }
    }
  }

  if (seed->nmem > 0 && !(seed->in_halos.empty())) { // si esta vacio de particulas no lo uso....
#pragma omp critical
    { pop_voids.push_back(*seed); }
    if (seed->npart != (int)seed->in_halos.size()) {
      cout << seed->npart << "   -------" << seed->id << "----------  " << (seed->in_halos).size() << endl;
      cout.flush();
    }
  }
  seed->membs.clear();
  seed->in_halos.clear();
  delete seed;
}

int popcorn_start(float min_radius, double eps, float limit, int rank, int numrank) {
  // int Ngrid=512;
  int Ngrid = find_ngrid(Tracer.size(), min_radius, boxsize.x);
  int nshell = 30;

  int nvds = Sphvd.size();
  // Cleaning small voids
  float radmax = -1.0;
  for (int i = 0; i < nvds; i++) {
    if (radmax < Sphvd[i].radius) {
      radmax = Sphvd[i].radius;
    }
  }

  float maxscale = radmax * 1.1;

  for (long unsigned int ii = 0; ii < Tracer.size(); ii++) {
    int i = int(Tracer[ii].x / boxsize.x * Ngrid);
    if (unlikely(i > (Ngrid - 1))) {
      i = Ngrid - 1;
    }
    int j = int(Tracer[ii].y / boxsize.y * Ngrid);
    if (unlikely(j > (Ngrid - 1))) {
      j = Ngrid - 1;
    }
    int k = int(Tracer[ii].z / boxsize.z * Ngrid);
    if (unlikely(k > (Ngrid - 1))) {
      k = Ngrid - 1;
    }
    int ind = i + Ngrid * (j + Ngrid * k);
    Tracer[ii].gbin = ind;
  }

  cout << "sorting particles by grid..." << endl;
  cout.flush();
  sort(Tracer.begin(), Tracer.end(), small_grid_index_first());
  cout << "done." << endl;

  Grid *G = new Grid(Ngrid, Tracer.size(), boxsize, maxscale);
  cout << "placing particles into the grid... " << Ngrid << endl;
  cout.flush();
  G->build(&Tracer);
  cout << "done." << endl;
  cout.flush();

  for (int i = 0; i < nvds; i++) {
    if (Sphvd[i].radius < G->maxscale / (float)nshell) {
      Sphvd[i].erase = true;
    }
  }

  int cc = 0;
  for (int i = 0; i < nvds; i++) {
    if (!Sphvd[i].erase) {
      cc++;
    }
  }

  cout << LGREEN << "total de candidatos: " << cc << DEFA << endl;
  /// Main computation /////////////////////////////////////////////////////
  cout << CYAN << "procesando.... " << DEFA << endl;

  clock_t t = clock();
  // for(cnt=0;cnt<Sphvd.size(); cnt++)
  //
  int start = (int)((float)Sphvd.size() / (float)numrank * rank);
  int end = (int)((float)Sphvd.size() / (float)numrank * (rank + 1));
  if (rank == (numrank - 1)) {
    end = Sphvd.size();
  }

  for (int cnt = start; cnt < end; cnt++) {
    if (cnt % 10 == 0) {
      cout << cnt << " in "
           << "(" << start << "," << end << ")" << endl;
      cout.flush();
    }
    // cout << cnt << " in " << "("<< start << "," << end << ")" << endl; cout.flush();

    if (Sphvd[cnt].erase) {
      continue;
    }

    popcorn_work(cnt, nshell, min_radius, eps, limit, G);
  }

  Time(t, NCORES);

  return 0;
}
