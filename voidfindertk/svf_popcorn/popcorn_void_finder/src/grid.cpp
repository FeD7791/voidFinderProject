#include "grid.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "allvars.h"
#include "arvo_functions.h"
#include "constants.h"
#include "objects.h"
#include "template_popcorn.h"

using namespace std;
using namespace arvo;

namespace grid {
Grid::Grid(int nbin, unsigned long long int nt, float3 box, float maxsc) {
  maxscale = maxsc;
  Ngrid = nbin;
  Npart = nt;
  boxsize.x = box.x;
  boxsize.y = box.y;
  boxsize.z = box.z;
  /// OJO!!!!!!!!!!!!!!!!!!1
  gridsize = box.x / (float)nbin;

  long long int n3 = nbin * nbin * nbin;
  first = new long long int[n3];
  cnt = new unsigned long long int[n3];
  fill_n(first, n3, -1);
  fill_n(cnt, n3, 0);

  next = new long long int[nt];
  nsearch = (int)ceil(maxsc / gridsize) + 1;
  nsten = (nsearch + 2);
  // int ns2 = nsten * nsten;
  // int ns3 = ns2 * nsten;
}

// version_original
int indx_period_old(int ind, int ngrid) {
  int i = ind;

  if (i > (ngrid - 1))
    i = i - ngrid;
  if (i < 0)
    i = i + ngrid;
  return (i);
}

// mas simple en assembler pero mas lenta que _old
int indx_period_old2(int ind, int ngrid) {
  int i = ind % ngrid;
  i = i < 0 ? i + ngrid : i;

  return (i);
}

// mas simple en assembler pero mas lenta que _old
int indx_period_old3(int ind, int ngrid) {
  int i = (ind + ngrid) % ngrid;
  return (i);
}

// ngrid debe ser potencia de dos
int indx_period(int ind, int ngrid) {
  int i = ind + ngrid;
  return (i & (ngrid - 1));
}

cell_state is_there(int shifti, int shiftj, int shiftk, float rad) {
  // rad is in gridsize units
  int i = abs(shifti);
  int j = abs(shiftj);
  int k = abs(shiftk);
  int caso = (int)(i != 0) + (int)(j != 0) + (int)(k != 0);

  cell_state stat;
  stat.is_out = (((float)(i * (i - 1) + j * (j - 1) + k * (k - 1)) - rad * rad - sqrt((float)caso) * rad) > 0.0);
  stat.is_in = (((float)(i * (i + 1) + j * (j + 1) + k * (k + 1)) - rad * rad + sqrt((float)caso) * rad) < 0.0);

  return (stat);
}

float periodicity_delta(float coord, float boxsize) {
  float returned;
  returned = coord;

  if (returned > boxsize * 0.5) {
    returned = boxsize - returned;
  }
  return (returned);
};

//  cell_state **compute_steps(Grid* g,int nstep)
//  {
//  	 int ntotal = (*g).nsearch;
//  	 int dim=2*ntotal+1;
//  	 int ind;
//  	 int i,j,k,irad;
//  	 float rad;
//
//  	 cell_state **status_rad;
//
//  	 status_rad = new cell_state*[nstep];
//  	 for(i = 0; i < nstep; ++i)
//  		 status_rad[i] = new cell_state[dim*dim*dim];
//
//  	 for(irad=0;irad<nstep;irad++)
//  	 {
//  	    rad=(*g).maxscale*(float)(irad+1)/(float)nstep;
//  	    for(k=0; k< dim; k++)
//              {
//                   for(j=0; j < dim; j++)
//                   {
//                       for(i=0; i< dim; i++)
//                       {
//  	           	 ind=i+dim*(j+dim*k);
//  	    		 //status_rad[irad][ind]=
//  is_there(i-ntotal,j-ntotal,k-ntotal,rad, g);
//			 //verificando si tengo mucho cache misssing
//  	    		 status_rad[nstep-irad-1][ind]=
//  is_there(i-ntotal,j-ntotal,k-ntotal,rad, g);
//  		       }
//  	           }
//  	      }
//  	 }
//
//  	 return(status_rad);
//  }
//
float periodicity(float coord, float boxsize) {
  float returned;
  returned = coord;
  if (returned > boxsize) {
    returned = returned - boxsize;
  }
  if (returned < 0.0) {
    returned = returned + boxsize;
  }

  return (returned);
};

float distance(float x1, float y1, float z1, float x2, float y2, float z2, float3 boxsize) {
  float dx, dy, dz;
  float d;

  dx = fabs((x1 - x2));
  dx = periodicity_delta(dx, boxsize.x);
  dy = fabs((y1 - y2));
  dy = periodicity_delta(dy, boxsize.y);
  dz = fabs((z1 - z2));
  dz = periodicity_delta(dz, boxsize.z);

  d = sqrt(dx * dx + dy * dy + dz * dz);
  return (d);
}

ostream &error() { return cerr << "ERROR : "; }

int index3(int ix, int iy, int iz, int nside) {
  int ind;
  int indx, indy, indz;
  indx = ix;
  indy = iy;
  indz = iz;

  if (indx > (nside - 1))
    indx = 0;
  if (indx < 0)
    indx = nside - 1;
  if (indy > (nside - 1))
    indy = 0;
  if (indy < 0)
    indy = nside - 1;
  if (indz > (nside - 1))
    indz = 0;
  if (indz < 0)
    indz = nside - 1;
  ind = (indz * nside + indy) * nside + indx;
  return ind;
}

vector<int> revertindex(int ind, int nside) {
  vector<int> indx;
  int ix, iy, iz;

  indx.clear();
  ix = ind % nside;
  iy = (ind - ix) / nside % nside;
  iz = ((ind - ix) / nside - iy) / nside;

  indx.push_back(ix);
  indx.push_back(iy);
  indx.push_back(iz);

  return indx;
}

bool check_inside(Sphere cta, vector<Sphere> *pop, float3 boxsize, double eps) {
  vector<Sphere>::iterator it;
  bool retbool;
  float d;

  retbool = false;

  for (it = (*pop).begin(); it != (*pop).end(); ++it) {
    if (cta.id == (*it).id)
      continue;
    d = distance(cta.x, cta.y, cta.z, (*it).x, (*it).y, (*it).z, boxsize);
    retbool = (d < (((*it).radius) + eps));
    if (retbool)
      break;
  }
  return (retbool);
}

void fibonacci_sphere(int i, int samples, double *vec) {
  int rnd = 1;

  double offset = 2.0 / ((double)samples);
  double increment = PI * (3.0 - sqrt(5.0));
  double r;
  double phi;

  vec[1] = ((i * offset) - 1) + (offset / 2);
  r = sqrt(1 - vec[1] * vec[1]);
  phi = (double)((i + rnd) % samples) * increment;
  vec[0] = cos(phi) * r;
  vec[2] = sin(phi) * r;
}

vector<Sphere> pochoclea(Sphere root, vector<Sphere> *pop, float3 boxsize, double eps) {
  int lim = (int)round(4 * PI * root.radius * root.radius); // un punto por mpc cuadrado
                                                            // de area de la esfera
  if (lim > 120) {
    lim = 120;
  }

  vector<Sphere> collar;
  for (int i = 0; i < lim; i++) {
    double vec[3];
    fibonacci_sphere(i, lim, vec);

    Sphere cta;
    cta.x = ((float)vec[0]) * root.radius + root.x;
    cta.y = ((float)vec[1]) * root.radius + root.y;
    cta.z = ((float)vec[2]) * root.radius + root.z;

    cta.x = periodicity(cta.x, boxsize.x);
    cta.y = periodicity(cta.y, boxsize.y);
    cta.z = periodicity(cta.z, boxsize.z);
    cta.id = root.id;

    bool inner = check_inside(cta, pop, boxsize, eps);

    if (!inner) {
      collar.push_back(cta);
    }
  }
  return (collar);
}

void clean_candidates(Sphere *sph, int nside, vector<class cell> *grid, float3 boxsize, vector<Sphere> rts) {
  int ind;
  int ip;
  float d;
  vector<class cell>::iterator ucell;

  ind = which_cell((*sph), nside, boxsize);

  for (int j = 0; j < 27; j++) {
    ucell = (*grid).begin() + (*grid)[ind].next[j];
    for (long unsigned int i = 0; i < (*ucell).in.size(); i++) {
      ip = (*ucell).in[i];
      if (!rts[ip].erase) // si no se cumple quiere decir que el candidato ya es
                          // interno a algun void
      {
        d = distance((*sph).x, (*sph).y, (*sph).z, rts[ip].x, rts[ip].y, rts[ip].z, boxsize);
        if (d < (*sph).radius) {
          rts[ip].erase = true; // esta adentro, lo descartamos como candidato
        }
      }
    }
  }
};

void show_progress(float progress) {
  int barWidth = 120;

  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      std::cout << "=";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

double unidim_split(double c0, double c1, float boxsize) {
  double returned;

  returned = c1;
  if ((c1 - c0) > boxsize * 0.5)
    returned = c1 - boxsize;
  if ((c1 - c0) < (-boxsize * 0.5))
    returned = c1 + boxsize;
  return (returned);
}

void check_boundary_splits(double *test_pop, int sizes, float3 boxsize) {
  int i;
  double xc, yc, zc;
  double xm, ym, zm;
  vector<int> to_erase;

  xc = test_pop[0];
  yc = test_pop[1];
  zc = test_pop[2];

  xm = 0.0;
  ym = 0.0;
  zm = 0.0;
  for (i = 0; i < sizes; i++) {
    test_pop[i * 4] = unidim_split(xc, test_pop[i * 4], boxsize.x);         /// .x
    test_pop[i * 4 + 1] = unidim_split(yc, test_pop[i * 4 + 1], boxsize.y); /// .y
    test_pop[i * 4 + 2] = unidim_split(zc, test_pop[i * 4 + 2], boxsize.z); /// .z

    xm = xm + test_pop[i * 4];
    ym = ym + test_pop[i * 4 + 1];
    zm = zm + test_pop[i * 4 + 2];
  }

  xm = xm / (double)sizes;
  ym = ym / (double)sizes;
  zm = zm / (double)sizes;

  for (i = 0; i < sizes; i++) {
    test_pop[i * 4] = test_pop[i * 4] - xm;         /// .x
    test_pop[i * 4 + 1] = test_pop[i * 4 + 1] - ym; /// .y
    test_pop[i * 4 + 2] = test_pop[i * 4 + 2] - zm; /// .z
  }
}

void check_inner_spheres(double **test_pop, int *sizes, double eps) {
  int i, j, newsize;
  double dx, dy, dz, d;
  bool *erase;
  double *buffer;

  buffer = new double[(*sizes) * 4];
  erase = new bool[(*sizes)];
  for (i = 0; i < (*sizes); i++)
    erase[i] = false;

  for (i = 0; i < (*sizes); i++) {
    buffer[i * 4] = (*test_pop)[i * 4];
    buffer[i * 4 + 1] = (*test_pop)[i * 4 + 1];
    buffer[i * 4 + 2] = (*test_pop)[i * 4 + 2];
    buffer[i * 4 + 3] = (*test_pop)[i * 4 + 3];

    for (j = 0; j < (*sizes); j++) {
      if ((*test_pop)[i * 4 + 3] < (*test_pop)[j * 4 + 3])
        continue;
      // ((*test_pop)[i * 4 + 3] >= (*test_pop)[j * 4 + 3])
      if (i == j)
        continue;

      dx = (*test_pop)[i * 4] - (*test_pop)[j * 4];
      dy = (*test_pop)[i * 4 + 1] - (*test_pop)[j * 4 + 1];
      dz = (*test_pop)[i * 4 + 2] - (*test_pop)[j * 4 + 2];

      d = sqrt(dx * dx + dy * dy + dz * dz);
      if (d < ((*test_pop)[i * 4 + 3] - (*test_pop)[j * 4 + 3] + eps)) {

        if ((*test_pop)[i * 4 + 3] == (*test_pop)[j * 4 + 3]) {
          if (j > i)
            erase[j] = true;
        } else {
          erase[i] = true;
        }
      }
    }
  }

  newsize = 0;
  for (i = 0; i < (*sizes); i++) {
    if (!(erase[i]))
      newsize++;
  }

  delete[] *test_pop;
  *test_pop = new double[newsize * 4];

  j = 0;
  for (i = 0; i < (*sizes); i++) {
    if (!(erase[i])) {
      (*test_pop)[j * 4] = buffer[i * 4];
      (*test_pop)[j * 4 + 1] = buffer[i * 4 + 1];
      (*test_pop)[j * 4 + 2] = buffer[i * 4 + 2];
      (*test_pop)[j * 4 + 3] = buffer[i * 4 + 3];
      j++;
    }
  }

  *sizes = newsize;
  delete[] buffer;
  delete[] erase;
}

void check_inner_spheres_v2(double **test_pop, int *sizes) {
  int i, j, newsize;
  bool *erase;
  double *buffer;
  double area1, vol, volprev;

  buffer = new double[(*sizes) * 4];
  erase = new bool[(*sizes)];
  for (i = 0; i < (*sizes); i++) {
    erase[i] = false;
    buffer[i * 4] = (*test_pop)[i * 4];
    buffer[i * 4 + 1] = (*test_pop)[i * 4 + 1];
    buffer[i * 4 + 2] = (*test_pop)[i * 4 + 2];
    buffer[i * 4 + 3] = (*test_pop)[i * 4 + 3];
  }

  volprev = 4.0 / 3.0 * PI * (*test_pop)[3] * (*test_pop)[3] * (*test_pop)[3];

  for (i = 1; i < (*sizes); i++) {
    Arvo_class *arvo_obj = new Arvo_class();
    (*arvo_obj).compute_volume_arvo(&area1, &vol, (*test_pop), i + 1);
    arvo_obj->Cleanup();
    delete (arvo_obj);

    if ((vol <= volprev)) {
      erase[i] = true;
    } else {
      volprev = vol;
    }
  }

  newsize = 0;
  for (i = 0; i < (*sizes); i++) {
    if (!(erase[i]))
      newsize++;
  }

  delete[] *test_pop;
  *test_pop = new double[newsize * 4];

  j = 0;
  for (i = 0; i < (*sizes); i++) {
    if (!(erase[i])) {
      (*test_pop)[j * 4] = buffer[i * 4];
      (*test_pop)[j * 4 + 1] = buffer[i * 4 + 1];
      (*test_pop)[j * 4 + 2] = buffer[i * 4 + 2];
      (*test_pop)[j * 4 + 3] = buffer[i * 4 + 3];
      j++;
    }
  }

  *sizes = newsize;
  delete[] buffer;
  delete[] erase;
}
bool inflate_fast(Sphere *sph, class popcorn *pope,
                  // Grid* g,cell_state **status_rad,int nstep,
                  Grid *g, int nstep, vector<bool> *visit, vector<int> *in_halos, double *volret, int *nhret,
                  float limit, float density, vector<halo> &Tracer, float min_radius) {
  int ind, i, j, k, indx;
  long long int it, jt, kt, fst, nxt;
  vector<Sphere>::iterator itpop;
  vector<double> test_pop;
  int ns;
  unsigned long long int npx, npx2;
  float dd, rad, radx;
  float deltax, delta_inf;
  double area, vol;
  int npart;
  float vol_new = -1.0;
  vector<class cell>::iterator ucell;
  int sphsize;
  double *test_arr;
  float_in res, d;

  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  halo cent;
  // int ntotal = (*g).nsearch;
  // int dim = 2 * ntotal + 1;
  vector<int> visitados;
  vector<float_in> dist;
  cell_state stat_out, stat_in, stat;

  cent.x = (*sph).x;
  cent.y = (*sph).y;
  cent.z = (*sph).z;

  npart = (*pope).npart; // numero de particulas que ya habia antes

  sphsize = ((int)(((*pope).membs).size()));
  test_arr = new double[4 * (sphsize + 1)];
  i = 0;
  for (itpop = ((*pope).membs).begin(); itpop != ((*pope).membs).end(); ++itpop) {
    test_arr[4 * i] = (double)(*itpop).x;
    test_arr[4 * i + 1] = (double)(*itpop).y;
    test_arr[4 * i + 2] = (double)(*itpop).z;
    test_arr[4 * i + 3] = (double)(*itpop).radius;
    i++;
  }
  /// esfera a pruena
  test_arr[4 * i] = (double)(*sph).x;
  test_arr[4 * i + 1] = (double)(*sph).y;
  test_arr[4 * i + 2] = (double)(*sph).z;
  test_arr[4 * i + 3] = 1e-5;

  check_boundary_splits(test_arr, sphsize + 1, boxsize);

  // calculo de la densidad en esferas crecientes primero de manera cruda para
  // identificar una escala aproximada de inflado
  int ic = int(cent.x / boxsize.x * Ngrid);
  int jc = int(cent.y / boxsize.y * Ngrid);
  int kc = int(cent.z / boxsize.z * Ngrid);

  res.d = -1.0;
  res.cnt = 0;
  res.id = -1;

  int *counter = new int[nstep];
  for (int irad = 0; irad < nstep; irad++)
    counter[irad] = 0;

  // ns=(int)ceil(rad/(*g).gridsize)+1;
  ns = (int)ceil((*g).maxscale / (*g).gridsize) + 1;

  for (k = (-ns); k < (ns + 1); k++) // k1
  {
    kt = indx_period(k + kc, Ngrid);
    for (j = (-ns); j < (ns + 1); j++) // j1
    {
      jt = indx_period(j + jc, Ngrid);
      for (i = (-ns); i < (ns + 1); i++) // i1
      {
        it = indx_period(i + ic, Ngrid);
        // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
        // stat=status_rad[nstep-irad-1][indi];

        for (int irad = 0; irad < nstep; irad++) {
          rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
          float radg = rad / (*g).gridsize;
          stat = is_there(i, j, k, radg);

          if (stat.is_in) // para la couta inferior basta con los grides
                          // totalmente contenidos
          {
            ind = it + Ngrid * (jt + Ngrid * kt);
            fst = (*g).first[ind];
            if (fst >= 0) // grid no vacio
            {
              nxt = fst;
              do {
                if (!(*visit)[nxt]) // que no sea ya interior al pop
                {
                  counter[irad]++;
                }
                nxt = (*g).next[nxt];
                if (nxt == fst)
                  break;
              } while (1);
            }
          }
        }
      } // i1
    }   // j1
  }     // k1

  for (int irad = nstep - 1; irad >= 0; irad--) // loop1
  {
    rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
    test_arr[4 * (sphsize) + 3] = rad;
    Arvo_class *arvo_obj = new Arvo_class();
    try {
      // calculo el volumen
      arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
    } catch (const std::bad_alloc &e) {
      arvo_obj->Cleanup();
      exit(-1);
    }
    delete (arvo_obj);

    // cota inferior de la densidad
    delta_inf = ((float)(counter[irad] + npart) / vol) / density - 1.0;
    if (delta_inf < limit) {
      indx = irad;
      // calculamos cuanto es la delta exacta al radio correspondiente a indx,
      // si es mayor diminuimos y volvemos a calcular, el bucle termina cuando
      // tenemos una cascara donde la delta es menor al limite (y por lo tanto
      // cascara inmediata sup es mayor) o cuando llegamos al ultimo cascaron y
      // la delta sigue siendo superior

      do // loop 2
      {
        rad = (*g).maxscale * (float)(indx + 1) / (float)nstep;
        float radg = rad / (*g).gridsize;
        ns = (int)ceil(rad / (*g).gridsize) + 1;
        visitados.clear();

        npx = 0;
        for (k = (-ns); k < (ns + 1); k++) // k2
        {
          kt = indx_period(k + kc, (*g).Ngrid);
          for (j = (-ns); j < (ns + 1); j++) // j2
          {
            jt = indx_period(j + jc, (*g).Ngrid);
            for (i = (-ns); i < (ns + 1); i++) // i2
            {
              it = indx_period(i + ic, (*g).Ngrid);
              stat = is_there(i, j, k, radg);
              if (!stat.is_out) {
                ind = it + (*g).Ngrid * (jt + (*g).Ngrid * kt);
                fst = (*g).first[ind];
                if (fst >= 0) {
                  nxt = fst;
                  do {
                    if (!(*visit)[nxt]) {
                      if (stat.is_in) {
                        npx++;
                        visitados.push_back(nxt);
                      } else {
                        dd = distance(Tracer[nxt], cent, boxsize);
                        if (dd < rad) {
                          npx++;
                          visitados.push_back(nxt);
                        }
                      }
                    }
                    nxt = (*g).next[nxt];
                    if (nxt == fst)
                      break;
                  } while (1);
                }
              }

            } // i2
          }   // j2
        }     // k2

        test_arr[4 * sphsize + 3] = rad;
        Arvo_class *arvo_obj = new Arvo_class();
        try {
          // calculo el volumen
          arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
        } catch (const std::bad_alloc &e) {
          arvo_obj->Cleanup();
          exit(-1);
        }
        delete (arvo_obj);

        // delta exacta al radio indx
        deltax = ((float)(npx + npart) / vol) / density - 1.0;

        if (deltax < limit) {
          // dado que estamos seguros que la cota minima de indx+1 es superior
          // al limite (asi funciona el bucle anterior) y sabemos que en este
          // indx estamos debajo, el limite se cruza en un radio intermedio
          // busca entre los radios correspondientes a los cascarones de indx y
          // indx+1, el indice de la particula que da el cambio de densidad para
          // eso utiliza solo los boxcell que debe abrir segun la desigualdad
          // del triangulo (como fast_counter hace)
          // res=shell_finder(cent,limit,density,indx,nstep,g,&Tracer,status_rad);

          // shell_finder es la version esferica, debemos tener en cuenta la
          // geometria del pop
          dist.clear();
          float radsup = (*g).maxscale * (float)(indx + 2) / (float)nstep;
          float radinf = (*g).maxscale * (float)(indx + 1) / (float)nstep;
          float radsupg = radsup / (*g).gridsize;
          float radinfg = radinf / (*g).gridsize;
          for (k = (-ns); k < (ns + 1); k++) // k3
          {
            kt = indx_period(k + kc, Ngrid);
            for (j = (-ns); j < (ns + 1); j++) // j3
            {
              jt = indx_period(j + jc, Ngrid);
              for (i = (-ns); i < (ns + 1); i++) // i3
              {
                it = indx_period(i + ic, Ngrid);

                ind = it + Ngrid * (jt + Ngrid * kt);
                fst = (*g).first[ind];
                if (fst >= 0) {
                  // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
                  // en el máximo radio de inflado ya estoy por debajo de la
                  // delta
                  if (indx == nstep - 1) {
                    // stat_in =status_rad [nstep- indx   -1][indi];
                    stat_in = is_there(i, j, k, radinfg);
                    // opero los logicos para asegurarme que se abra el grid
                    // para todas las particulas internas y asi se carguen
                    stat_in.is_in = !stat_in.is_in;
                    stat_out.is_out = false;
                  } else {
                    // stat_out=status_rad [nstep-(indx+1)-1][indi];
                    stat_out = is_there(i, j, k, radsupg);
                    // stat_in =status_rad [nstep- indx   -1][indi];
                    stat_in = is_there(i, j, k, radinfg);
                  }
                  if (!(stat_out.is_out) && !(stat_in.is_in)) {
                    // open grid
                    nxt = fst;
                    do {
                      if (!(*visit)[nxt]) {
                        d.d = distance(Tracer[nxt], cent, boxsize);
                        d.id = nxt;
                        if ((d.d < radsup) && (d.d >= radinf)) {
                          dist.push_back(d);
                        }
                      }
                      nxt = (*g).next[nxt];
                      if (nxt == fst)
                        break;
                    } while (1);
                  }
                }

              } // i3
            }   // j3
          }     // k3

          // sort (dist.begin(), dist.end());
          // for(i=(dist.size()-1);i>=0;i--)
          sort(dist.begin(), dist.end(), greater<float_in>());
          long unsigned int ii;
          for (ii = 0; ii < dist.size(); ii++) {
            npx2 = (dist.size() - ii) + npx + npart;
            d = dist[ii];
            radx = d.d;

            test_arr[4 * sphsize + 3] = radx;
            Arvo_class *arvo_obj = new Arvo_class();
            try {
              // calculo el volumen
              arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
            } catch (const std::bad_alloc &e) {
              arvo_obj->Cleanup();
              exit(-1);
            }
            delete (arvo_obj);
            deltax = ((float)npx2 / vol) / density - 1.0;

            if (deltax < limit) {
              // happy ending 1
              res.d = radx;
              res.cnt = npx2;
              res.id = d.id;
              vol_new = vol;

              if (radx < min_radius) {
                delete[] test_arr;
                return (false); // no sirve el centro la esfera es menor al
                                // radio minimo
              }
              break;
            }
          }
          // cargando las nuevas particulas internas
          // cargando la cascara
          for (long unsigned int jj = ii; jj < dist.size(); jj++) {
            (*in_halos).push_back(dist[jj].id);
          }

          // cargando dentro del radio inferior
          for (long unsigned int jj = 0; jj < visitados.size(); jj++) {
            (*in_halos).push_back(visitados[jj]);
          }

          break; // loop2 //happy ending
        } else {
          indx--;
          if ((*g).maxscale * (float)(indx + 2) / (float)nstep < min_radius) {
            break;
          }
        }
      } while (indx >= 0); // loop2

      if (indx == (-1)) {
        rad = (*g).maxscale / (float)nstep;
        float radg = rad / (*g).gridsize;
        if (rad < min_radius)
          break;
        // res=ball_finder(cent,limit,density,0,nshell,g,&Tracer,status_rad);
        // hacemos lo mismo pero teniendo en cuenta la geometria del pop
        ns = (int)ceil(rad / (*g).gridsize) + 1;
        dist.clear();
        for (k = (-ns); k < (ns + 1); k++) // k4
        {
          kt = indx_period(k + kc, Ngrid);
          for (j = (-ns); j < (ns + 1); j++) // j4
          {
            jt = indx_period(j + jc, Ngrid);
            for (i = (-ns); i < (ns + 1); i++) // i4
            {
              it = indx_period(i + ic, Ngrid);
              ind = it + Ngrid * (jt + Ngrid * kt);
              fst = (*g).first[ind];
              if (fst >= 0) {
                // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
                // stat=status_rad[nstep-irad-1][indi]; //is_there(i,j,k,rad,
                // g);
                stat = is_there(i, j, k, radg);
                if (!(stat.is_out)) {
                  // open grid
                  nxt = fst;
                  do {
                    if (!(*visit)[nxt]) {
                      d.d = distance(Tracer[nxt], cent, boxsize);
                      d.id = nxt;
                      if (d.d < rad) {
                        dist.push_back(d);
                      }
                    }
                    nxt = (*g).next[nxt];
                    if (nxt == fst)
                      break;
                  } while (1);
                }
              }
            }
          }
        }
        // sort (dist.begin(), dist.end());
        // for(i=(dist.size()-1);i>=0;i--)
        sort(dist.begin(), dist.end(), greater<float_in>());
        long unsigned int ii;
        for (ii = 0; ii < dist.size(); ii++) {
          npx2 = dist.size() - ii + npart;
          d = dist[ii];
          radx = d.d;

          test_arr[4 * sphsize + 3] = radx;
          Arvo_class *arvo_obj = new Arvo_class();
          try {
            // calculo el volumen
            arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
          } catch (const std::bad_alloc &e) {
            arvo_obj->Cleanup();
            exit(-1);
          }
          delete (arvo_obj);
          deltax = ((float)npx2 / vol) / density - 1.0;

          if (deltax < limit) {
            if (radx < min_radius) {
              delete[] test_arr;
              return (false); // no sirve el centro la esfera es menor al radio
                              // minimo
            }
            res.d = radx;
            res.cnt = npx2;
            res.id = d.id;
            vol_new = vol;
            break; // happy ending 2
          }
        }
        // cargando las nuevas particulas internas
        for (long unsigned int jj = ii; jj < dist.size(); jj++) {
          (*in_halos).push_back(dist[jj].id);
        }
      }
      break; // loop1
    }
  } // loop1

  ////////////////////work in progress//////////////////////////////////

  if (res.id < 0) {
    delete[] test_arr;
    return (false); // no sirve el centro
  } else {
    if (res.d <= 0.0) {
      cout << "error in inflate..." << endl;
      cout.flush();
      exit(-1);
    }

    (*volret) = vol_new;   // actualizo  el volumen
    (*sph).radius = res.d; // inflando hasta la ultima particula antes de superar el limite
    (*nhret) = res.cnt;    // actualizando numero de particulas

    delete[] test_arr;
    return (true);
  }
}

bool inflate_fast_old(Sphere *sph, class popcorn *pope,
                      // Grid* g,cell_state **status_rad,int nstep,
                      Grid *g, int nstep, vector<bool> *visit, vector<int> *in_halos, double *volret, int *nhret,
                      float limit, float density, vector<halo> &Tracer, float min_radius) {
  int ind, i, j, k, indx;
  long long int it, jt, kt, fst, nxt;
  vector<Sphere>::iterator itpop;
  vector<double> test_pop;
  int ns;
  unsigned long long int npx, npx2;
  float dd, rad, radx;
  float deltax, delta_inf;
  double area, vol;
  int npart, counter;
  float vol_new = -1.0;
  vector<class cell>::iterator ucell;
  int sphsize;
  double *test_arr;
  float_in res, d;

  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  halo cent;
  // int ntotal = (*g).nsearch;
  // int dim = 2 * ntotal + 1;
  vector<int> visitados;
  vector<float_in> dist;
  cell_state stat_out, stat_in, stat;

  cent.x = (*sph).x;
  cent.y = (*sph).y;
  cent.z = (*sph).z;

  npart = (*pope).npart; // numero de particulas que ya habia antes

  sphsize = ((int)(((*pope).membs).size()));
  test_arr = new double[4 * (sphsize + 1)];
  i = 0;
  for (itpop = ((*pope).membs).begin(); itpop != ((*pope).membs).end(); ++itpop) {
    test_arr[4 * i] = (double)(*itpop).x;
    test_arr[4 * i + 1] = (double)(*itpop).y;
    test_arr[4 * i + 2] = (double)(*itpop).z;
    test_arr[4 * i + 3] = (double)(*itpop).radius;
    i++;
  }
  /// esfera a pruena
  test_arr[4 * i] = (double)(*sph).x;
  test_arr[4 * i + 1] = (double)(*sph).y;
  test_arr[4 * i + 2] = (double)(*sph).z;
  test_arr[4 * i + 3] = 1e-5;

  check_boundary_splits(test_arr, sphsize + 1, boxsize);

  // calculo de la densidad en esferas crecientes primero de manera cruda para
  // identificar una escala aproximada de inflado
  int ic = int(cent.x / boxsize.x * Ngrid);
  int jc = int(cent.y / boxsize.y * Ngrid);
  int kc = int(cent.z / boxsize.z * Ngrid);

  res.d = -1.0;
  res.cnt = -1;
  res.id = -1;
  for (int irad = nstep - 1; irad >= 0; irad--) // loop1
  {
    counter = 0;

    rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
    float radg = rad / (*g).gridsize;
    ns = (int)ceil(rad / (*g).gridsize) + 1;

    for (k = (-ns); k < (ns + 1); k++) // k1
    {
      for (j = (-ns); j < (ns + 1); j++) // j1
      {
        for (i = (-ns); i < (ns + 1); i++) // i1
        {
          // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
          // stat=status_rad[nstep-irad-1][indi];
          stat = is_there(i, j, k, radg);

          it = indx_period(i + ic, Ngrid);
          jt = indx_period(j + jc, Ngrid);
          kt = indx_period(k + kc, Ngrid);

          ind = it + Ngrid * (jt + Ngrid * kt);
          if (stat.is_in) // para la couta inferior basta con los grides
                          // totalmente contenidos
          {
            fst = (*g).first[ind];
            if (fst >= 0) // grid no vacio
            {
              nxt = fst;
              do {
                if (!(*visit)[nxt]) // que no sea ya interior al pop
                {
                  counter++;
                }
                nxt = (*g).next[nxt];
                if (nxt == fst)
                  break;
              } while (1);
            }
          }
        } // i1
      }   // j1
    }     // k1

    test_arr[4 * (sphsize) + 3] = rad;
    Arvo_class *arvo_obj = new Arvo_class();
    try {
      // calculo el volumen
      arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
    } catch (const std::bad_alloc &e) {
      arvo_obj->Cleanup();
      exit(-1);
    }
    delete (arvo_obj);

    // cota inferior de la densidad
    delta_inf = ((float)(counter + npart) / vol) / density - 1.0;
    if (delta_inf < limit) {
      indx = irad;
      // calculamos cuanto es la delta exacta al radio correspondiente a indx,
      // si es mayor diminuimos y volvemos a calcular, el bucle termina cuando
      // tenemos una cascara donde la delta es menor al limite (y por lo tanto
      // cascara inmediata sup es mayor) o cuando llegamos al ultimo cascaron y
      // la delta sigue siendo superior

      do // loop 2
      {
        rad = (*g).maxscale * (float)(indx + 1) / (float)nstep;
        float radg = rad / (*g).gridsize;
        ns = (int)ceil(rad / (*g).gridsize) + 1;
        visitados.clear();

        npx = 0;
        for (k = (-ns); k < (ns + 1); k++) // k2
        {
          kt = indx_period(k + kc, (*g).Ngrid);
          for (j = (-ns); j < (ns + 1); j++) // j2
          {
            jt = indx_period(j + jc, (*g).Ngrid);
            for (i = (-ns); i < (ns + 1); i++) // i2
            {
              it = indx_period(i + ic, (*g).Ngrid);
              ind = it + (*g).Ngrid * (jt + (*g).Ngrid * kt);
              /// indi= i  + ntotal +dim*((j+ntotal)+dim*(k+ntotal));
              // stat=status_rad[nstep-indx-1][indi]; //is_there(i,j,k,rad, g);
              stat = is_there(i, j, k, radg);
              if (!stat.is_out) {
                fst = (*g).first[ind];
                if (fst >= 0) {
                  nxt = fst;
                  do {
                    if (!(*visit)[nxt]) {
                      if (stat.is_in) {
                        npx++;
                        visitados.push_back(nxt);
                      } else {
                        dd = distance(Tracer[nxt], cent, boxsize);
                        if (dd < rad) {
                          npx++;
                          visitados.push_back(nxt);
                        }
                      }
                    }
                    nxt = (*g).next[nxt];
                    if (nxt == fst)
                      break;
                  } while (1);
                }
              }

            } // i2
          }   // j2
        }     // k2

        test_arr[4 * sphsize + 3] = rad;
        Arvo_class *arvo_obj = new Arvo_class();
        try {
          // calculo el volumen
          arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
        } catch (const std::bad_alloc &e) {
          arvo_obj->Cleanup();
          exit(-1);
        }
        delete (arvo_obj);

        // delta exacta al radio indx
        deltax = ((float)(npx + npart) / vol) / density - 1.0;

        if (deltax < limit) {
          // dado que estamos seguros que la cota minima de indx+1 es superior
          // al limite (asi funciona el bucle anterior) y sabemos que en este
          // indx estamos debajo, el limite se cruza en un radio intermedio
          // busca entre los radios correspondientes a los cascarones de indx y
          // indx+1, el indice de la particula que da el cambio de densidad para
          // eso utiliza solo los boxcell que debe abrir segun la desigualdad
          // del triangulo (como fast_counter hace)
          // res=shell_finder(cent,limit,density,indx,nstep,g,&Tracer,status_rad);

          // shell_finder es la version esferica, debemos tener en cuenta la
          // geometria del pop
          dist.clear();
          float radsup = (*g).maxscale * (float)(indx + 2) / (float)nstep;
          float radinf = (*g).maxscale * (float)(indx + 1) / (float)nstep;
          float radsupg = radsup / (*g).gridsize;
          float radinfg = radinf / (*g).gridsize;
          for (k = (-ns); k < (ns + 1); k++) // k3
          {
            for (j = (-ns); j < (ns + 1); j++) // j3
            {
              for (i = (-ns); i < (ns + 1); i++) // i3
              {
                it = indx_period(i + ic, Ngrid);
                jt = indx_period(j + jc, Ngrid);
                kt = indx_period(k + kc, Ngrid);

                ind = it + Ngrid * (jt + Ngrid * kt);

                fst = (*g).first[ind];

                if (fst >= 0) {
                  // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
                  // en el máximo radio de inflado ya estoy por debajo de la
                  // delta
                  if (indx == nstep - 1) {
                    // stat_in =status_rad [nstep- indx   -1][indi];
                    stat_in = is_there(i, j, k, radinfg);
                    // opero los logicos para asegurarme que se abra el grid
                    // para todas las particulas internas y asi se carguen
                    stat_in.is_in = !stat_in.is_in;
                    stat_out.is_out = false;
                  } else {
                    // stat_out=status_rad [nstep-(indx+1)-1][indi];
                    stat_out = is_there(i, j, k, radsupg);
                    // stat_in =status_rad [nstep- indx   -1][indi];
                    stat_in = is_there(i, j, k, radinfg);
                  }
                  if (!(stat_out.is_out) && !(stat_in.is_in)) {
                    // open grid
                    nxt = fst;
                    do {
                      if (!(*visit)[nxt]) {
                        d.d = distance(Tracer[nxt], cent, boxsize);
                        d.id = nxt;
                        if ((d.d < radsup) && (d.d >= radinf)) {
                          dist.push_back(d);
                        }
                      }
                      nxt = (*g).next[nxt];
                      if (nxt == fst)
                        break;
                    } while (1);
                  }
                }

              } // i3
            }   // j3
          }     // k3

          // sort (dist.begin(), dist.end());
          // for(i=(dist.size()-1);i>=0;i--)
          sort(dist.begin(), dist.end(), greater<float_in>());
          long unsigned int ii;
          for (ii = 0; ii < dist.size(); ii++) {
            npx2 = (dist.size() - ii) + npx + npart;
            d = dist[ii];
            radx = d.d;

            test_arr[4 * sphsize + 3] = radx;
            Arvo_class *arvo_obj = new Arvo_class();
            try {
              // calculo el volumen
              arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
            } catch (const std::bad_alloc &e) {
              arvo_obj->Cleanup();
              exit(-1);
            }
            delete (arvo_obj);
            deltax = ((float)npx2 / vol) / density - 1.0;

            if (deltax < limit) {
              // happy ending 1
              res.d = radx;
              res.cnt = npx2;
              res.id = d.id;
              vol_new = vol;

              if (radx < min_radius) {
                delete[] test_arr;
                return (false); // no sirve el centro la esfera es menor al
                                // radio minimo
              }
              break;
            }
          }
          // cargando las nuevas particulas internas
          // cargando la cascara
          for (long unsigned int jj = ii; jj < dist.size(); jj++) {
            (*in_halos).push_back(dist[jj].id);
          }

          // cargando dentro del radio inferior
          for (long unsigned int jj = 0; jj < visitados.size(); jj++) {
            (*in_halos).push_back(visitados[jj]);
          }

          break; // loop2 //happy ending
        } else {
          indx--;
          if ((*g).maxscale * (float)(indx + 2) / (float)nstep < min_radius) {
            break;
          }
        }
      } while (indx >= 0); // loop2

      if (indx == (-1)) {
        rad = (*g).maxscale / (float)nstep;
        float radg = rad / (*g).gridsize;
        if (rad < min_radius)
          break;
        // res=ball_finder(cent,limit,density,0,nshell,g,&Tracer,status_rad);
        // hacemos lo mismo pero teniendo en cuenta la geometria del pop
        ns = (int)ceil(rad / (*g).gridsize) + 1;
        dist.clear();
        for (k = (-ns); k < (ns + 1); k++) // k4
        {
          for (j = (-ns); j < (ns + 1); j++) // j4
          {
            for (i = (-ns); i < (ns + 1); i++) // i4
            {
              it = indx_period(i + ic, Ngrid);
              jt = indx_period(j + jc, Ngrid);
              kt = indx_period(k + kc, Ngrid);

              ind = it + Ngrid * (jt + Ngrid * kt);

              fst = (*g).first[ind];

              if (fst >= 0) {
                // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
                // stat=status_rad[nstep-irad-1][indi]; //is_there(i,j,k,rad,
                // g);
                stat = is_there(i, j, k, radg);
                if (!(stat.is_out)) {
                  // open grid
                  nxt = fst;
                  do {
                    if (!(*visit)[nxt]) {
                      d.d = distance(Tracer[nxt], cent, boxsize);
                      d.id = nxt;
                      if (d.d < rad) {
                        dist.push_back(d);
                      }
                    }
                    nxt = (*g).next[nxt];
                    if (nxt == fst)
                      break;
                  } while (1);
                }
              }
            }
          }
        }
        // sort (dist.begin(), dist.end());
        // for(i=(dist.size()-1);i>=0;i--)
        sort(dist.begin(), dist.end(), greater<float_in>());
        long unsigned int ii;
        for (ii = 0; ii < dist.size(); ii++) {
          npx2 = dist.size() - ii + npart;
          d = dist[ii];
          radx = d.d;

          test_arr[4 * sphsize + 3] = radx;
          Arvo_class *arvo_obj = new Arvo_class();
          try {
            // calculo el volumen
            arvo_obj->compute_volume_arvo(&area, &vol, test_arr, sphsize + 1);
          } catch (const std::bad_alloc &e) {
            arvo_obj->Cleanup();
            exit(-1);
          }
          delete (arvo_obj);
          deltax = ((float)npx2 / vol) / density - 1.0;

          if (deltax < limit) {
            if (radx < min_radius) {
              delete[] test_arr;
              return (false); // no sirve el centro la esfera es menor al radio
                              // minimo
            }
            res.d = radx;
            res.cnt = npx2;
            res.id = d.id;
            vol_new = vol;
            break; // happy ending 2
          }
        }
        // cargando las nuevas particulas internas
        for (long unsigned int jj = ii; jj < dist.size(); jj++) {
          (*in_halos).push_back(dist[jj].id);
        }
      }
      break; // loop1
    }
  } // loop1

  ////////////////////work in progress//////////////////////////////////

  if (res.id < 0) {
    delete[] test_arr;
    return (false); // no sirve el centro
  } else {
    if (res.d <= 0.0 || vol_new < 0.0) {
      cout << "error in inflate..." << endl;
      cout.flush();
      exit(-1);
    }

    (*volret) = vol_new;   // actualizo  el volumen
    (*sph).radius = res.d; // inflando hasta la ultima particula antes de superar el limite
    (*nhret) = res.cnt;    // actualizando numero de particulas

    delete[] test_arr;
    return (true);
  }
}

// int init_flags(Sphere *centre, int nstep, Grid* g, vector<halo>
// *pnts,cell_state **status_rad,vector <int> * in_halos)
int init_flags(Sphere *centre, Grid *g, vector<halo> *pnts, vector<int> *in_halos) {
  cell_state stat;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int((*centre).x / boxsize.x * Ngrid);
  int jc = int((*centre).y / boxsize.y * Ngrid);
  int kc = int((*centre).z / boxsize.z * Ngrid);
  float rad = (*centre).radius;
  float radg = (*centre).radius / (*g).gridsize;
  // int ntotal = (*g).nsearch;

  // int irad = (int)(rad / ((*g).maxscale) * nstep) - 1;
  int ns = (int)ceil(rad / gridsize) + 1;

  int it, jt, kt, ind;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  long long int fst, nxt;

  float d;
  vector<float_in> dist;
  halo cent;
  cent.x = (*centre).x;
  cent.y = (*centre).y;
  cent.z = (*centre).z;

  int counter = 0;

  for (k = (-ns); k < (ns + 1); k++) {
    for (j = (-ns); j < (ns + 1); j++) {
      for (i = (-ns); i < (ns + 1); i++) {
        it = indx_period(i + ic, Ngrid);
        jt = indx_period(j + jc, Ngrid);
        kt = indx_period(k + kc, Ngrid);

        ind = it + Ngrid * (jt + Ngrid * kt);

        fst = (*g).first[ind];

        if (fst >= 0) {
          // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
          // stat=status_rad[nstep-irad-1][indi]; //is_there(i,j,k,rad, g);
          stat = is_there(i, j, k, radg);
          if (!(stat.is_out)) {
            // open grid
            nxt = fst;
            do {
              if (stat.is_in) {
                (*in_halos).push_back(nxt);
                counter++;
              } else {
                d = distance((*pnts)[nxt], cent, boxsize);
                if (d < rad) {
                  (*in_halos).push_back(nxt);
                  counter++;
                }
              }
              nxt = (*g).next[nxt];
              if (nxt == fst)
                break;
            } while (1);
          }
        }
      }
    }
  }

  return (counter);
}
} // namespace grid
int iper(int i, int N) {
  int ip;

  if (i >= N) {
    ip = i - N;
  } else if (i < 0) {
    ip = i + N;
  } else {
    ip = i;
  }

  return ip;
}
