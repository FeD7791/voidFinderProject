#include "finder.h"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "allvars.h"
#include "grid.h"
#include "objects.h"
#include "t.h"

using namespace std;
using namespace grid;

void svf_work(int wrk, float limit, float rhobar, int nshell, Grid *G, float min_radius) {
  // fusion de
  // FindCenters_fast(limit, rhobar, nshell, G, min_radius);
  // FindSphvds_fastV2(limit, rawfname, nshell, G, min_radius);
  //
  int Ngrid = (*G).Ngrid;
  float gridsize = (*G).gridsize;
  vector<int> ii;
  Sphere cndt;
  float Radius;
  float xr[NRAND][3];
  float Volume, Delta;
  vector<float_in> resran;
  int indx, ns;

  ii = revertindex(wrk, Ngrid);
  int res_irad = raw_finder(ii[0], ii[1], ii[2], limit, rhobar, nshell, G, min_radius);

  if (res_irad >= 0) {
    if (res_irad >= (nshell - 1))
      res_irad = nshell - 2;
    cndt.x = ((float)ii[0] + 0.5) * gridsize;
    cndt.y = ((float)ii[1] + 0.5) * gridsize;
    cndt.z = ((float)ii[2] + 0.5) * gridsize;
    cndt.radius = 0.0;
    cndt.irad = res_irad;
    cndt.ToF = false;
  }

  for (int iran = 0; iran < NRAND; iran++) {
    xr[iran][0] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
    xr[iran][1] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
    xr[iran][2] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
  }

  float_in res;
  res.d = -1.0;
  res.cnt = -1;
  res.id = -1;

  resran.clear();
  for (int iran = 0; iran < NRAND; iran++)
    resran.push_back(res);

  // ahora indx es calculado a la hora de buscar el cndt
  indx = cndt.irad;
  // nos da la mayor cáscara cuya cota inferior esta por debajo del limite

  // ahora indx es calculado a la hora de buscar el cndt
  // indx=raw_finder(centre,limit,NMEAN,nshell,g,pnts,status_rad);

  if (indx > 0) {
    // calculamos cuanto es la delta exacta al radio correspondiente a indx, si
    // es mayor diminuimos y volvemos a calcular, el bucle termina cuando
    // tenemos una cascara donde la delta es menor al limite (y por lo tanto
    // cascara inmediata sup es mayor) o cuando llegamos al ultimo cascaron y la
    // delta sigue siendo superior
    do {
      int ic = int(cndt.x / boxsize.x * (*G).Ngrid);
      int jc = int(cndt.y / boxsize.y * (*G).Ngrid);
      int kc = int(cndt.z / boxsize.z * (*G).Ngrid);

      float rad = (*G).maxscale * (float)(indx + 1) / (float)nshell;
      float vol = (4.0 / 3.0 * M_PI * rad * rad * rad);
      float radg = rad / (*G).gridsize;
      ns = (int)ceil(rad / (*G).gridsize) + 1;

#pragma omp parallel for default(none) schedule(dynamic)                                                               \
    shared(resran, ns, G, rad, vol, radg, ic, jc, kc, limit, boxsize, cndt, xr, Tracer, NMEAN, nshell, indx)
      for (int iran = 0; iran < NRAND; iran++) {
        halo center;
        center.x = xper(cndt.x + xr[iran][0], boxsize.x);
        center.y = xper(cndt.y + xr[iran][1], boxsize.y);
        center.z = xper(cndt.z + xr[iran][2], boxsize.z);
        int npx = 0;

        for (int k = (-ns); k < (ns + 1); k++) {
          for (int j = (-ns); j < (ns + 1); j++) {
            for (int i = (-ns); i < (ns + 1); i++) {
              int it = indx_period(i + ic, (*G).Ngrid);
              int jt = indx_period(j + jc, (*G).Ngrid);
              int kt = indx_period(k + kc, (*G).Ngrid);

              int ind = it + (*G).Ngrid * (jt + (*G).Ngrid * kt);

              int fst = (*G).first[ind];

              if (fst >= 0) {
                cell_state stat = is_there(i, j, k, radg);
                if (stat.is_in) {
                  npx = npx + (*G).cnt[ind];
                } else if (!(stat.is_out)) {
                  // open grid
                  int nxt = fst;
                  do {

                    float d = distance(Tracer[nxt], center, boxsize);
                    if (d < rad)
                      npx = npx + 1;
                    //}
                    nxt = (*G).next[nxt];
                    if (nxt == fst)
                      break;
                  } while (1);
                }
              }
            }
          }
        }

        float deltax = ((float)npx / vol) / NMEAN - 1.0;
        if (deltax < limit) {
          resran[iran] = shell_finder(center, limit, NMEAN, indx, nshell, G, &Tracer);
          resran[iran].id = iran;
        }
      }

      sort(resran.begin(), resran.end());
      if ((resran.back()).d > -1.0) {
        res = resran.back();
        break;
      } else {
        indx--;
      }
    } while (indx >= 0);

    if (indx == (-1)) {
      halo center;
      center.x = cndt.x;
      center.y = cndt.y;
      center.z = cndt.z;
      // llegamos al ultimo cascaron sin caer debajo del limite, vamos a
      // calcular dentro de la bola y vemos si encontramos algo o lo descartamos
      res = ball_finder(center, limit, NMEAN, 0, nshell, G, &Tracer);
      res.id = -1;
    }

    if (res.cnt > 0.0 && res.d > min_radius) {
      Radius = res.d;
      Volume = (4.0 / 3.0) * PI * Radius * Radius * Radius;
      Delta = (float)(res.cnt) / Volume / NMEAN - 1.0;
      int iran = res.id;

      Sphere cndt2;
      cndt2.radius = Radius;
      cndt2.delta = Delta;
      cndt2.np = res.cnt;
      cndt2.ToF = true;

      if (iran == -1) {
        cndt2.x = cndt.x;
        cndt2.y = cndt.y;
        cndt2.z = cndt.z;
      } else {
        cndt2.x = xper(cndt.x + xr[iran][0], boxsize.x);
        cndt2.y = xper(cndt.y + xr[iran][1], boxsize.y);
        cndt2.z = xper(cndt.z + xr[iran][2], boxsize.z);
      }
#pragma omp critical
      { Sphvd.push_back(cndt2); }
    }
  }
}

// int raw_finder(int ic,int jc,int kc, float limit,float rhobar,int nstep,
// Grid* g, cell_state **status_rad,float rad_min)
int raw_finder(int ic, int jc, int kc, float limit, float rhobar, int nstep, Grid *g, float rad_min) {
  int Ngrid = (*g).Ngrid;
  float rad;
  int ntotal = (*g).nsearch;
  float delta_inf, vol;
  int indx = -1;

  bool flag = true;

  int *counter = new int[nstep];
  for (int irad = 0; irad < nstep; irad++)
    counter[irad] = 0;

  float fct = (*g).maxscale / (float)nstep / (*g).gridsize;

#pragma omp parallel for default(none) shared(nstep, ntotal, Ngrid, kc, jc, ic, fct, g) reduction(+ : counter[ : nstep])
  for (int k = (-ntotal); k < (ntotal + 1); k++) {
    int kt = indx_period(k + kc, Ngrid);
    int kk = abs(k);
    for (int j = (-ntotal); j < (ntotal + 1); j++) {
      int jt = indx_period(j + jc, Ngrid);
      int jj = abs(j);
      for (int i = (-ntotal); i < (ntotal + 1); i++) {
        int it = indx_period(i + ic, Ngrid);
        int ii = abs(i);
        int ind = it + Ngrid * (jt + Ngrid * kt);
        int nn = (*g).cnt[ind];
        if (nn > 0) {
          int caso = (int)(ii != 0) + (int)(jj != 0) + (int)(kk != 0);
          float val = (float)(ii * (ii + 1) + jj * (jj + 1) + kk * (kk + 1));
          for (int irad = 0; irad < nstep; irad++) {
            float radg = fct * ((float)(irad + 1));
            if ((val + radg * (sqrt((float)caso) - radg)) < 0.0) {
              counter[irad] = counter[irad] + nn;
            }
          }
        }
      }
    }
  }

  for (int irad = nstep - 1; irad >= 0; irad--) {
    rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
    vol = (4.0 / 3.0 * M_PI * rad * rad * rad);
    delta_inf = ((float)counter[irad]) / vol;
    delta_inf = delta_inf / rhobar - 1.0;
    if (delta_inf < limit && rad >= rad_min) {
      indx = irad;
      flag = false;
      break;
    }
  }
  delete[] counter;

  // if(indx<=0 && (delta_inf > limit))indx=-1;
  if (flag)
    indx = -1;
  return (indx);
}

// void FindCenters_fast(float limseed,float rhobar,int nstep, Grid *G,
// vector<halo> *pnts, cell_state **status_rad,float min_radius)
void FindCenters_fast(float limseed, float rhobar, int nstep, Grid *G, float min_radius) {
  clock_t t;
  int Ngrid = (*G).Ngrid;
  float gridsize = (*G).gridsize;
  int res;
  int indx;
  vector<int> ii;
  Sphere cndt;

  fprintf(stdout, "\n Busqueda de centros subdensos\n");
  fflush(stdout);
  t = clock();

  int Ngrid3 = Ngrid * Ngrid * Ngrid;

  fprintf(stdout, " | Ngrid**3 = %d \n", Ngrid3);

#pragma omp declare reduction(merge : std::vector<Sphere> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  int cnt = 0;
#pragma omp parallel for default(none) reduction(merge : seed)                                                         \
    shared(Ngrid, Ngrid3, gridsize, limseed, rhobar, nstep, G, min_radius, cnt) private(indx, ii, cndt, res)
  for (indx = 0; indx < Ngrid3; indx++) {

    if (indx % 1000 == 0) {
#pragma omp critical
      {
        cnt++;
        Progress(cnt, Ngrid3 / 1000);
      }
    }

    ii = revertindex(indx, Ngrid);
    res = raw_finder(ii[0], ii[1], ii[2], limseed, rhobar, nstep, G, min_radius);

    if (res >= 0) {
      if (res >= (nstep - 1))
        res = nstep - 2;
      cndt.x = ((float)ii[0] + 0.5) * gridsize;
      cndt.y = ((float)ii[1] + 0.5) * gridsize;
      cndt.z = ((float)ii[2] + 0.5) * gridsize;
      // cndt.Ref = 0.0;
      cndt.radius = 0.0;
      cndt.irad = res;
      cndt.ToF = false;
      seed.push_back(cndt);
    }
  }

  NVOID = seed.size();

  fprintf(stdout, " | Candidatos a void = %d \n", NVOID);
  fflush(stdout);
  Time(t, 1);
}

void FindCenters(float limseed) {
  int i;
  clock_t t;

  fprintf(stdout, "\n Busqueda de centros subdensos\n");
  fflush(stdout);
  t = clock();

  NVOID = 0;
  for (i = 0; i < NTRAC; i++) {
    if (Tessell[i].Delta <= limseed) {
      Sphvd.push_back(Sphere());
      Sphvd[NVOID].x = Tessell[i].Cen[0];
      Sphvd[NVOID].y = Tessell[i].Cen[1];
      Sphvd[NVOID].z = Tessell[i].Cen[2];
      /// Sphvd[NVOID].Ref = pow(0.75*Tessell[i].Volume/PI,1.0/3.0);
      Sphvd[NVOID].radius = 0.0;
      Sphvd[NVOID].ToF = false;
      NVOID++;
    }
  }

  fprintf(stdout, " | Candidatos a void = %d \n", NVOID);
  fflush(stdout);
  Time(t, 1);
}

// void FindSphvds_fast(float limit,string fname,int nshell, Grid* G, cell_state
// **status_rad)
void FindSphvds_fast(float limit, string fname, int nshell, Grid *G) {
  int iran;
  int iv, check;
  float Radius, RadiusMax;
  float xr[3], xc[3];
  float Volume, Delta;

  clock_t t;
  FILE *fd;
  halo center;
  float_in res;
  float gridsize = (*G).gridsize;
  Sphere cndt;

  fprintf(stdout, "\n Identificacion de voids \n");
  fflush(stdout);
  t = clock();

  srand(time(NULL));

  // Selecciono vecinos
#pragma omp declare reduction(merge : std::vector<Sphere> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  fprintf(stdout, " | NVOID = %d \n", NVOID);

  int cnt = 0;
#pragma omp parallel for default(none) schedule(dynamic) reduction(merge : Sphvd)                                      \
    shared(NVOID, NMEAN, seed, Tracer, stdout, radmax, boxsize, cout, limit, G, nshell,                                \
               gridsize) private(center, iv, xc, xr, check, Radius, RadiusMax, iran, cnt, res, Volume, Delta, cndt)
  for (iv = 0; iv < NVOID; iv++) {

    if (omp_get_thread_num() == 0)
      Progress(iv, NVOID);

    cnt++;

    cndt.radius = -1.0;
    cndt.delta = 1000.0;
    cndt.np = -1;
    cndt.ToF = false;

    RadiusMax = 0.0;
    check = 0;
    for (iran = 0; iran <= NRAND; iran++) {

      check++;
      if (check == NRAND / 4)
        break;

      if (iran == 0) {
        for (int k = 0; k < 3; k++)
          xr[k] = 0.0;

      } else {
        xr[0] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
        xr[1] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
        xr[2] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
      }

      xc[0] = xper(seed[iv].x + xr[0], boxsize.x);
      xc[1] = xper(seed[iv].y + xr[1], boxsize.y);
      xc[2] = xper(seed[iv].z + xr[2], boxsize.z);

      center.x = xc[0];
      center.y = xc[1];
      center.z = xc[2];
      res = fast_finder(center, limit, NMEAN, nshell, G, &Tracer, seed[iv].irad);
      if (res.cnt > 0.0) {
        Radius = res.d;
        Volume = (4.0 / 3.0) * PI * Radius * Radius * Radius;
        Delta = (float)(res.cnt) / Volume / NMEAN - 1.0;
        if (iran == 0) {
          RadiusMax = Radius;
          cndt.radius = Radius;
          cndt.delta = Delta;
          cndt.np = res.cnt;
          cndt.x = xc[0];
          cndt.y = xc[1];
          cndt.z = xc[2];
          cndt.ToF = true;
          check = 0;
        } else {
          if (Radius > RadiusMax) {
            cndt.radius = Radius;
            cndt.delta = Delta;
            cndt.np = res.cnt;
            cndt.x = xc[0];
            cndt.y = xc[1];
            cndt.z = xc[2];
            cndt.ToF = true;
            RadiusMax = Radius;
            check = 0;
          }
        }
      }

    } /* Fin lazo random */
    Sphvd.push_back(cndt);

    // if (Sphvd[iv].ToF) {
    //   lambda = (4.0/3.0)*PI*pow(Sphvd[iv].Rad,3)*NMEAN;
    //}

  } /* Fin lazo voids */

  fd = fopen(fname.c_str(), "wb");
  for (long unsigned int p = 0; p < Sphvd.size(); p++) {
    if (Sphvd[p].ToF) {
      fwrite(&(Sphvd[p].id), sizeof(int), 1, fd);
      fwrite(&(Sphvd[p].radius), sizeof(float), 1, fd);
      // fwrite(&(Sphvd[p].Ref    ),sizeof(float),1,fd);
      fwrite(&(Sphvd[p].x), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].y), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].z), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].delta), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].ToF), sizeof(bool), 1, fd);
      fwrite(&(Sphvd[p].np), sizeof(int), 1, fd);
      fwrite(&(Sphvd[p].irad), sizeof(int), 1, fd);
    }
  }
  fclose(fd);

  cnt = 0;
  for (long unsigned int p = 0; p < Sphvd.size(); p++) {
    if (Sphvd[p].ToF)
      cnt++;
  }
  cout << "N ToF " << cnt << endl;

  Time(t, NCORES);
}

// void FindSphvds_fastV2(float limit,string fname,int nshell, Grid* G,
// cell_state **status_rad,float min_radius)
void FindSphvds_fastV2(float limit, string fname, int nshell, Grid *G, float min_radius) {
  int iran;
  int i, j, k, iv;
  float Radius;
  float xr[NRAND][3];
  int npran[NRAND];
  float Volume, Delta;
  float rad;
  int it, jt, kt, ind;
  int fst, nxt;
  float d;

  clock_t t;
  FILE *fd;
  halo center;
  float_in res;
  vector<float_in> resran;
  float gridsize = (*G).gridsize;
  Sphere cndt;
  int indx, ns;
  int ic, jc, kc;

  int npx;
  float radx, vol;
  float deltax;
  cell_state stat;

  fprintf(stdout, "\n Identificacion de voids \n");
  fflush(stdout);
  t = clock();

  srand(time(NULL));

  for (iran = 0; iran < NRAND; iran++) {
    xr[iran][0] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
    xr[iran][1] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
    xr[iran][2] = ((float)(rand()) / (float)RAND_MAX - 0.5) * gridsize;
  }

// Selecciono vecinos
#pragma omp declare reduction(merge : std::vector<Sphere> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  fprintf(stdout, " | NVOID = %d \n", NVOID);

  int cnt = 0;
#pragma omp parallel for default(none) schedule(dynamic) reduction(merge : Sphvd)                                      \
    shared(NVOID, NMEAN, seed, Tracer, stdout, radmax, boxsize, cout, limit, G, nshell, gridsize, xr,                  \
               min_radius) private(center, iv, Radius, iran, cnt, res, Volume, Delta, cndt, resran, indx, ic, jc, kc,  \
                                       rad, ns, npx, npran, i, j, k, it, jt, kt, ind, stat, fst, nxt, d, radx, vol,    \
                                       deltax)
  for (iv = 0; iv < NVOID; iv++) {
    if (omp_get_thread_num() == 0)
      Progress(iv, NVOID);

    res.d = -1.0;
    res.cnt = -1;
    res.id = -1;

    cnt++;

    cndt.radius = -1.0;
    cndt.delta = 1000.0;
    cndt.np = -1;
    cndt.ToF = false;

    resran.clear();
    for (iran = 0; iran < NRAND; iran++)
      resran.push_back(res);

    // ahora indx es calculado a la hora de buscar el seed
    indx = seed[iv].irad;
    // nos da la mayor cáscara cuya cota inferior esta por debajo del limite

    // ahora indx es calculado a la hora de buscar el seed
    // indx=raw_finder(centre,limit,NMEAN,nshell,g,pnts,status_rad);

    if (indx < 0)
      continue;

    // calculamos cuanto es la delta exacta al radio correspondiente a indx, si
    // es mayor diminuimos y volvemos a calcular, el bucle termina cuando
    // tenemos una cascara donde la delta es menor al limite (y por lo tanto
    // cascara inmediata sup es mayor) o cuando llegamos al ultimo cascaron y la
    // delta sigue siendo superior
    do {
      ic = int(seed[iv].x / boxsize.x * (*G).Ngrid);
      jc = int(seed[iv].y / boxsize.y * (*G).Ngrid);
      kc = int(seed[iv].z / boxsize.z * (*G).Ngrid);

      rad = (*G).maxscale * (float)(indx + 1) / (float)nshell;
      float radg = rad / (*G).gridsize;
      ns = (int)ceil(rad / (*G).gridsize) + 1;

      npx = 0;
      for (iran = 0; iran < NRAND; iran++)
        npran[iran] = 0;

      for (k = (-ns); k < (ns + 1); k++) {
        for (j = (-ns); j < (ns + 1); j++) {
          for (i = (-ns); i < (ns + 1); i++) {
            it = indx_period(i + ic, (*G).Ngrid);
            jt = indx_period(j + jc, (*G).Ngrid);
            kt = indx_period(k + kc, (*G).Ngrid);

            ind = it + (*G).Ngrid * (jt + (*G).Ngrid * kt);

            fst = (*G).first[ind];

            if (fst >= 0) {
              // indi=(i+(*G).nsearch)+dim*((j+(*G).nsearch)+dim*(k+(*G).nsearch));
              // stat=status_rad[nshell-indx-1][indi]; //is_there(i,j,k,rad, g);
              stat = is_there(i, j, k, radg);
              if (stat.is_in) {
                npx = npx + (*G).cnt[ind];
              } else if (!(stat.is_out)) {
                // open grid
                nxt = fst;
                do {
                  for (iran = 0; iran < NRAND; iran++) {
                    center.x = xper(seed[iv].x + xr[iran][0], boxsize.x);
                    center.y = xper(seed[iv].y + xr[iran][1], boxsize.y);
                    center.z = xper(seed[iv].z + xr[iran][2], boxsize.z);

                    d = distance(Tracer[nxt], center, boxsize);

                    if (d < rad)
                      npran[iran] = npran[iran] + 1;
                  }
                  nxt = (*G).next[nxt];
                  if (nxt == fst)
                    break;
                } while (1);
              }
            }
          }
        }
      }

      radx = (*G).maxscale * (float)(indx + 1) / (float)nshell;
      vol = (4.0 / 3.0 * M_PI * radx * radx * radx);

      for (iran = 0; iran < NRAND; iran++)

      {
        deltax = ((float)(npx + npran[iran]) / vol) / NMEAN - 1.0;

        if (deltax < limit) {
          center.x = xper(seed[iv].x + xr[iran][0], boxsize.x);
          center.y = xper(seed[iv].y + xr[iran][1], boxsize.y);
          center.z = xper(seed[iv].z + xr[iran][2], boxsize.z);

          // dado que estamos seguros que la cota minima de indx+1 es superior
          // al limite (asi funciona raw_finder y el bucle anterior) y sabemos
          // que en este indx estamos debajo, el limite se cruza en un radio
          // intermedio
          resran[iran] = shell_finder(center, limit, NMEAN, indx, nshell, G, &Tracer);
          resran[iran].id = iran; // para conservar el indice en el sort posterior
          // busca entre los radios correspondientes a los cascarones de indx y
          // indx+1, el indice de la particula que da el cambio de densidad para
          // eso utiliza solo los boxcell que debe abrir segun la desigualdad
          // del triangulo (como fast_counter hace)
        }
      }
      sort(resran.begin(), resran.end());
      if ((resran.back()).d > -1.0) {
        res = resran.back();
        break;
      } else {
        indx--;
      }
    } while (indx >= 0);

    if (indx == (-1)) {
      center.x = seed[iv].x;
      center.y = seed[iv].y;
      center.z = seed[iv].z;
      // llegamos al ultimo cascaron sin caer debajo del limite, vamos a
      // calcular dentro de la bola y vemos si encontramos algo o lo descartamos
      res = ball_finder(center, limit, NMEAN, 0, nshell, G, &Tracer);
      res.id = -1;
    }

    if (res.cnt > 0.0 && res.d > min_radius) {
      Radius = res.d;
      Volume = (4.0 / 3.0) * PI * Radius * Radius * Radius;
      Delta = (float)(res.cnt) / Volume / NMEAN - 1.0;
      iran = res.id;

      cndt.radius = Radius;
      cndt.delta = Delta;
      cndt.np = res.cnt;
      cndt.ToF = true;

      if (iran == -1) {
        cndt.x = seed[iv].x;
        cndt.y = seed[iv].y;
        cndt.z = seed[iv].z;
      } else {
        cndt.x = xper(seed[iv].x + xr[iran][0], boxsize.x);
        cndt.y = xper(seed[iv].y + xr[iran][1], boxsize.y);
        cndt.z = xper(seed[iv].z + xr[iran][2], boxsize.z);
      }

      Sphvd.push_back(cndt);
    }

  } /* Fin lazo voids */

  fd = fopen(fname.c_str(), "wb");
  for (long unsigned int p = 0; p < Sphvd.size(); p++) {
    if (Sphvd[p].ToF) {
      fwrite(&(Sphvd[p].id), sizeof(int), 1, fd);
      fwrite(&(Sphvd[p].radius), sizeof(float), 1, fd);
      // fwrite(&(Sphvd[p].Ref    ),sizeof(float),1,fd);
      fwrite(&(Sphvd[p].x), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].y), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].z), sizeof(float), 1, fd);
      // fwrite(&(Sphvd[p].Vel[0] ),sizeof(float),1,fd);
      // fwrite(&(Sphvd[p].Vel[1] ),sizeof(float),1,fd);
      // fwrite(&(Sphvd[p].Vel[2] ),sizeof(float),1,fd);
      fwrite(&(Sphvd[p].delta), sizeof(float), 1, fd);
      fwrite(&(Sphvd[p].ToF), sizeof(bool), 1, fd);
      fwrite(&(Sphvd[p].np), sizeof(int), 1, fd);
      fwrite(&(Sphvd[p].irad), sizeof(int), 1, fd);
    }
  }
  fclose(fd);

  cnt = 0;
  for (long unsigned int p = 0; p < Sphvd.size(); p++) {
    if (Sphvd[p].ToF)
      cnt++;
  }
  cout << "N ToF " << cnt << endl;

  Time(t, NCORES);
}

// void CleanSphvds_Tol0(Grid *G,float maxscale,cell_state **status_rad,int
// nstep)
std::vector<int> svf_touch_list(Sphere centre, Grid *g, std::vector<Sphere> *others, float maxscale) {
  cell_state stat;
  cell_state stat2;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  float radsearch = centre.radius + maxscale;
  int ns = (int)ceil(radsearch / gridsize) + 1;

  int it, jt, kt, ind;
  int i, j, k;
  int fst, nxt;
  std::vector<float_in> dist;
  float d;
  std::vector<int> res;

  float radsearchg = radsearch / (*g).gridsize;

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
          // stat=status_rad[nstep-irad-1][indi]; //stat=is_there(i,j,k,radsearch, g);
          stat = is_there(i, j, k, radsearchg);
          if (!(stat.is_out)) {
            // open grid
            nxt = fst;
            do {
              if (((*others)[nxt].radius > centre.radius)) {
                float jradg = ((*others)[nxt].radius + centre.radius) / ((*g).gridsize);
                // stat2=status_rad[nstep-jrad-1][indi]; //stat2=is_there(i,j,k,(*others)[nxt].Rad+centre.Rad, g);
                stat2 = is_there(i, j, k, jradg);
                if (!stat2.is_out) {
                  if (stat2.is_in) {
                    res.push_back(nxt);
                    break;
                  } else {
                    d = distance((*others)[nxt], centre, boxsize);
                    if (d < ((*others)[nxt].radius + centre.radius)) {
                      res.push_back(nxt);
                      break;
                    }
                  }
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
  return (res);
}

void CleanSphvds(Grid *G, float maxscale) {
  int ID;
  clock_t t;
  float_in dv;

  struct small_first {
    bool operator()(par const &a, par const &b) const { return a.vol1 < b.vol1; }
  };

  fprintf(stdout, "\n Limpiando voids por superposicion\n");
  fflush(stdout);
  fprintf(stdout, " | Tolerancia de interseccion = 0\n");
  fflush(stdout);
  t = clock();

  fprintf(stdout, " | Computing intersections and cleaning\n");
  fflush(stdout);

  std::vector<par> svf_pairs;
  par para;
#pragma omp declare reduction(merge : std::vector<par> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

#pragma omp parallel for default(none) reduction(merge : svf_pairs) shared(Sphvd, G, maxscale, cout) private(para)
  for (long unsigned int ii = 0; ii < Sphvd.size(); ii++) {
    if (Sphvd.size() > 1000) {
      Progress(ii, Sphvd.size());
    } else {
      cout << ii << " / " << Sphvd.size() << endl;
      cout.flush();
    }

    std::vector<int> list;
    list = svf_touch_list(Sphvd[ii], G, &Sphvd, maxscale);
    for (long unsigned int jj = 0; jj < list.size(); jj++) {
      int cur = list[jj];
      int idmenor = Sphvd[ii].radius < Sphvd[cur].radius ? ii : cur;
      int idmayor = Sphvd[ii].radius < Sphvd[cur].radius ? cur : ii;

      para.idmin = idmenor;
      para.idmax = idmayor;
      para.vol1 = Sphvd[idmenor].radius;
      para.vol2 = Sphvd[idmayor].radius;
      svf_pairs.push_back(para);
    }
    list.clear();
  }

  fprintf(stdout, " | Sorting...\n");
  fflush(stdout);
  std::sort(svf_pairs.begin(), svf_pairs.end(), small_first());

  for (long unsigned int i = 0; i < svf_pairs.size(); i++) {
    if (!Sphvd[svf_pairs[i].idmin].ToF || !Sphvd[svf_pairs[i].idmax].ToF)
      continue;
    Sphvd[svf_pairs[i].idmin].ToF = false;
  }
  ID = 1;
  for (long unsigned int i = 0; i < Sphvd.size(); i++) {
    if (!Sphvd[i].ToF)
      continue;
    Sphvd[i].id = ID;
    ID++;
  }

  Time(t, 1);
}

void CleanSphvds_Tol0(Grid *G, float maxscale) {
  int ID;
  clock_t t;
  float_in dv;
  vector<float_in> szs;

  fprintf(stdout, "\n Limpiando voids por superposicion\n");
  fflush(stdout);
  fprintf(stdout, " | Tolerancia de interseccion = 0\n");
  fflush(stdout);
  t = clock();

  // sortin by sizes
  fprintf(stdout, " | Sorting by size\n");
  fflush(stdout);
  for (long unsigned int i = 0; i < Sphvd.size(); i++) {
    dv.d = Sphvd[i].radius;
    dv.id = i;
    szs.push_back(dv);
  }
  sort(szs.begin(), szs.end());

  fprintf(stdout, " | Computing intersections and cleaning\n");
  fflush(stdout);

  // if there are too many voids, lets make a clean litle bit inaccurate (5% in radius)
  if (Sphvd.size() > 500000) {
    int nstep = 1000;
    long unsigned int stepsize = (long unsigned int)((double)Sphvd.size() / (double)nstep);

    // bin in sizes
    for (int st = 0; st < nstep; st++) {
      Progress(st, nstep);
      long unsigned int from = stepsize * st;
      long unsigned int to = stepsize * (st + 1);
      long unsigned int it1 = szs[szs.size() - 1 - from].id;
      long unsigned int it2 = szs[szs.size() - 1 - to].id;
      double diff = (Sphvd[it1].radius - Sphvd[it2].radius) / Sphvd[it1].radius;

      // in this parallel loop we have a race condition, hopefully, only it will produce
      // a relative error of 5% in radius
      if (diff > 0.05)
        continue;
#pragma omp parallel for default(none) shared(szs, Sphvd, G, maxscale, it1, from, to)
      for (long unsigned int i = from; i < to; i++) {
        long unsigned int it = szs[szs.size() - 1 - i].id;
        if (Sphvd[it].ToF) {
          if (is_touched(Sphvd[it], G, &Sphvd, Sphvd[it1].radius)) {
#pragma omp critical
            { Sphvd[it].ToF = false; }
          }
        }
      }
    }
  }

  for (long unsigned int i = 0; i < szs.size(); i++) {
    if (Sphvd.size() > 1000) {
      if (i % 100 == 0) {
        Progress(i, szs.size());
      }
    } else {
      cout << i << " / " << szs.size() << endl;
      cout.flush();
    }
    long unsigned int it = szs[szs.size() - 1 - i].id;
    if (Sphvd[it].ToF) {
      if (is_touched(Sphvd[it], G, &Sphvd, maxscale)) {
        Sphvd[it].ToF = false;
      }
    }
  }

  ID = 1;
  for (long unsigned int i = 0; i < Sphvd.size(); i++) {
    if (!Sphvd[i].ToF)
      continue;
    Sphvd[i].id = ID;
    ID++;
  }

  Time(t, 1);
}

void CleanSphvds_Tol0_opt(Grid *G, float maxscale) {
  unsigned long int i, j, it, ic, ID;
  clock_t t;
  float_in dv;
  vector<float_in> szs;

  fprintf(stdout, "\n Limpiando voids por superposicion\n");
  fflush(stdout);
  fprintf(stdout, " | Tolerancia de interseccion = 0 \n");
  fflush(stdout);
  t = clock();

  // sortin by sizes
  fprintf(stdout, " | Sorting by size\n");
  fflush(stdout);
  for (i = 0; i < Sphvd.size(); i++) {
    dv.d = Sphvd[i].radius;
    dv.id = i;
    szs.push_back(dv);
  }
  sort(szs.begin(), szs.end());

  fprintf(stdout, " | Computing intersections and cleaning\n");
  fflush(stdout);
#pragma omp parallel for default(none) shared(szs, Sphvd, G, maxscale, cout) private(i, it)
  for (i = 0; i < szs.size(); i++) {
    it = szs.size() - 1 - i;
    szs[i].vec = touch_list(Sphvd[szs[it].id], G, &Sphvd, maxscale);
  }

  for (i = 0; i < szs.size(); i++)
    Sphvd[szs[i].id].ToF = true;

  for (i = 0; i < szs.size(); i++) {
    for (j = 0; j < szs[i].vec.size(); j++) {
      ic = szs[i].vec[j];
      if (Sphvd[ic].ToF) {
        Sphvd[szs[i].id].ToF = false;
        break;
      }
    }
  }

  ID = 1;
  for (i = 0; i < Sphvd.size(); i++) {
    if (!Sphvd[i].ToF)
      continue;
    Sphvd[i].id = ID;
    ID++;
  }

  Time(t, 1);
}

float xper(float x, float boxsize) {
  float xp;

  if (x >= boxsize) {
    xp = x - boxsize;
  } else if (x < 0.0) {
    xp = x + boxsize;
  } else {
    xp = x;
  }

  return xp;
}

// float logfactorial(int n)
//{
//  int   i;
//  float f;
//
//  if (n > 1) {
//     f = 0.0;
//     for (i=2; i<=n; i++)
//         f += log((double)i);
//  } else {
//     f = 0.0;
//  }
//
//  return f;
//}
