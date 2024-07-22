#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#define likely(x) __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)

namespace grid {
template <typename T> void Grid::build(T *rts) {
  long long int i, j, k;
  long long int ind;

  assert(Npart == (*rts).size());
  for (unsigned long int ii = 0; ii < Npart; ii++) {
    i = (long long int)((*rts)[ii].x / boxsize.x * Ngrid);
    if (unlikely(i > (Ngrid - 1)))
      i = Ngrid - 1;
    j = (long long int)((*rts)[ii].y / boxsize.y * Ngrid);
    if (unlikely(j > (Ngrid - 1)))
      j = Ngrid - 1;
    k = (long long int)((*rts)[ii].z / boxsize.z * Ngrid);
    if (unlikely(k > (Ngrid - 1)))
      k = Ngrid - 1;
    ind = i + Ngrid * (j + Ngrid * k);
    first[ind] = ii;
    cnt[ind]++;
  }

  for (unsigned long int ii = 0; ii < Npart; ii++) {
    i = (long long int)((*rts)[ii].x / boxsize.x * Ngrid);
    if (unlikely(i > (Ngrid - 1)))
      i = Ngrid - 1;
    j = (long long int)((*rts)[ii].y / boxsize.y * Ngrid);
    if (unlikely(j > (Ngrid - 1)))
      j = Ngrid - 1;
    k = (long long int)((*rts)[ii].z / boxsize.z * Ngrid);
    if (unlikely(k > (Ngrid - 1)))
      k = Ngrid - 1;
    ind = i + Ngrid * (j + Ngrid * k);
    next[first[ind]] = ii;
    first[ind] = ii;
  }
};

template <typename T> float distance(T p1, T p2, float3 boxsize) {
  float dx, dy, dz;
  float d;

  dx = fabs((p1.x - p2.x));
  dx = periodicity_delta(dx, boxsize.x);
  dy = fabs((p1.y - p2.y));
  dy = periodicity_delta(dy, boxsize.y);
  dz = fabs((p1.z - p2.z));
  dz = periodicity_delta(dz, boxsize.z);

  d = sqrt(dx * dx + dy * dy + dz * dz);
  return (d);
}

template <typename T> int fast_counter(T centre, int irad, int nstep, Grid *g, std::vector<T> *pnts) {
  cell_state stat;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  // int ntotal = (*g).nsearch;
  float rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
  float radg = rad / (*g).gridsize;
  int ns = (int)ceil(rad / gridsize) + 1;

  long long int it, jt, kt;
  long long int ind;
  // int indi;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  long long int fst, nxt;
  float d;

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
          if (stat.is_in) {
            counter = counter + (*g).cnt[ind];
          } else if (!(stat.is_out)) {
            // open grid
            nxt = fst;
            do {
              d = distance((*pnts)[nxt], centre, boxsize);
              if (d < rad)
                counter = counter + 1;
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

template <typename T> int slow_counter(T centre, float rad, Grid *g, std::vector<T> *pnts) {
  cell_state stat;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);
  // int ns=(*g).nsearch;
  int ns = (int)ceil(rad / gridsize) + 1;
  long long int it, jt, kt, ind;
  int i, j, k;

  int fst, nxt;
  float d;

  int counter = 0;

  for (k = (-ns); k < (ns + 1); k++) {
    for (j = (-ns); j < (ns + 1); j++) {
      for (i = (-ns); i < (ns + 1); i++) {
        it = indx_period(i + ic, Ngrid);
        jt = indx_period(j + jc, Ngrid);
        kt = indx_period(k + kc, Ngrid);

        ind = it + Ngrid * (jt + Ngrid * kt);

        // open grid
        fst = (*g).first[ind];
        nxt = fst;
        if (fst >= 0) {
          do {
            d = distance((*pnts)[nxt], centre, boxsize);
            if (d < rad)
              counter = counter + 1;
            nxt = (*g).next[nxt];
            if (nxt == fst)
              break;
          } while (1);
        }
      }
    }
  }
  return (counter);
}

template <typename T>
unsigned long long int *raw_counter(T centre, int irad, int nstep, Grid *g, std::vector<T> *pnts) {
  cell_state stat;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);
  // int ns=(*g).nsearch;

  float rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
  float radg = rad / (*g).gridsize;
  int ns = (int)ceil(rad / gridsize) + 1;
  long long int it, jt, kt, ind;
  int i, j, k;
  int ntotal = (*g).nsearch;
  int dim = 2 * ntotal + 1;

  int fst, nxt;
  float d;

  unsigned long long int counter = 0;
  unsigned long long int counter2 = 0;

  for (k = (-ns); k < (ns + 1); k++) {
    for (j = (-ns); j < (ns + 1); j++) {
      for (i = (-ns); i < (ns + 1); i++) {
        // indi=(i+ntotal)+dim*((j+ntotal)+dim*(k+ntotal));
        // stat=status_rad[nstep-irad-1][indi]; //is_there(i,j,k,rad, g);
        stat = is_there(i, j, k, radg);

        it = indx_period(i + ic, Ngrid);
        jt = indx_period(j + jc, Ngrid);
        kt = indx_period(k + kc, Ngrid);

        ind = it + Ngrid * (jt + Ngrid * kt);

        if (stat.is_in) {
          counter = counter + (*g).cnt[ind];
        } else if (!(stat.is_out)) {
          counter2 = counter2 + (*g).cnt[ind];
        }
      }
    }
  }
  unsigned long long int *res;
  res = new unsigned long long int[2];
  res[0] = counter;
  res[1] = counter2;

  return (res);
}

template <typename T>
float_in fast_finder(T centre, float limit, float rhobar, int nstep, Grid *g, std::vector<T> *pnts, int indx) {
  int npx;
  float radx, vol;
  float deltax;
  float_in res;

  // nos da la mayor cáscara cuya cota inferior esta por debajo del limite

  res.d = -1.0;
  res.cnt = -1;
  res.id = -1;

  if (indx < 0)
    return (res);

  // calculamos cuanto es la delta exacta al radio correspondiente a indx, si
  // es mayor diminuimos y volvemos a calcular, el bucle termina cuando tenemos
  // una cascara donde la delta es menor al limite (y por lo tanto cascara inmediata sup es mayor)
  // o cuando llegamos al ultimo cascaron y la delta sigue siendo superior
  do {
    npx = fast_counter(centre, indx, nstep, g, pnts);
    radx = (*g).maxscale * (float)(indx + 1) / (float)nstep;
    vol = (4.0 / 3.0 * M_PI * radx * radx * radx);
    deltax = ((float)npx / vol) / rhobar - 1.0;

    if (deltax < limit) {
      break;
    } else {
      indx--;
    }

  } while (indx >= 0);

  if (indx == (-1)) {
    // llegamos al ultimo cascaron sin caer debajo del limite, vamos a calcular
    // dentro de la bola y vemos si encontramos algo o lo descartamos
    radx = (*g).maxscale / (float)nstep;
    res = ball_finder(centre, limit, rhobar, 0, nstep, g, pnts);
  } else {
    // dado que estamos seguros que la cota minima de indx+1 es superior al limite (asi funciona raw_finder y el bucle
    // anterior) y sabemos que en este indx estamos debajo, el limite se cruza en un radio intermedio
    res = shell_finder(centre, limit, rhobar, indx, nstep, g, pnts);

    // busca entre los radios correspondientes a los cascarones de indx y indx+1, el indice
    // de la particula que da el cambio de densidad para eso utiliza solo los boxcell que debe abrir
    // segun la desigualdad del triangulo (como fast_counter hace)
  }
  return (res);
}

template <typename T>
float_in ball_finder(T centre, float limit, float rhobar, int irad, int nstep, Grid *g, std::vector<T> *pnts) {
  cell_state stat;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  float radx;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  // int ntotal = (*g).nsearch;
  float rad = (*g).maxscale * (float)(irad + 1) / (float)nstep;
  float radg = rad / (*g).gridsize;
  int ns = (int)ceil(rad / gridsize) + 1;

  long long int it, jt, kt, ind;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  long long int fst, nxt;
  unsigned long long int npx;
  std::vector<float_in> dist;
  float_in d, res;
  float vol, deltax;

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
              d.d = distance((*pnts)[nxt], centre, boxsize);
              d.id = nxt;
              if (d.d < rad) {
                dist.push_back(d);
              };
              nxt = (*g).next[nxt];
              if (nxt == fst)
                break;
            } while (1);
          }
        }
      }
    }
  }
  sort(dist.begin(), dist.end());

  res.d = -1.0;
  res.cnt = -1;
  res.id = -1;
  for (i = (dist.size() - 1); i >= 0; i--) {
    npx = (unsigned long long int)(i + 1);
    d = dist[i];
    radx = d.d;
    vol = (4.0 / 3.0 * M_PI * radx * radx * radx);
    deltax = ((float)npx / vol) / rhobar - 1.0;

    if (deltax < limit) {
      res.d = radx;
      res.cnt = npx;
      res.id = d.id;
      break;
    }
  }

  return (res);
}

template <typename T>
float_in shell_finder(T centre, float limit, float rhobar, int irad, int nstep, Grid *g, std::vector<T> *pnts) {
  cell_state stat_out, stat_in;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  // int ntotal = (*g).nsearch;
  float radsup = (*g).maxscale * (float)(irad + 2) / (float)nstep;
  float radinf = (*g).maxscale * (float)(irad + 1) / (float)nstep;
  float radsupg = radsup / (*g).gridsize;
  float radinfg = radinf / (*g).gridsize;
  float radx, vol, deltax;
  int ns = (int)ceil(radsup / gridsize) + 1;

  long long int it, jt, kt, ind;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  long long int fst, nxt;
  unsigned long long int npx;
  std::vector<float_in> dist;
  float_in d, res;

  unsigned long long int counter = 0;

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
          // stat_out=status_rad [nstep-(irad+1)-1] [indi]; //is_there(i,j,k,rad, g);
          stat_out = is_there(i, j, k, radsupg);
          // stat_in =status_rad [nstep-irad-1 ] [indi]; //is_there(i,j,k,rad, g);
          stat_in = is_there(i, j, k, radinfg);
          if (stat_in.is_in) {
            counter = counter + (*g).cnt[ind];
          } else if (!(stat_out.is_out)) {
            // open grid
            nxt = fst;
            do {
              d.d = distance((*pnts)[nxt], centre, boxsize);
              d.id = nxt;
              if ((d.d < radsup) && (d.d >= radinf)) {
                dist.push_back(d);
              } else if (d.d < radinf) {
                counter++;
              };
              nxt = (*g).next[nxt];
              if (nxt == fst)
                break;
            } while (1);
          }
        }
      }
    }
  }

  sort(dist.begin(), dist.end());

  res.d = -1.0;
  res.cnt = -1;
  res.id = -1;
  for (i = (dist.size() - 1); i >= 0; i--) {
    npx = (unsigned long long int)(i + 1 + counter);
    d = dist[i];
    radx = d.d;
    vol = (4.0 / 3.0 * M_PI * radx * radx * radx);
    deltax = ((float)npx / vol) / rhobar - 1.0;

    if (deltax < limit) {
      res.d = radx;
      res.cnt = npx;
      res.id = d.id;
      break;
    }
  }

  return (res);
}

template <typename T> bool is_touched(T centre, Grid *g, std::vector<T> *others, float maxscale) {
  cell_state stat;
  cell_state stat2;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  // int ntotal = (*g).nsearch;
  float radsearch = centre.radius + maxscale;
  int ns = (int)ceil(radsearch / gridsize) + 1;

  long long int it, jt, kt, ind;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  long long int fst, nxt;
  std::vector<float_in> dist;
  float d;
  bool res;

  // int jrad;
  // int irad = (int)(radsearch / ((*g).maxscale) * nstep) - 1;
  float radsearchg = radsearch / (*g).gridsize;

  res = false;

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
              if ((*others)[nxt].ToF && ((*others)[nxt].radius > centre.radius)) {
                float jradg = ((*others)[nxt].radius + centre.radius) / ((*g).gridsize);
                // stat2=status_rad[nstep-jrad-1][indi]; //stat2=is_there(i,j,k,(*others)[nxt].Rad+centre.Rad, g);
                stat2 = is_there(i, j, k, jradg);
                if (!stat2.is_out) {
                  if (stat2.is_in) {
                    res = true;
                    break;
                  } else {
                    d = distance((*others)[nxt], centre, boxsize);
                    if (d < ((*others)[nxt].radius + centre.radius)) {
                      res = true;
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
        if (res)
          break;
      }
      if (res)
        break;
    }
    if (res)
      break;
  }
  return (res);
}

template <typename T> std::vector<int> touch_list(T centre, Grid *g, std::vector<T> *others, float maxscale) {
  cell_state stat;
  cell_state stat2;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  // int ntotal = (*g).nsearch;
  float radsearch = centre.radius + maxscale;
  int ns = (int)ceil(radsearch / gridsize) + 1;

  long long int it, jt, kt, ind;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  long long int fst, nxt;
  std::vector<float_in> dist;
  float d;
  std::vector<int> res;

  for (k = (-ns); k < (ns + 1); k++) {
    for (j = (-ns); j < (ns + 1); j++) {
      for (i = (-ns); i < (ns + 1); i++) {
        it = indx_period(i + ic, Ngrid);
        jt = indx_period(j + jc, Ngrid);
        kt = indx_period(k + kc, Ngrid);

        ind = it + Ngrid * (jt + Ngrid * kt);
        fst = (*g).first[ind];
        if (fst >= 0) {
          stat = is_there(i, j, k, radsearch / (*g).gridsize);
          if (!(stat.is_out)) {
            // open grid
            nxt = fst;
            do {
              if ((*others)[nxt].radius > centre.radius) {
                stat2 = is_there(i, j, k, ((*others)[nxt].radius + centre.radius) / (*g).gridsize);
                if (!stat2.is_out) {
                  if (stat2.is_in) {
                    res.push_back(nxt);
                  } else {
                    d = distance((*others)[nxt], centre, boxsize);
                    if (d < ((*others)[nxt].radius + centre.radius)) {
                      res.push_back(nxt);
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

template <typename T>
std::vector<float_in> fast_ball(T centre, float rad, int nstep, Grid *g, std::vector<T> *pnts, std::vector<bool> *visit)
// std::vector <float_in> fast_ball(T centre, float rad,int nstep, Grid* g, std::vector<T> *pnts,cell_state
// **status_rad,std::vector<bool>* visit)
{
  cell_state stat;
  float3 boxsize = (*g).boxsize;
  int Ngrid = (*g).Ngrid;
  float gridsize = (*g).gridsize;
  int ic = int(centre.x / boxsize.x * Ngrid);
  int jc = int(centre.y / boxsize.y * Ngrid);
  int kc = int(centre.z / boxsize.z * Ngrid);

  // int ntotal = (*g).nsearch;

  int irad = (int)(rad / ((*g).maxscale) * nstep) - 1;
  int ns = (int)ceil(rad / gridsize) + 1;

  long long int it, jt, kt, ind;
  int i, j, k;
  // int dim = 2 * ntotal + 1;

  float radg = rad / (*g).gridsize;
  long long int fst, nxt;

  float_in d;
  std::vector<float_in> dist;

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
              if (!(*visit)[nxt]) {
                d.d = distance((*pnts)[nxt], centre, boxsize);
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

  sort(dist.begin(), dist.end());

  return (dist);
}

} // namespace grid
