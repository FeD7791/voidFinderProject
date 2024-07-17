#include "compute_intersecs.h"

#include "allvars.h"
#include "arvo_functions.h"
#include "colors.h"
#include "grid.h"
#include "objects.h"
#include "t.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

using namespace std;
using namespace arvo;
using namespace grid;

int compute_intersecs_start() {

  double *spheres;
  int i1, i2;
  int idmenor, idmayor;
  long unsigned int ip1, ip2;
  bool touch;
  float r1, r2, d;
  float joinvol, volint, minvol;
  float tdist;
  int isCave;
  double area1, vol1;
  int id1, id2, s1, s2, common;
  float pop_dens, pop_Delta;
  int sphsize;
  float density;
  float maxscale;
  int id_rad_menor, id_rad_mayor;
  par p;

  maxscale = -1.0;
  for (long unsigned int i = 0; i < pop_voids.size(); i++) {
    for (i1 = 0; i1 < pop_voids[i].nmem; i1++) {
      r1 = pop_voids[i].membs[i1].radius;
      if (maxscale < r1) {
        maxscale = r1;
      }
    }
  }

  density = NTRAC / boxsize.x / boxsize.y / boxsize.z;

#pragma omp declare reduction(merge : std::vector<par> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  cout << CYAN << "procesando.... " << DEFA << endl;

#pragma omp parallel for default(none) reduction(merge : pairs)                                                        \
    shared(pop_voids, boxsize, density, pop_Delta, cout) private(                                                      \
            ip1, ip2, i1, i2, r1, r2, d, touch, tdist, sphsize, spheres, area1, vol1, joinvol, volint, minvol,         \
                idmenor, idmayor, id_rad_menor, id_rad_mayor, common, s1, s2, id1, id2, pop_dens, p, isCave)
  for (ip1 = 0; ip1 < pop_voids.size(); ip1++) {
    if (pop_voids.size() > 1000) {
      Progress(ip1, pop_voids.size());
    } else {
      cout << ip1 << " / " << pop_voids.size() << endl;
      cout.flush();
    }
    for (ip2 = 0; ip2 < pop_voids.size(); ip2++) {
      if (ip2 == ip1)
        continue;
      touch = false;

      for (i1 = 0; i1 < pop_voids[ip1].nmem; i1++) {
        r1 = pop_voids[ip1].membs[i1].radius;
        for (i2 = 0; i2 < (pop_voids[ip2]).nmem; i2++) {
          d = distance(pop_voids[ip1].membs[i1].x, pop_voids[ip1].membs[i1].y, pop_voids[ip1].membs[i1].z,
                       (pop_voids[ip2]).membs[i2].x, (pop_voids[ip2]).membs[i2].y, (pop_voids[ip2]).membs[i2].z,
                       boxsize);
          r2 = pop_voids[ip2].membs[i2].radius;
          if (d < (r1 + r2) * 0.99) {
            touch = true;
            tdist = d / (r1 + r2);
            break;
          }
        }
        if (touch)
          break;
      }

      if (touch) {
        sphsize = pop_voids[ip1].membs.size() + pop_voids[ip2].membs.size();
        spheres = new double[4 * sphsize];

        i2 = 0;
        for (i1 = 0; i1 < (int)pop_voids[ip1].membs.size(); i1++) {
          spheres[4 * i2] = (double)(pop_voids[ip1].membs[i1]).x;
          spheres[4 * i2 + 1] = (double)(pop_voids[ip1].membs[i1]).y;
          spheres[4 * i2 + 2] = (double)(pop_voids[ip1].membs[i1]).z;
          spheres[4 * i2 + 3] = (double)(pop_voids[ip1].membs[i1]).radius;
          i2++;
        }

        for (i1 = 0; i1 < (int)pop_voids[ip2].membs.size(); i1++) {
          spheres[4 * i2] = (double)(pop_voids[ip2].membs[i1]).x;
          spheres[4 * i2 + 1] = (double)(pop_voids[ip2].membs[i1]).y;
          spheres[4 * i2 + 2] = (double)(pop_voids[ip2].membs[i1]).z;
          spheres[4 * i2 + 3] = (double)(pop_voids[ip2].membs[i1]).radius;
          i2++;
        }

        vol1 = 0;
        check_boundary_splits(spheres, sphsize, boxsize);
        check_inner_spheres_v2(&spheres, &sphsize);

        Arvo_class *arvo_obj0 = new Arvo_class();
        try {
          isCave = (*arvo_obj0).compute_volume_arvo(&area1, &vol1, spheres, sphsize);
        } catch (const std::bad_alloc &e) {
          arvo_obj0->Cleanup();
          exit(-1);
        }
        arvo_obj0->Cleanup();
        delete (arvo_obj0);

        joinvol = vol1;
        // double area = area1;

        volint = pop_voids[ip1].vol + pop_voids[ip2].vol - joinvol;

        if (volint < 0.0) {
          if (abs(volint) / joinvol < 1e-3) {
            volint = -volint;
          } else {
            cout << tdist << endl;
            cout << pop_voids[ip1].membs.size() << " + " << pop_voids[ip2].membs.size() << " = " << sphsize << endl;
            cout << pop_voids[ip1].nmem << " + " << pop_voids[ip2].nmem << " = " << sphsize << endl;
            cout << pop_voids[ip1].vol << "  " << pop_voids[ip2].vol << " | " << joinvol << " " << volint << endl;
            cout << ip1 << " " << ip2 << endl;
            cout << "pasa algo feo pepe" << endl;
            cout.flush();
            volint = 1.0e15;
            // exit(-1);
          }
        }
        delete[] spheres;

        if (pop_voids[ip1].vol < pop_voids[ip2].vol) {
          minvol = pop_voids[ip1].vol;
          idmenor = ip1;
          idmayor = ip2;
        } else {
          minvol = pop_voids[ip2].vol;
          idmenor = ip2;
          idmayor = ip1;
        }

        /// cambio de estretegia... ordenamos no por volumen total, si no por el
        /// tamaño del void esférico de semilla (comenzar con un void chico
        /// puede ser muy de puente)
        if (pop_voids[ip1].membs[0].radius < pop_voids[ip2].membs[0].radius) {
          id_rad_menor = ip1;
          id_rad_mayor = ip2;
        } else {
          id_rad_menor = ip2;
          id_rad_mayor = ip1;
        }

        i1 = 0;
        i2 = 0;
        common = 0;
        s1 = pop_voids[ip1].in_halos.size();
        s2 = pop_voids[ip2].in_halos.size();
        do {
          id1 = pop_voids[ip1].in_halos[i1];
          id2 = pop_voids[ip2].in_halos[i2];

          if (id1 == id2) {
            common++;
            i1++;
            i2++;
          } else {
            if (id1 < id2) {
              i1++;
            } else {
              i2++;
            }
          }
        } while (i1 < s1 && i2 < s2);

        pop_dens = (float)(pop_voids[ip1].npart + pop_voids[ip2].npart - common) / joinvol;
        pop_Delta = pop_dens / density - 1.0;

        p.idmin = pop_voids[idmenor].id;
        p.idmax = pop_voids[idmayor].id;
        p.root_rad_min = pop_voids[id_rad_menor].membs[0].radius;
        p.root_rad_max = pop_voids[id_rad_mayor].membs[0].radius;
        p.id_rad_min = pop_voids[id_rad_menor].id;
        p.id_rad_max = pop_voids[id_rad_mayor].id;
        p.volint = volint;
        p.minvol = minvol;
        p.pop_Delta = pop_Delta;
        p.joinvol = joinvol;
        p.npart1 = pop_voids[idmenor].npart;
        p.npart2 = pop_voids[idmayor].npart;
        p.common = common;
        p.vol1 = pop_voids[idmenor].vol;
        p.vol2 = pop_voids[idmayor].vol;
        p.isCave = isCave;
        pairs.push_back(p);
      }
    }
  }
  return 0;
}
