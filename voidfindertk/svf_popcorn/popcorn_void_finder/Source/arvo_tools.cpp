#include "arvo_tools.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include "arvo_functions.h"
#include "constants.h"

#define NORTH_POLE_REDUCE 0.9999
#define EPS_NORTH_POLE 0.0001
#define EPS_TWO_PI 1e-12

using namespace std;
namespace arvo {
// fraction evaluation for integral
double Fract(double A, double B, double C, double sinphi, double cosphi, double k) {
  return (-B * sinphi + C * cosphi) / pow((A + B * cosphi + C * sinphi), k);
}

void PrintUsage() {
  printf("ARVO 2.0\n\nUsage: arvo protein=file1 [log=file2]\n\n");
  printf("Mandatory:\n  protein - name of input file as created by pdbread ");
  printf("(you can also use short form 'p=')\n");
  printf("Optional:\n");
  printf("      log - name of file for results and log messages\n");
}

/* unused function*/
int QuickSpheresEstimate(FILE *fpAts) {
  char buffer[85];
  int numSpheres = 0;
  rewind(fpAts);
  while (fgets(buffer, sizeof(buffer), fpAts)) {
    if (buffer[0] != '#')
      numSpheres++;
  }
  return numSpheres;
}

void NorthPoleFix(const int numAtoms, double *sphereLocal) {
  /*
   * Here we check local spheres and if two north poles are on same level,
   * we change radius of not 0 atom by factor of NORTH_POLE_REDUCE
   *
   * dmin - square of minimal distance of the North Pole to neighbor sphere
   */
  // here will come north pole test on local spheres!!!
  // all except atom 0
  // Test precision - MAY BE CHANGED
  for (int idx = 1; idx < numAtoms; idx++) {
    double d = fabs(sqrt((sphereLocal[0] - sphereLocal[idx * 4]) * (sphereLocal[0] - sphereLocal[idx * 4]) +
                         (sphereLocal[1] - sphereLocal[idx * 4 + 1]) * (sphereLocal[1] - sphereLocal[idx * 4 + 1]) +
                         (sphereLocal[2] + sphereLocal[3] - sphereLocal[idx * 4 + 2]) *
                             (sphereLocal[2] + sphereLocal[3] - sphereLocal[idx * 4 + 2])) -
                    sphereLocal[idx * 4 + 3]);
    if (d < EPS_NORTH_POLE) {
      sphereLocal[idx * 4 + 3] = sphereLocal[idx * 4 + 3] * NORTH_POLE_REDUCE;
      idx--;
      //  printf("Reduced sphere radius\n");
    }
  }
}

// computing integrals over arcs given in arc structure
// according to paper Hayrian, Dzurina, Plavka, Busa
void AvIntegral(const int nArcs, double *pVolume, double *pArea, const double r1, const double z1,
                const double *circles, const double *arcs) {
  int idx;
  double t, s, r, A, B, C, S, rr, ca, sa, cb, sb, be, al;
  double vIone, vItwo, vIthree, vJone, vJtwo, vJthree, delta_vint, delta_aint;

  *pVolume = 0.0;
  *pArea = 0.0;

  for (idx = 0; idx < nArcs; idx++) { // cycle over all arcs
    t = circles[((int)(arcs[idx * 3])) * 4];
    s = circles[((int)(arcs[idx * 3])) * 4 + 1];
    r = circles[((int)(arcs[idx * 3])) * 4 + 2];
    A = (4.0 * r1 * r1 + t * t + s * s + r * r) / 2.0;
    B = t * r;
    C = s * r;
    S = sqrt(A * A - B * B - C * C);
    rr = r * r - A;
    if (fabs(fabs(arcs[idx * 3 + 2]) - 2.0 * PI) < EPS_TWO_PI) { // full circle arc
      vIone = 2.0 * PI / S;
      vItwo = 2.0 * PI * A / (pow(S, 3.));
      vIthree = PI * (2.0 * A * A + B * B + C * C) / (pow(S, 5.));
      vJone = PI + rr / 2.0 * vIone;
      vJtwo = (vIone + rr * vItwo) / 4.0;
      vJthree = (vItwo + rr * vIthree) / 8.0;
      delta_vint = (128.0 * vJthree * pow(r1, 7.) + 8.0 * vJtwo * pow(r1, 5.) + 2.0 * vJone * pow(r1, 3.)) / 3.0 -
                   8.0 * pow(r1, 4.) * vJtwo * (z1 + r1);
      delta_aint = 2.0 * vJone * r1 * r1;
      if (arcs[idx * 3 + 2] < 0) {
        delta_vint = -delta_vint;
        delta_aint = -delta_aint;
      }
      *pVolume = *pVolume + delta_vint;
      *pArea = *pArea + delta_aint;
    } else { // integration over arcs
      if (arcs[idx * 3 + 2] < 0) {
        al = arcs[idx * 3 + 1] + arcs[idx * 3 + 2];
        be = arcs[idx * 3 + 1];
      } else {
        be = arcs[idx * 3 + 1] + arcs[idx * 3 + 2];
        al = arcs[idx * 3 + 1];
      }
      vIone = 2.0 *
              (PI / 2.0 - atan((A * cos((be - al) / 2.0) + B * cos((al + be) / 2.0) + C * sin((al + be) / 2.0)) /
                               (S * sin((be - al) / 2.0)))) /
              S;

      sb = sin(be);
      cb = cos(be);
      sa = sin(al);
      ca = cos(al);
      vItwo = (Fract(A, B, C, sb, cb, 1) - Fract(A, B, C, sa, ca, 1) + A * vIone) / (S * S);
      vIthree =
          (Fract(A, B, C, sb, cb, 2) - Fract(A, B, C, sa, ca, 2) +
           (Fract(A, B, C, sb, cb, 1) - Fract(A, B, C, sa, ca, 1)) / A + (2.0 * A * A + B * B + C * C) * vItwo / A) /
          (2.0 * S * S);
      vJone = ((be - al) + rr * vIone) / 2.0;

      vJtwo = (vIone + rr * vItwo) / 4.0;
      vJthree = (vItwo + rr * vIthree) / 8.0;
      //      delta_vint=(128.0*vJthree*pow(r1,7.)+8.0*vJtwo*pow(r1,5.)+ 2.0*vJone*pow(r1,3.))/3.0-8.0*pow(r1,4.)*vJtwo*(z1+r1);
      delta_vint = (64 * vJthree * pow(r1, 4) - 4 * vJtwo * r1 * (3 * z1 + 2 * r1) + vJone);
      delta_vint *= 2 * pow(r1, 3) / 3.0;
      delta_aint = 2.0 * vJone * r1 * r1;
      if (arcs[idx * 3 + 2] < 0) {
        delta_vint = -delta_vint;
        delta_aint = -delta_aint;
      }
      *pVolume = *pVolume + delta_vint;
      *pArea = *pArea + delta_aint;
    }
  }
}

} // namespace arvo
