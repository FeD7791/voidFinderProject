#include "arvo_functions.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>

#include "arvo_tools.h"
#include "arvo_variables_tools.h"
#include "constants.h"

// #define ARRAY_INCREMENT 10 // by this size are arrays resized when necessary

#define EPS_DELTAT 1e-12
#define EPS_ANGLE 1e-12
// #define EPS_TWO_PI 1e-12

// #define RANDOM_SA 0.324 // "random" sin value for spheres rotation
// #define RANDOM_CA 0.946057081 // \sqrt{1-RANDOM_SA^2}
using namespace arvo;
using namespace std;

void Arvo_class::Initialization() {
  //  maxNeighbors = 0;
  ind = NULL;
  spheres = NULL;
  neighborsNumber = NULL;
  neighborsIndices = NULL;
  indexStart = NULL;
  sizeNeighborsIndices = 0;
  sizeInd = 0;
}

static inline void CirclesIntersection(const double *Circ1, const double *Circ2, double *a1, double *a2) {
  /*
    Function returns angles of two intersection points of circles with
    indices ic1 and ic2 in circles structure circles (we will use it ONLY
    IN CASE, WHEN 2 INTERSECTION POINTS EXIST)
    a1 and a2 are corresponding angles with respect to the center of 1st circ
    b1 and b2 are corresponding angles with respect to the center of 2nd circ
  */

  //     (t,s) - circle center, r - circle radius
  double A, B, C, D;
  *a1 = 0;
  *a2 = 0;
  double b1 = 0;
  double b2 = 0;

  double t1 = Circ1[0];
  double s1 = Circ1[1];
  double r1 = Circ1[2];
  double t2 = Circ2[0];
  double s2 = Circ2[1];
  double r2 = Circ2[2];

  if (fabs(t2 - t1) < EPS_DELTAT) // t2 == t1
  {
    B = ((r1 * r1 - r2 * r2) / (s2 - s1) - (s2 - s1)) / 2.0;
    A = sqrt(r2 * r2 - B * B);
    if (B == 0) {
      b1 = 0.0;
      b2 = PI;
    } else if (B > 0) {
      b1 = atan(fabs(B / A));
      b2 = PI - b1;
    } else {
      b1 = PI + atan(fabs(B / A));
      b2 = 3.0 * PI - b1;
    }
    B = B + s2 - s1;
    if (B == 0) {
      *a1 = 0.0;
      *a2 = PI;
    } else if (B > 0) {
      *a1 = atan(fabs(B / A));
      *a2 = PI - *a1;
    } else {
      *a1 = PI + atan(fabs(B / A));
      *a2 = 3.0 * PI - *a1;
    }
  } else // t2 != t1
  {
    C = ((r1 * r1 - r2 * r2 - (s2 - s1) * (s2 - s1)) / (t2 - t1) - (t2 - t1)) / 2.0;
    D = (s1 - s2) / (t2 - t1);
    B = (-C * D + sqrt((D * D + 1.0) * r2 * r2 - C * C)) / (D * D + 1.0);
    A = C + D * B;
    if (A == 0)
      b1 = (B > 0) ? PI / 2.0 : -PI / 2.0;
    else if (A > 0)
      b1 = atan(B / A);
    else
      b1 = PI + atan(B / A);
    B = B + s2 - s1;
    A = A + t2 - t1;
    if (A == 0)
      *a1 = (B > 0) ? PI / 2.0 : -PI / 2.0;
    else if (A > 0)
      *a1 = atan(B / A);
    else
      *a1 = PI + atan(B / A);
    B = (-C * D - sqrt((D * D + 1.0) * r2 * r2 - C * C)) / (D * D + 1.0);
    A = C + D * B;
    if (A == 0)
      b2 = (B > 0) ? PI / 2.0 : -PI / 2.0;
    else if (A > 0)
      b2 = atan(B / A);
    else
      b2 = PI + atan(B / A);
    B = B + s2 - s1;
    A = A + t2 - t1;
    if (A == 0)
      *a2 = (B > 0) ? PI / 2.0 : -PI / 2.0;
    else if (A > 0)
      *a2 = atan(B / A);
    else
      *a2 = PI + atan(B / A);
  }

  if (*a1 < 0)
    *a1 = *a1 + 2.0 * PI;
  if (*a2 < 0)
    *a2 = *a2 + 2.0 * PI;
  if (b1 < 0)
    b1 = b1 + 2.0 * PI;
  if (b2 < 0)
    b2 = b2 + 2.0 * PI;
}

// 1  if Circ1-th circle is inside Circ2-th positive circle or outside Circ2-th
// negative circle 0  otherwise WE KNOW, THAT CIRCLES HAVE LESS THAN 2
// INTERSECTION POINTS! i -> Circ1, k->Circ2
static inline int CircleInCircle(const double *Circ1, const double *Circ2) {
  int CirclesInCircle = 0;

  double d2 = (Circ1[0] + Circ1[2] - Circ2[0]) * (Circ1[0] + Circ1[2] - Circ2[0]) +
              (Circ1[1] - Circ2[1]) * (Circ1[1] - Circ2[1]);

  if (d2 < Circ2[2] * Circ2[2])
    CirclesInCircle = (Circ2[3] > 0) ? 1 : 0;
  else if (d2 > Circ2[2] * Circ2[2])
    CirclesInCircle = (Circ2[3] > 0) ? 0 : 1;
  else {
    d2 = (Circ1[0] - Circ2[0]) * (Circ1[0] - Circ2[0]) + (Circ1[1] - Circ2[1]) * (Circ1[1] - Circ2[1]);
    if (d2 < Circ2[2] * Circ2[2])
      CirclesInCircle = (Circ2[3] > 0) ? 1 : 0;
    else
      CirclesInCircle = (Circ2[3] > 0) ? 0 : 1;
  }
  return CirclesInCircle;
}

// 1  if point (t,s) is inside k-th positive circle or outside k-th negative
// circle 0  otherwise WE KNOW, THAT POINT IS NOT ON THE CIRCLE!
static inline int PointInCircle(const double t, const double s, const double *circle) {
  int PointInCircle = 0;

  double d2 = (t - circle[0]) * (t - circle[0]) + (s - circle[1]) * (s - circle[1]);
  if (d2 < circle[2] * circle[2]) {
    PointInCircle = (circle[3] > 0) ? 1 : 0;
  } else {
    PointInCircle = (circle[3] > 0) ? 0 : 1;
  }
  return PointInCircle;
}

int Arvo_class::NewArcs(const int Circle,
                        const int NumAtoms) //, double **circlesTemp, double **angles_temp, int
                                            // sizeAngles_temp, double **newArcs_temp,  int
                                            // sizeNewArcs_temp)
{

  // Function prepares arcs, which are part of i-th circle in circle structure
  // circles. Interesting are these arcs, which are inside other positive
  // circles or outside other negative circles
  //
  // Matrix arcsnew in each row has elements
  //
  // arcsnew(i,1)=ic - ic is the index of arc-circle in circle
  // arcsnew(i,2)=sigma - sigma is the starting angle of arc
  // arcsnew(i,3)=delta - delta is oriented arc angle

  //  double *circlesAux= *circlesTemp;
  //  double *anglesAux= *angles;
  // double *newArcsAux = *newArcs_temp;

  //  if(anglesAux != NULL){
  //    cout<<"angles_   \n";
  //   cout<<"angles+   "<< angles;
  //   cout<<"angles*   \n"<<* angles;
  //    cout<<"angles&   \n"<<&angles;
  // }

  // if(anglesAux != NULL){
  //   cout << ": angles :"<<angles<<": * :"<<*angles<<": & :"<<&angles<<"
  //   anglesTemp :"<< angles_temp <<": * :"<< *angles_temp<<": & :"<<
  //   &angles_temp << ": anglesAux  :"<<anglesAux<<": * : "<<*anglesAux<<": &
  //   :"<<&anglesAux<<"\n" ;
  //}

  //      integer circle_in_circle,point_in_circle,delete_equal
  double a1, a2;
  int numArc = 0;
  int numCond = 0;
  int numAngle = 0;
  angles.clear();

  double ti = circles[Circle * 4 + 0];
  double si = circles[Circle * 4 + 1];
  double ri = circles[Circle * 4 + 2];

  for (int idx = 0; idx < NumAtoms; idx++) // composition of angles vector, consisting of intersection points
  {
    //  for (idx = 28; idx<29; idx++) { // composition of angles vector,
    //  consisting of intersection points
    if (idx != Circle) {
      double t = circles[idx * 4];
      double s = circles[idx * 4 + 1];
      double r = circles[idx * 4 + 2];
      double d2 = (ti - t) * (ti - t) + (si - s) * (si - s);
      if ((d2 < (r + ri) * (r + ri)) && (fabs(r - ri) * fabs(r - ri) < d2)) // 2 intersection points exist
      {
        CirclesIntersection(&circles[Circle * 4], &circles[idx * 4], &a1, &a2);

        angles.emplace_back(a1);
        angles.emplace_back(a2);
        numAngle += 2;
      }
    }
  }

  if (!numAngle) // there are no double intersections of idx-th circles with
                 // others
  {
    numCond = 0;
    for (int idx = 0; idx < NumAtoms; idx++) // if i-th circle is inside of all other positive and outside of
                                             // all other negative circles, it will be new arc
    {
      if (idx != Circle) {
        int cic = CircleInCircle(&circles[Circle * 4], &circles[idx * 4]);
        if (cic == 0) {
          break;
        }
        numCond = numCond + cic;
      }
    }
    if (numCond == (NumAtoms - 1)) // all conditions hold
    {
      newArcs.emplace_back(Circle);
      newArcs.emplace_back(0.0);
      newArcs.emplace_back(2.0 * PI * circles[Circle * 4 + 3]);
      numArc++;
    }
  } else // there are double intersection points
  {
    if (circles[Circle * 4 + 3] > 0) {
      std::sort(angles.begin(), angles.end());
    } else {
      std::sort(angles.begin(), angles.end(), std::greater<double>());
    }
    auto last = std::unique(angles.begin(), angles.end(), [](double a, double b) { return fabs(a - b) <= EPS_ANGLE; });
    numAngle = std::distance(angles.begin(), last);
    for (int idx = 0; idx < numAngle - 1; idx++) {
      numCond = 0;
      double t = ti + ri * cos((angles[idx] + angles[idx + 1]) / 2.0);
      double s = si + ri * sin((angles[idx] + angles[idx + 1]) / 2.0);
      for (int jdx = 0; jdx < NumAtoms; jdx++) {
        if (jdx != Circle) {
          int pic = PointInCircle(t, s, &circles[jdx * 4]);
          if (pic == 0) {
            break;
          }
          numCond = numCond + pic;
        }
      }
      if (numCond == (NumAtoms - 1)) // all conditions hold
      {
        newArcs.emplace_back(Circle);
        newArcs.emplace_back(angles[idx]);
        newArcs.emplace_back(angles[idx + 1] - angles[idx]);
        numArc++; //  zero based indices!
      }
    }

    numCond = 0;
    double t = ti + ri * cos((angles[0] + 2.0 * PI + angles[numAngle - 1]) / 2.0);
    double s = si + ri * sin((angles[0] + 2.0 * PI + angles[numAngle - 1]) / 2.0);
    for (int idx = 0; idx < NumAtoms; idx++) {
      if (idx != Circle) {
        int pic = PointInCircle(t, s, &circles[idx * 4]);
        if (pic == 0) {
          break;
        }
        numCond = numCond + pic;
      }
    }
    if (numCond == (NumAtoms - 1)) // all conditions hold
    {
      newArcs.emplace_back(Circle);
      newArcs.emplace_back(angles[numAngle - 1]);
      newArcs.emplace_back(angles[0] + circles[Circle * 4 + 3] * 2.0 * PI - angles[numAngle - 1]);
      numArc++;
    }
  }
  return numArc;
}

int Arvo_class::CirclesToArcs(const int NumAtoms) { // circles,arcs,kl,nls,ka){

  // Computing integration arcs

  // arcs(i,1)=ci        - corresponding circle index
  // arcs(i,2)=sigma     - starting arc angle
  // arcs(i,3)=delta     - oriented arc angle

  // Arcs (with their orientation) are parts of circles, which
  // bounds are circles intersection points. If the center of
  // arc lies inside all other positive and outside all other
  // negative circles, then we will put it inside arcs structure

  int idx, nna;
  int nArcs = 0;
  arcs.clear();
  if (NumAtoms == 1) // we have only 1 circle
  {
    nArcs = 1;
    arcs.emplace_back(0.0);
    arcs.emplace_back(0.0);
    arcs.emplace_back(2.0 * PI * circles[3]);
  } else // more than 1 circle
  {
    for (idx = 0; idx < NumAtoms; idx++) {
      newArcs.clear();
      nna = NewArcs(idx, NumAtoms); //,&circles,&angles,sizeAngles,&newArcs,sizeNewArcs );
      if (nna) {
        std::copy(newArcs.begin(), newArcs.end(), std::back_inserter(arcs));
        nArcs += nna;
      }
    }
  }
  return nArcs;
}

void Arvo_class::MakeTsCircles(const int NumAtoms) // sphere_local,circles,kl,nls)
{

  // Preparing circles structure for 1st sphere in array      circles
  // according to the paper Hayrjan, Dzurina, Plavka, Busa
  //
  // circles(i,1)=ti
  // circles(i,2)=si     - ith circle's center coordinates
  // circles(i,3)=ri    - ith circle's radius
  // circles(i,4)=+1/-1 - circle orientation

  //  dimension circles(kl,4),sphere_local(kl,4)
  double r1 = sphereLocal[3];
  circles.clear();
  for (int idx = 0; idx < NumAtoms; idx++) {
    double dx = sphereLocal[0] - sphereLocal[(idx + 1) * 4];
    double dy = sphereLocal[1] - sphereLocal[(idx + 1) * 4 + 1];
    double a = dx * dx + dy * dy +
               (sphereLocal[2] + r1 - sphereLocal[(idx + 1) * 4 + 2]) *
                   (sphereLocal[2] + r1 - sphereLocal[(idx + 1) * 4 + 2]) -
               sphereLocal[(idx + 1) * 4 + 3] * sphereLocal[(idx + 1) * 4 + 3];
    double b = 8.0 * r1 * r1 * dx;
    double c = 8.0 * r1 * r1 * dy;
    double d = 4.0 * r1 * r1 *
               (dx * dx + dy * dy +
                (sphereLocal[2] - r1 - sphereLocal[(idx + 1) * 4 + 2]) *
                    (sphereLocal[2] - r1 - sphereLocal[(idx + 1) * 4 + 2]) -
                sphereLocal[(idx + 1) * 4 + 3] * sphereLocal[(idx + 1) * 4 + 3]);
    circles.emplace_back(-b / (2.0 * a));
    circles.emplace_back(-c / (2.0 * a));
    circles.emplace_back(sqrt((b * b + c * c - 4.0 * a * d) / (4.0 * a * a)));
    if (a > 0) {
      circles.emplace_back(-1);
    } else {
      circles.emplace_back(1);
    }
  }

  // for (int idx = 0; idx < NumAtoms; idx++)
  //   printf("%d %f %f %f
  //   %f\n",idx,circles[idx*4],circles[idx*4+1],circles[idx*4+2],circles[idx*4+3]);
}

// take sphere_local out of the main array spheres
void Arvo_class::LocalSpheres(const int numAtoms) {
  sphereLocal.clear();
  for (int idx = 0; idx < numAtoms; idx++) {
    for (int jdx = 0; jdx < 4; jdx++) {
      sphereLocal.emplace_back(spheres[ind[idx] * 4 + jdx]);
    }
  }
  arvo::NorthPoleFix(numAtoms, sphereLocal.data());
}

// Function computes Atom-th part of the whole volume - the volume
// of domain inside Atom-th and outside of all other spheres
void Arvo_class::AreaVolume(const int Atom, double *pVolume, double *pArea) {
  // nls was originaly neighborsNumber[Atom]+1, shifted to
  // neighborsNumber[Atom]!!!!
  *pVolume = 0;
  *pArea = 0;
  // dimension spheres(ks,4),circles(kl,4),arcs(ka,3),
  // sphere_local(kl,4),ind(kl),neighbors_number(ks),
  // index_start(ks),neighbors_indices(ki),av(2),avi(2)

  // circles, arcs, sphere_local are described below
  // integer circles_to_arcs

  // Determination of i-th sphere's neighbors (row indices in matrix spheres)
  if (neighborsNumber[Atom] < 0) {
    // ith sphere is subset of other sphere, sphere(i,4) will be done negative
    *pVolume = 0.0;
    *pArea = 0.0;
  } else if (neighborsNumber[Atom] == 0) {
    // there are no neighbors (nls - number of local spheres = ith sphere +
    // neigh
    *pVolume = 4.0 * PI * spheres[Atom * 4 + 3] * spheres[Atom * 4 + 3] * spheres[Atom * 4 + 3] / 3.0;
    *pArea = 4.0 * PI * spheres[Atom * 4 + 3] * spheres[Atom * 4 + 3];
  } else {
    int idx, nArcs, nPos;
    double z1, r1, partVol, partArea;
    // there are neighbors

    arvo::SetInd(0, Atom, &ind, &sizeInd);
    // SetInd(0,Atom);

    for (idx = 0; idx < neighborsNumber[Atom]; idx++) {
      // SetInd(idx+1,neighborsIndices[indexStart[Atom]+idx]);
      arvo::SetInd(idx + 1, neighborsIndices[indexStart[Atom] + idx], &ind, &sizeInd);
    }

    // SetInd(idx+1,neighborsIndices[indexStart[Atom]+idx]);
    // we will work only with ith and neighbors spheres
    LocalSpheres(neighborsNumber[Atom] + 1);

    *pVolume = 0.0;
    *pArea = 0.0;

    MakeTsCircles(neighborsNumber[Atom]);         // sphere_local,circles,kl,neighborsNumber[Atom])
    nArcs = CirclesToArcs(neighborsNumber[Atom]); // circles,arcs,kl,neighborsNumber[Atom],ka)

    nPos = 0;
    for (idx = 0; idx < neighborsNumber[Atom]; idx++) {
      if (circles[idx * 4 + 3] > 0)
        nPos++;
    }

    z1 = sphereLocal[2];
    r1 = sphereLocal[3];

    AvIntegral(nArcs, &partVol, &partArea, r1, z1, circles.data(), arcs.data());
    if (nPos) // there exists positive oriented circle
    {
      *pVolume = *pVolume + partVol;
      *pArea = *pArea + partArea;
    } else // all circles are negative oriented - we have computed complement
    {
      *pVolume = *pVolume + partVol + 4.0 * PI * sphereLocal[3] * sphereLocal[3] * sphereLocal[3] / 3.0;
      *pArea = *pArea + partArea + 4.0 * PI * sphereLocal[3] * sphereLocal[3];
    }
  }
}

int Arvo_class::Neighbors(int Atom) {

  // If ith sphere is a subset of other sphere, index_number(i)=-1 and we change
  // radius in matrix spheres to -radius!
  // If some other sphere is subset of ith sphere, than we change its radius to

  int NeighborsNum = 0, idx; // number of neighbors
  double xi, yi, zi, ri, dd, rk;

  xi = spheres[Atom * 4];
  yi = spheres[Atom * 4 + 1];
  zi = spheres[Atom * 4 + 2];
  ri = spheres[Atom * 4 + 3];
  for (idx = 0; idx < spheresNumber; idx++) {
    if (idx == Atom)
      continue;
    if (fabs(xi - spheres[idx * 4]) < ri + spheres[idx * 4 + 3]) {
      dd = sqrt((xi - spheres[idx * 4]) * (xi - spheres[idx * 4]) +
                (yi - spheres[idx * 4 + 1]) * (yi - spheres[idx * 4 + 1]) +
                (zi - spheres[idx * 4 + 2]) * (zi - spheres[idx * 4 + 2]));
      rk = spheres[idx * 4 + 3];
      if (dd < ri + rk) {
        if (dd + ri <= rk) // ith sphere is inside of other sphere
        {
          NeighborsNum = -1;
          return NeighborsNum;
        } else if (dd + rk > ri) // kth sphere is neighbor
        {
          // SetInd(NeighborsNum++,idx);
          arvo::SetInd(NeighborsNum++, idx, &ind, &sizeInd);
        }
      }
    }
  }
  return NeighborsNum;
}

void Arvo_class::MakeNeighbors(void) {

  // Determination of neighbors for all atoms. We construct next structure:
  // neighbors_number(i)=neighbors number for ith atom
  // index_start(i)=start of neighbors indices for ith atom in array
  // neighbors_indices neighbors_indices = array of neighbors indices for each
  // atom
  // neighbors_indices(index_start(i)):neighbors_ind(index_start(i)+neighbors_number
  //
  // For example: 1. atom has neighbors with indices 2, 4, 7
  // 2. atom has neighbors with indices 1, 3
  // 3. atom has neighbors with indices 2, 4
  // 4. atom has neighbors with indices 1, 3
  // 5. atom is subset of some atom
  // 6. atom has no neighbors
  // 7. atom has neighbors with index 1
  // then we have:
  // neighbors_number=(3,2,2,2,-1,0,1)
  // index_start=(1,4,6,8,10,10,10,11)
  // neighbors_indices(2,4,7,1,3,2,4,1,3,1)

  int idx, jdx;
  indexStart[0] = 0;
  for (idx = 0; idx < spheresNumber; idx++) {
    neighborsNumber[idx] = Neighbors(idx);
    // it is not used anywhere   if (neighborsNumber[idx] > maxNeighbors)
    //      maxNeighbors = neighborsNumber[idx];

    if (neighborsNumber[idx] <= 0) // sphere is subset or there are no neighbors
      indexStart[idx + 1] = indexStart[idx];
    else { //  there are neighbors
      indexStart[idx + 1] = indexStart[idx] + neighborsNumber[idx];
      for (jdx = 0; jdx < neighborsNumber[idx]; jdx++)
        arvo::SetNI(indexStart[idx] + jdx, ind[jdx], &neighborsIndices, &sizeNeighborsIndices);

      // SetNI(indexStart[idx]+jdx,ind[jdx]);
    }
  }
}

void Arvo_class::Cleanup() {
  if (ind) {
    free(ind);
    ind = NULL;
  }
  // if (newArcs) {free(newArcs); newArcs=NULL;}
  // if (spheres) {free(spheres); spheres=NULL;}
  if (neighborsNumber) {
    free(neighborsNumber);
    neighborsNumber = NULL;
  }
  if (neighborsIndices) {
    free(neighborsIndices);
    neighborsIndices = NULL;
  }
  if (indexStart) {
    free(indexStart);
    indexStart = NULL;
  }
}

// Computing surface area and volume of the overlapping spheres

int Arvo_class::compute_volume_arvo(double *area1, double *vol1, double *test_arr, int sphsize) {
  // FILE *fpLog, *fpAts;
  // char *atsFile, *logFile; // Name of the *.pdb and *.log file
  int idx; // parameters set
  double Volume, Area, retVolume, retArea;

  // ks - maximal spheres' number
  // kl - maximal neighbors' number of one sphere (local spheres' nu
  // ka - maximal angles' or arcs' number
  // ki - maximal neighbors' relations' number = cca. spheres' number * maximal
  // neighbors' number eps_north_pole - accuracy level in the function
  // North_Pole_test eps_deltat - accuracy level in the subroutine
  // circles_intersection eps_angle - accuracy level in the subroutine
  // delete_equal (angles)

  // probe radius -> rWater
  // radii set
  // protein name
  // pdbread(inputFile, ivdwrSet);

  Initialization();
  spheresNumber = sphsize;
  spheres = test_arr;
  indexStart = (int *)malloc((spheresNumber + 1) * sizeof(int)); // bug por eso el +1? no se entiende
                                                                 // la explicacion de makeneighbors
  neighborsNumber = (int *)malloc(spheresNumber * sizeof(int));

  if (spheresNumber <= 0) {
    printf("No spheres to calculate volume. Exiting.\n");
    Cleanup();
    return -1;
  }

  // shift the molecule so, that it is centered around its center of gravity
  double avgX, avgY, avgZ;
  avgX = avgY = avgZ = 0.0;
  for (int idx = 0; idx < spheresNumber; ++idx) {
    avgX += spheres[idx * 4];
    avgY += spheres[idx * 4 + 1];
    avgZ += spheres[idx * 4 + 2];
  }
  avgX /= (spheresNumber * 1.0);
  avgY /= (spheresNumber * 1.0);
  avgZ /= (spheresNumber * 1.0);
  for (int idx = 0; idx < spheresNumber; ++idx) {
    spheres[idx * 4] -= avgX;
    spheres[idx * 4 + 1] -= avgY;
    spheres[idx * 4 + 2] -= avgZ;
  }

  // Study the neighborhood relations
  MakeNeighbors();

  // Computation of area and volume as a sum of surface integrals
  // clock_t start, end;
  // start = clock();
  Volume = 0;
  Area = 0;
  for (idx = 0; idx < spheresNumber; idx++) {
    AreaVolume(idx, &retVolume, &retArea);
    //    printf("Element %d\tvolume %f\tarea %f\n",idx,retVolume,retArea);
    Volume += retVolume;
    Area += retArea;
  }

  *area1 = Area;
  *vol1 = Volume;
  // printf("Volume: %f\tArea: %f\tSpheres num: %d\n", Volume, Area,
  // spheresNumber);

  // end = clock();
  // double cas = (end-start)/(1.0*CLOCKS_PER_SEC);
  // printf ("It took %f seconds\n", cas);

  Cleanup();
  return (0);
}

arvo::Arvo_class::Arvo_class() { Initialization(); }
arvo::Arvo_class::~Arvo_class() { Cleanup(); }
