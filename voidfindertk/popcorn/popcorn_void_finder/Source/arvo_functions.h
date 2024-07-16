#ifndef ARVO_FUNCTIONS_H_
#define ARVO_FUNCTIONS_H_

#include <vector>

namespace arvo {

class Arvo_class {
  int sizeInd;
  int *ind; // there shouldn't be more neighbors

  int spheresNumber;
  double *spheres;      // 4*count
  int *neighborsNumber; // same size as spheres
  int *indexStart;      // same size as spheres
  // rest will have dynamic size
  int sizeNeighborsIndices;
  int *neighborsIndices;
  std::vector<double> sphereLocal;
  std::vector<double> circles;
  std::vector<double> arcs;
  std::vector<double> newArcs;
  std::vector<double> angles;

public:
  void Initialization();
  /*
  Function prepares arcs, which are part of i-th circle in circle structure
  circles. Interesting are these arcs, which are inside other positive circles
  or outside other negative circles
  /Matrix arcsnew in each row has elements

  arcsnew(i,1)=ic - ic is the index of arc-circle in circle
  arcsnew(i,2)=sigma - sigma is the starting angle of arc
  arcsnew(i,3)=delta - delta is oriented arc angle
  */
  int NewArcs(const int Circle, const int NumAtoms);

  /*Computing integration arcs
  arcs(i,1)=ci        - corresponding circle index
  arcs(i,2)=sigma     - starting arc angle
  arcs(i,3)=delta     - oriented arc angle
  Arcs (with their orientation) are parts of circles, which
  bounds are circles intersection points. If the center of
  arc lies inside all other positive and outside all other
  negative circles, then we will put it inside arcs structure
  */
  int CirclesToArcs(const int NumAtoms);

  /*
  Preparing circles structure for 1st sphere in array      circles
  according to the paper Hayrjan, Dzurina, Plavka, Busa

  circles(i,1)=ti
  circles(i,2)=si     - ith circle's center coordinates
  circles(i,3)=ri    - ith circle's radius
  circles(i,4)=+1/-1 - circle orientation
  */
  void MakeTsCircles(const int NumAtoms); // sphere_local,circles,kl,nls)

  /* take sphere_local out of the main array spheres*/
  void LocalSpheres(const int numAtoms);

  /*
  Function computes Atom-th part of the whole volume - the volume
  of domain inside Atom-th and outside of all other spheres
  */
  void AreaVolume(const int Atom, double *pVolume, double *pArea);

  /*
  If ith sphere is a subset of other sphere, index_number(i)=-1 and we change
  radius in matrix spheres to -radius!
  If some other sphere is subset of ith sphere, than we change its radius to
  */
  int Neighbors(int Atom);

  /*
  Determination of neighbors for all atoms. We construct next structure:
  neighbors_number(i)=neighbors number for ith atom
  index_start(i)=start of neighbors indices for ith atom in array
  neighbors_indices neighbors_indices = array of neighbors indices for each atom
  neighbors_indices(index_start(i)):neighbors_ind(index_start(i)+neighbors_number

  For example: 1. atom has neighbors with indices 2, 4, 7
  2. atom has neighbors with indices 1, 3
  3. atom has neighbors with indices 2, 4
  4. atom has neighbors with indices 1, 3
  5. atom is subset of some atom
  6. atom has no neighbors
  7. atom has neighbors with index 1
  then we have:
  neighbors_number=(3,2,2,2,-1,0,1)
  index_start=(1,4,6,8,10,10,10,11)
  neighbors_indices(2,4,7,1,3,2,4,1,3,1)*/
  void MakeNeighbors(void);

  void Cleanup();

  // Computing surface area and volume of the overlapping spheres
  int compute_volume_arvo(double *area1, double *vol1, double *test_arr, int sphsize);
  Arvo_class();
  ~Arvo_class();
};

} // namespace arvo

#endif // ARVO_FUNCTIONS
