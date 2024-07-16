#include "clean_duplicates.h"

#include <algorithm>
#include <vector>

#include "allvars.h"
#include "colors.h"
#include "grid.h"
#include "objects.h"

using namespace grid;

struct small_first {
  bool operator()(par const &a, par const &b) const { return a.minvol < b.minvol; }
};

struct small_root_first {
  bool operator()(par const &a, par const &b) const { return a.root_rad_min < b.root_rad_min; }
};
int clean_duplicates_start(float overlap, std::vector<bool> *erase) {
  int id;
  int cnt = -1;

  for (long unsigned int i = 0; i < pairs.size(); i++) {
    id = pairs[i].idmin;
    cnt = (id > cnt) ? id : cnt;
    id = pairs[i].idmax;
    cnt = (id > cnt) ? id : cnt;
  }
  cnt++; // quiero el size, para contener hasta el maximo indice

  (*erase).assign(cnt, false);
  // hacemos un sort de volumenes minimos y comencemos borrando por el menor
  // volumen, como en svf
  // sort(pairs.begin(),pairs.end(),small_first());

  /// cambio de estretegia... ordenamos no por volumen total, si no por el
  /// tamaño del void esférico de semilla (comenzar con un void chico puede ser
  /// muy de puente)
  std::sort(pairs.begin(), pairs.end(), small_root_first());
  std::cout << "numero de popcorns antes de la limpieza...  " << pop_voids.size() << std::endl;

  std::cout << CYAN << "procesando.... " << DEFA << std::endl;
  for (long unsigned int i = 0; i < pairs.size(); i++) {
    // if( (*erase)[pairs[i].idmin] ||
    //    (*erase)[pairs[i].idmax]) continue;
    if ((*erase)[pairs[i].id_rad_min] || (*erase)[pairs[i].id_rad_max])
      continue;

    // if(pairs[i].volint/pairs[i].minvol > overlap)
    if (pairs[i].volint > overlap) {
      //  (*erase)[pairs[i].idmin]=1;

      /// cambio de estretegia... ordenamos no por volumen total, si no por el
      /// tamaño del void esférico de semilla (comenzar con un void chico puede
      /// ser muy de puente)
      (*erase)[pairs[i].id_rad_min] = true;
    }
  }

  return 0;
}
