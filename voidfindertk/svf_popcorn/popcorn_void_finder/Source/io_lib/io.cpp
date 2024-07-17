#include "io.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "hdf5.h"

#include "io_tracers_ascii.h"
#include "io_tracers_gadget.h"
#include "io_tracers_hdf5.h"
#include "io_tracers_stream.h"

#include "../allvars.h"
#include "../colors.h"
#include "../grid.h"
#include "../t.h"

using namespace grid;
using std::cout;
using std::endl;

void ReadTracers(std::string fmt, std::string fname, int nfile) {

  std::cout << fmt << std::endl;
  if (!fmt.compare("ASCII")) {
    ReadTracers_ascii(fname);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("STREAM")) {
    ReadTracers_stream(fname);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("HDF5")) {
    ReadTracers_HDF5(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("SWIFT")) {
    ReadTracers_HDF5_swift(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("HDF5_SUBFIND_GROUPS")) {
    ReadTracers_HDF5_subfind_groups(fname, nfile, masslim);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("HDF5_SUBFIND_SUBHALOS")) {
    ReadTracers_HDF5_subfind_subhalos(fname, nfile, masslim);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("GADGET1")) {
    ReadTracers_Gadget1(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("GADGET4_TYPE1")) {
    ReadTracers_Gadget4_type1(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("GADGET2")) {
    ReadTracers_Gadget2(fname, nfile);
    std::cout << fmt << std::endl;
  } else {
    std::cout << "invalid specification of file format" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void ReadHeaderTracers(std::string fmt, std::string fname, int nfile) {

  std::cout << fmt << std::endl;
  if (!fmt.compare("ASCII")) {
    ReadHeaderTracers_ascii(fname);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("STREAM")) {
    ReadHeaderTracers_stream(fname);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("HDF5")) {
    ReadHeaderTracers_HDF5(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("SWIFT")) {
    ReadHeaderTracers_HDF5_swift(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("HDF5_SUBFIND_GROUPS")) {
    ReadHeaderTracers_HDF5_subfind_groups(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("HDF5_SUBFIND_SUBHALOS")) {
    ReadHeaderTracers_HDF5_subfind_subhalos(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("GADGET1")) {
    ReadHeaderTracers_Gadget1(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("GADGET4_TYPE1")) {
    ReadHeaderTracers_Gadget4_type1(fname, nfile);
    std::cout << fmt << std::endl;
  } else if (!fmt.compare("GADGET2")) {
    ReadHeaderTracers_Gadget2(fname, nfile);
    std::cout << fmt << std::endl;
  } else {
    std::cout << "invalid specification of file format" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Write_sv(std::string fname) {

  FILE *fd;
  clock_t t;

  cout << "Writing... " << fname << endl;
  fd = fopen(fname.c_str(), "w");
  t = clock();

  for (long unsigned int i = 0; i < Sphvd.size(); i++) {
    if (Sphvd[i].ToF) {
      fprintf(fd, "%6d %25.15e %25.15e %25.15e %25.15e %25.15e \n", Sphvd[i].id, Sphvd[i].radius, Sphvd[i].x,
              Sphvd[i].y, Sphvd[i].z, Sphvd[i].delta);
    }
  }
  fclose(fd);

  Time(t, 1);
  fprintf(stdout, "\n");
  fflush(stdout);
}

void Read_sv(std::string fname_rootvds, float3 boxsize) {

  // reading spherical voids ////////////////////////////////////////////////
  float rad, d1; // delete

  int id;
  int i;
  Sphere vv;
  std::ifstream file_roots;
  file_roots.open(fname_rootvds.c_str());
  int numcand = 0;
  while (file_roots >> id >> rad >> vv.x >> vv.y >> vv.z >> d1) {
    numcand++;
  }

  file_roots.clear();
  file_roots.seekg(0);

  i = 0;

  while (file_roots >> id >> vv.radius >> vv.x >> vv.y >> vv.z >> d1) {
    vv.x = periodicity(vv.x, boxsize.x);
    vv.y = periodicity(vv.y, boxsize.y);
    vv.z = periodicity(vv.z, boxsize.z);
    vv.id = id;
    vv.erase = false;
    Sphvd.push_back(vv);
    i++;
  }
  file_roots.close();
}

void Write_raw_popcorn(std::string outFile) {
  // output file
  cout << LBLUE << "escribiendo data.... " << YELLOW << pop_voids.size() << " " << Sphvd.size() << DEFA << endl;
  std::ofstream file_out;
  file_out.open(outFile.c_str());
  file_out << pop_voids.size() << endl;

  for (long unsigned int i = 0; i < pop_voids.size(); i++) {
    sort(pop_voids[i].in_halos.begin(), pop_voids[i].in_halos.end());
    file_out << pop_voids[i];
  }

  file_out.close();

  // Restore text format for standard output
  cout << RESTORE << ".";
}

void Read_raw_popcorn(std::string pop_file) {
  std::ifstream file_pops;
  int nvds, i;
  popcorn *seed;

  cout << "Leyendo... " << pop_file << endl;
  file_pops.open(pop_file.c_str());
  file_pops >> nvds;
  for (i = 0; i < nvds; i++) {
    seed = new popcorn;
    file_pops >> (*seed);
    pop_voids.push_back(*seed);
    (*seed).membs.clear();
    (*seed).in_halos.clear();
    delete seed;
  }
  file_pops.close();
}

void Write_pairs(std::string outFile) {
  std::ofstream file_out1;
  cout << "Writting file with the volume intersection of raw popcorn pairs... " << outFile << endl;
  file_out1.open(outFile.c_str());

  for (long unsigned int i = 0; i < pairs.size(); i++) {
    file_out1 << pairs[i].idmin << " " << pairs[i].idmax << " " << pairs[i].root_rad_min << " " << pairs[i].root_rad_max
              << " " << pairs[i].id_rad_min << " " << pairs[i].id_rad_max << " " << pairs[i].volint << " "
              << pairs[i].minvol << " " << pairs[i].pop_Delta << " " << pairs[i].joinvol << " " << pairs[i].npart1
              << " " << pairs[i].npart2 << " " << pairs[i].common << " " << pairs[i].vol1 << " " << pairs[i].vol2 << " "
              << pairs[i].isCave << endl;
  }

  file_out1.close();
}

void Read_pairs(std::string inFile) {
  std::ifstream file_pairs;
  cout << "Reading file with the volume intersection of raw popcorn pairs... " << inFile << endl;
  file_pairs.open(inFile.c_str());
  par p;
  while (file_pairs >> p.idmin >> p.idmax >> p.root_rad_min >> p.root_rad_max >> p.id_rad_min >> p.id_rad_max >>
         p.volint >> p.minvol >> p.pop_Delta >> p.joinvol >> p.npart1 >> p.npart2 >> p.common >> p.vol1 >> p.vol2 >>
         p.isCave) {
    pairs.push_back(p);
  }

  file_pairs.close();
}

void Write_clean_popcorn(std::string pop_file, std::vector<bool> *erase) {
  std::ofstream file_out;

  file_out.open(pop_file.c_str());
  int nsob = 0;
  for (long unsigned int i = 0; i < pop_voids.size(); i++) {
    if ((*erase)[pop_voids[i].id])
      continue;
    if (pop_voids[i].nmem == 0)
      continue;
    nsob++;
  }

  file_out << nsob << endl;
  cout << "Number of popcorns after cleaning overlaping...  " << nsob << endl;

  for (long unsigned int i = 0; i < pop_voids.size(); i++) {
    if ((*erase)[pop_voids[i].id])
      continue;
    if (pop_voids[i].nmem == 0)
      continue;
    sort(pop_voids[i].in_halos.begin(), pop_voids[i].in_halos.end());
    file_out << pop_voids[i];
  }
  file_out.close();

  cout << "Writing output: " << pop_file.c_str() << endl;
}
