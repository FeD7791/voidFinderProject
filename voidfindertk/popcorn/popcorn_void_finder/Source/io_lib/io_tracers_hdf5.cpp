#include "io_tracers_hdf5.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "hdf5.h"

#include "../allvars.h"
#include "../grid.h"
#include "../t.h"

using namespace grid;
using std::cout;
using std::endl;

void ReadTracers_HDF5_swift(std::string fname, int nfile) {

  struct GadgetHeader {
    int Npart[7];
    double Mass[7];
    double Time;
    double Redshift;
    int NpartTotal[7];
    int NumFiles;
    double BoxSize[3];
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
  } Header;

  int i = 0, j, Npart;
  char snapshot[500];
  float *pos, *vel;
  int *id;
  hid_t file_id, dataset_id, group_id, attribute_id;
  halo TT;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "NumPart_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, Header.NpartTotal);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, Header.BoxSize);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Redshift");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.Redshift);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Time");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.Time);
  H5Aclose(attribute_id);

  H5Gclose(group_id);

  if ((float)Header.BoxSize[0] != (float)boxsize.x && (float)Header.BoxSize[0] / 1000.0 != (float)boxsize.x) {
    fprintf(stdout, " | ERROR!! boxsize = %f - BoxSize = %f \n", boxsize.x, Header.BoxSize[0]);
    exit(EXIT_FAILURE);
  }
  /// TEST en x y z
  NTRAC = Header.NpartTotal[1];

  // Tracer = (struct tracers *) malloc(NTRAC*sizeof(struct tracers));
  Tracer.reserve(NTRAC);

  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);
    }

    file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

    attribute_id = H5Aopen_name(group_id, "NumPart_ThisFile");
    H5Aread(attribute_id, H5T_NATIVE_INT, Header.Npart);
    H5Aclose(attribute_id);

    Npart = Header.Npart[1];

    H5Gclose(group_id);

    pos = (float *)malloc(3 * Npart * sizeof(float));
    vel = (float *)malloc(3 * Npart * sizeof(float));
    id = (int *)malloc(Npart * sizeof(int));

    dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/PartType1/Velocities", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/PartType1/ParticleIDs", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
    H5Dclose(dataset_id);

    H5Fclose(file_id);

    for (j = 0; j < Npart; j++) {
      TT.x = pos[3 * j];
      TT.y = pos[3 * j + 1];
      TT.z = pos[3 * j + 2];

      TT.vx = vel[3 * j] * sqrt(Header.Time);
      TT.vy = vel[3 * j + 1] * sqrt(Header.Time);
      TT.vz = vel[3 * j + 2] * sqrt(Header.Time);
      Tracer.push_back(TT);
    }

    free(pos);
    free(vel);
    free(id);
  }
}
void ReadTracers_HDF5(std::string fname, int nfile) {

  struct GadgetHeader {
    int Npart[6];
    double Mass[6];
    double Time;
    double Redshift;
    int NpartTotal[6];
    int NumFiles;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
  } Header;

  int i = 0, j, Npart, pid;
  char snapshot[500];
  float *pos, *vel;
  int *id;
  hid_t file_id, dataset_id, group_id, attribute_id;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "NumPart_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, Header.NpartTotal);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.BoxSize);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Omega0");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.Omega0);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "OmegaLambda");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.OmegaLambda);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Redshift");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.Redshift);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Time");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Header.Time);
  H5Aclose(attribute_id);

  H5Gclose(group_id);

  if (Header.BoxSize != boxsize.x && Header.BoxSize / 1000.0 != boxsize.x) {
    fprintf(stdout, " | ERROR!! boxsize = %f - BoxSize = %f \n", boxsize.x, Header.BoxSize);
    exit(EXIT_FAILURE);
  }

  NTRAC = Header.NpartTotal[1];
  // Tracer = (struct tracers *) malloc(NTRAC*sizeof(struct tracers));
  Tracer.reserve(NTRAC);

  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);
    }

    file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

    attribute_id = H5Aopen_name(group_id, "NumPart_ThisFile");
    H5Aread(attribute_id, H5T_NATIVE_INT, Header.Npart);
    H5Aclose(attribute_id);

    Npart = Header.Npart[1];

    H5Gclose(group_id);

    pos = (float *)malloc(3 * Npart * sizeof(float));
    vel = (float *)malloc(3 * Npart * sizeof(float));
    id = (int *)malloc(Npart * sizeof(int));

    dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/PartType1/Velocities", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/PartType1/ParticleIDs", H5P_DEFAULT);
    H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
    H5Dclose(dataset_id);

    H5Fclose(file_id);

    for (j = 0; j < Npart; j++) {
      pid = id[j] - 1;
      Tracer[pid].x = pos[3 * j];
      Tracer[pid].y = pos[3 * j + 1];
      Tracer[pid].z = pos[3 * j + 2];

      Tracer[pid].vx = vel[3 * j] * sqrt(Header.Time);
      Tracer[pid].vy = vel[3 * j + 1] * sqrt(Header.Time);
      Tracer[pid].vz = vel[3 * j + 2] * sqrt(Header.Time);
    }

    free(pos);
    free(vel);
    free(id);
  }
}

void ReadTracers_HDF5_subfind_groups(std::string fname, int nfile, float mass_lim) {
  int i = 0, j, Npart;
  char snapshot[500];
  float *pos, *vel, *mass;
  hid_t file_id, dataset_id, group_id, attribute_id;
  double Redshift;
  double Time;
  double BoxSize;
  halo TT;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "Ngroups_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, &Npart);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Redshift");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Redshift);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Time");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Time);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &BoxSize);
  H5Aclose(attribute_id);

  if (BoxSize != boxsize.x && BoxSize / 1000.0 != boxsize.x) {
    fprintf(stdout, " | ERROR!! boxsize = %f - BoxSize = %f \n", boxsize.x, BoxSize);
    exit(EXIT_FAILURE);
  }

  H5Gclose(group_id);
  H5Fclose(file_id);

  NTRAC = Npart;
  // Tracer = (struct tracers *) malloc(NTRAC*sizeof(struct tracers));
  // Tracer.reserve(NTRAC);
  int npart_eff = 0;

  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);
    }

    file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

    attribute_id = H5Aopen_name(group_id, "Ngroups_ThisFile");
    H5Aread(attribute_id, H5T_NATIVE_INT, &Npart);
    H5Aclose(attribute_id);
    H5Gclose(group_id);
    if (Npart > 0) {
      pos = (float *)malloc(3 * Npart * sizeof(float));
      vel = (float *)malloc(3 * Npart * sizeof(float));
      mass = (float *)malloc(Npart * sizeof(float));

      dataset_id = H5Dopen2(file_id, "/Group/GroupPos", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
      H5Dclose(dataset_id);

      dataset_id = H5Dopen2(file_id, "/Group/GroupVel", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel);
      H5Dclose(dataset_id);

      dataset_id = H5Dopen2(file_id, "/Group/GroupMass", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
      H5Dclose(dataset_id);
      H5Fclose(file_id);

      for (j = 0; j < Npart; j++) {
        if (mass[j] > mass_lim) {
          TT.x = pos[3 * j];
          TT.y = pos[3 * j + 1];
          TT.z = pos[3 * j + 2];

          TT.vx = vel[3 * j] * sqrt(Time);
          TT.vy = vel[3 * j + 1] * sqrt(Time);
          TT.vz = vel[3 * j + 2] * sqrt(Time);
          Tracer.push_back(TT);
          npart_eff++;
        }
      }

      free(mass);
      free(pos);
      free(vel);
    }
  }

  NTRAC = npart_eff;
}

void ReadTracers_HDF5_subfind_subhalos(std::string fname, int nfile, float mass_lim) {
  int i = 0, j, Npart;
  char snapshot[500];
  float *pos, *vel, *mass;
  hid_t file_id, dataset_id, group_id, attribute_id;
  double Redshift;
  double Time;
  double BoxSize;
  halo TT;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "Nsubgroups_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, &Npart);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Redshift");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Redshift);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Time");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Time);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &BoxSize);
  H5Aclose(attribute_id);

  if (BoxSize != boxsize.x && BoxSize / 1000.0 != boxsize.x) {
    fprintf(stdout, " | ERROR!! boxsize = %f - BoxSize = %f \n", boxsize.x, BoxSize);
    exit(EXIT_FAILURE);
  }

  H5Gclose(group_id);
  H5Fclose(file_id);

  NTRAC = Npart;
  // Tracer = (struct tracers *) malloc(NTRAC*sizeof(struct tracers));
  // Tracer.reserve(NTRAC);
  int npart_eff = 0;

  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);
    }

    file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

    attribute_id = H5Aopen_name(group_id, "Nsubgroups_ThisFile");
    H5Aread(attribute_id, H5T_NATIVE_INT, &Npart);
    H5Aclose(attribute_id);
    H5Gclose(group_id);
    if (Npart > 0) {
      pos = (float *)malloc(3 * Npart * sizeof(float));
      vel = (float *)malloc(3 * Npart * sizeof(float));
      mass = (float *)malloc(Npart * sizeof(float));

      dataset_id = H5Dopen2(file_id, "/Subhalo/SubhaloPos", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
      H5Dclose(dataset_id);

      dataset_id = H5Dopen2(file_id, "/Subhalo/SubhaloVel", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel);
      H5Dclose(dataset_id);

      dataset_id = H5Dopen2(file_id, "/Subhalo/SubhaloMass", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
      H5Dclose(dataset_id);
      H5Fclose(file_id);

      for (j = 0; j < Npart; j++) {
        if (mass[j] > mass_lim) {
          TT.x = pos[3 * j];
          TT.y = pos[3 * j + 1];
          TT.z = pos[3 * j + 2];

          TT.vx = vel[3 * j] * sqrt(Time);
          TT.vy = vel[3 * j + 1] * sqrt(Time);
          TT.vz = vel[3 * j + 2] * sqrt(Time);
          Tracer.push_back(TT);
          npart_eff++;
        }
      }

      free(mass);
      free(pos);
      free(vel);
    }
  }

  NTRAC = npart_eff;
}

void ReadHeaderTracers_HDF5_subfind_groups(std::string fname, int nfile) {
  int i = 0, Npart;
  char snapshot[500];
  hid_t file_id, group_id, attribute_id;
  double Redshift;
  double Time;
  double BoxSize;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  cout << fname.c_str() << endl;

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "Ngroups_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, &Npart);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Redshift");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Redshift);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Time");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Time);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &BoxSize);
  H5Aclose(attribute_id);

  boxsize.x = boxsize.y = boxsize.z = BoxSize;

  H5Gclose(group_id);
  NTRAC = Npart;
  H5Fclose(file_id);
}

void ReadHeaderTracers_HDF5_subfind_subhalos(std::string fname, int nfile) {
  int i = 0, Npart;
  char snapshot[500];
  hid_t file_id, group_id, attribute_id;
  double Redshift;
  double Time;
  double BoxSize;

  i = 0;
  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  cout << fname.c_str() << endl;

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "Nsubgroups_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, &Npart);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Redshift");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Redshift);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "Time");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &Time);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &BoxSize);
  H5Aclose(attribute_id);

  boxsize.x = boxsize.y = boxsize.z = BoxSize;

  H5Gclose(group_id);
  NTRAC = Npart;
  H5Fclose(file_id);
}

void ReadHeaderTracers_HDF5(std::string fname, int nfile) {

  int NpartTotal[6];
  double BoxSize;
  int i = 0;
  char snapshot[500];
  hid_t file_id, group_id, attribute_id;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "NumPart_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, NpartTotal);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &BoxSize);
  H5Aclose(attribute_id);

  H5Gclose(group_id);
  boxsize.x = BoxSize;
  boxsize.y = BoxSize;
  boxsize.z = BoxSize;
  NTRAC = NpartTotal[1];
}
void ReadHeaderTracers_HDF5_swift(std::string fname, int nfile) {

  int NpartTotal[7];
  double BoxSize[3];
  int i = 0;
  char snapshot[500];
  hid_t file_id, group_id, attribute_id;

  if (nfile == 1)
    sprintf(snapshot, "%s", fname.c_str());
  else
    sprintf(snapshot, "%s.%d.hdf5", fname.c_str(), i);

  file_id = H5Fopen(snapshot, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attribute_id = H5Aopen_name(group_id, "NumPart_Total");
  H5Aread(attribute_id, H5T_NATIVE_INT, NpartTotal);
  H5Aclose(attribute_id);

  attribute_id = H5Aopen_name(group_id, "BoxSize");
  H5Aread(attribute_id, H5T_NATIVE_DOUBLE, BoxSize);
  H5Aclose(attribute_id);

  H5Gclose(group_id);
  boxsize.x = BoxSize[0];
  boxsize.y = BoxSize[1];
  boxsize.z = BoxSize[2];
  NTRAC = NpartTotal[1];
  printf("111111111111 %lli %f\n", NTRAC, boxsize.x);
}
