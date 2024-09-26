#include "io_tracers_gadget.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "../allvars.h"
#include "../grid.h"
#include "../t.h"

using namespace grid;
using std::cout;
using std::endl;

void ReadTracers_Gadget1(std::string fname, int nfile) {

  struct GadgetHeader {
    int Npart[6];
    double Mass[6];
    double Time;
    double Redshift;
    int Flag_1;
    int Flag_2;
    int NpartTotal[6];
    int Flag_3;
    int NumFiles;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[96];
  } Header;

  int j, SkipBlock, dummy, Np, id, NC;
  float pos[3], vel[3];
  FILE *f1, *f2;
  char snapshot[500];

  int *nstart;
  int *npartfile;
  halo dumm;

  nstart = (int *)malloc(nfile * sizeof(int));
  npartfile = (int *)malloc(nfile * sizeof(int));

  if (nfile == 1) {
    sprintf(snapshot, "%s", fname.c_str());
  } else {
    sprintf(snapshot, "%s.0", fname.c_str());
    cout << "Nfiles= " << nfile << endl;
  }

  f1 = fopen(snapshot, "r");

  fread(&dummy, sizeof(int), 1, f1);
  fread(&Header, sizeof(struct GadgetHeader), 1, f1);
  fclose(f1);

  NTRAC = Header.NpartTotal[1];
  boxsize.x = (float)Header.BoxSize;
  boxsize.y = (float)Header.BoxSize;
  boxsize.z = (float)Header.BoxSize;

  cout << "Boxsize " << boxsize.x << endl;

  dumm.x = (-100.0);
  dumm.y = (-100.0);
  dumm.z = (-100.0);
  for (long long int i = 0; i < NTRAC; i++)
    Tracer.push_back(dumm);

  if (nfile < NCORES)
    NC = nfile;
  else
    NC = NCORES;

  for (long long int i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%lld", fname.c_str(), i);
    }

    f1 = fopen(snapshot, "r"); // Pos

    fread(&dummy, sizeof(int), 1, f1);
    fread(&Header, sizeof(struct GadgetHeader), 1, f1);
    npartfile[i] = Header.Npart[1];
    fclose(f1);
  }

  for (long long int i = 0; i < nfile; i++) {
    nstart[i] = 0;
    for (j = 0; j < i; j++)
      nstart[i] = nstart[i] + npartfile[j];
  }

#pragma omp parallel for default(none) schedule(static)                                                                \
    num_threads(NC) private(snapshot, f1, f2, Np, Header, SkipBlock, pos, vel, id, j, dummy)                           \
    shared(Tracer, stdout, fname, nfile, nstart, cout, boxsize)
  for (int i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d", fname.c_str(), i);
    }

    f1 = fopen(snapshot, "r"); // Pos
    f2 = fopen(snapshot, "r"); // Vel

    fread(&dummy, sizeof(int), 1, f1);
    fread(&Header, sizeof(struct GadgetHeader), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);

    Np = Header.Npart[1];

    SkipBlock = sizeof(int) + sizeof(struct GadgetHeader) + sizeof(int) + sizeof(int) + 3 * sizeof(float) * Np +
                2 * sizeof(int);
    fseek(f2, SkipBlock, SEEK_CUR);

    id = nstart[i];
    for (j = 0; j < Np; j++) {
      fread(pos, sizeof(float), 3, f1);
      fread(vel, sizeof(float), 3, f2);

      Tracer[id].x = pos[0];
      Tracer[id].y = pos[1];
      Tracer[id].z = pos[2];

      Tracer[id].vx = vel[0] * sqrt(Header.Time);
      Tracer[id].vy = vel[1] * sqrt(Header.Time);
      Tracer[id].vz = vel[2] * sqrt(Header.Time);
      id++;
    }
    fclose(f1);
    fclose(f2);
  }

  free(nstart);
  free(npartfile);
}

void ReadTracers_Gadget4_type1(std::string fname, int nfile) {

  struct GadgetHeader {
    unsigned long int Npart[2];
    unsigned long long int NpartTotal[2];
    double Mass[2];
    double Time;
    double Redshift;
    double BoxSize;
    int NumFiles;
  } Header;

  int i, j, SkipBlock, dummy, Np, id, NC;
  float pos[3], vel[3];
  FILE *f1, *f2;
  char snapshot[500];
  int *nstart;
  int *npartfile;
  halo dumm;
  long int head_skyp;

  nstart = (int *)malloc(nfile * sizeof(int));
  npartfile = (int *)malloc(nfile * sizeof(int));

  if (nfile == 1) {
    sprintf(snapshot, "%s", fname.c_str());
  } else {
    sprintf(snapshot, "%s.0", fname.c_str());
    cout << "Nfiles= " << nfile << endl;
  }

  f1 = fopen(snapshot, "r");
  fread(&dummy, sizeof(int), 1, f1);
  fread(&Header, sizeof(struct GadgetHeader), 1, f1);
  fclose(f1);

  head_skyp = (long int)dummy - (long int)sizeof(struct GadgetHeader);

  cout << "Npart: ";
  for (i = 0; i < 2; ++i) {
    cout << Header.Npart[i] << " ";
  }
  cout << endl;

  cout << "NpartTotal: ";
  for (i = 0; i < 2; ++i) {
    cout << Header.NpartTotal[i] << " ";
  }
  cout << endl;

  cout << "Mass: ";
  for (i = 0; i < 2; ++i) {
    cout << Header.Mass[i] << " ";
    cout << "NpartTotal: ";
    for (i = 0; i < 2; ++i) {
      cout << Header.NpartTotal[i] << " ";
    }
    cout << endl;
  }
  cout << endl;

  NTRAC = Header.NpartTotal[1];
  boxsize.x = (float)Header.BoxSize;
  boxsize.y = (float)Header.BoxSize;
  boxsize.z = (float)Header.BoxSize;

  cout << "Boxsize X " << boxsize.x << endl;
  cout << "Boxsize Y " << boxsize.y << endl;
  cout << "Boxsize Z " << boxsize.z << endl;

  dumm.x = (-100.0);
  dumm.y = (-100.0);
  dumm.z = (-100.0);
  for (i = 0; i < NTRAC; i++)
    Tracer.push_back(dumm);

  if (nfile < NCORES)
    NC = nfile;
  else
    NC = NCORES;

  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d", fname.c_str(), i);
    }

    f1 = fopen(snapshot, "r"); // Pos

    fread(&dummy, sizeof(int), 1, f1);
    fread(&Header, sizeof(struct GadgetHeader), 1, f1);
    npartfile[i] = Header.Npart[1];
    fclose(f1);
  }

  for (i = 0; i < nfile; i++) {
    nstart[i] = 0;
    for (j = 0; j < i; j++)
      nstart[i] = nstart[i] + npartfile[j];
  }

#pragma omp parallel for default(none) schedule(static)                                                                \
    num_threads(NC) private(i, snapshot, f1, f2, Np, Header, SkipBlock, pos, vel, id, j, dummy)                        \
    shared(Tracer, stdout, fname, nfile, nstart, cout, boxsize, head_skyp)
  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d", fname.c_str(), i);
    }

    f1 = fopen(snapshot, "r"); // Pos
    f2 = fopen(snapshot, "r"); // Vel

    fread(&dummy, sizeof(int), 1, f1);
    fread(&Header, sizeof(struct GadgetHeader), 1, f1);
    fseek(f1, head_skyp, SEEK_CUR);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);

    Np = Header.Npart[1];

    SkipBlock = sizeof(int) + sizeof(struct GadgetHeader) + sizeof(int) + sizeof(int) + 3 * sizeof(float) * Np +
                2 * sizeof(int);
    fseek(f2, SkipBlock, SEEK_CUR);

    id = nstart[i];
    for (j = 0; j < Np; j++) {
      fread(pos, sizeof(float), 3, f1);
      fread(vel, sizeof(float), 3, f2);

      Tracer[id].x = pos[0];
      Tracer[id].y = pos[1];
      Tracer[id].z = pos[2];

      Tracer[id].vx = vel[0] * sqrt(Header.Time);
      Tracer[id].vy = vel[1] * sqrt(Header.Time);
      Tracer[id].vz = vel[2] * sqrt(Header.Time);
      id++;
    }
    fclose(f1);
    fclose(f2);
  }

  free(nstart);
  free(npartfile);
}

void ReadTracers_Gadget2(std::string fname, int nfile) {

  struct GadgetHeader {
    int Npart[6];
    double Mass[6];
    double Time;
    double Redshift;
    int Flag_1;
    int Flag_2;
    int NpartTotal[6];
    int Flag_3;
    int NumFiles;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[96];
  } Header;

  int i, j, dummy, Np, id, NC;
  float pos[3], vel[3];
  FILE *f1, *f2;
  char snapshot[500];

  char buffer[5];
  char key[5];
  int sizeblck;
  int *npart_file;
  int *nstart;

  npart_file = (int *)malloc(nfile * sizeof(int));
  nstart = (int *)malloc(nfile * sizeof(int));

  if (nfile == 1) {
    sprintf(snapshot, "%s", fname.c_str());
  } else {
    sprintf(snapshot, "%s.0", fname.c_str());
  }

  f1 = fopen(snapshot, "r");
  strcpy(key, "HEAD\0");
  do {
    fread(&dummy, sizeof(int), 1, f1);
    fread(buffer, 4 * sizeof(char), 1, f1);
    buffer[4] = '\0';
    fread(&sizeblck, sizeof(int), 1, f1);
    if (strcmp(key, buffer) == 0)
      break;
    fseek(f1, sizeblck + sizeof(int), SEEK_CUR);
  } while (1);
  fread(&dummy, sizeof(int), 1, f1);
  fread(&dummy, sizeof(int), 1, f1);
  fread(&Header, sizeof(struct GadgetHeader), 1, f1);
  fread(&dummy, sizeof(int), 1, f1);
  fclose(f1);

  NTRAC = Header.NpartTotal[1];
  boxsize.x = (float)Header.BoxSize;
  boxsize.y = (float)Header.BoxSize;
  boxsize.z = (float)Header.BoxSize;
  cout << "Boxsize: " << Header.BoxSize << endl;
  // Tracer = (struct tracers *) malloc(NTRAC*sizeof(struct tracers));
  Tracer.reserve(NTRAC);

  if (nfile < NCORES)
    NC = nfile;
  else
    NC = NCORES;

  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d", fname.c_str(), i);
    }

    f1 = fopen(snapshot, "r"); // Pos

    strcpy(key, "HEAD\0");
    do {
      fread(&dummy, sizeof(int), 1, f1);
      fread(buffer, 4 * sizeof(char), 1, f1);
      buffer[4] = '\0';
      fread(&sizeblck, sizeof(int), 1, f1);
      if (strcmp(key, buffer) == 0)
        break;
      fseek(f1, sizeblck + sizeof(int), SEEK_CUR);
    } while (1);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&Header, sizeof(struct GadgetHeader), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);

    npart_file[i] = Header.Npart[1];
  }

  for (i = 0; i < nfile; i++) {
    nstart[i] = 0;
    for (j = 0; j < i; j++)
      nstart[i] = nstart[i] + npart_file[j];
  }

#pragma omp parallel for default(none) schedule(static)                                                                \
    num_threads(NC) private(i, snapshot, f1, f2, Np, Header, pos, vel, id, j, dummy, key, buffer, sizeblck)            \
    shared(Tracer, stdout, fname, nfile, cout, nstart)
  for (i = 0; i < nfile; i++) {

    if (nfile == 1) {
      sprintf(snapshot, "%s", fname.c_str());
    } else {
      sprintf(snapshot, "%s.%d", fname.c_str(), i);
    }

    f1 = fopen(snapshot, "r"); // Pos
    f2 = fopen(snapshot, "r"); // Vel

    strcpy(key, "HEAD\0");
    do {
      fread(&dummy, sizeof(int), 1, f1);
      fread(buffer, 4 * sizeof(char), 1, f1);
      buffer[4] = '\0';
      fread(&sizeblck, sizeof(int), 1, f1);
      if (strcmp(key, buffer) == 0)
        break;
      fseek(f1, sizeblck + sizeof(int), SEEK_CUR);
    } while (1);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&Header, sizeof(struct GadgetHeader), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);

    Np = Header.Npart[1];

    strcpy(key, "POS \0");
    do {
      fread(&dummy, sizeof(int), 1, f1);
      fread(buffer, 4 * sizeof(char), 1, f1);
      buffer[4] = '\0';
      fread(&sizeblck, sizeof(int), 1, f1);
      if (strcmp(key, buffer) == 0)
        break;
      fseek(f1, sizeblck + sizeof(int), SEEK_CUR);
    } while (1);
    fread(&dummy, sizeof(int), 1, f1);
    fread(&dummy, sizeof(int), 1, f1);

    strcpy(key, "VEL \0");
    do {
      fread(&dummy, sizeof(int), 1, f2);
      fread(buffer, 4 * sizeof(char), 1, f2);
      buffer[4] = '\0';
      fread(&sizeblck, sizeof(int), 1, f2);
      if (strcmp(key, buffer) == 0)
        break;
      fseek(f2, sizeblck + sizeof(int), SEEK_CUR);
    } while (1);
    fread(&dummy, sizeof(int), 1, f2);
    fread(&dummy, sizeof(int), 1, f2);

    fseek(f1, sizeof(float) * Header.Npart[1] * 3, SEEK_CUR);
    fseek(f2, sizeof(float) * Header.Npart[1] * 3, SEEK_CUR);

    id = nstart[i];
    for (j = 0; j < Np; j++) {
      fread(pos, sizeof(float), 3, f1);
      fread(vel, sizeof(float), 3, f2);
      Tracer[id].x = pos[0];
      Tracer[id].y = pos[1];
      Tracer[id].z = pos[2];

      Tracer[id].vx = vel[0] * sqrt(Header.Time);
      Tracer[id].vy = vel[1] * sqrt(Header.Time);
      Tracer[id].vz = vel[2] * sqrt(Header.Time);
      id++;
    }
    fclose(f1);
    fclose(f2);
  }

  free(npart_file);
  free(nstart);
}

void ReadHeaderTracers_Gadget1(std::string fname, int nfile) {

  struct GadgetHeader {
    int Npart[6];
    double Mass[6];
    double Time;
    double Redshift;
    int Flag_1;
    int Flag_2;
    int NpartTotal[6];
    int Flag_3;
    int NumFiles;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[96];
  } Header;

  int dummy;
  FILE *f1;
  char snapshot[500];
  int *npart_file;

  npart_file = (int *)malloc(nfile * sizeof(int));

  if (nfile == 1) {
    sprintf(snapshot, "%s", fname.c_str());
  } else {
    sprintf(snapshot, "%s.0", fname.c_str());
  }

  f1 = fopen(snapshot, "r");

  fread(&dummy, sizeof(int), 1, f1);
  fread(&Header, sizeof(struct GadgetHeader), 1, f1);
  fclose(f1);

  boxsize.x = Header.BoxSize;
  boxsize.y = Header.BoxSize;
  boxsize.z = Header.BoxSize;
  NTRAC = Header.NpartTotal[1];

  free(npart_file);
}

void ReadHeaderTracers_Gadget4_type1(std::string fname, int nfile) {
  struct GadgetHeader {
    unsigned long int Npart[2];
    unsigned long long int NpartTotal[2];
    double Mass[2];
    double Time;
    double Redshift;
    double BoxSize;
    int NumFiles;
  } Header;

  int dummy;
  FILE *f1;
  char snapshot[500];
  int *npart_file;

  npart_file = (int *)malloc(nfile * sizeof(int));

  if (nfile == 1) {
    sprintf(snapshot, "%s", fname.c_str());
  } else {
    sprintf(snapshot, "%s.0", fname.c_str());
  }

  f1 = fopen(snapshot, "r");

  fread(&dummy, sizeof(int), 1, f1);
  fread(&Header, sizeof(struct GadgetHeader), 1, f1);
  fclose(f1);

  boxsize.x = Header.BoxSize;
  boxsize.y = Header.BoxSize;
  boxsize.z = Header.BoxSize;
  NTRAC = Header.NpartTotal[1];

  free(npart_file);
}
void ReadHeaderTracers_Gadget2(std::string fname, int nfile) {

  struct GadgetHeader {
    int Npart[6];
    double Mass[6];
    double Time;
    double Redshift;
    int Flag_1;
    int Flag_2;
    int NpartTotal[6];
    int Flag_3;
    int NumFiles;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[96];
  } Header;

  int dummy;
  FILE *f1;
  char snapshot[500];

  char buffer[5];
  char key[5];
  int sizeblck;
  int *npart_file;

  npart_file = (int *)malloc(nfile * sizeof(int));

  if (nfile == 1) {
    sprintf(snapshot, "%s", fname.c_str());
  } else {
    sprintf(snapshot, "%s.0", fname.c_str());
  }

  f1 = fopen(snapshot, "r");
  strcpy(key, "HEAD\0");
  do {
    fread(&dummy, sizeof(int), 1, f1);
    fread(buffer, 4 * sizeof(char), 1, f1);
    buffer[4] = '\0';
    fread(&sizeblck, sizeof(int), 1, f1);
    if (strcmp(key, buffer) == 0)
      break;
    fseek(f1, sizeblck + sizeof(int), SEEK_CUR);
  } while (1);
  fread(&dummy, sizeof(int), 1, f1);
  fread(&dummy, sizeof(int), 1, f1);
  fread(&Header, sizeof(struct GadgetHeader), 1, f1);
  fread(&dummy, sizeof(int), 1, f1);
  fclose(f1);

  NTRAC = Header.NpartTotal[1];
  boxsize.x = Header.BoxSize;
  boxsize.y = Header.BoxSize;
  boxsize.z = Header.BoxSize;
  cout << "Boxsize: " << Header.BoxSize << endl;

  free(npart_file);
}
