#ifndef IO_TRACERS_HDF5_H
#define IO_TRACERS_HDF5_H

#include <string>

void ReadTracers_HDF5(std::string fname, int nfile);
void ReadTracers_HDF5_swift(std::string fname, int nfile);
void ReadTracers_HDF5_subfind_groups(std::string fname, int nfile, float mass_lim);
void ReadTracers_HDF5_subfind_subhalos(std::string fname, int nfile, float mass_lim);

void ReadHeaderTracers_HDF5_subfind_groups(std::string fname, int nfile);
void ReadHeaderTracers_HDF5_subfind_subhalos(std::string fname, int nfile);
void ReadHeaderTracers_HDF5(std::string fname, int nfile);
void ReadHeaderTracers_HDF5_swift(std::string fname, int nfile);

#endif
