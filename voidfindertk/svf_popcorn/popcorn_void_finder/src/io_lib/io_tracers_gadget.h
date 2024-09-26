#ifndef IO_TRACERS_GADGET_H
#define IO_TRACERS_GADGET_H

#include <string>

void ReadTracers_Gadget1(std::string fname, int nfile);
void ReadTracers_Gadget4_type1(std::string fname, int nfile);
void ReadTracers_Gadget2(std::string fname, int nfile);

void ReadHeaderTracers_Gadget1(std::string fname, int nfile);
void ReadHeaderTracers_Gadget4_type1(std::string fname, int nfile);
void ReadHeaderTracers_Gadget2(std::string fname, int nfile);

#endif