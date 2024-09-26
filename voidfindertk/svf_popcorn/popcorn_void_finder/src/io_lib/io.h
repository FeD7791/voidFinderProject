#ifndef IO_H
#define IO_H

#include <string>
#include <vector>

#include "../objects.h"

void ReadTracers(std::string fmt, std::string fname, int nfile);
void ReadHeaderTracers(std::string fmt, std::string fname, int nfile);
void Write_sv(std::string fname);
void Write_raw_popcorn(std::string outFile);
void Read_raw_popcorn(std::string pop_file);
void Read_sv(std::string fname_rootvds, float3 box_size);
void Write_pairs(std::string outFile);
void Read_pairs(std::string inFile);
void Write_clean_popcorn(std::string pop_file, std::vector<bool> *erase);

#endif