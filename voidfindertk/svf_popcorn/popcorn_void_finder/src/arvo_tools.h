#ifndef ARVO_TOOL_H_
#define ARVO_TOOL_H_

namespace arvo {
// fraction evaluation for integral
double Fract(double A, double B, double C, double sinphi, double cosphi, double k);

void PrintUsage();

void NorthPoleFix(const int numAtoms, double *sphereLocal);

// computing integrals over arcs given in arc structure
// according to paper Hayrian, Dzurina, Plavka, Busa
void AvIntegral(const int nArcs, double *pVolume, double *pArea, const double r1, const double z1,
                const double *circles, const double *arcs);
} // namespace arvo

#endif
