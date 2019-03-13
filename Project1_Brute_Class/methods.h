#ifndef METHODS_H
#define METHODS_H
#include "system.h"


double mc_sampling( System pOld, System pNew, double wfOld, int numMCCycles, int ind, double alpha);
void mc_sampling_IMS( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha, int thermalization);
//void writeToFile( string x, double *E1, double *E2, int numVar, double alpha, double deltaAlpha);
void gradiendescent_brute( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha, int Thermalization);

#endif // METHODS_H
