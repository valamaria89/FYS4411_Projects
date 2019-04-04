#ifndef METHODS_H
#define METHODS_H
#include "system.h"


void mc_sampling(int numVar, int numOfPart, int dim, int numMCCycles, int ind, double alpha, double step, double delAl, double *Etot, double *Etot2, double be);

#endif // METHODS_H
