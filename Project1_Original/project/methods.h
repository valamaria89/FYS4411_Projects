#ifndef METHODS_H
#define METHODS_H
#include <random>
#include <iostream>
#include "matrix.h"

void mc_sampling( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha);
double wavefunction(Matrix r, double alpha, int dim, int NumOfPart);
double local_energy_num(Matrix r, double alpha, double wfold, int dim, int numOfPart);
double local_energy_analytic(Matrix r, double alpha, int dim, int numOfPart);
double deriv_num(Matrix r, double alpha, double wfold, int dim, int numOfPart);
double deriv_analytic(Matrix r, double alpha, int dim, int numOfPart);
void mc_sampling_IMS( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha);
void QuantumForce(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha);

class methods
{
public:
    methods();
};

#endif // METHODS_H
