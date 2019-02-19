#ifndef METHODS_H
#define METHODS_H
#include <random>
#include <iostream>
#include "matrix.h"
#include <armadillo>
using namespace arma;

void mc_sampling( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha, int Thermalization);
double wavefunction(mat r, double alpha, int dim, int NumOfPart);
double local_energy_num(mat r, double alpha, double wfold, int dim, int numOfPart);
double local_energy_analytic(mat r, double alpha, int dim, int numOfPart);
double deriv_num(mat r, double alpha, double wfold, int dim, int numOfPart);
double deriv_analytic(mat r, double alpha, int dim, int numOfPart);
void mc_sampling_IMS( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha, int thermalization);
void QuantumForce(mat r, mat QForce, int dim, int numOfPart, double alpha);
void writeToFile( string x, double *E1, double *E2, int numVar, double alpha, double deltaAlpha);
class methods
{
public:
    methods();
};

#endif // METHODS_H
