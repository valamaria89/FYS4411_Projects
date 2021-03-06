#ifndef METHODS_INTERACTION_H
#define METHODS_INTERACTION_H

#endif // METHODS_INTERACTION_H

#include <random>
#include <iostream>
#include <fstream>
#include<math.h>
#include "matrix.h"

void mc_sampling_INT(double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization, double beta);
void mc_sampling_IMS_INT(double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization, double beta);
double wavefunction(Matrix r, double alpha, double beta, int dim, int NumOfPart);
double local_energy_num(Matrix r, double alpha, double beta, double wfold, int dim, int numOfPart);
double local_energy_analytic(Matrix r, double alpha, double beta, int dim, int numOfPart);
void QuantumForce_num(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha, double beta);
void QuantumForce(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha, double beta);
double derWF(Matrix r, double alpha, double beta, int dim, int NumOfPart);
void gradiendescent_brute_INT( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind, double alpha, double deltaAlpha, int Thermalization, double beta);
