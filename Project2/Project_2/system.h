#ifndef SYSTEM_H
#define SYSTEM_H
#include "matrix.h"

class System
{
public:
    //System(){}
    System(int N, int M, double sigma, double omega);
    ~System(){}
    double waveFunction(Matrix x, Matrix a, Matrix W, double sigma);
    double localEnergy_IMS(Matrix x, Matrix a, Matrix W, double sigma);
    double localEnergy_Gibbs(Matrix x, Matrix a, Matrix W, double sigma);

    void QuantumForce_IMS(Matrix x, Matrix a, Matrix W, Matrix qForce);
    void QuantumForce_Gibbs(Matrix x, Matrix a, Matrix W, Matrix qForce);

    void grad_a(Matrix x, Matrix a, Matrix W, Matrix psiDa);
    void grad_b(Matrix x, Matrix a, Matrix W, Matrix psiDb);
    void grad_W(Matrix x, Matrix a, Matrix W, Matrix psiDW);
    void grad_sigma(Matrix x, Matrix a, Matrix W, Matrix psiDW);

    double u(Matrix x, Matrix a, Matrix b, Matrix W, double sigma, int k);

private:
  int N;
  int M;
  double sigma;
  double omega;

};


#endif // SYSTEM_H
