#ifndef SYSTEM_H
#define SYSTEM_H
#include "matrix.h"
#include <Eigen/Dense>
using namespace Eigen;

class System
{
public:
    //System(){}
    System(int inp_N, int inp_M, double inp_sigma, double inp_omega);
    ~System(){}

    double localEnergy(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma, int samp, bool interact, int P);
    VectorXd grad_a(VectorXd x, VectorXd a, double sigma, int samp);
    VectorXd grad_b(VectorXd x, VectorXd b, MatrixXd W, double sigma,int samp);
    MatrixXd grad_W(VectorXd x, VectorXd b, MatrixXd W, double sigma,int samp);
    double grad_sigma(VectorXd x,VectorXd a, VectorXd b, MatrixXd W, double sigma);
    double QuantumForce(VectorXd x,VectorXd a, VectorXd b, MatrixXd W, double sigma, int i);
    double waveFunction(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma);
    VectorXd logSigm(VectorXd u);
    double GreenFunction(VectorXd xold,VectorXd xnew, VectorXd a,VectorXd b, MatrixXd W, double sigma, int par, int P,double difConst, double timeStep);

private:
  int s_N;
  int s_M;
  double s_sigma;
  double s_omega;

};


#endif // SYSTEM_H
