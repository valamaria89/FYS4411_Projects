#include "system.h"
#include "matrix.h"
#include "gradientdescent.h"
#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;

int main()
{
    int N = 2;

    int P =  1;
    int D = 1;

    double sigma = 1;
    double omega = 1;
    int NumofMC = 10000;
    int numVar = 200;

   GradientDescent(P, D, N, NumofMC, numVar, sigma, omega, 0.45, 0.1);

    RBM R(N, P, D, sigma);
    System Hamiltonian(N,P, sigma, omega );

/* cout << "Hei!" << endl;
 VectorXd a =VectorXd::Random(5);
 for(int i=0;i<5;i++){
     cout<<a(i)<<" "<<endl;

 }*/
}


