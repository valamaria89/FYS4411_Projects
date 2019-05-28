#include "system.h"
#include "matrix.h"
#include "gradientdescent.h"
#include <iostream>
#include <Eigen/Dense>
#include "mpi.h"


using namespace Eigen;

using namespace std;

int main()
{ MPI_Init(nullptr, nullptr);
    int N = 1;

    int P =  1;
    int D = 1;

    double sigma = 1;
    double omega = 1;
    int NumofMC = 1000000;
    int numVar = 200;
    int samp=0;
    bool interact = false;
    double stepsize = 0.1;
    double learnRate = 0.1;
    double diffConst = 0.5;
    double timestep = 0.01;

   GradientDescent(P, D, N, NumofMC, numVar, sigma, omega, stepsize, learnRate ,samp,diffConst,timestep, interact);

    RBM R(N, P, D, sigma);
    System Hamiltonian(N,P, sigma, omega );


    MPI_Finalize();
}


