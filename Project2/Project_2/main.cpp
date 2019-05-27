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
    int N = 2;

    int P =  2;
    int D = 2;

    double sigma = 1;
    double omega = 1;
    int NumofMC = 2000;
    int numVar = 50000;
    int samp=1;
    bool interact = true;
    cout << "Hei" << endl;
   GradientDescent(P, D, N, NumofMC, numVar, sigma, omega, 0.45, 0.01,samp,0.5,0.1, interact);

    RBM R(N, P, D, sigma);
    System Hamiltonian(N,P, sigma, omega );

/* cout << "Hei!" << endl;
 VectorXd a =VectorXd::Random(5);
 for(int i=0;i<5;i++){
     cout<<a(i)<<" "<<endl;

 }*/
    MPI_Finalize();
}


