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
    int NumofMC = 262144;
    //262144
    int numVar = 500;
    int samp=0;
    bool interact = false;
    double stepsize = 1;
    double learnRate = 0.2;
    double diffConst = 0.5;
    double timestep = 0.01;

    clock_t time_start = clock();
    //for (int i = 1; i < N; i++ ){
    GradientDescent(P, D, N, NumofMC, numVar, sigma, omega, stepsize, learnRate ,samp,diffConst,timestep, interact);



    // Time start


     cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;


    MPI_Finalize();
}


