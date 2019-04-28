#include "system.h"
#include "matrix.h"
#include "gradientdescent.h"
#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;

using namespace std;

int main()
{
    int N = 2;

    int P =  1;
    int D = 1;

    double sigma = 1;
    double omega = 1;
    int NumofMC = 1000;
    int numVar = 50;

    GradientDescent(P, D, N, NumofMC, numVar, sigma, omega, 1, 0.5);



}


