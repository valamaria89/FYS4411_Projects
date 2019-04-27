#include <iostream>
#include "system.h"
#include "matrix.h"
#include <iostream>
#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;

using namespace std;

int main()
{
    int N = 2;
    int M = 2;
    double sigma = 3;
    double omega = 5;

System h(N,M,sigma,omega);

}


