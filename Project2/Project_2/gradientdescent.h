#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H
#include "system.h"
#include "rbm.h"
#include <iostream>
#include <Eigen/Dense>
#include <random>

void GradientDescent( int P, int dim, int N, int MCcycles, int numOfvar, double sigma, double omega, double step, double learnRate );

#endif // GRADIENTDESCENT_H
