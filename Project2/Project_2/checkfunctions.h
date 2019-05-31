#ifndef CHECKFUNCTIONS_H
#define CHECKFUNCTIONS_H
#include "system.h"
#include "rbm.h"
#include <iostream>
#include <Eigen/Dense>
#include <random>

void check_convergence(double eng, int M, double omega, bool inter, int var);

void check_energy();

void check_units();

#endif // CHECKFUNCTIONS_H
