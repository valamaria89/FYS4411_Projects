#ifndef RBM_H
#define RBM_H
#include <iostream>
#include <Eigen/Dense>
#include <random>

using namespace Eigen;

class RBM
{
private:
    std::mt19937_64 rbm_randomEngine;
    void setup_RBM(int h, int x, int dim, double sigma);
    void setup_rweights();
    void setup_rposition();
public:
    RBM(int h, int x, int dim, double sigma);
    ~RBM(){}
    int rbm_M;
    int rbm_N;
    int rbm_dim;
    double rbm_sigma;
    VectorXd rbm_x;
    VectorXd rbm_h;
    VectorXd rbm_a;
    VectorXd rbm_b;
    MatrixXd rbm_W;

};

#endif // RBM_H



