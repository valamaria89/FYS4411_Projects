#include "rbm.h"
double factor = 0.05;
//constructor
RBM::RBM(int h, int p, int dim, double sigma){
    std::random_device rd;
    rbm_randomEngine = std::mt19937_64(rd());
    setup_RBM(h, p, dim, sigma);
}

//set up for RBM Machine
void RBM::setup_RBM(int h, int p, int dim, double sigma){
    rbm_M = dim*p;
    rbm_N = h;
    rbm_dim = dim;
    rbm_sigma = sigma;
    rbm_x.resize(rbm_M);
    rbm_h.resize(rbm_N);
    rbm_a.resize(rbm_M);
    rbm_b.resize(rbm_N);
    rbm_W.resize(rbm_M,rbm_N);
    setup_rweights();
    setup_rposition();
}

//sets random numbers for weights and biases
void RBM::setup_rweights(){
    rbm_x =VectorXd::Random(rbm_M);
    rbm_a = VectorXd::Random(rbm_M)*factor;
    rbm_b = VectorXd::Random(rbm_N)*factor;
    rbm_W = MatrixXd::Random(rbm_M, rbm_N)*factor;
}

//sets random position
void RBM::setup_rposition(){
    std::uniform_real_distribution<double> distribution_initX(-0.5,0.5);
        for(int i=0; i<rbm_M; i++){
            rbm_x(i)=distribution_initX(rbm_randomEngine);
        }

}
