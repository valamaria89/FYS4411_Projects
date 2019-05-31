#include "system.h"
#include "matrix.h"
#include "gradientdescent.h"
#include "checkfunctions.h"
#include <iostream>
#include <Eigen/Dense>
#include "mpi.h"


void check_convergence(double eng, int M, double omega, bool inter, int var) {

    if (var>400){
        if(inter){
            if(M == 4){
                if(fabs(eng-3*omega)>0.15){
                    cout<<"Warning! Energy fluctuations are out of range!"<<endl;
            }
        }
        else{
                if(fabs(eng-0.5*M*omega)>0.1){
                     cout<<"Warning! Energy fluctuations are out of range!"<<endl;
                  }
             }
        }
    }
}

void check_units(){

    int N=2;
    int P=2;
    int dim=2;
    bool inter=false;
    int sampling=0;

    double omega=1;
    double sigma=1;

    RBM R(N, P, dim, sigma);
    System Hamiltonian(R.rbm_N,R.rbm_M, R.rbm_sigma, omega );

    double eng=0;
    eng = Hamiltonian.localEnergy(R.rbm_x, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma,sampling, inter, P);

    cout<<"Result of the energy test = "<<eng<<endl;

    cout<<"Gradient check"<<endl;

    cout<<"grad_a  grad_b  grad_W  grad_sigma"<<endl;

    for(int i=0;i<P*dim;i++){
        cout<<Hamiltonian.grad_a(R.rbm_x, R.rbm_a, R.rbm_sigma,sampling)<<" ";
        cout<<Hamiltonian.grad_b(R.rbm_x, R.rbm_b, R.rbm_W, R.rbm_sigma,sampling)<<" ";
        cout<<Hamiltonian.grad_W(R.rbm_x, R.rbm_b,R.rbm_W, R.rbm_sigma,sampling)<<" ";
        cout<<Hamiltonian.grad_sigma(R.rbm_x, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma)<<endl;
    }
}
