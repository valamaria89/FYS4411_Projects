#include "gradientdescent.h"


using namespace Eigen;
using namespace std;

// RNG


void GradientDescent( int P, int dim, int N, int MCcycles, int numOfvar, double sigma, double omega, double step, double learnRate ){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    uniform_int_distribution<> hrand(0, 1);

    RBM R(N, P, dim, sigma);
    VectorXd Xold = R.rbm_x;
    VectorXd Xnew = VectorXd::Zero(R.rbm_M);

    System Hamiltonian(R.rbm_N,R.rbm_M, R.rbm_sigma, omega );

    uniform_int_distribution<> mrand(0, R.rbm_M-1);

    for (int var = 0; var < numOfvar; var++){
        //energies:
        double eng = 0, eng2= 0, delta_e =0;

        //Derivatives of variables
        VectorXd Grad_a = VectorXd::Zero(R.rbm_M);
        VectorXd Grad_b = VectorXd::Zero(R.rbm_N);
        MatrixXd Grad_W = MatrixXd::Zero(R.rbm_M, R.rbm_N);
        double Grad_sigma = 0;

        //Average of energy and derivatives
        VectorXd Eng_Grad_a = VectorXd::Zero(R.rbm_M);
        VectorXd Eng_Grad_b = VectorXd::Zero(R.rbm_N);
        MatrixXd Eng_Grad_W = MatrixXd::Zero(R.rbm_M, R.rbm_N);
        double Eng_Grad_sigma = 0;


        //acceptance rate
        int accept = 0;

        for (int cycl = 0; cycl < MCcycles; cycl ++){
            Xnew = Xold;
            int M_random = mrand(gen);
            Xnew(M_random) = Xold(M_random) + (2*dis(gen)-1)*step;
            double wfold = Hamiltonian.waveFunction(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);
            double wfnew = Hamiltonian.waveFunction(Xnew, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);

            if (wfnew*wfnew/wfold*wfold >= dis(gen)){
                accept++;
                Xold = Xnew;

            }
            delta_e = Hamiltonian.localEnergy(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);

            //Summation of energies
            eng += delta_e;
            eng2 += pow(delta_e,2);

          //Derivatives of variables
            Grad_a += Hamiltonian.grad_a(Xold, R.rbm_a, R.rbm_sigma);
            Grad_b += Hamiltonian.grad_b(Xold, R.rbm_b, R.rbm_W, R.rbm_sigma);
            Grad_W += Hamiltonian.grad_W(Xold, R.rbm_b,R.rbm_W, R.rbm_sigma);
            Grad_sigma += Hamiltonian.grad_sigma(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);

             //Multiplied with energy
            Eng_Grad_a += Hamiltonian.grad_a(Xold, R.rbm_a, R.rbm_sigma)*delta_e;
            Eng_Grad_b += Hamiltonian.grad_b(Xold, R.rbm_b, R.rbm_W, R.rbm_sigma)*delta_e;
            Eng_Grad_W += Hamiltonian.grad_W(Xold, R.rbm_b,R.rbm_W, R.rbm_sigma)*delta_e;
            Eng_Grad_sigma += Hamiltonian.grad_sigma(Xold, R.rbm_a, R.rbm_b,R.rbm_W, R.rbm_sigma)*delta_e;

        }

        eng /= MCcycles;
        eng2 /= MCcycles;
        Grad_a /= MCcycles;
        Grad_b /= MCcycles;
        Grad_W /= MCcycles;
        Grad_sigma /= MCcycles;

        R.rbm_a -= 2*learnRate*(Eng_Grad_a/MCcycles - eng*Grad_a);
        R.rbm_b -= 2*learnRate*(Eng_Grad_b/MCcycles - eng*Grad_b);
        R.rbm_W -= 2*learnRate*(Eng_Grad_W/MCcycles - eng*Grad_W);
        R.rbm_sigma -= 2*learnRate*(Eng_Grad_sigma/MCcycles - eng*Grad_sigma);




    }
cout << "Hei!" << endl;
}
