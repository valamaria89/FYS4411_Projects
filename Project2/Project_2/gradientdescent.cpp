#include "gradientdescent.h"
#include "mpi.h"


using namespace Eigen;
using namespace std;

// RNG


void GradientDescent( int P, int dim, int N, int MCcycles, int numOfvar, double sigma, double omega, double step, double learnRate, int samp,double difConst,double timeStep, bool interact){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    uniform_int_distribution<> hrand(0, 1);
    normal_distribution<double> gauss(0,1);
    cout<<"Random = "<<dis(gen)<<endl;

    //RBM R(N, P, dim, sigma);
    //VectorXd Xold = R.rbm_x;
    //VectorXd Xnew = VectorXd::Zero(R.rbm_M);
    //cout<<"X initial="<<Xold(0)<<endl;
    //System Hamiltonian(R.rbm_N,R.rbm_M, R.rbm_sigma, omega );
    //cout<<Hamiltonian.localEnergy(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma)<<endl;
    //uniform_int_distribution<> mrand(0, R.rbm_M-1);

    int counter = 0;
    RBM R(N, P, dim, sigma);
    System Hamiltonian(R.rbm_N,R.rbm_M, R.rbm_sigma, omega );
    //VectorXd Xold = R.rbm_x;
    VectorXd Xold = VectorXd::Zero(R.rbm_M);


    for (int var = 0; var < numOfvar; var++){
        //energies:
        double eng = 0, eng2= 0, delta_e =0;
//RBM R(N, P, dim, sigma);
//System Hamiltonian(R.rbm_N,R.rbm_M, R.rbm_sigma, omega );
//VectorXd Xold = R.rbm_x;
VectorXd Xnew = VectorXd::Zero(R.rbm_M);
Xold =VectorXd::Random(R.rbm_M);
uniform_int_distribution<> mrand(0, R.rbm_M-1);
uniform_int_distribution<> nrand(0, R.rbm_N-1);


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
        double condition=0;
        int n_procs ;
        MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

   //cout<<"Xold= "<<Xold(0)<<endl;
//#pragma ompi
        for (int cycl = 0; cycl < MCcycles / n_procs; cycl ++){
            Xnew = Xold;
            //cout<<"Xnew= "<<Xnew(0)<<endl;
            int M_random = mrand(gen);
            int N_random = nrand(gen);
            //Xnew(M_random) = Xold(M_random) + (2*dis(gen)-1)*step;

            //Metropolis Brute
            if(samp==0){
            Xnew(M_random) = Xold(M_random) + (dis(gen)-0.5)*step;
            double wfold = Hamiltonian.waveFunction(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);
            double wfnew = Hamiltonian.waveFunction(Xnew, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);
            condition = wfnew*wfnew/wfold*wfold;
            }

            // Metropolis Hastings
            double GreenFunc=0;
            if(samp==1){
            Xnew(M_random) = Xold(M_random) + difConst*Hamiltonian.QuantumForce(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma,M_random)*timeStep + gauss(gen)*sqrt(timeStep);
            double wfold = Hamiltonian.waveFunction(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);
            double wfnew = Hamiltonian.waveFunction(Xnew, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);
            GreenFunc=Hamiltonian.GreenFunction(Xold,Xnew,R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma,0,P,difConst, timeStep);
            condition = GreenFunc*wfnew*wfnew/wfold*wfold;
            }

            if((samp==0)||(samp==1)){
            if (condition >= dis(gen)){
                accept++;
                Xold = Xnew;
               }
            }
            // Gibbs sampling
            if(samp==2){
                Xold(M_random) = R.pick_x (M_random);
                R.rbm_h(N_random) = R.pick_h(N_random);
            }
            delta_e = Hamiltonian.localEnergy(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma,samp, interact, P);

            //Summation of energies
            eng += delta_e;
            eng2 += pow(delta_e,2);

          //Derivatives of variables
            Grad_a += Hamiltonian.grad_a(Xold, R.rbm_a, R.rbm_sigma,samp);
            Grad_b += Hamiltonian.grad_b(Xold, R.rbm_b, R.rbm_W, R.rbm_sigma,samp);
            Grad_W += Hamiltonian.grad_W(Xold, R.rbm_b,R.rbm_W, R.rbm_sigma,samp);
            Grad_sigma += Hamiltonian.grad_sigma(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma);

          //Multiplied with energy
            Eng_Grad_a += Hamiltonian.grad_a(Xold, R.rbm_a, R.rbm_sigma,samp)*delta_e;
            Eng_Grad_b += Hamiltonian.grad_b(Xold, R.rbm_b, R.rbm_W, R.rbm_sigma,samp)*delta_e;
            Eng_Grad_W += Hamiltonian.grad_W(Xold, R.rbm_b,R.rbm_W, R.rbm_sigma,samp)*delta_e;
            Eng_Grad_sigma += Hamiltonian.grad_sigma(Xold, R.rbm_a, R.rbm_b,R.rbm_W, R.rbm_sigma)*delta_e;

        }

        double global_eng, global_eng2;
        VectorXd global_Grad_a = VectorXd::Zero(R.rbm_M);
        VectorXd global_Grad_b = VectorXd::Zero(R.rbm_N);
        MatrixXd global_Grad_W = MatrixXd::Zero(R.rbm_M, R.rbm_N);

        VectorXd global_Eng_Grad_a = VectorXd::Zero(R.rbm_M);
        VectorXd global_Eng_Grad_b = VectorXd::Zero(R.rbm_N);
        MatrixXd global_Eng_Grad_W = MatrixXd::Zero(R.rbm_M, R.rbm_N);
        int global_accept;
        MPI_Allreduce(&accept, &global_accept, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&eng2, &global_eng2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&eng, &global_eng, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Grad_a.data(), global_Grad_a.data(), Grad_a.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Grad_b.data(), global_Grad_b.data(), Grad_b.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Grad_W.data(), global_Grad_W.data(), Grad_W.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Eng_Grad_a.data(), global_Eng_Grad_a.data(), Eng_Grad_a.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Eng_Grad_b.data(), global_Eng_Grad_b.data(), Eng_Grad_b.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(Eng_Grad_W.data(), global_Eng_Grad_W.data(), Eng_Grad_W.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   //cout<<Hamiltonian.localEnergy(Xold, R.rbm_a, R.rbm_b, R.rbm_W, R.rbm_sigma)<<endl;
       // for(int i=0;i<R.rbm_M;i++){
            //cout<<R.rbm_a(i)<<" ";
       // }
        //cout <<endl;
        //cout<<accept;
        //cout<<"X="<<Xnew(0)<<endl;

        global_eng /= MCcycles;
        global_eng2 /= MCcycles;
        global_Grad_a /= MCcycles;
        global_Grad_b /= MCcycles;
        global_Grad_W /= MCcycles;
        Grad_sigma /= MCcycles;

        R.rbm_a -= 2*learnRate*(global_Eng_Grad_a/MCcycles - global_eng*global_Grad_a);
        R.rbm_b -= 2*learnRate*(global_Eng_Grad_b/MCcycles - global_eng*global_Grad_b);
        R.rbm_W -= 2*learnRate*(global_Eng_Grad_W/MCcycles - global_eng*global_Grad_W);
        R.rbm_sigma -= 0; //2*learnRate*(Eng_Grad_sigma/MCcycles - eng*Grad_sigma);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank ==0){
        cout <<"Energy = " <<global_eng<< " Variance = "<<global_eng2-global_eng*global_eng<< " Acceptance = " << global_accept << endl;
        }
/*if((global_eng<0.6)&&(eng>0.45)){
counter++;
}
    }
    cout<< "Counter="<<counter<<endl;
*/
}
}
