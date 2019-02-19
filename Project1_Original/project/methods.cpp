#include "methods.h"
#include "matrix.h"
#include <fstream>
#include <random>
#include <cassert>
using namespace std;

// Uniform distribution for s in [0, 1]
//std::uniform_real_distribution<double> RNG(0.0,1.0);
// Constants
#define h 0.0001
#define h2 100000000
#define mass 1
#define omega 1

void mc_sampling( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization)
{
//The variational parameter
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;
        mat rOld;
        mat rNew;

//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        rOld = zeros<mat>(numOfPart, dim);
        rNew = zeros<mat>(numOfPart, dim);

        ofstream fileblock;
        fileblock.open("E_block.txt");

        random_device rd;
        mt19937_64 gen(rd());
        uniform_real_distribution<double> RNG(0.0,1.0);

//The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 //  initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rOld(i,j) = step*(RNG(gen)-0.5);
                 }
           }

 // initial trial WF
    wfold = wavefunction(rOld, alpha, dim, numOfPart);
    //cout<< wfold<<endl;

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position
      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rNew(i,j) = rOld(i,j) + step*(RNG(gen)-0.5);
               }
         }
    wfnew = wavefunction(rNew, alpha, dim, numOfPart);

 // Metropolis test
    if(RNG(gen) < wfnew*wfnew/(wfold*wfold)){
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rOld(i,j) = rNew(i,j);
                 }
           }
  wfold=wfnew;
  accept++;
    }
    if (cycl > thermalization){
 // compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rOld, alpha, dim, numOfPart);
    }
        else{
            delta_e = local_energy_num(rOld, alpha, wfold, dim, numOfPart);
    }
// update energies
    eng  += delta_e;
    eng2 += delta_e*delta_e;

    if(alpha==0.5){
    if(!((cycl) % 1000)){

        fileblock  <<eng/cycl<< endl;
        }
       ;};
    }

        }
// end of loop over MC trials
    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
    cout << " " << accept <<endl;

//Increase the variational parameter
        alpha += deltaAlpha;
 }

// end of loop over variational  steps
        fileblock.close();

}
//The wave function without interaction
double wavefunction(mat r, double alpha, int dim, int NumOfPart){

      double wf=0, r_sqr=0;

      for (int i = 0; i < NumOfPart; i++) {
        for (int j = 0; j < dim; j++) {
          r_sqr  += r(i,j)*r(i,j);
        }
      }
      wf = exp(-alpha*r_sqr) ;
      return wf;
    }

//The numerical calculation of local energy
double local_energy_num(mat r, double alpha, double wfold, int dim, int numOfPart){
      double  wfminus, wfplus, ekin=0, epot=0, eloc=0, radsq = 0;
      mat rplus(numOfPart, dim), rminus(numOfPart, dim);


      rplus = zeros<mat>(numOfPart, dim);
      rminus = zeros<mat>(numOfPart, dim);
      rplus = r;
      rminus = r;

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 radsq  += r(i,j)*r(i, j);
             }
         }
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus(i,j) += h;
                 rminus(i,j) -= h;
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
                 ekin += (wfminus+wfplus-2*wfold);
                 rplus (i,j) = r(i,j);
                 rminus(i,j)=  r(i,j);
             }
         }
      ekin = -0.5*ekin*h2/wfold;

      epot = radsq/2;
      eloc = epot+ekin;

      return eloc;
  }


//The alytical calculation of local energy
double local_energy_analytic(mat r, double alpha, int dim, int numOfPart){
    double  ekin = 0,epot=0, eloc=0, radsq = 0;

    for(int i = 0; i < numOfPart; i++){
        for(int j = 0; j < dim; j++){
            radsq += r(i,j)*r(i,j);
        }
    }
    ekin = -2*alpha*alpha*radsq+alpha*dim*numOfPart;

    epot = radsq/2;
    eloc = epot+ekin;
    return eloc;
}

void QuantumForce(mat r, mat QForce, int dim, int numOfPart, double alpha){
    mat rplus, rminus;
        rplus = zeros<mat>(numOfPart, dim);
        rminus = zeros<mat>(numOfPart, dim);
        rplus = r;
        rminus = r;

    double wfold ;
    wfold=wavefunction(r, alpha, dim, numOfPart);
    double  wfminus = 0, wfplus = 0;

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus(i,j) += h;
                 rminus(i,j) -= h;
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
                 QForce(i,j)= (wfplus-wfminus)/(wfold*h);
                 rplus (i,j) = r(i,j);
                 rminus(i,j)=  r(i,j);
          }
     }
}



void mc_sampling_IMS( double timestep, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization)
{
//The variational parameter
       // double alpha=0.50, deltaAlpha=0.01;
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;

//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        mat rOld;
        mat rNew;
        mat QFOld;
        mat QFNew;


        rOld = zeros<mat>(numOfPart, dim);
        rNew = zeros<mat>(numOfPart, dim);

        QFOld = zeros<mat>(numOfPart, dim);
        QFNew = zeros<mat>(numOfPart, dim);

      random_device rd;
      mt19937_64 gen(rd());
      uniform_real_distribution<double> RNG(0.0,1.0);
      normal_distribution<double> ND(0.0,1.0);

//The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 //  initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rOld(i,j) = sqrt(timestep)*(ND(gen));
                 }
           }

 // initial trial WF
    wfold = wavefunction(rOld, alpha, dim, numOfPart);
    QuantumForce(rOld, QFOld, dim, numOfPart, alpha);


 //The MC сycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position

      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rNew(i,j) = rOld(i,j) + sqrt(timestep)*ND(gen)+0.5*QFOld(i,j)*timestep;
               }
         }

    wfnew = wavefunction(rNew, alpha, dim, numOfPart);
    QuantumForce(rNew, QFNew, dim, numOfPart, alpha);
    double GreenFunction = 0;
     for (int i = 0; i < numOfPart; i++) {
         for (int j = 0; j < dim; j++){

             GreenFunction += 0.5*(QFOld(i,j) + QFNew(i,j))*(0.25*timestep*(QFOld(i, j)-QFNew(i, j))+rOld(i, j) - rNew(i,j));
         }
     }
     GreenFunction = exp(GreenFunction);

 // Metropolis-Hastings test

    if(RNG(gen) <= GreenFunction*wfnew*wfnew/(wfold*wfold) ) {
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                 rOld(i,j) = rNew(i,j);
                 QFOld(i,j) = QFNew(i,j);
                 }
           }
  wfold=wfnew;
  accept++;
    }
    if (cycl > thermalization){
 // compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rOld, alpha, dim, numOfPart);
        }
        else{
            delta_e = local_energy_num(rOld, alpha, wfold, dim, numOfPart);
        }
// update energies
        eng  += delta_e;
        eng2 += delta_e*delta_e;

        }
        }
// end of loop over MC trials
    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
    cout << " " << accept <<endl;

  //Increase the variational parameter
        alpha += deltaAlpha;
 }    // end of loop over variational  steps

 }

void writeToFile( string x, double *E1, double *E2, int numVar, double alpha, double deltaAlpha)
{

    string a;
    //cout << "What are we looking at?" << endl;
    ofstream myfile;
    myfile.open(x);
    //cin >> a;
    //myfile  << a << endl;

    for (int E = 0; E <= numVar; E++){
        myfile << alpha+deltaAlpha*E <<" "<< E1[E] << endl ;
        myfile << "                variance "<<  E2[E]-E1[E]*E1[E] << endl ;
    }
    myfile.close();
}




