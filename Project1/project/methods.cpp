#include "methods.h"
#include "matrix.h"
#include <random>

using namespace std;

// Uniform distribution for s in [0, 1]
std::uniform_real_distribution<double> RNG(0.0,1.0);
// Constants
#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1

void mc_sampling(std::mt19937_64& seed, double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind)
{
//The variational parameter
        double alpha=0.001, deltaAlpha=0.001;
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;

//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rold(numOfPart, dim);
        Matrix rnew(numOfPart, dim);

      for (int i = 0; i < numOfPart; i++){
               for (int j = 0; j < dim; j++){
                   rold.set_Elem(i, j, 0);
                   rnew.set_Elem(i ,j, 0);
               }
           }

//The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 //  initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, step*RNG(seed));
                 }
           }

  //  initial trial WF
    wfold = wavefunction(rold, alpha, dim, numOfPart);

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles; cycl++) {

 // We set the new position
      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
                  rnew.set_Elem(i, j, rold.get_Elem(i,j)+step*RNG(seed));
               }
         }
    wfnew = wavefunction(rnew, alpha, dim, numOfPart);
 // Metropolis test
//cout << wfnew<<" " << wfold <<endl;
    if(RNG(seed) <= wfnew*wfnew/wfold/wfold ) {
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rold.set_Elem(i, j, rnew.get_Elem(i,j));
                 }
           }
  wfold = wfnew;
  accept++;
    }
 // compute local energy
    if (ind==1)
  delta_e = local_energy_analytic(rold, alpha, dim, numOfPart);
    else
  delta_e = local_energy_num(rold, alpha, wfold, dim, numOfPart);
// update energies
    eng = eng + delta_e;
    eng2 = eng2+(delta_e)*(delta_e);
         }
// end of loop over MC trials
    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
//cout << " " << accept <<endl;
  //Increase the variational parameter
        alpha+= var*deltaAlpha;
 }    // end of loop over variational  steps


 }

//The wave function without interaction
double wavefunction(Matrix r, double alpha, int dim, int NumOfPart){

      double wf=0, r_sqr=0;

      for (int i = 0; i < NumOfPart; i++) {
        for (int j = 0; j < dim; j++) {
          r_sqr  += r.get_Elem(i,j)*r.get_Elem(i,j);
        }
      }
      wf = exp(-alpha*r_sqr) ;
      return wf;
    }

//The numerical calculation of local energy
double local_energy_num(Matrix r, double alpha, double wfold, int dim, int numOfPart){
      double  wfminus, wfplus, ekin=0,  radsq = 0;
      Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
                 radsq  += r.get_Elem(i,j)*r.get_Elem(i, j);
             }
         }
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
                 ekin += (wfminus+wfplus-2*wfold);
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
             }
         }
      ekin = ekin*h2/wfold;

    return ekin;
}

//The nanalytical calculation of local energy
double local_energy_analytic(Matrix r, double alpha, int dim, int numOfPart){
    double  ekin = 0, radsq = 0;

    for(int i = 0; i < numOfPart; i++){
        for(int j = 0; j < dim; j++){
            radsq = r.get_Elem(i,j)*r.get_Elem(i,j);
        }
    }
    ekin = 4*alpha*alpha*radsq-2*alpha*dim*numOfPart;

    return ekin;
}

double deriv_num(Matrix r, double alpha, double wfold, int dim, int numOfPart){
      double  wfminus, wfplus, ekin=0,  radsq = 0;
      Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
                 radsq  += r.get_Elem(i,j)*r.get_Elem(i, j);
             }
         }
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
                 ekin += (wfminus+wfplus-2*wfold);
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
             }
         }
      ekin = ekin*h2/wfold;

    return ekin;
}

//The nanalytical calculation of local energy
double deriv_analytic(Matrix r, double alpha, int dim, int numOfPart){
    double eloc = 0, ekin = 0, epot = 0, radsq = 0;

    for(int i = 0; i < numOfPart; i++){
        for(int j = 0; j < dim; j++){
            radsq = r.get_Elem(i,j)*r.get_Elem(i,j);
         };
    };
    ekin = (4*alpha*alpha*radsq-2*alpha*dim*numOfPart);
    return ekin;
}
