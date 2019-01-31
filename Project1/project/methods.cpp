#include "methods.h"
#include "matrix.h"
#include <random>

using namespace std;

// Uniform distribution for s in [0, 1]
std::uniform_real_distribution<double> RNG(0.0,1.0);
#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1




void mc_sampling(std::mt19937_64& seed, double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2)
{
//The variational parameter
        double alpha=1.0, deltaAlpha=0.01;
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
        for (int i=0; i<=numVar; i++) {
            eng = eng2 = 0; accept =0; delta_e=0;

 //  initial trial position

        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, RNG(seed));
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

  delta_e = local_energy_num(rold, alpha, wfold, dim, numOfPart);
          // update energies
              eng += delta_e;
              eng2 += delta_e*delta_e;
         }

         // end of loop over MC trials
        Etot[i] = eng/numMCCycles;
        Etot2[i] = eng2/numMCCycles;

  //Increase the variational parameter
        alpha+= i*deltaAlpha;
 }    // end of loop over variational  steps


 }

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


double local_energy_num(Matrix r, double alpha, double wfold, int dim, int numOfPart){
      double eloc, wfminus, wfplus, ekin, epot, radsq = 0;
      Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 radsq  += r.get_Elem(i,j)*r.get_Elem(i, j);
             }
         }
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
                 ekin += (wfminus+wfplus-2*wfold)*h2;
             }
         }
      ekin = ekin/wfold;
      epot = mass*radsq*omega*omega/2;
      eloc = epot+ekin;
    return eloc;
}

double local_energy_analytic(Matrix r, double alpha, int dim, int numOfPart){
    double eloc = 0, ekin = 0, epot = 0, radsq = 0;

    for(int i = 0; i < numOfPart; i++){
        for(int j = 0; j < dim; j++){
            radsq = r.get_Elem(i,j)*r.get_Elem(i,j);
        }
    }
    ekin = -alpha*alpha*radsq*2/mass-3*alpha/2;
    epot = mass*radsq*omega*omega/2;
    eloc = ekin + epot;
    return eloc;
}


