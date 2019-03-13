#include <fstream>
#include <random>
#include <cassert>
#include <iostream>
#include "methods.h"
#include "system.h"
#include "position.h"
using namespace std;

// Uniform distribution for s in [0, 1]
//std::uniform_real_distribution<double> RNG(0.0,1.0);
// Constants
#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1
#define a 0.0043

random_device cd;
mt19937_64 ge(cd());
uniform_real_distribution<double> GNR(0.0,1.0);

double mc_sampling( System pOld, System pNew, double wfold, int numMCCycles, int ind, double alpha)
{

//The variational parameter
        double wfnew, eng, eng2, delta_e;
        int accept = 0;
        int counter=0;
        eng = 0;
        eng2 = 0;
        delta_e=0;


//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])


       // ofstream fileblock;
        //fileblock.open("E_block.txt");



 // initial trial WF


 //The MC Cycle
            for (int cycl = 1; cycl <= numMCCycles; cycl++) {

// We set the new position
            pNew.setRandomPosition();
            wfnew = pNew.wavefunction();

 // Metropolis test
            if(GNR(ge) < wfnew*wfnew/(wfold*wfold)){
                pOld.setNewToOld(pNew.get_Matrix());
                //p.setNewToOld();


                wfold=wfnew;
                accept++;
    }
// compute local energy
              if (ind==1){
                 delta_e = pOld.local_energy_analytic();
                  }
              else{
                 delta_e = pOld.local_energy_num();
    }
// update energies
            eng  += delta_e;
            eng2 += delta_e*delta_e;

            if(alpha==0.5){
            if(!((cycl) % 1000)){
            if(counter<512){
                //fileblock <<eng/cycl<< endl;
        }
                counter++;
      }
    }
}
// End of loop over MC trials

       // fileblock.close();
        return eng;

}
/*

void mc_sampling_IMS( double timestep, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization)
{
// The variational parameter
       // double alpha=0.50, deltaAlpha=0.01;
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;
        int counter=0;

// Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rold(numOfPart, dim);
        Matrix rnew(numOfPart, dim);
        Matrix QForceNew(numOfPart, dim);
        Matrix QForceOld(numOfPart, dim);
        ofstream fileblock;
        fileblock.open("E_block.txt");



      for (int i = 0; i < numOfPart; i++){
               for (int j = 0; j < dim; j++){
                   rold.set_Elem(i,j, 0);
                   rnew.set_Elem(i,j, 0);
                   QForceNew.set_Elem(i,j,0);
                   QForceOld.set_Elem(i,j, 0);
               }
           }

      random_device rd;
      mt19937_64 gen(rd());
      uniform_real_distribution<double> RNG(0.0,1.0);
      normal_distribution<double> ND(0.0,1.0);

// The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 // Initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, sqrt(timestep)*(ND(gen)));
                 }
           }

 // Initial trial WF
    wfold = wavefunction(rold, alpha, dim, numOfPart);
    QuantumForce(rold, QForceOld, dim, numOfPart, alpha);


 // The MC cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position

      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+sqrt(timestep)*ND(gen)+0.5*QForceOld.get_Elem(i,j)*timestep);
               }
         }

    wfnew = wavefunction(rnew, alpha, dim, numOfPart);
    QuantumForce(rnew, QForceNew, dim, numOfPart, alpha);
    double GreenFunction = 0;
     for (int i = 0; i < numOfPart; i++) {
         for (int j = 0; j < dim; j++){

             GreenFunction += 0.5*(QForceOld.get_Elem(i,j) + QForceNew.get_Elem(i,j))*(0.25*timestep*(QForceOld.get_Elem(i, j)-QForceNew.get_Elem(i, j))+rold.get_Elem(i, j) - rnew.get_Elem(i,j));
         }
     }
     GreenFunction = exp(GreenFunction);

 // Metropolis-Hastings test

    if(RNG(gen) <= GreenFunction*wfnew*wfnew/(wfold*wfold) ) {
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rold.set_Elem(i, j, rnew.get_Elem(i,j));
                    QForceOld.set_Elem(i, j, QForceNew.get_Elem(i,j));
                 }
           }
  wfold=wfnew;
  accept++;
    }
    if (cycl > thermalization){
 // compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rold, alpha, dim, numOfPart);
        }
        else{
            delta_e = local_energy_num(rold, alpha, wfold, dim, numOfPart);
        }
// update energies
        eng  += delta_e;
        eng2 += delta_e*delta_e;

        if(alpha==0.5){
        if(!((cycl) % 1000)){
            if(counter<512){
            fileblock <<eng/cycl<< endl;
            }
            counter++;
          }
        }
      }
    }
// end of loop over MC trials
    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
    cout << " " << accept <<endl;

  //Increase the variational parameter
        alpha += deltaAlpha;
 }    // end of loop over variational  steps
    fileblock.close();

 }

void writeToFile( string x, double *E1, double *E2, int numVar, double alpha, double deltaAlpha){
    string a;

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

*/
