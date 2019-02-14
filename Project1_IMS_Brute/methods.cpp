#include "methods.h"
//#include "matrix.h"
#include <random>
#include <cassert>
#include <fstream>
using namespace std;

// Constants for Hamiltonian and numerical calculation
#define h 0.0001
#define h2 100000000
#define mass 1
#define omega 1

//Monte Carlo cycles
void Methods::mc_sampling()
{
  cout<< dim <<" "<<numOfPart<<" "<<numMCCycles<< " "<<numVar<<" "<<step<<" "<<ISvsBrute<< " " << alpha << " " << deltaAlpha << endl;
//The variational parameter
        double wfnew = 0, wfold = 0, eng, eng2, delta_e;
        int accept;

//Stting of the old and the new positions to zero (2d matrices [number of particles, dimension])
        rOld = zeros<mat>(numOfPart, dim);
        rNew = zeros<mat>(numOfPart, dim);

        QFOld = zeros<mat>(numOfPart, dim);
        QFNew = zeros<mat>(numOfPart, dim);

//Initialize the random device
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

            if(ISvsBrute){
            for (int i = 0; i < numOfPart; i++) {
                for (int j = 0; j < dim; j++){
                    rOld(i,j) = step*(RNG(gen)-0.5);
                 }
               }
             }
            else{
                for (int i = 0; i < numOfPart; i++) {
                    for (int j = 0; j < dim; j++){
                        rOld(i,j) = sqrt(step)*(ND(gen));
                 }
               }
             }
 // initial trial WF and quantum force
       wfold = wavefunction(rOld);

       if(!ISvsBrute){
       QuantumForce(rOld, QFOld);}

 //The MC Cycle
            for (int cycl = 0; cycl <= numMCCycles; cycl++) {
 // We set the new position

            if(ISvsBrute){
            for (int i = 0; i < numOfPart; i++) {
                for (int j = 0; j < dim; j++){
                    rNew(i,j) = rOld(i,j) + step*(RNG(gen)-0.5);
                }
              }
            }
            else{
            for (int i = 0; i < numOfPart; i++) {
                for (int j = 0; j < dim; j++){
                    rNew(i,j) = rOld(i,j) + sqrt(step)*ND(gen)+0.5*QFOld(i,j)*step;
                }
              }
            }
 // The new trial WF
        wfnew = wavefunction(rNew);
        double GreenFunction = 0;
        QuantumForce(rNew, QFNew);

             for (int i = 0; i < numOfPart; i++) {
                 for (int j = 0; j < dim; j++){
                     GreenFunction += 0.5*(QFOld(i,j) + QFNew(i,j))*(0.25*step*(QFOld(i, j)-QFNew(i, j))+rOld(i, j) - rNew(i,j));
                 }
             }
             GreenFunction = exp(GreenFunction);



 // Metropolis test
         if(ISvsBrute){
             if(RNG(gen) <= wfnew*wfnew/(wfold*wfold) ) {
                  for (int i = 0; i < numOfPart; i++) {
                      for (int j = 0; j < dim; j++){
                            rOld(i,j) = rNew(i,j);
                            wfold=wfnew;
                            accept++;
                       }
                    }
                 }
              }
         else{
             if(RNG(gen) <= GreenFunction*wfnew*wfnew/(wfold*wfold) ) {
                  for (int i = 0; i < numOfPart; i++) {
                      for (int j = 0; j < dim; j++){
                            rOld(i,j) = rNew(i,j);
                            QFOld(i,j) = QFNew(i,j);
                            wfold=wfnew;
                            accept++;
                       }
                    }
                 }

         }
 // Computing of local energy
              if (ind){
                  delta_e = local_energy_analytic(rNew);
                 }
              else{
                  delta_e = local_energy_num(rNew);
                   }
// Updating of energies
                eng  += delta_e;
                eng2 += delta_e*delta_e;

         }
// The end of loop over MC trials

// Averaging over the number of MC trials
            Etot[var] = eng/numMCCycles;
            Etot2[var] = eng2/numMCCycles;
    cout<< alpha<< " " << Etot[var]<<" " <<Etot2[var]<<" " << accept<<endl;

//Increase the variational parameter
        alpha+=deltaAlpha;

 }
// The end of loop over variational  steps

 }
// The end of the MCsimpling function

//The wave function without interaction
double Methods::wavefunction(const mat &r){

      double wf=0, r_sqr=0;

      for (int i = 0; i < numOfPart; i++) {
        for (int j = 0; j < dim; j++) {
          r_sqr  += r(i,j)*r(i,j);
             }
           }

      wf = exp(-alpha*r_sqr) ;
      return wf;
    }

// The numerical calculation of local energy
double Methods::local_energy_num(const mat &r){
      double  wfminus, wfplus, wfold,ekin=0, epot=0, eloc=0, radsq = 0;
      mat rplus, rminus;

      rplus = zeros<mat>(numOfPart, dim);
      rminus = zeros<mat>(numOfPart, dim);
      rplus = r;
      rminus = r;

        for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 radsq  += r(i,j)*r(i, j);
             }
        }

    wfold=wavefunction(r);
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus(i,j) += h;
                 rminus(i,j) -= h;
                 wfminus = wavefunction(rminus);
                 wfplus  = wavefunction(rplus);
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


// The alytical calculation of local energy
double Methods::local_energy_analytic(const mat &r){
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

// The alytical calculation of quantum force
void Methods::QuantumForce(const mat &r, mat &QForce){
    mat rplus, rminus;
    rplus = zeros<mat>(numOfPart, dim);
    rminus = zeros<mat>(numOfPart, dim);
    rplus = r;
    rminus = r;

    double wfold ;
    wfold=wavefunction(r);
    double  wfminus = 0, wfplus = 0;

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus(i,j) += h;
                 rminus(i,j) -= h;
                 wfminus = wavefunction(rminus);
                 wfplus  = wavefunction(rplus);
                 QForce(i,j)= (wfplus-wfminus)/(wfold*h);
                 rplus (i,j) = r(i,j);
                 rminus(i,j)=  r(i,j);
          }
     }
}



/*void Methods::mc_sampling_IMS()
{
//The variational parameter
       // double alpha=0.50, deltaAlpha=0.01;
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;

//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rold(numOfPart, dim);
        Matrix rnew(numOfPart, dim);
        Matrix QForceNew(numOfPart, dim);
        Matrix QForceOld(numOfPart, dim);


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

//The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 //  initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, sqrt(step)*(ND(gen)));
                 }
           }

 // initial trial WF
    wfold = wavefunction(rold, al);
    QuantumForce(rold, QForceOld);

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles; cycl++) {

 // We set the new position

      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+sqrt(step)*ND(gen)+0.5*QForceOld.get_Elem(i,j)*step);
               }
         }

    wfnew = wavefunction(rnew, al);
    QuantumForce(rnew, QForceNew);
    double GreenFunction = 0;
     for (int i = 0; i < numOfPart; i++) {
         for (int j = 0; j < dim; j++){

             GreenFunction += 0.5*(QForceOld.get_Elem(i,j) + QForceNew.get_Elem(i,j))*(0.25*step*(QForceOld.get_Elem(i, j)-QForceNew.get_Elem(i, j))+rold.get_Elem(i, j) - rnew.get_Elem(i,j));
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

 // compute local energy
    if (ind==1){
        delta_e = local_energy_analytic(rold, al);
    }
    else{
        delta_e = local_energy_num(rold, al);
    }
// update energies
    eng  += delta_e;
    eng2 += delta_e*delta_e;

         }
// end of loop over MC trials
    Etot[var] = eng/numMCCycles;
   // cout<< "enrg= "<<eng<<endl;
    Etot2[var] = eng2/numMCCycles;
cout << " " << accept <<endl;

  //Increase the variational parameter
        //Methods al(alpha);
        //al.setalpha();
        //alpha += al.getalpha();
       // cout << "alpha" << alpha << endl;
 }    // end of loop over variational  steps

 }*/

void Methods:: writeToFile()
{
    ofstream myfile;
    myfile.open("E_average.txt");
    myfile  <<"Variational parameter E E2 variance " << endl ;
    for (int E = 0; E <= numVar; E++){
        myfile << 0.1+0.1*E <<" "<< Etot[E] << endl ;
        myfile << "                  variance "<< Etot2[E]- Etot[E]*Etot[E]<< endl ;
    }
    myfile.close();
}

