#include <iostream>
#include <random>
#include <fstream>
#include <cmath>
#include "time.h"
#include "matrix.h"
#include "methods.h"
#include "position.h"
#include "system.h"


using namespace std;


//normal_distribution<double> ND(0.0,1.0);



void initialise(int&, int&, int&,int&,double&) ;


int main(){

    //The number of dimensions, particles, MC cycles and choise parameter
           int dim, numOfPart, numMCCycles,ind;

    //The number of variational parameters and variational parameter
           int numVar=4;
           double stepSize;
           double alpha=0.3, deltaAlpha=0.1;
           //double beta = 2.82843;
          double beta = 1.00;

    //Initialise the number of dimensions, particles, MC cycles and energy arrays
           initialise(dim, numOfPart, numMCCycles, ind, stepSize);
           double *Etot, *Etot2, wfold, energy;
           Etot = new double[numVar+1];
           Etot2 = new double[numVar+1];


    //Time start
           clock_t time_start = clock();

    //MC cycles
           for (int var=0; var<=numVar; var++) {
               System pOld(numOfPart, dim, alpha, beta, stepSize);
               System pNew(numOfPart,dim, alpha, beta, stepSize);
               pOld.setPosition_initial();
               pNew.setPosition_initial();
               pOld.setRandomPosition();
               wfold = pOld.wavefunction();
               energy = mc_sampling(pOld, pNew, wfold, numMCCycles, ind, alpha);
               Etot[var] = energy/numMCCycles;
               //Etot[var] = eng/numMCCycles;
               //Etot2[var] = eng2/numMCCycles;
               //cout << " " << accept <<endl;
               cout << "Energy: " << Etot[var] << "  alpha: " << alpha << endl;
               //Increase the variational parameter
               alpha += deltaAlpha;

           }
           //mc_sampling_IMS(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind,alpha, deltaAlpha, thermailization);
           // writeToFile("E_average_LA.txt", Etot, Etot2, numVar, alpha, deltaAlpha);
           //gradiendescent_brute(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind, alpha, deltaAlpha, thermailization);
           delete [] Etot;
           delete [] Etot2;

           cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;

    return 0;
}

// Beginning of function initialise
void initialise(int& dim, int& numOfPart, int& numMCCycles, int& ind, double& stepSize){

  cout << "Insert the number of dimensions = ";
  cin >> dim;
  cout << endl;
  cout << "Insert the number of particles = ";
  cin >> numOfPart;
  cout << endl;
  cout << "Insert the number of MC cycles = ";
  cin >> numMCCycles;
  cout <<  endl;
  cout << "Insert the stepsize = ";
  cin >> stepSize;
  cout <<  endl;
  cout << "1.Analytical - type 1"<<endl;
  cout << "1.Numerical - type 2"<<endl;
  cin >> ind;
}
// end of function initialise
