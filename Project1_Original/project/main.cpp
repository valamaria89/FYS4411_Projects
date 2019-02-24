#include <iostream>
#include <random>
#include <fstream>
#include "time.h"
#include "exercises.h"
#include "matrix.h"


using namespace std;



void initialise(int&, int&, int&,int&,double&, int& ) ;


int main(){

//The number of dimensions, particles, MC cycles and choise parameter
       int dim, numOfPart, numMCCycles,ind, thermailization;

//The number of variational parameters and variational parameter
       int numVar=10;
       double stepSize;
       double alpha=0.1, deltaAlpha=0.1;

//Initialise the number of dimensions, particles, MC cycles and energy arrays
       initialise(dim, numOfPart, numMCCycles, ind, stepSize, thermailization);
       double *Etot, *Etot2;
       Etot = new double[numVar+1];
       Etot2 = new double[numVar+1];


//Time start
       clock_t time_start = clock();

//MC cycles
       mc_sampling(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind, alpha, deltaAlpha, thermailization);
       //mc_sampling_IMS(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind,alpha, deltaAlpha, thermailization);
       writeToFile("E_average_LA.txt", Etot, Etot2, numVar, alpha, deltaAlpha);
       //gradiendescent_brute(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind, alpha, deltaAlpha, thermailization);
       delete [] Etot;
       delete [] Etot2;

       cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;


    return 0;
}

// Beginning of function initialise
void initialise(int& dim, int& numOfPart, int& numMCCycles, int& ind, double& stepSize, int& thermalization){

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
  cout << "Insert Thermalization =  ";
  cin >> thermalization;
  cout << endl;
  cout << "1.Analytical - type 1"<<endl;
  cout << "1.Numerical - type 2"<<endl;
  cin >> ind;
}
// end of function initialise
