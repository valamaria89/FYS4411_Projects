#include <iostream>
#include <random>
#include <fstream>
#include "time.h"
#include "exercises.h"
#include "matrix.h"
#include "methods.h"

using namespace std;



void initialise(int&, int&, int&,int&,double& ) ;


int main(){

//The number of dimensions, particles, MC cycles and choise parameter
       //int dim, numOfPart, numMCCycles,ind;

//The number of variational parameters and variational parameter
       //int numVar=20;
       //double stepSize;
       //double alpha=0.40, deltaAlpha=0.01;

//Initialise the number of dimensions, particles, MC cycles and energy arrays
       //initialise(dim, numOfPart, numMCCycles, ind, stepSize);

//Time start
       clock_t time_start = clock();

       Methods ourMC(1,1,1000000,10,1, 1, 0.1, 0.1);
       ourMC.mc_sampling();
       ourMC.writeToFile();



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
