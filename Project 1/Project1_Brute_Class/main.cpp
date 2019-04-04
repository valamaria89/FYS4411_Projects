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
//random_device cd;
//mt19937_64 ge(cd());
//uniform_real_distribution<double> GNR(0.0,1.0);

void initialise(int&, int&, int&,int&,double&) ;


int main(){

    //The number of dimensions, particles, MC cycles and choise parameter
           int dim, numOfPart, numMCCycles,ind;

    //The number of variational parameters and variational parameter
           int numVar=4;
           double stepSize = 1;
           double alpha=0.3, deltaAlpha=0.1;
           //double beta = 2.82843;
           double beta = 1.00;
           dim = 1;
           numOfPart = 20;
           numMCCycles = 1000000;
           ind =1;


    //Initialise the number of dimensions, particles, MC cycles and energy arrays
          // initialise(dim, numOfPart, numMCCycles, ind, stepSize);
           double *Etot, *Etot2;
           Etot = new double[numVar+1];
           Etot2 = new double[numVar+1];




    //Time start
           clock_t time_start = clock();

    //MC cycles


           ofstream file;
           file.open("time.txt");

           for (int num = 1; num < numOfPart+1; num ++){
           mc_sampling(numVar, num, dim, numMCCycles, ind, alpha, stepSize, deltaAlpha, Etot, Etot2, beta);
           file <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << ", ";
           //cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;
}
           file.close();
           delete [] Etot;
           delete [] Etot2;



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
