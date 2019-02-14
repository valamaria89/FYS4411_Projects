#include <iostream>
#include <random>
#include <fstream>
#include <cmath>
#include "time.h"
#include "exercises.h"
#include "matrix.h"
#include "methods.h"


using namespace std;

random_device rd;
mt19937_64 gen(rd());
uniform_real_distribution<double> RNG(0.0,1.0);
normal_distribution<double> ND(0.0,1.0);



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

       Methods ourMC(1,1,1000000,10,0.1, 2, 0.1, 0.1);
       ourMC.mc_sampling();
       //ourMC.writeToFile();
       /*for (int i = 0; i < 1000; i++){
            mat r= zeros<mat>(1, 1);
            for (int i = 0; i < 1; i++) {
                for (int j = 0; j < 1; j++){
                    r(i,j) = 0.1*(RNG(gen)-0.5);
            }
       }
        double eloc = ourMC.getFunc(r);
        if (abs(eloc - 0.5) < 1e-20){
           cout << "fail" << endl;
           break;
       }
       }
      // cout << "Alpha: " << alpha << "    " <<"Energy: " << eloc << endl;
*/


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
