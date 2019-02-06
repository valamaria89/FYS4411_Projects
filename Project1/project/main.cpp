#include <iostream>
#include <random>
#include <fstream>
#include "time.h"
#include "exercises.h"
#include "matrix.h"


using namespace std;



void initialise(int&, int&, int&,int&,double& ) ;


int main(){



    //The number of dimensions, particles, MC cycles and choise parameter
       int dim, numOfPart, numMCCycles,ind;

    //The number of variational parameters
       int numVar=20;
       double stepSize;

   //Initialise the number of dimensions, particles, MC cycles and energy arrays
       initialise(dim, numOfPart, numMCCycles, ind, stepSize);
       double *Etot, *Etot2;
       Etot = new double[numVar+1];
       Etot2 = new double[numVar+1];
clock_t time_start = clock();
   //MC cycles
double alpha=0.40, deltaAlpha=0.01;
       mc_sampling(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind, alpha, deltaAlpha);

       ofstream myfile;
       myfile.open("E_average.txt");
       myfile  <<"Variational parameter E E2 variance " << endl ;
       for (int E = 0; E <= numVar; E++){
           myfile << alpha+deltaAlpha*E <<" "<< Etot[E] << endl ;
           myfile << "                variance "<< Etot2[E]-Etot[E]*Etot[E]  << endl ;
       }
       myfile.close();
       delete [] Etot;
       delete [] Etot2;

     cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;


    /* double testen;
     testen=local_energy_analytic(rtest,0.5,dim,numOfPart);
for(int k=0;k<=100;k++){
    for (int i = 0; i < numOfPart; i++) {
         for (int j = 0; j < dim; j++){
               rtest.set_Elem( i, j, 0.0+0.01*k);
               cout<< " Energy= "<< local_energy_analytic(rtest,0.75,dim,numOfPart)<<endl;
             }
       }
     }*/

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
