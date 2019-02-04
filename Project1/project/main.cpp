#include <iostream>
#include <random>
#include <fstream>
#include "time.h"
#include "exercises.h"
#include "matrix.h"


using namespace std;



void initialise(int&, int&, int&,int& ) ;


int main(){

    std::mt19937_64 seed(1234);

    //The number of dimensions, particles, MC cycles and choise parameter
       int dim, numOfPart, numMCCycles,ind;

    //The number of variational parameters
       int numVar=50;
       double stepSize=0.00001;

   //Initialise the number of dimensions, particles, MC cycles and energy arrays
       initialise(dim, numOfPart, numMCCycles, ind);
       double *Etot, *Etot2;
       Etot = new double[numVar+1];
       Etot2 = new double[numVar+1];
clock_t time_start = clock();
   //MC cycles
       mc_sampling(seed, stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind);

       ofstream myfile;
       myfile.open("E_average.txt");
       myfile  <<"Variational parameter E E2 variance " << endl ;
       for (int E = 0; E <= numVar; E++){
           myfile << 0.001+0.001*E <<" "<< Etot[E] << endl ;
           //myfile << " variance"<< Etot2[E]  << endl ;
       }
       myfile.close();
       delete [] Etot;
       delete [] Etot2;


       Matrix rnew(numOfPart, dim);

     for (int i = 0; i < numOfPart; i++){
              for (int j = 0; j < dim; j++){
                  rnew.set_Elem(i ,j, 1);
              }
          }
     double wf=wavefunction(rnew,0.1,dim,numOfPart);
     double ander= deriv_analytic(rnew, 0.1, dim, numOfPart);
      double numder= deriv_num(rnew, 0.1, wf, dim, numOfPart);
     cout << ander <<" and "<<numder<< endl;
     cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;
    return 0;
}

// Beginning of function initialise
void initialise(int& dim, int& numOfPart, int& numMCCycles, int& ind){

  cout << "Insert the number of dimensions = ";
  cin >> dim;
  cout << endl;
  cout << "Insert the number of particles = ";
  cin >> numOfPart;
  cout << endl;
  cout << "Insert the number of MC cycles = ";
  cin >> numMCCycles;
  cout <<  endl;
  cout << "1.Analytical - type 1"<<endl;
  cout << "1.Numerical - type 2"<<endl;
  cin >> ind;
}
// end of function initialise
