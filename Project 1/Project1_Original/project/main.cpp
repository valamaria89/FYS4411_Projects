#include <iostream>
#include <random>
#include <fstream>
#include "time.h"
#include "exercises.h"
#include "matrix.h"
#include "methods_interaction.h"

using namespace std;



void initialise(int&, int&, int&,int&,int&,double&, int& ) ;


int main(){

// The number of dimensions, particles, MC cycles and choise parameter
       int dim, numOfPart, numMCCycles, ind_1,ind_2, thermailization;

// The number of variational parameters and variational parameter
       int numVar=4;
       double stepSize = 1;
       double alpha=0.3, deltaAlpha=0.1;
       double beta = 2.82843;



// Initialise the number of dimensions, particles, MC cycles and energy arrays
       initialise(dim, numOfPart, numMCCycles, ind_1, ind_2, stepSize, thermailization);
       double *Etot, *Etot2;
       Etot = new double[numVar+1];
       Etot2 = new double[numVar+1];


// Time start

       clock_t time_start = clock();




        cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;



// Choose function
       switch (ind_2)
       {
       case 1:
         mc_sampling(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind_1, alpha, deltaAlpha, thermailization);
         writeToFile("E_average_LA.txt", Etot, Etot2, numVar, alpha, deltaAlpha);
         break;
       case 2:
         mc_sampling_IMS(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind_1, alpha, deltaAlpha, thermailization);
         writeToFile("E_average_LA.txt", Etot, Etot2, numVar, alpha, deltaAlpha);
         break;
       case 3:
         mc_sampling_INT(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind_1, alpha, deltaAlpha, thermailization, beta);
         writeToFile("E_average_LA.txt", Etot, Etot2, numVar, alpha, deltaAlpha);
         break;
       case 4:
         mc_sampling_IMS_INT(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind_1, alpha, deltaAlpha, thermailization, beta);
         writeToFile("E_average_LA.txt", Etot, Etot2, numVar, alpha, deltaAlpha);
         break;
       case 5:
         gradiendescent_brute(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind_1, alpha, deltaAlpha, thermailization);
         break;
       case 6:
         gradiendescent_brute_INT(stepSize, dim, numOfPart, numMCCycles, numVar, Etot, Etot2, ind_1, alpha, deltaAlpha, thermailization, beta);
         break;
       default:
       cout << "Coose any option next time!";
       }
       cout << "Main running for " << " " <<  double((clock()-time_start)/double(CLOCKS_PER_SEC)) << " seconds" << endl;

       alpha=alpha-numVar*deltaAlpha;
       if((ind_2==1)||(ind_2==2)){
       for(int i=0; i<numVar; i++){
           double e_GP=0;
           e_GP=numOfPart*(1/(4*alpha)+alpha)*(1+beta/2);
           if(fabs(e_GP-Etot[i])>5){
               cout<<"Results do not correspond to the GP equation!Please, check the parameters!"<<endl;
           }
         }
       }

       delete [] Etot;
       delete [] Etot2;
    return 0;
}

// Beginning of function initialise
void initialise(int& dim, int& numOfPart, int& numMCCycles, int& ind_1, int& ind_2, double& stepSize, int& thermalization){

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
  cin >> ind_1;
  cout << "1.MC with brute force, non-interacting - 1"<<endl;
  cout << "2.MC with importance sampling, non-interacting - 2"<<endl;
  cout << "3.MC with brute force, non-interacting - 3"<<endl;
  cout << "4.MC with importance sampling, interacting - 4"<<endl;
  cout << "5.MC with importance sampling, interacting - 5"<<endl;
  cout << "6.GD, non-interacting - 6"<<endl;
  cout << "7.GD, interacting - 7"<<endl;
  cin >> ind_2;
}
// end of function initialise
