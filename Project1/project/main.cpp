#include <iostream>
#include <random>
#include "time.h"
#include "exercises.h"
#include "matrix.h"


using namespace std;



void initialise(int&, int&, int&) ;

int main(){
   /*Matrix a(5,2);
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 2; j++){
            a.set_Elem(i, j, i);
            cout << a.get_Elem(i, j) << ", ";
        }

         cout << endl;
    }*/


    std::mt19937_64 seed(1234);

    //The number of dimensions, particles, MC cycles
       int dim, numOfPart, numMCCycles;

    //The number of variational parameters
       int numVar=5;

   //Initialise the number of dimensions, particles, MC cycles
       initialise(dim, numOfPart, numMCCycles);
          cout << dim << endl;
          cout << numOfPart<< endl;
          cout << numMCCycles << endl;


            cout << endl;
      Matrix r(numOfPart, dim);
       for (int i = 0; i < numOfPart; i++){
           for (int j = 0; j < dim ; j++){
               r.set_Elem(i, j, 2);
           }}
       double alpha = 1;
       double wf1 = wavefunction(r, alpha, dim, numOfPart);
       double E_l = local_energy_num(r, alpha, wf1, dim, numOfPart);
       double E_la = local_energy_analytic(r, alpha, dim, numOfPart);
       double *Etot, *Etot2;
       Etot = new double[numVar+1];
       Etot2 = new double[numVar+1];
       mc_sampling(seed, 0.001, dim, numOfPart, numMCCycles, numVar, Etot, Etot2);
       cout << E_la << endl;

       for (int i = 0; i <= numVar; i++){
           cout << Etot[i] << endl;

       }



/////////////////////////////////////////



    return 0;
}

// Beginning of function initialise
void initialise(int& dim, int& numOfPart, int& numMCCycles)
{
  cout << "Insert the number of dimensions ";
  cin >> dim;
   cout << endl;
  cout << "Insert the number of particles = ";
  cin >> numOfPart;
    cout << endl;
  cout << "Insert the number of MC cycles = ";
  cin >> numMCCycles;
    cout <<  endl;

}
// end of function initialise
