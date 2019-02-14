#ifndef METHODS_H
#define METHODS_H
#include <random>
#include <iostream>
#include <armadillo>
#include "matrix.h"

using namespace arma;
using namespace std;
//double alpha;
//double deltaAlpha;

class Methods
{
private:
    int dim;
    int numOfPart;
    int numMCCycles;
    double *Etot, *Etot2;
    int numVar;
    double step;
    int ind;
    double alpha;
    double deltaAlpha;
    mat rOld;
    mat rNew;
    //mat QFOld;
    //mat QFNew;
    double wavefunction(const mat &r);
    double local_energy_analytic(const mat &r);
    double local_energy_num(const mat &r);
    void QuantumForce(const mat &r, mat &QForce);

public:


    Methods();
    Methods(int di, int numPart, int MCCycles, int numV,double stepsz, int id, double al, double dal){
       dim = di;
       numOfPart = numPart;
       numMCCycles = MCCycles;
       numVar = numV;
       Etot = new double[numV+1];
       Etot2 = new double[numV+1];
       step = stepsz;
       ind = id;
       alpha = al;
       deltaAlpha = dal;
    }
   // Methods(double alp, double dalp){
     //   alpha = alp;
       // deltaAlpha = dalp;
    //}

    //Methods(int var){
      //  Etot = new double[var+1];
    //}


    ~Methods(){
        delete [] Etot;
        delete [] Etot2;
        cout << "woooooorking" << endl;
    }
    //void setdeltaAlpha(double deltaAlpha){deltaAlpha = deltaAlpha; };
    //double get_deltaAlpha(){return deltaAlpha;};
    //void setalpha(double alp, double deltaAl){ alpha = alp+deltaAl;};
    //double getalpha(){return alpha;};
    //void setEtot(int var, double value){Etot[var] = value;};
    //double getEtot(int var){return Etot[var];};
    void mc_sampling();


    //void mc_sampling_IMS();
    //void QuantumForce(Matrix r, Matrix QForce);
    void writeToFile();

};
#endif // METHODS_H
