#include <fstream>
#include <random>
#include <cassert>
#include <iostream>
#include "methods.h"
#include "system.h"
#include "position.h"
using namespace std;

// Constants
#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1
#define a 0.0043


void mc_sampling(int numVar, int numOfPart, int dim, int numMCCycles, int ind, double alpha, double step, double delAl, double *Etot, double *Etot2, double be)
{


//The variational parameter
        double wfold, wfnew, eng, eng2, delta_e;
        int accept;


        random_device rd;
        mt19937_64 gen(rd());
        uniform_real_distribution<double> RNG(0.0,1.0);


//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rOld(numOfPart, dim);
        Matrix rNew(numOfPart, dim);
        Position* pos = new Position();

        pos->setPosition_initial(rOld);
        pos->setPosition_initial(rNew);

        for (int var=0; var<=numVar; var++) {
             accept=0;
             eng = 0;
             eng2 = 0;
             delta_e=0;
             System s(alpha, be);
             pos->setRandomPosition(rOld, step, RNG(gen));
 // initial trial WF

             wfold = s.wavefunction(rOld);

 //The MC Cycle
            for (int cycl = 1; cycl <= numMCCycles; cycl++) {

// We set the new position
                for (int i = 0; i < numOfPart; i++) {
                     for (int j = 0; j < dim; j++){
                        rNew.set_Elem(i, j, rOld.get_Elem(i,j)+step*(RNG(gen)-0.5));
                         }
                   }



            //  pos->setRandomStep(rOld, rNew, step, RNG(gen)); //this can be called set position


             wfnew = s.wavefunction(rNew);

 // Metropolis test
            if(RNG(gen) < wfnew*wfnew/(wfold*wfold)){
                pos->setNewToOld(rOld, rNew);
                wfold=wfnew;
                accept++;
    }
// compute local energy
              if (ind==1){
                 delta_e = s.local_energy_analytic(rOld);

                  }
              else{
                 delta_e = s.local_energy_num(rOld);
    }
// update energies
            eng  += delta_e;
            eng2 += delta_e*delta_e;


}
// End of loop over MC trials

            Etot[var] = eng/numMCCycles;
            Etot2[var] = eng2/numMCCycles;
            //cout << "Energy: " << Etot[var] << "  alpha: " << alpha << endl;
            //Increase the variational parameter
            alpha += delAl;

}
//End of variational steps



}
