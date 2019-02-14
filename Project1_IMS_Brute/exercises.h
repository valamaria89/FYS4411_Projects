
#include "methods.h"
#include <random>

using namespace std;

/*void exe_b(){

    double alpha = 1, deltalpha = 0.1;
    for (int i = 0; i < 2; i++){



        int MC = 100000000;
        double a_ho = 1.; //CHECK THIS UNIT

    //double alpha = 1./(2*a_ho*a_ho);


        double R = 5;

        double step = 10E-3;

        double energy = 0; double energy2 = 0; double dE; double wf;

        std::mt19937_64 seed(1234); //Mersienne twister

        for (int cycle = 0; cycle < MC; cycle++){

            dE = 0; wf = 0;
            wf = wavefunction(alpha, R);
            Metropolis(R, seed, step, alpha, wf);
            dE = local_energy(step, alpha, R, wf);
            energy += dE;
            energy2 += dE*dE;


    }
        double energy_norm = energy/MC; double energy2_norm = energy2/MC;
        double variance = energy2_norm-energy_norm*energy_norm;
        cout << variance << endl;
        alpha += i*deltalpha;

}
}*/
