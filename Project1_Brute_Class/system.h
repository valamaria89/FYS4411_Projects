#ifndef SYSTEM_H
#define SYSTEM_H

#include "position.h"


class System : public Position
{
public:
    System(){}
    System(int the_Np, int the_dim, double the_alpha, double b, double st);
    ~System(){}
    double wavefunction();
    double local_energy_num();
    double local_energy_analytic();
    double wavefunction_int();
    double local_energy_analytical_int();
    void QuantumForce_Analytical();
    void QuantumForce_num();
    double derWF();
   // void set_alpha( double alp){alpha = alp;}
   // void set_beta(double b){beta = b;}
    double get_alpha() {return alpha;}


protected:
    double alpha;


};

#endif // SYSTEM_H
