#ifndef SYSTEM_H
#define SYSTEM_H

#include "position.h"


class System
{
public:
    //System(){}
    System(double the_alpha, double b);
    ~System(){}
    double wavefunction(Matrix r);
    double local_energy_num(Matrix r);
    double local_energy_analytic(Matrix r);
    double wavefunction_int(Matrix r, double b);
    double local_energy_analytical_int();
    void QuantumForce_Analytical();
    void QuantumForce_num();
    double derWF();
   // void set_alpha( double alp){alpha = alp;}
   // void set_beta(double b){beta = b;}
    double get_alpha() {return alpha;}
    double get_beta() {return beta;}
    //double get_step() {return step;}


private:
    //int m_numberofpart;
    //int m_dimensions;
    double alpha;
    double beta;
   // Position m_p;

};

#endif // SYSTEM_H
