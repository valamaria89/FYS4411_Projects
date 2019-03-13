#include <iostream>
#include <random>
#include "position.h"
#include "system.h"
#include "matrix.h"
#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1
#define a 0.0043


System::System(int the_Np, int the_dim, double the_alpha, double b, double st): Position(the_Np, the_dim, b, st), alpha(the_alpha){

}

double System::wavefunction(){

    double wf=0, r_sqr=0;
    r_sqr = getrsqrt();

    wf = exp(-get_alpha()*r_sqr) ;
    return wf;

}


double System::local_energy_num(){
    double  wfminus, wfplus, ekin=0, epot=0, eloc=0, radsq = 0;
    System p_plus(r.get_NumOfPart(), r.get_dim(), get_alpha(), get_beta(), get_step()), p_minus(r.get_NumOfPart(), r.get_dim(), get_alpha(), get_beta(), get_step());

    p_plus.setPosition(getPosition(r.get_NumOfPart(), r.get_dim()));
    p_minus.setPosition(getPosition(r.get_NumOfPart(), r.get_dim()));
    radsq = getrsqrt();


  for (int i = 0; i < r.get_NumOfPart(); i++){
           for (int j = 0; j <r.get_dim(); j++){
               p_plus.set_Elem(i,j, r.get_Elem(i,j)+h);
               p_minus.set_Elem(i, j,r.get_Elem(i,j)-h);
               wfminus = p_minus.wavefunction();
               wfplus  = p_plus.wavefunction();
               ekin += (wfminus+wfplus-2*wavefunction());
               p_plus.set_Elem(i,j, get_Elem(i,j));
               p_minus.set_Elem(i, j, get_Elem(i,j));
           }
       }
    ekin = -0.5*ekin*h2/wavefunction();

    epot = radsq/2;
    eloc = epot+ekin;

    return eloc;
}


double System::local_energy_analytic(){
    double  ekin = 0,epot=0, eloc=0, radsq = 0;

    radsq = getrsqrt();
    ekin = -2*get_alpha()*get_alpha()*radsq+get_alpha()*r.get_dim()*r.get_NumOfPart();
    epot = radsq/2;
    eloc = epot+ekin;
    return eloc;
}

double System::wavefunction_int(){
    double wf=0, r_sqr=0;

    r_sqr = getrsqrt();
    wf = exp(-get_alpha()*r_sqr);

          return wf;
}

double System::local_energy_analytical_int(){
    double  ekin = 0,epot=0, eloc=0, radsq = 0, r_ij=0,r_ik=0;
    double ekin_1=0, ekin_2=0, ekin_3=0, ekin_4=0;

    //This is the first term of kinetic energy which works fine!
        for(int i=0; i<r.get_NumOfPart(); i++){
            ekin_1+=-4*get_alpha()-2*get_alpha()*get_beta()+4*get_alpha()*get_alpha()*(r.get_Elem(i,0)*r.get_Elem(i,0)+r.get_Elem(i,1)*r.get_Elem(i,1)+get_beta()*get_beta()*r.get_Elem(i,2)*r.get_Elem(i,2));
        }
    //This is the second term
         for(int i=0; i<r.get_NumOfPart(); i++){
             for(int j=0; j<r.get_NumOfPart(); j++){
                 if(j!=i){
                     for ( int di = 0; di < r.get_dim(); di++){
                         r_ij += (r.get_Elem(i, di) - r.get_Elem(j, di))*(r.get_Elem(i, di) - r.get_Elem(j, di));
                      }
                     r_ij=sqrt(r_ij);
                     ekin_2+= 2*a /(r_ij*r_ij+a*r_ij)/r_ij*(-2*get_alpha()*r.get_Elem(i,0)*(r.get_Elem(i,0)-r.get_Elem(j,0))-2*get_alpha()*r.get_Elem(i,1)*(r.get_Elem(i,1)-r.get_Elem(j,1))-2*get_alpha()*get_beta()*r.get_Elem(i,2)*(r.get_Elem(i,2)-r.get_Elem(j,2)));
                     r_ij=0;
                 }
             }
         }
    //This is the third term
         for(int i=0; i<r.get_NumOfPart(); i++){
             for(int j=0; j<r.get_NumOfPart(); j++){
                 if(j!=i){
                     for(int k=0; k<r.get_NumOfPart(); k++){
                         if(k!=i){
                             for ( int di = 0; di < r.get_dim(); di++){
                                 r_ij += (r.get_Elem(i, di) - r.get_Elem(j, di))*(r.get_Elem(i, di) - r.get_Elem(j, di));
                                 r_ik += (r.get_Elem(i, di) - r.get_Elem(k, di))*(r.get_Elem(i, di) - r.get_Elem(k, di));
                              }
                      r_ij=sqrt(r_ij);
                      r_ik=sqrt(r_ik);
                      ekin_3+=a*a/((r_ij*r_ij+a*r_ij)*(r_ik*r_ik+a*r_ik)*r_ij*r_ik)*((r.get_Elem(i,0)-r.get_Elem(j,0))*(r.get_Elem(i,0)-r.get_Elem(k,0))+(r.get_Elem(i,1)-r.get_Elem(j,1))*(r.get_Elem(i,1)-r.get_Elem(k,1))+(r.get_Elem(i,2)-r.get_Elem(j,2))*(r.get_Elem(i,2)-r.get_Elem(k,2)));
                      r_ik=0;
                      r_ij=0;
                      }
                    }
                  }
                }
              }

    //This is the fourth term
         for(int i=0; i<r.get_NumOfPart(); i++){
             for(int j=0; j<r.get_NumOfPart(); j++){
                 if(j!=i){
                     for ( int di = 0; di < get_dim(); di++){
                         r_ij += (r.get_Elem(i, di) - r.get_Elem(j, di))*(r.get_Elem(i, di) - r.get_Elem(j, di));
                      }
                     r_ij=sqrt(r_ij);
                     ekin_4+=-a*(2*r_ij+a)/((r_ij*r_ij+a*r_ij)*(r_ij*r_ij+a*r_ij))+2/r_ij*a/(r_ij*r_ij+a*r_ij);
                     r_ij=0;
                 }
             }
         }

        radsq = getrsqrt();
        ekin = -(ekin_1+ekin_2+ekin_3+ekin_4)/2;
        epot = radsq/2;
        eloc = epot+ekin;
        return eloc;
}

void System::QuantumForce_Analytical(){
        for (int i = 0; i < r.get_NumOfPart(); i++){
             for (int j = 0; j < r.get_dim(); j++){
                 if (j == 2){
                     r.set_Elem(i, j, -get_beta()*get_alpha()*r.get_Elem(i,j));}
                 else {r.set_Elem(i, j, -get_alpha()*r.get_Elem(i,j));}
                 }
            }
       }





void System::QuantumForce_num(){
    System p_plus(r.get_NumOfPart(), r.get_dim(), get_alpha(), get_beta(), get_step()), p_minus(r.get_NumOfPart(), r.get_dim(), get_alpha(), get_beta(), get_step());

    p_plus.setPosition(getPosition(r.get_NumOfPart(), r.get_dim()));
    p_minus.setPosition(getPosition(r.get_NumOfPart(), r.get_dim()));

    double wfold ;
    wfold= wavefunction();
    double  wfminus = 0, wfplus = 0;

        for (int i = 0; i < r.get_NumOfPart(); i++){
                 for (int j = 0; j < r.get_dim(); j++){
                     p_plus.set_Elem(i,j, r.get_Elem(i,j)+h);
                     p_minus.set_Elem(i, j, r.get_Elem(i,j)-h);
                     wfminus = p_minus.wavefunction();
                     wfplus  = p_plus.wavefunction();
                     r.set_Elem(i, j, (wfplus-wfminus)/(wfold*h));
                     p_plus.set_Elem(i,j, get_Elem(i,j));
                     p_minus.set_Elem(i, j, get_Elem(i,j));
              }
         }
}


double System::derWF(){
    double wfder=0, r_sqr=0;
    r_sqr = getrsqrt();

    wfder = -r_sqr;
    return wfder;
}




