#include <iostream>
#include <random>
#include "position.h"
#include "matrix.h"
#define a 0.0043


using namespace std;

random_device rd;
mt19937_64 gen(rd());
uniform_real_distribution<double> RNG(0.0,1.0);

Position::Position(int Np, int Di, double b, double step):r(Np,Di), beta(b), sd(step){

}

void Position::setPosition_initial(){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
         for (int j = 0; j < r.get_dim(); j++){
              r.set_Elem(i, j, 0);

             }
      }


}

void Position::setPosition(double value){

    for (int i = 0; i < r.get_NumOfPart(); i++) {
         for (int j = 0; j < r.get_dim(); j++){
             r.set_Elem(i,j, value);


         }
         }
      }

void Position::setRandomPosition(){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
         for (int j = 0; j < r.get_dim(); j++){
              r.set_Elem(i, j, get_step()*(RNG(gen)-0.5));
         }
    }


}



void Position::setNewToOld(Matrix r){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
        for (int j = 0; j < r.get_dim(); j++){
            r.set_Elem(i, j, r.get_Elem(i,j));
                     }
}
}

double Position::getrsqrt(){
    double r_sqr = 0, r_ik = 0, arg = 0;
    int v = 1;
    for(int i = 0; i < r.get_NumOfPart(); i++){
        for(int j = 0; j < r.get_dim(); j++){
            r_sqr  += r.get_Elem(i,j)*r.get_Elem(i,j);
        }
            if (j == 2){ //beta squared??

                r_sqr  += get_beta()*r.get_Elem(i,j)*r.get_Elem(i,j);

        }
            else r_sqr += r.get_Elem(i,j)*r.get_Elem(i,j);

        if(get_beta() != 1){
            for (int k = i+1; k < r.get_NumOfPart(); k++){
              for ( int di = 0; di < r.get_dim(); di++){
                   r_ik += (r.get_Elem(i, di) - r.get_Elem(k, di))*(r.get_Elem(i, di) - r.get_Elem(k, di));
              }
              r_ik = sqrt(r_ik);
              if( r_ik <= a){
                 v = 0;
             }
              else arg += log(1-a/(r_ik));

              r_ik=0;

        }
            else arg = 0;
}

}

    return r_sqr;
   }



ostream& operator<<(ostream& os, const Matrix & r)
{

        for (int i = 0; i < r.get_NumOfPart(); i++) {
            os << endl;
            for (int j = 0; j < r.get_dim(); j++) {
                os << r.get_Elem(i,j) << " ";
            }
        }



    return os;
}

