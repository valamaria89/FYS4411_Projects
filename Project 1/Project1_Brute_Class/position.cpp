#include <iostream>
#include <random>
#include "position.h"
#include "matrix.h"
#define a 0.0043


using namespace std;

//random_device rd;
//mt19937_64 gen(rd());
//uniform_real_distribution<double> RNG(0.0,1.0);



void Position::setPosition_initial(Matrix r){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
         for (int j = 0; j < r.get_dim(); j++){
              r.set_Elem(i, j, 0);
              //cout <<"Position initial: " << r.get_Elem(i,j) << endl;

             }
      }

}

void Position::setPosition_(Matrix r1, Matrix r){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
         for (int j = 0; j < r.get_dim(); j++){
               r1.set_Elem(i, j, r.get_Elem(i,j));
      }
    }
   // cout <<"D: " << r.get_dim() << endl;
//cout <<"N: " << r.get_NumOfPart() << endl;

}


void Position::setRandomPosition(Matrix r, double step, double random){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
         for (int j = 0; j < r.get_dim(); j++){
              r.set_Elem(i, j, step*(random-0.5));
              //cout <<"Position random: " << r.get_Elem(i,j) << endl;
         }
    }



}

void Position::setRandomStep(Matrix rold, Matrix rnew, double step, double random){
    for (int i = 0; i < rnew.get_NumOfPart(); i++) {
        for (int j = 0; j < rnew.get_dim(); j++){
            rnew.set_Elem(i, j,  rold.get_Elem(i,j)+step*(random-0.5));
                     }
}

}




void Position::setNewToOld(Matrix rold, Matrix rnew){
    for (int i = 0; i < rold.get_NumOfPart(); i++) {
        for (int j = 0; j < rold.get_dim(); j++){
            rold.set_Elem(i, j, rnew.get_Elem(i,j));
                     }
}
}

void Position::set_h(Matrix r, double h){
    for (int i = 0; i < r.get_NumOfPart(); i++) {
        for (int j = 0; j < r.get_dim(); j++){
            r.set_Elem(i, j, r.get_Elem(i,j)+h);
                     }
}
}

double Position::getrsqrt(Matrix r){
    double r_sqr = 0;

    for(int i = 0; i < r.get_NumOfPart(); i++){
        for(int j = 0; j < r.get_dim(); j++){
            r_sqr  += r.get_Elem(i,j)*r.get_Elem(i,j);

           }}

    return r_sqr;
   }





