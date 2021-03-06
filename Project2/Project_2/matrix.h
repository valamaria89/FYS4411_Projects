#ifndef MATRIX_H
#define MATRIX_H

#pragma once
#include <iostream>

using namespace std;

class Matrix {

private:

    int NumOfPart;
    int dim;
    double ** r;

public:

    Matrix();
    Matrix(int NP, int di){
        r = new double * [NP];
        for (int i = 0; i < NP; i++ ){
            r[i] = new double [di];
        }
        dim = di;
        NumOfPart = NP;
    }

    ~Matrix(){
        }

    void set_Elem(int N, int M, double value){
        r[N][M] =  value;
    }

    int& get_dim(){
        return dim;
    }
    int& get_NumOfPart(){
       return NumOfPart;
    }

    double get_Elem(int N, int M){
        return r[N][M];
    }


};

#endif // MATRIX_H
