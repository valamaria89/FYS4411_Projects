#include <iostream>
#include <random>
#include "system.h"
#include "matrix.h"

System::System(int inp_N, int inp_M, double inp_sigma, double inp_omega):N(inp_N), M(inp_N),sigma(inp_sigma),omega(inp_omega){
}

double System::waveFunction(Matrix x, Matrix a, Matrix W, double sigma){
    x.set_Elem(N-1,M-1,67);
    return x.get_Elem(N-1,M-1);
}

double System::localEnergy_IMS(Matrix x, Matrix a, Matrix W, double sigma){
    return 15;
}

double System::localEnergy_Gibbs(Matrix x, Matrix a, Matrix W, double sigma){
    return 15;
}

void System::QuantumForce_IMS(Matrix x, Matrix a, Matrix W, Matrix qForce){

}

void System::QuantumForce_Gibbs(Matrix x, Matrix a, Matrix W, Matrix qForce){

}
void System::grad_a(Matrix x, Matrix a, Matrix W, Matrix psiDa){

}

void System::grad_b(Matrix x, Matrix a, Matrix W, Matrix psiDb){

}

void System::grad_W(Matrix x, Matrix a, Matrix W, Matrix psiDW){

}

void System::grad_sigma(Matrix x, Matrix a, Matrix W, Matrix psiDW){

}

double System::logSigm(Matrix x, Matrix a, Matrix b, Matrix W, double sigma, int k){
 double res=0;
 for (int i = 0; i < M ; i++) {
        res+=x.get_Elem(i,0)*W.get_Elem(i,k-1);
    }
  res/=pow(sigma,2);
  res+=b.get_Elem(k-1,0);
  return exp(res);
}
