#include <iostream>
#include <random>
#include "system.h"
#include "matrix.h"


System::System(int inp_N, int inp_M, double inp_sigma, double inp_omega):s_N(inp_N), s_M(inp_M),s_sigma(inp_sigma),s_omega(inp_omega){
}

double System::localEnergy(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma){

   VectorXd u = b + (W.transpose() * x)/pow(sigma,2);

   double E_kin = 0;
   double E_pot = 0;
   double E_loc = 0;

   for(int i =0;i<s_M;i++){
       double sum_1=0;
       double sum_2=0;
         for(int j =0;j<s_N;j++){
             sum_1+=W(i,j)*logSigm(u)(j);
             sum_2+=W(i,j)*W(i,j)*logSigm(-u)(j)*logSigm(u)(j);
         }
         double d1lnPsi = -(x(i)-a(i))/pow(sigma,2)+sum_1/pow(sigma,2);
         double d2lnPsi = -1/pow(sigma,2)+sum_2/pow(sigma,4);

         E_kin += -1/2*(d2lnPsi+d1lnPsi*d1lnPsi);
         E_pot += 1/2*x(i)*x(i)*s_omega*s_omega;
   }

   E_loc = E_kin+E_pot;

   return E_loc;
}


VectorXd System::grad_a(VectorXd x, VectorXd a, double sigma){

   return -(x-a)/pow(sigma,2);

}

VectorXd System::grad_b(VectorXd x, VectorXd b, MatrixXd W, double sigma){

    VectorXd u = b + (W.transpose() * x)/pow(sigma,2);
    return logSigm(u);

}


MatrixXd System::grad_W(VectorXd x, VectorXd b, MatrixXd W, double sigma){

       VectorXd u = b + (W.transpose() * x)/pow(sigma,2);
       return x*logSigm(u).transpose()/pow(sigma,2);

}

double System::grad_sigma(VectorXd x,VectorXd a, VectorXd b, MatrixXd W, double sigma){

       VectorXd u = b + (W.transpose() * x)/pow(sigma,2);
       VectorXd xa = x-a;
       double term = (W.transpose()*x).transpose()*logSigm(u);
       return (xa.transpose()*xa+term)/pow(sigma,3);
}

VectorXd System::QuantumForce(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma){

    VectorXd u = b + (W.transpose() * x)/pow(sigma,2);
    return (-2*(x-a)+2*logSigm(u).transpose()*W)/pow(sigma,2);
}





/*double System::waveFunction(Matrix x, Matrix a, Matrix W, double sigma){

}*/

VectorXd System::logSigm(VectorXd u){
 VectorXd res = VectorXd::Zero(s_N);

 for (int i = 0; i < s_N ; i++) {
        res(i) = 1/(1+exp(-u(i)));
    }

  return res;
}
