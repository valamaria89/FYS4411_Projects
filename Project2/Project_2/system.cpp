#include <iostream>
#include <random>
#include "system.h"
#include "matrix.h"


System::System(int inp_N, int inp_M, double inp_sigma, double inp_omega):s_N(inp_N), s_M(inp_M),s_sigma(inp_sigma),s_omega(inp_omega){
}

double System::localEnergy(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma, int samp){
   double coeff=1;
   VectorXd u = b + (x.transpose()*W).transpose()/pow(sigma,2);
   if(samp==2)
   {
     coeff=0.5;
   }
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
         double d1lnPsi = -(x(i)-a(i))/pow(sigma,2)*coeff+sum_1/pow(sigma,2)*coeff;
         double d2lnPsi = -1/pow(sigma,2)*coeff+sum_2/pow(sigma,4)*coeff;

         E_kin += -(d2lnPsi+d1lnPsi*d1lnPsi)/2;
         E_pot += x(i)*x(i)*s_omega*s_omega/2;
   }

   E_loc = E_kin+E_pot;

   return E_loc;
}


VectorXd System::grad_a(VectorXd x, VectorXd a, double sigma,int samp){
   double coeff=1;
   if(samp==2)
   {
     coeff=0.5;
   }
   return (x-a)/pow(sigma,2)*coeff;


}

VectorXd System::grad_b(VectorXd x, VectorXd b, MatrixXd W, double sigma, int samp){
     double coeff=1;
     if(samp==2)
     {
       coeff=0.5;
     }
     VectorXd u = b + (x.transpose()*W).transpose()/pow(sigma,2);
    return logSigm(u)*coeff;

}


MatrixXd System::grad_W(VectorXd x, VectorXd b, MatrixXd W, double sigma, int samp){
        double coeff=1;
        if(samp==2)
        {
          coeff=0.5;
        }
        VectorXd u = b + (x.transpose()*W).transpose()/pow(sigma,2);

       return x*logSigm(u).transpose()/pow(sigma,2)*coeff;

}

double System::grad_sigma(VectorXd x,VectorXd a, VectorXd b, MatrixXd W, double sigma){

        VectorXd u = b + (x.transpose()*W).transpose()/pow(sigma,2);
       VectorXd xa = x-a;
       double term = (W.transpose()*x).transpose()*logSigm(u);
       return (xa.transpose()*xa+term)/pow(sigma,3);
}

double System::QuantumForce(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma, int i){

     VectorXd u = b + (x.transpose()*W).transpose()/pow(sigma,2);
    return (-2*(x(i)-a(i))+2*W.row(i)*logSigm(u))/pow(sigma,2);
}

double System::waveFunction(VectorXd x, VectorXd a,VectorXd b, MatrixXd W, double sigma){
     VectorXd u = b + (x.transpose()*W).transpose()/pow(sigma,2);
    double prod1 = 1;
    VectorXd Xa = x-a;
    double prod2 = Xa.transpose()*Xa;
    for (int i = 0; i < s_N; i++){
        prod1 += (1+ exp(u(i)));

    }
    return exp(-prod2/(2*sigma*sigma))*prod1;
}

VectorXd System::logSigm(VectorXd u){
 VectorXd res = VectorXd::Zero(s_N);

 for (int i = 0; i < s_N ; i++) {
        res(i) = 1/(1+exp(-u(i)));
    }

  return res;
}

double System::GreenFunction(VectorXd xold,VectorXd xnew, VectorXd a,VectorXd b, MatrixXd W, double sigma, int par, int P,double difConst,double timeStep){
    double total_G  = 0;
    int dim = s_M/P;
    if(par==0){
    for(int i=0; i<P; i++) {
            double sum_G = 0;
            for(int j=0; j<dim; j++) {
                int index = dim*i+j;
                double QForceOld = QuantumForce(xold,a,b,W,sigma,index);
                double QForceNew = QuantumForce(xnew,a,b,W,sigma,index);
                sum_G += 0.5*(QForceOld + QForceNew) * (0.5*difConst*timeStep*(QForceOld - QForceNew)-xnew(index)+xold(index));
            }
            total_G += exp(sum_G);
        }
        return total_G;
    }
    if(par==1){
    for(int i=0; i<P; i++) {
            for(int j=0; j<dim; j++) {
                int index = dim*i+j;
                double QForceOld = QuantumForce(xold,a,b,W,sigma,index);
                double QForceNew = QuantumForce(xnew,a,b,W,sigma,index);
                total_G += 0.5*(QForceOld + QForceNew) * (0.5*difConst*timeStep*(QForceOld - QForceNew)-xnew(index)+xold(index));
            }
        }
        return exp(total_G);
    }


}
