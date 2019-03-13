#include "methods_interaction.h"

#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1
<<<<<<< HEAD
//#define a 0.43
#define a 0.43
=======
#define a 0.00

>>>>>>> 0490687d288bc3175cb820ac966b9180b7193c8e

using namespace std;

//*********************************BRUTE FORCE**************************************

void mc_sampling_INT( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization, double beta)
{
//The variational parameter
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;
        int counter=0;
//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rold(numOfPart, dim);
        Matrix rnew(numOfPart, dim);
        ofstream fileblock;
        fileblock.open("E_block.txt");

      for (int i = 0; i < numOfPart; i++){
               for (int j = 0; j < dim; j++){
                   rold.set_Elem(i, j, 0);
                   rnew.set_Elem(i ,j, 0);

               }
           }

      //cout<<local_energy_analytic(rold, 0.5, dim, numOfPart)<<endl;

      random_device rd;
      mt19937_64 gen(rd());
      uniform_real_distribution<double> RNG(0.0,1.0);

//The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 //  initial trial position
        for (int i = 0; i < numOfPart; i++) {

             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, step*(RNG(gen)-0.5));
                 }
           }


 // initial trial WF
    wfold = wavefunction(rold, alpha, beta, dim, numOfPart);
    //One body densities

           int numOfBins =120;
           double radMax =3;
           double radStep;
           radStep=radMax/numOfBins;
           double *volumes, *partInBins,*distances;
           volumes = new double[numOfBins];
           partInBins = new double[numOfBins];
           distances = new double[numOfBins];

           for(int i = 0; i < numOfBins; i++){
               distances[i]=i*radStep;
               partInBins[i]=0;
           }
        // This part we don't use
           volumes[0]=(4*3.14/3)*pow(distances[0],3);
           for(int j = 1; j < numOfBins; j++) {
                  volumes[j] = (4*M_PI/3)*pow((distances[j]), 3) - volumes[j-1];
              }

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position
      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+step*(RNG(gen)-0.5));
               }
         }
    wfnew = wavefunction(rnew, alpha, beta, dim, numOfPart);
     //cout<< wfnew<<endl;
 // Metropolis test
    if(RNG(gen) < wfnew*wfnew/(wfold*wfold)){
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rold.set_Elem(i, j, rnew.get_Elem(i,j));
                 }
           }
  wfold=wfnew;
  accept++;
    }
    if (cycl > thermalization){
 // compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rold, alpha, beta, dim, numOfPart);
    }
        else{
            delta_e = local_energy_num(rold, alpha, beta, wfold, dim, numOfPart);
    }
// update energies
    eng  += delta_e;
    eng2 += delta_e*delta_e;
    for(int l=0; l< numOfPart; l++){
      double rdist = sqrt(rold.get_Elem(l,0)*rold.get_Elem(l,0) + rold.get_Elem(l,1)*rold.get_Elem(l,1) +rold.get_Elem(l,2)*rold.get_Elem(l,2));
      int bin = 0;
      //double err = 1000000;
    for(int k=1; k<numOfBins; k++) {
      if((distances[k-1]<=rdist)&&(rdist < distances[k]))
        //double e = fabs(volumes[k] - rdist);
        //if(e < err) {
             //err = e;
             bin = k;
                 //}
               }
    partInBins[bin] += 1;
             }
    if(alpha==0.5){
    if(!((cycl) % 1000)){
        if(counter<512){
        fileblock <<eng/cycl<< endl;
        }
        counter++;
      }
    }
  }
// End of thermalization loop
        }
// End of loop over MC trials
    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
    int ar=0;
    //cout << " " << accept <<endl;
    cout << "alpha: " << alpha << " Energy " << Etot[var] <<" Variance " << Etot2[var] - Etot[var]* Etot[var] << endl;
    for(int k=0; k<numOfBins; k++){
           cout<<distances[k]<<" "<< partInBins[k]/(numMCCycles)<< endl;
           ar+=partInBins[k];
           //cout<<distances[k]<<" "<< partInBins[k]/(numMCCycles*volumes[k])<< endl;
    }
    //for (int i = 0; i < numOfPart; i++) {
         //for (int j = 0; j < dim; j++){
           // cout<< rold.get_Elem(i,j)<<" ";
           //  }
         //cout<<endl;
       //}
    cout<<" Numb= "<< ar/(numMCCycles)<<endl;
// Increase the variational parameter
        alpha += deltaAlpha;
        delete [] volumes;
        delete [] partInBins;
        delete [] distances;
        //cout << wfold << endl;
 }

// End of loop over variational  steps
        fileblock.close();

}

//*********************************IMPORTANCE SAMPLING**************************************

void mc_sampling_IMS_INT( double timestep, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization,double beta)
{
// The variational parameter
       // double alpha=0.50, deltaAlpha=0.01;
        double wfnew, wfold, eng, eng2, delta_e;
        int accept;
        int counter=0;

// Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rold(numOfPart, dim);
        Matrix rnew(numOfPart, dim);
        Matrix QForceNew(numOfPart, dim);
        Matrix QForceOld(numOfPart, dim);
        ofstream fileblock;
        fileblock.open("E_block.txt");

      for (int i = 0; i < numOfPart; i++){
               for (int j = 0; j < dim; j++){
                   rold.set_Elem(i,j, 0);
                   rnew.set_Elem(i,j, 0);
                   QForceNew.set_Elem(i,j,0);
                   QForceOld.set_Elem(i,j, 0);
               }
           }

      random_device rd;
      mt19937_64 gen(rd());
      uniform_real_distribution<double> RNG(0.0,1.0);
      normal_distribution<double> ND(0.0,1.0);

// The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 // Initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, sqrt(timestep)*(ND(gen)));
                 }
           }

 // Initial trial WF
    wfold = wavefunction(rold, alpha, beta, dim, numOfPart);
    QuantumForce(rold, QForceOld, dim, numOfPart, alpha, beta);


 // The MC cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position

      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+sqrt(timestep)*ND(gen)+0.5*QForceOld.get_Elem(i,j)*timestep);
               }
         }

    wfnew = wavefunction(rnew, alpha, beta, dim, numOfPart);
    QuantumForce(rnew, QForceNew, dim, numOfPart, alpha, beta);
    double GreenFunction = 0;
     for (int i = 0; i < numOfPart; i++) {
         for (int j = 0; j < dim; j++){

             GreenFunction += 0.5*(QForceOld.get_Elem(i,j) + QForceNew.get_Elem(i,j))*(0.25*timestep*(QForceOld.get_Elem(i, j)-QForceNew.get_Elem(i, j))+rold.get_Elem(i, j) - rnew.get_Elem(i,j));
         }
     }
     GreenFunction = exp(GreenFunction);

 // Metropolis-Hastings test

    if(RNG(gen) <= GreenFunction*wfnew*wfnew/(wfold*wfold) ) {
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rold.set_Elem(i, j, rnew.get_Elem(i,j));
                    QForceOld.set_Elem(i, j, QForceNew.get_Elem(i,j));
                 }
           }
  wfold=wfnew;
  accept++;
    }
    if (cycl > thermalization){
 // Compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rold, alpha, beta, dim, numOfPart);
        }
        else{
            delta_e = local_energy_num(rold, alpha, beta, wfold, dim, numOfPart);
        }
// Update energies
        eng  += delta_e;
        eng2 += delta_e*delta_e;

        if(alpha==0.5){
        if(!((cycl) % 1000)){
            if(counter<512){
            fileblock <<eng/cycl<< endl;
            }
            counter++;
          }
        }
      }
    }
// End of loop over MC trials
    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
    //cout << " " << accept <<endl;
    //cout << "alpha: " << alpha << " Energy " << Etot[var] <<" Variance " << Etot2[var] - Etot[var]* Etot[var] << endl;

//Increase the variational parameter
        alpha += deltaAlpha;
 }
// End of loop over variational  steps
    fileblock.close();

 }


//******************************GRADIENT DESCENT BRUTE FORCE*****************************

void gradiendescent_brute_INT( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization, double beta)
{
//The variational parameter
        double wfnew, wfold, eng, eng2, delta_e,engDer,delta_wf,wfEngDeriv,deriv_wf;
        int accept;
        double alphagrad=0;
        double energy_average=0;
        double energy_average2=0,eng_derivative=0;
        double *count;
        count = new double[10];

//Here we set the old and the new positions to zero (2d matrices [number of particles, dimension])
        Matrix rold(numOfPart, dim);
        Matrix rnew(numOfPart, dim);

      for (int i = 0; i < numOfPart; i++){
               for (int j = 0; j < dim; j++){
                   rold.set_Elem(i, j, 0);
                   rnew.set_Elem(i ,j, 0);

               }
           }

      random_device rd;
      mt19937_64 gen(rd());
      uniform_real_distribution<double> RNG(0.0,1.0);

//The loop over the variational parameter
        for (int var=0; var<=numVar; var++) {
            eng = 0;
            eng2 = 0;
            accept =0;
            delta_e=0;

 //  initial trial position
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                   rold.set_Elem( i, j, step*(RNG(gen)-0.5));
                 }
           }

 // initial trial WF
    wfold = wavefunction(rold, alpha, beta, dim, numOfPart);
    //cout<< wfold<<endl;

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position
      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+step*(RNG(gen)-0.5));
               }
         }
    wfnew = wavefunction(rnew, alpha, beta, dim, numOfPart);

 // Metropolis test
    if(RNG(gen) < wfnew*wfnew/(wfold*wfold)){
        for (int i = 0; i < numOfPart; i++) {
             for (int j = 0; j < dim; j++){
                    rold.set_Elem(i, j, rnew.get_Elem(i,j));
                 }
           }
  wfold=wfnew;
  accept++;
    }

 // compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rold, alpha, beta, dim, numOfPart);
            delta_wf = derWF(rold, alpha, beta, dim, numOfPart);
    }
        else{
            delta_e = local_energy_num(rold, alpha, beta, wfold, dim, numOfPart);
            delta_wf = derWF(rold, alpha, beta, dim, numOfPart);
    }

// update energies
    eng  += delta_e;
    eng2 += delta_e*delta_e;
    deriv_wf +=delta_wf;
    wfEngDeriv +=delta_wf*delta_e;

        }
// end of loop over MC trials
    energy_average = eng/numMCCycles;
    energy_average2 = eng2/numMCCycles;
    wfEngDeriv=wfEngDeriv/numMCCycles;
    deriv_wf=deriv_wf/numMCCycles;

    eng_derivative=2*(wfEngDeriv-energy_average*deriv_wf);


//Increase the variational parameter
        alphagrad=eng_derivative/numOfPart;
        alpha-=deltaAlpha*alphagrad;

//The selection rules for the gradient descent stop

//If alpha decrease is less than epsilon
        cout<< "alpha="<<alpha<<" "<<deltaAlpha*alphagrad<<endl;
       if(abs(deltaAlpha*alphagrad)<10e-6){
           var=numVar;
        }

//If alpha over ten last steps fluctuates less than epsilon
       double aver=0,std=0,avers=0;
       if(var<9){
           count[var]=alphagrad;}
       else{
           for(int m=0; m<9; m++){
               count[m]=count[m+1];
               aver+=count[m];
           }
          count[9]=alphagrad;
          aver=(aver+alphagrad)/10.0;
       }

       if(var>9)
       {for(int m=0; m<10; m++){
               avers+=(count[m]-aver)*(count[m]-aver);
           }
           std=sqrt(avers/10.0);
           if(std<10e-8)
           {var=numVar;}
       }

 }

// end of loop over variational  steps

}



//The wave function with interaction
double wavefunction(Matrix r, double alpha, double beta, int dim, int NumOfPart){

      double wf=0, r_sqr=0, r_ik =0, arg = 0;
      int p = 1;
      for (int i = 0; i < NumOfPart; i++) {

        for (int j = 0; j < dim; j++) {
            if (j == 2){
                r_sqr  += beta*r.get_Elem(i,j)*r.get_Elem(i,j);
        }
            else r_sqr += r.get_Elem(i,j)*r.get_Elem(i,j);
        }
        for (int k = i+1; k < NumOfPart; k++){
            for ( int di = 0; di < dim; di++){
                r_ik += (r.get_Elem(i, di) - r.get_Elem(k, di))*(r.get_Elem(i, di) - r.get_Elem(k, di));
             }
            r_ik = sqrt(r_ik);
            if( r_ik <= a){
<<<<<<< HEAD
                p = 0;
            }
            else arg += log(1-a/(r_ik));
=======
                p = 0;      
            }
            else arg += log(1-a/(r_ik));

            r_ik=0;
        }
>>>>>>> 0490687d288bc3175cb820ac966b9180b7193c8e

            r_ik=0;
        }
      }

      wf = exp(-alpha*r_sqr+arg);

      return wf*p;
    }

// The numerical calculation of local energy
double local_energy_num(Matrix r, double alpha, double beta, double wfold, int dim, int numOfPart){
      double  wfminus, wfplus, ekin=0, epot=0, eloc=0, radsq = 0;
      Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 if (j == 2){
                     radsq  += beta*beta*r.get_Elem(i,j)*r.get_Elem(i,j);
             }
                 else {radsq += r.get_Elem(i,j)*r.get_Elem(i,j);
                 }

                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
             }
           }

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 wfminus = wavefunction(rminus, alpha, beta, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, beta, dim, numOfPart);
                 ekin += (wfminus+wfplus-2*wfold);
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
             }
         }
      ekin = -0.5*ekin*h2/wfold;

      epot = radsq/2;
      eloc = epot+ekin;

      return eloc;
  }


//The alytical calculation of local energy
double local_energy_analytic(Matrix r, double alpha, double beta, int dim, int numOfPart){
    double  ekin = 0,epot=0, eloc=0, radsq = 0, r_ij=0,r_ik=0;
    double ekin_1=0, ekin_2=0, ekin_3=0, ekin_4=0;

//This is the first term of kinetic energy which works fine!
    for(int i=0; i<numOfPart; i++){
        ekin_1+=-4*alpha-2*alpha*beta+4*alpha*alpha*(r.get_Elem(i,0)*r.get_Elem(i,0)+r.get_Elem(i,1)*r.get_Elem(i,1)+beta*beta*r.get_Elem(i,2)*r.get_Elem(i,2));
    }
//This is the second term
     for(int i=0; i<numOfPart; i++){
         for(int j=0; j<numOfPart; j++){
             if(j!=i){
                 for ( int di = 0; di < dim; di++){
                     r_ij += (r.get_Elem(i, di) - r.get_Elem(j, di))*(r.get_Elem(i, di) - r.get_Elem(j, di));
                  }
                 r_ij=sqrt(r_ij);
                 ekin_2+= 2*a /(r_ij*r_ij+a*r_ij)/r_ij*(-2*alpha*r.get_Elem(i,0)*(r.get_Elem(i,0)-r.get_Elem(j,0))-2*alpha*r.get_Elem(i,1)*(r.get_Elem(i,1)-r.get_Elem(j,1))-2*alpha*beta*r.get_Elem(i,2)*(r.get_Elem(i,2)-r.get_Elem(j,2)));
                 r_ij=0;
             }
         }
     }
//This is the third term
     for(int i=0; i<numOfPart; i++){
         for(int j=0; j<numOfPart; j++){
             if(j!=i){
                 for(int k=0; k<numOfPart; k++){
                     if(k!=i){
                         for ( int di = 0; di < dim; di++){
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
     for(int i=0; i<numOfPart; i++){
         for(int j=0; j<numOfPart; j++){
             if(j!=i){
                 for ( int di = 0; di < dim; di++){
                     r_ij += (r.get_Elem(i, di) - r.get_Elem(j, di))*(r.get_Elem(i, di) - r.get_Elem(j, di));
                  }
                 r_ij=sqrt(r_ij);
                 ekin_4+=-a*(2*r_ij+a)/((r_ij*r_ij+a*r_ij)*(r_ij*r_ij+a*r_ij))+2/r_ij*a/(r_ij*r_ij+a*r_ij);
                 r_ij=0;
             }
         }
     }
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 if (j == 2){
                     radsq  += beta*beta*r.get_Elem(i,j)*r.get_Elem(i,j);
                   }
                 else {radsq += r.get_Elem(i,j)*r.get_Elem(i,j);
             }
          }
      }
    ekin = -(ekin_1+ekin_2+ekin_3+ekin_4)/2;
    epot = radsq/2;
    eloc = epot+ekin;
    return eloc;
}
// The analytical quantum force calculation
void QuantumForce(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha, double beta){
    Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                     if (j == 2){
                 QForce.set_Elem(i, j, -beta*alpha*r.get_Elem(i,j));}
                    else {QForce.set_Elem(i, j, -alpha*r.get_Elem(i,j));}
          }
     }
}
// The numerical quantum force calculation
void QuantumForce_num(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha,double beta){
    Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
             }
    }
    double wfold ;
    wfold=wavefunction(r, alpha, beta, dim, numOfPart);
    double  wfminus = 0, wfplus = 0;

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 wfminus = wavefunction(rminus, alpha, beta, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, beta, dim, numOfPart);
                 QForce.set_Elem(i, j, (wfplus-wfminus)/(wfold*h));
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
          }
     }
}

double derWF(Matrix r, double alpha, double beta, int dim, int NumOfPart){

      double wfder=0, r_sqr=0;

      for (int i = 0; i < NumOfPart; i++) {
        for (int j = 0; j < dim; j++) {
          if (j == 2){
              r_sqr   += beta*r.get_Elem(i,j)*r.get_Elem(i,j);
            }
          else { r_sqr  += r.get_Elem(i,j)*r.get_Elem(i,j);
      }
   }
}
      wfder = -r_sqr;
      return wfder;
    }
