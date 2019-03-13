#include "methods.h"
#include "matrix.h"
#include <fstream>
#include <random>
#include <cassert>
using namespace std;

// Uniform distribution for s in [0, 1]
//std::uniform_real_distribution<double> RNG(0.0,1.0);
// Constants
#define h 0.1
#define h2 100
#define mass 1
#define omega 1

//*********************************BRUTE FORCE**************************************

void mc_sampling( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization)
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
    wfold = wavefunction(rold, alpha, dim, numOfPart);
    /* One body densities

       int numOfBins =200;
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
    */// ////////////////////

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position
      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+step*(RNG(gen)-0.5));
               }
         }
    wfnew = wavefunction(rnew, alpha, dim, numOfPart);

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

 // Compute local energy
        if (ind==1){
            delta_e = local_energy_analytic(rold, alpha, dim, numOfPart);
    }
        else{
            delta_e = local_energy_num(rold, alpha, wfold, dim, numOfPart);
    }
        /* /////////////////////////
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


        // ///////////////////////*/
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
// End of thermalization loop
        }
// End of loop over MC trials

    Etot[var] = eng/numMCCycles;
    Etot2[var] = eng2/numMCCycles;
    cout << " " << accept <<endl;
   /* for(int k=1; k<numOfBins; k++){
        cout<<distances[k]<<" "<< partInBins[k]/(numMCCycles)<< endl;
       //cout<<distances[k]<<" "<< partInBins[k]/(numMCCycles*volumes[k])<< endl;
    }*/
    for (int i = 0; i < numOfPart; i++) {
         for (int j = 0; j < dim; j++){
            cout<< rold.get_Elem(i,j)<<" ";
             }
         cout<<endl;
       }
    cout << " " << accept <<endl;
   /* delete [] volumes;
    delete [] partInBins;
    delete [] distances;*/
// Increase the variational parameter
        alpha += deltaAlpha;
 }
// End of loop over variational  steps

        fileblock.close();

}


//*********************************IMPORTANCE SAMPLING**************************************

void mc_sampling_IMS( double timestep, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization)
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
    wfold = wavefunction(rold, alpha, dim, numOfPart);
    QuantumForce(rold, QForceOld, dim, numOfPart, alpha);

 /*/ One body densities

    int numOfBins =200;
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
 // ////////////////////*/


 // The MC cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position

      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+sqrt(timestep)*ND(gen)+0.5*QForceOld.get_Elem(i,j)*timestep);
               }
         }

    wfnew = wavefunction(rnew, alpha, dim, numOfPart);
    QuantumForce(rnew, QForceNew, dim, numOfPart, alpha);
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
            delta_e = local_energy_analytic(rold, alpha, dim, numOfPart);
        }
        else{
            delta_e = local_energy_num(rold, alpha, wfold, dim, numOfPart);
        }
/*/ /////////////////////////
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


// ///////////////////////*/
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
    /*for(int k=1; k<numOfBins; k++){
        cout<<distances[k]<<" "<< partInBins[k]/(numMCCycles)<< endl;
       //cout<<distances[k]<<" "<< partInBins[k]/(numMCCycles*volumes[k])<< endl;
    }
    cout << " " << accept <<endl;
    delete [] volumes;
    delete [] partInBins;
    delete [] distances;*/
cout << " " << accept <<endl;
//Increase the variational parameter
        alpha += deltaAlpha;
 }
// End of loop over variational  steps
    fileblock.close();

 }

//******************************GRADIENT DESCENT BRUTE FORCE*****************************

void gradiendescent_brute( double step, int dim, int numOfPart, int numMCCycles, int numVar,  double *Etot, double *Etot2, int ind,double alpha, double deltaAlpha, int thermalization)
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
    wfold = wavefunction(rold, alpha, dim, numOfPart);
    //cout<< wfold<<endl;

 //The MC Cycle
    for (int cycl = 1; cycl <= numMCCycles+thermalization; cycl++) {

 // We set the new position
      for (int i = 0; i < numOfPart; i++) {
           for (int j = 0; j < dim; j++){
               rnew.set_Elem(i, j, rold.get_Elem(i,j)+step*(RNG(gen)-0.5));
               }
         }
    wfnew = wavefunction(rnew, alpha, dim, numOfPart);

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
            delta_e = local_energy_analytic(rold, alpha, dim, numOfPart);
            delta_wf = derWF(rold, alpha, dim, numOfPart);
    }
        else{
            delta_e = local_energy_num(rold, alpha, wfold, dim, numOfPart);
            delta_wf = derWF(rold, alpha, dim, numOfPart);
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
       if(abs(deltaAlpha*alphagrad)<10e-8){
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

//**********************************FUNCTIONS**************************************

// The wave function without interaction
double wavefunction(Matrix r, double alpha, int dim, int NumOfPart){

      double wf=0, r_sqr=0;

      for (int i = 0; i < NumOfPart; i++) {
        for (int j = 0; j < dim; j++) {
          r_sqr  += r.get_Elem(i,j)*r.get_Elem(i,j);
        }
      }
      wf = exp(-alpha*r_sqr) ;
      return wf;
    }

// The numerical calculation of local energy
double local_energy_num(Matrix r, double alpha, double wfold, int dim, int numOfPart){
      double  wfminus, wfplus, ekin=0, epot=0, eloc=0, radsq = 0;
      Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
                 radsq  += r.get_Elem(i,j)*r.get_Elem(i, j);
             }
         }
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
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


// The alytical calculation of local energy
double local_energy_analytic(Matrix r, double alpha, int dim, int numOfPart){
    double  ekin = 0,epot=0, eloc=0, radsq = 0;

    for(int i = 0; i < numOfPart; i++){
        for(int j = 0; j < dim; j++){
            radsq += r.get_Elem(i,j)*r.get_Elem(i,j);
        }
    }
    ekin = -2*alpha*alpha*radsq+alpha*dim*numOfPart;

    epot = radsq/2;
    eloc = epot+ekin;
    return eloc;
}

// The analytical quantum force calculation
void QuantumForce(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha){
    Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){

                 QForce.set_Elem(i, j, -alpha*r.get_Elem(i,j));
          }
     }
}

// The numerical quantum force calculation (NOT USED)
void QuantumForce_num(Matrix r, Matrix QForce, int dim, int numOfPart, double alpha){
    Matrix rplus(numOfPart, dim), rminus(numOfPart, dim);
    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
             }
    }
    double wfold ;
    wfold=wavefunction(r, alpha, dim, numOfPart);
    double  wfminus = 0, wfplus = 0;

    for (int i = 0; i < numOfPart; i++){
             for (int j = 0; j < dim; j++){
                 rplus.set_Elem(i,j, r.get_Elem(i,j)+h);
                 rminus.set_Elem(i, j, r.get_Elem(i,j)-h);
                 wfminus = wavefunction(rminus, alpha, dim, numOfPart);
                 wfplus  = wavefunction(rplus, alpha, dim, numOfPart);
                 QForce.set_Elem(i, j, (wfplus-wfminus)/(wfold*h));
                 rplus.set_Elem(i,j, r.get_Elem(i,j));
                 rminus.set_Elem(i, j, r.get_Elem(i,j));
          }
     }
}

void writeToFile( string x, double *E1, double *E2, int numVar, double alpha, double deltaAlpha){
    string a;

    ofstream myfile;
    myfile.open(x);
    //cin >> a;
    //myfile  << a << endl;

    for (int E = 0; E <= numVar; E++){
        myfile << alpha+deltaAlpha*E <<" "<< E1[E] << endl ;
        myfile << "                variance "<<  E2[E]-E1[E]*E1[E] << endl ;
    }
    myfile.close();
}


double derWF(Matrix r, double alpha, int dim, int NumOfPart){

      double wfder=0, r_sqr=0;

      for (int i = 0; i < NumOfPart; i++) {
        for (int j = 0; j < dim; j++) {
          r_sqr  += r.get_Elem(i,j)*r.get_Elem(i,j);
        }
      }
      wfder = -r_sqr;
      return wfder;
    }
