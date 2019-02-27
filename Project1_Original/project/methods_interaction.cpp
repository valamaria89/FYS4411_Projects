#include "methods_interaction.h"

#define h 0.001
#define h2 1000000
#define mass 1
#define omega 1
#define a 0.0043
//#define a 0.43

using namespace std;

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
    //cout << " " << accept <<endl;
    cout << "alpha: " << alpha << " Energy " << Etot[var] << endl;

//Increase the variational parameter
        alpha += deltaAlpha;
        //cout << wfold << endl;
 }

// end of loop over variational  steps
        fileblock.close();

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
                r_ik = sqrt(r_ik);

             }
            if( r_ik <= a){
                p = 0;
                //return 0;
            }
            else arg += log(1-a/(r_ik));
                //cout << arg << endl;
        }

      }



      wf = exp(-alpha*r_sqr+0*arg);

      return wf*p;
    }

//The numerical calculation of local energy
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
    double  ekin = 0,epot=0, eloc=0, radsq = 0;

    for(int i = 0; i < numOfPart; i++){
        for(int j = 0; j < dim; j++){
            if (j == 2){
                radsq  += beta*beta*r.get_Elem(i,j)*r.get_Elem(i,j);
        }
            else radsq += r.get_Elem(i,j)*r.get_Elem(i,j);


        }
    }
    ekin = -2*alpha*alpha*radsq+alpha*dim*numOfPart;

    epot = radsq/2;
    eloc = epot+ekin;
    return eloc;
}
