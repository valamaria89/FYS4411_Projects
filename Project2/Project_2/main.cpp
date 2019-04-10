#include <iostream>
#include "system.h"
#include "matrix.h"

using namespace std;

int main()
{
    int N = 2;
    int M = 2;
    double sigma = 3;
    double omega = 5;

System h(N,M,sigma,omega);

Matrix x(M,N);
Matrix a(M,N);
Matrix b(M,N);
Matrix W(M,N);

for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++){
        x.set_Elem(i, j, 0);
        a.set_Elem(i, j, 0);
        b.set_Elem(i, j, 1);
        W.set_Elem(i, j, 0);
    }
};
for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++){
       cout<<""<< x.get_Elem(i, j);
    };
    cout<<endl;
};

cout<<"This is our result "<<h.waveFunction(x,a,W,sigma)<<endl;
cout<<"This is also our result "<<h.localEnergy_IMS(x,a,W,sigma)<<endl;
cout<<"And this is u "<<h.u(x,a,b,W,sigma,1)<<endl;
}
