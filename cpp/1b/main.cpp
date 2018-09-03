#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
using namespace std;

int main(int argc, char *argv[])
{
    int n = 1000000;
    double h = 1./(n+1);

    double *a,*b,*c,*d,*x,*f,*u;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    d = new double[n+2];
    x = new double[n+2];
    f = new double[n+2];
    u = new double[n+2];
    double start, stop;

    for (int i = 0; i<n+2;i++){
        x[i] = i*h;
        f[i] = 100*exp(-10*x[i]);
        u[i] = 1- (1-exp(-10))*x[i] - exp(-10*x[i]);
        d[i] = h*h*f[i];
    }


    d[0]  = 0;
    d[n+1] = 0;
    for (int i = 0; i < n; i++){
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }

    a[0]  = 0;
    c[-1] = 0;
    double cof;

    start = clock();
    for (int i=1;i<n;i++){
        cof     = a[i]/b[i-1];
        b[i]   -= cof*c[i-1];
        d[i+1] -= cof*d[i];
    }


    for (int i = n-1; i>0;i--){
        d[i] -= (c[i-1]/b[i])*d[i+1];
    }

    stop = clock();
    for (int i = 0; i<n;i++){
        d[i+1] = d[i+1]/b[i];
    }
    printf("%f\n", (stop-start)/CLOCKS_PER_SEC);
    ofstream outfile;
    outfile.open("outfile.txt");//argv[1])
    for (int i = 0; i<n+2;i++){
        outfile<< x[i]<<" "<<u[i]<<" "<< d[i]<<"\n";
    }
    outfile.close();

    return 0;
}
