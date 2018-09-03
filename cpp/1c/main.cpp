#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

using namespace std;

int main()
{
    int n = 10;
    double h = 1./(n+1);

    double *b,*d,*x,*f,*u;
    b = new double[n];
    d = new double[n+2];
    x = new double[n+2];
    f = new double[n+2];
    u = new double[n+2];
    double a = -1;
    double c = -1;

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
        b[i] = 2;
    }

    double cof;

    start = clock();
    for (int i=1;i<n;i++){
        cof     = a/b[i-1];
        b[i]   += cof;
        d[i+1] -= cof*d[i];
    }


    for (int i = n-1; i>0;i--){
        d[i] += d[i+1]/b[i];
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
