#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <algorithm>

using namespace std;

int main()
{
    int N = 1+100000000;
    double *err;
    err = new double[int(floor(log10(N)))];
    for (int n = 10;n<N;n*=10){
        double h = 1./(n+1);

        double *b,*d,*x,*f,*u,*eps;
        b = new double[n];
        d = new double[n+2];
        x = new double[n+2];
        f = new double[n+2];
        u = new double[n+2];
        eps = new double[n+2];
        double a = -1;

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
        for(int i = 1; i < n+1;i++){
            eps[i] = log10(abs(d[i]-u[i])/u[i]);


        }
        err[int(floor(log10(n)))-1] = *min_element(eps,eps+n+2);
        cout<<err[int(floor(log10(n)))-1]<<endl;


    //    printf("%f\n", (stop-start)/CLOCKS_PER_SEC);
    //   ofstream outfile;
    //    outfile.open("outfile.txt");//argv[1])
    //    for (int i = 0; i<n+2;i++){
    //        outfile<< x[i]<<" "<<u[i]<<" "<< d[i]<<"\n";
    //    }
    //    outfile.close();

    }
    return 0;
}
