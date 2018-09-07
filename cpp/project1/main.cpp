#include <iostream>
#include <cmath>
#include <armadillo>
#include <algorithm>

using namespace std;
using namespace arma;

vec f(vec x);
vec analytical(vec x);
vec general_elim_solve(vec a, vec b,vec c, vec b_vec,int n);
vec specialized_elim_solve(vec b, vec b_vec,int n, int s);
void errors(int max_n);
void compare(int max_n);

int main(){
    int n = 1000;
    vec x,a,b,c,d;
    double h = 1./(n+1);
    x = linspace(0,1,n+2);
    a.ones(n); a *= -1;   // diagonals in the dridiagonal matrix A. a=lower diagonal, c=upper b=middle
    b.ones(n); b *= 2;
    c.ones(n); c *= -1;
    d = h*h*f(x);
    d[0] = 0; d[n+1] = 0; a[0] = 0; c[n-1] = 0; // randbetingelser

    vec solution = general_elim_solve(a,b,c,d, n);     //1b
    vec fast_solution = specialized_elim_solve(b,d,n,1); //c
    cout<<endl;
    //errors(int(1e7)); //d output to file
    compare(1000);      //e output to file


    //cout << solution << endl;
    //cout << fast_solution <<endl;
    return 0;
}

vec f(vec x){
    vec f = 100*exp(-10*x);
    return f;
}

vec analytical(vec x){
    vec u = 1- (1-exp(-10))*x - exp(-10*x);
    return u;
}

vec general_elim_solve(vec a, vec b,vec c, vec b_vec, int n){
    double cof,start,stop;

    start = clock();
    for (int i=1;i<n;i++){
        cof     = a[i]/b[i-1];
        b[i]   -= cof*c[i-1];
        b_vec[i+1] -= cof*b_vec[i];
    }


    for (int i = n-1; i>0;i--){
        b_vec[i] -= (c[i-1]/b[i])*b_vec[i+1];
    }

    stop = clock();
    for (int i = 0; i<n;i++){
        b_vec[i+1] = b_vec[i+1]/b[i];
    }
    printf("General algoritm %f secs with n= %i \n", (stop-start)/CLOCKS_PER_SEC,n);
    return b_vec;
}

vec specialized_elim_solve(vec b,vec b_vec,int n,int s){
    double cof,start,stop;
    start = clock();
    for (int i=1;i<n;i++){
        cof     = -1/b[i-1];
        b[i]   += cof;
        b_vec[i+1] -= cof*b_vec[i];
    }


    for (int i = n-1; i>0;i--){
        b_vec[i] += b_vec[i+1]/b[i];
    }

    stop = clock();
    for (int i = 0; i<n;i++){
        b_vec[i+1] = b_vec[i+1]/b[i];
    }
    if(s==1){
        printf("Specialized algoritm: %f secs with n= %i \n", (stop-start)/CLOCKS_PER_SEC,n);
    }
    return b_vec;
}

void errors(int max_n){
    vec err;
    err.zeros(int(log10(max_n)));
    for (int n = 10;n<=max_n;n = n*10){
        vec x,b,d,exact,approx,eps;
        x = linspace(0,1,n+2);
        b.ones(n); b *= 2;
        d = (1./(n+1))*(1./(n+1))*f(x);
        exact  = analytical(x);
        approx = specialized_elim_solve(b,d,n,0);
        eps.zeros(n+2);
        for(int i = 0; i < n+2;i++){
            eps[i] = log10(abs(exact[i]-approx[i])/exact[i]);

        }
        eps[0] = -200; // unngår randbetingelsene
        eps[n+1] = -200;
        err[int(log10(n))-1] = eps.max();

    }
    //cout << err << endl;
}

void compare(int max_n){
    //Comparing LU decomposition with the specialiced gaussian
    for(int n = 10; n<=max_n; n = n*10){

        //LU
        mat A(n,n);
        mat L,U;
        A.zeros();
        A(0,0) = 2;
        for(int i = 1;i<n;i++){
            A(i,i)   = 2;
            A(i,i-1) = -1;
            A(i-1,i) = -1;

            }
        double start,stop;
        start = clock();
        lu(L,U,A);
        stop = clock();
        vec arma_solution;
        printf("Armadillos LU decomp %f secs with n= %i \n", (stop-start)/CLOCKS_PER_SEC,n);

        //Gauss
        vec x,b,d,gaus_solution;
        x = linspace(0,1,n+2);
        b.ones(n); b *= 2;
        d = (1./(n+1))*(1./(n+1))*f(x);
        gaus_solution = specialized_elim_solve(b,d,n,1);
        cout<<endl;


    }
    //cout << L << endl << endl << U << endl;
}


