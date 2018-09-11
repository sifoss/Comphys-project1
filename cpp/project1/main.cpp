#include <iostream>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include <cstdlib>

using namespace std;
using namespace arma;

vec f(vec x);
vec analytical(vec x);
vec general_elim_solve(vec a, vec b,vec c, vec b_vec,int n);
vec specialized_elim_solve(vec b, vec b_vec,int n, int s);
void errors(int max_n);
void compare(int max_n);


int main(int argc, char* argv[]){

    int n = 1e+6;
    vec x,a,b,c,d,exact;
    double h = 1./(n+1);
    x = linspace(0,1,n+2);
    a.ones(n); a *= -1;   // diagonals in the tridiagonal matrix A. a=lower diagonal, c=upper b=middle
    b.ones(n); b *= 2;
    c.ones(n); c *= -1;
    d = h*h*f(x);
    exact = analytical(x);
    d[0] = 0; d[n+1] = 0; a[0] = 0; c[n-1] = 0; // randbetingelser
    vec solution =  general_elim_solve(a,b,c,d, n);
    vec fast_solution = specialized_elim_solve(b,d,n,1); //c

    /*
    ostringstream iss;
    iss << "u_exact_approx_" << n << ".txt";
    string filename = iss.str();

    ofstream outfile;     //write to file for exercise b
    outfile.open(filename);
    for (int i = 0; i<n+2;i++){
        outfile<< x[i]<<" "<<exact[i]<<" "<< solution[i]<<"\n";
    }
    outfile.close();
    */

    cout<<endl;
    //errors(n); //d output to file
    //compare(n);      //e output to file

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
    double cof,start,stop, time;
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
    printf("General algoritm: %f secs with n= %i \n", (stop-start)/CLOCKS_PER_SEC,n);
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
    vec err,n_vec;
    int m = int(log10(max_n));
    err.zeros(m);
    n_vec.zeros(m);
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
        n_vec[int(log10(n))-1] = n;
    }
    ofstream outfile;
    outfile.open("errors.txt");//argv[1])
    for (int i = 0; i<m;i++){
        outfile<< n_vec[i]<<" "<<err[i]<<" "<<"\n";
    }
    outfile.close();

}

void compare(int max_n){
    //Comparing LU decomposition with the specialized gaussian
    ofstream outfile;
    outfile.open("arma_errors.txt");
    for(int i = 10; i<=max_n; i = i*10){
        vec x,b,d,gaus_solution, analytic_sol, errors;
        x = linspace(0,1,i+2);
        d = (1./(i+1))*(1./(i+1))*f(x);
        double h = 1./(i+1);

        // Initializing arrays and solving with LU
        vec d_for_matrix(i);
        for(int j = 0; j<i;j++){
            d_for_matrix[j] = d[j+1];
        }

        mat A(i,i);
        mat L,U;
        A.zeros();
        A(0,0) = 2;
        for(int j = 1;j<i;j++){
            A(j,j)   = 2;
            A(j,j-1) = -1;
            A(j-1,j) = -1;

            }
        double start,stop, time;
        for(int j = 0; j <50; j++){
            start = clock(); //timing of LU algo
            lu(L,U,A);
            stop = clock();
            time += (stop -start)/CLOCKS_PER_SEC;
        }


        cout<<"time :"<< time/50 <<endl;
        vec y = solve(L,d_for_matrix);
        cout << i  << endl;
        vec arma_solution = solve(U,y);

        vec x2 = linspace(h, 1-h, i);
        analytic_sol = analytical(x2);
        errors = log10(abs((arma_solution - analytic_sol)/analytic_sol));  // calculating errors in LU method

        outfile<<max(errors) <<endl;
        cout<<max(errors) <<endl;

        //Gauss
        b.ones(i); b *= 2;
        gaus_solution = specialized_elim_solve(b,d,i,1);
        cout<<endl;

    outfile.close();

    }
    //cout << L << endl << endl << U << endl;
}


