#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

void csv_make(double* x, double *f, int n, string name){
    string csv_file_name = name;
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "x,f\n";
    for(int i=0; i<=n; i++){
        csv_file << x[i] << "," << f[i] << "\n";
    }
    csv_file.close();
}

double dif1(double* x, double* f, int n, int i, double h){
    if(i>0 && i<n){
        return (f[i+1]-f[i-1])/2/h;
    }
    else if(i==0){
        return (4*f[1]-3*f[0]-f[2])/2/h;
    }
    else if(i==n){
        return (f[n-2]-4*f[n-1]+3*f[n])/2/h;
    }
    return 0;
}

double dif2_2(double* x, double* f, int n, int i, double h){
    if(i>0 && i<n){
        return (f[i+1]+f[i-1]-2*f[i])/h/h;
    }
    else if(i==0){
        return (2*f[0]-5*f[1]+4*f[2]-f[3])/h/h;
    }
    else if(i==n){
        return (2*f[n]-5*f[n-1]+4*f[n-2]-f[n-3])/h/h;
    }
    return 0;
}

double func(double x){
    return atan(x);
}

double real_dif(double x){
    return 1/(1+x*x);
}

double real_dif2(double x){
    return -2*x/pow((1+x*x),2);
}

int main(){
    double a = -3.;
    double b = 3.;
    int n = 50;
    double h=(b-a)/n;
    double x[n+1];
    double f[n+1];
    double r_dif[n+1];
    double r_dif2[n+1];

    for(int i=0; i<=n; i++){
        x[i]=a+h*i;
        f[i]=func(x[i]);
        r_dif[i]=real_dif(x[i]);
        r_dif2[i]=real_dif2(x[i]);
    }
    csv_make(x, r_dif, n, "r_dif.csv");
    csv_make(x, r_dif2, n, "r_dif2.csv");

    double dif1_f[n+1];
    for(int i=0; i<=n; i++){
        dif1_f[i]=dif1(x, f, n, i, h);
    }
    csv_make(x, dif1_f, n, "dif1.csv");

    double dif2_f[n+1];
    for(int i=0; i<=n; i++){
        dif2_f[i]=dif2_2(x, f, n, i, h);
    }
    csv_make(x, dif2_f, n, "dif2_2.csv");

    return 0;
}