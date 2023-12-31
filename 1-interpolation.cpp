#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

double func(double x){
    return (2+x*x+10*cos(x))/(10+x);
}

double Lagrange(double x, int n, double* xarr, double* farr){
    double sum=0;
    double mul=1;
    double *psum=&sum;
    double *pmul=&mul;
    for(int i=0; i<n; i++){
        *pmul=1;
        for(int j=0; j<n; j++){
            if(j==i){
                *pmul=mul*1;
            }
            else{
                *pmul=mul*(x-xarr[j])/(xarr[i]-xarr[j]);
            }
        }
        *psum+=farr[i]*(*pmul);
    }
    return sum;
}

double dif(double* xarr, double* farr, double k){
    double sum=0;
    double mul=1;
    double *psum=&sum;
    double *pmul=&mul;
    for(int i=0; i<k+1;i++){
        *pmul=1;
        for(int j=0; j<k+1; j++){
            if(j!=i){
                *pmul=mul*(xarr[i]-xarr[j]);
            }
        }
        *psum+=farr[i]/(mul);
    }
    return sum;
}

double Newton(double x, int n, double* xarr, double* farr){
    double s=dif(xarr, farr, n-1);
    double *ps=&s;
    for(int i=n-2; i>=0; i--){
        *ps=s*(x-xarr[i])+dif(xarr, farr, i);
    }
    return s;
}

double Chebyshev(double a, double b, int i, int n){
    return ((a+b)/2+(b-a)/2*cos((2*i+1)*M_PI/(2*n)));
}

int main(){
    double a = 0.;
    double b = 10.;
    int n = 10;
    double h=(b-a)/n;
    double x[n];
    double f[n];
    double L[n];
    double N[n];

    for(int i=0; i<n; i++){
        x[i]=a+h*i;
        f[i]=func(x[i]);
    }
    //cout << "x f(x) Lagrange Newton\n";
    for(int i=0; i<n; i++){
        L[i]=Lagrange(x[i], n, x, f);
        N[i]=Newton(x[i], n, x, f);
        //cout << x[i] << " " << f[i] << " " << L[i] << " " << N[i] <<"\n";
    }

    int n_new=100;
    double h_new=(b-a)/n_new;
    double x_new[n_new];
    double f_new[n_new];
    double L_new[n_new];
    double N_new[n_new];

    const char csv_file_name1[64] = "data_fixed.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name1);
    csv_file << "x,f,Lagrange,Newton\n";

    //cout << "Fixxed step\nx f L Newton\n";
    for(int i=0; i<n_new; i++){
        x_new[i]=a+h_new*i;
        f_new[i]=func(x_new[i]);
    }
    for(int i=0; i<n_new; i++){
        L_new[i]=Lagrange(x_new[i], n, x, f);
        N_new[i]=Newton(x_new[i], n, x, f);
        csv_file << x_new[i] << "," << f_new[i] << "," << L_new[i] << "," << N_new[i] <<"\n";
        //cout << x_new[i] << " " << f_new[i] << " " << L_new[i] << " " << N_new[i] <<"\n";
    }
    csv_file.close();

    double x_ch[n];
    double f_ch[n];
    double L_ch[n];
    double N_ch[n];

    //cout << "Chebyshev step\nx f(x) Lagrange Newton\n";
    for(int i=0; i<n; i++){
        x_ch[i]=Chebyshev(a, b, i, n);
        f_ch[i]=func(x_ch[i]);
        //cout << x_ch[i] << " " << f_ch[i] <<"\n";
    }
    for(int i=0; i<n; i++){
        L_ch[i]=Lagrange(x_ch[i], n, x_ch, f_ch);
        N_ch[i]=Newton(x_ch[i], n, x_ch, f_ch);
        //cout << x_ch[i] << " " << f_ch[i] << " " << L_ch[i] << " " << N_ch[i] <<"\n";
    }

    double L_ch_new[n_new];
    double N_ch_new[n_new];

    const char csv_file_name2[64] = "data_Chebyshev.csv";
    csv_file.open(csv_file_name2);
    csv_file << "x,f,Lagrange,Newton\n";

    //cout << "Chebyshev step\nx f(x) Lagrange Newton\n";
    for(int i=0; i<n_new; i++){
        L_ch_new[i]=Lagrange(x_new[i], n, x_ch, f_ch);
        N_ch_new[i]=Newton(x_new[i], n, x_ch, f_ch);
        csv_file << x_new[i] << "," << f_new[i] << "," << L_ch_new[i] << "," << N_ch_new[i] <<"\n";
        //cout << x_new[i] << " " << f_new[i] << " " << L_ch_new[i] << " " << N_ch_new[i] <<"\n";
    }
    csv_file.close();
    return 0;
}