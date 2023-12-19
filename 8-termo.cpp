#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double F(double t, double x){
    return (2.-8./3.*x*x)*exp(t-x*x);
}

double Phi(double x){
    return 2./3.*exp(-x*x);
}

double Gamma0(double t){
    return 2./3.*exp(t);
}

double Gammal(double t){
    return -2./3.*exp(t-1);
}

double U0(double t, double x){
    return 2./3.*exp(t-x*x);
}

vector <double> progon (vector<double> X, vector<double> A, vector<double> B, vector<double> C, vector<double> d, int N){
    vector<double> u(N-1);
    vector<double> a;
    vector<double> b;

    a.push_back(-C[0]/B[0]);
    b.push_back(d[0]/B[0]);

    for(int i = 1; i < N; i++){
        a.push_back(-C[i]/(A[i]*a[i-1]+B[i]));
        b.push_back((d[i]-A[i]*b[i-1])/(A[i]*a[i-1]+B[i]));
    }
    /*
    for(int i = 0; i < a.size(); i++){
        cout << i << " " << a[i] << " " << b[i] << " " << d[i]  << endl;   
    }
    */
    u.push_back((d[N-1]-A[N-1]*b[N-2])/(B[N-1]+A[N-1]*a[N-2]));

    for(int i = N-2; i >= 0; i--){
        u[i] = a[i] * u[i+1] + b[i];
    }

    return u;
}

vector<vector<double>> Sol1(vector<double> x, vector<double> t, vector<vector<double>> f, double tau, double h, const double a, const double alpha0, const double alphal, const double beta0, const double betal, const double sigma){
    const int N_x = x.size();
    const int N_t = t.size();
    vector<vector<double>> u(N_t, vector<double>(N_x));
    for (int i = 0; i < N_x; i++){
        u[0][i] = Phi(x[i]);    
    }
    
    vector<double> A(N_x);
    vector<double> B(N_x);
    vector<double> C(N_x);
    const double AA = a*a*(1.-sigma)*tau/h/h;
    const double BB = -2.*a*a*(1.-sigma)*tau/h/h-1.;
    const double CC = a*a*(1.-sigma)*tau/h/h;
    A[0] = 0.;
    B[0] = alpha0*h - beta0;
    C[0] = beta0;
    for(int i = 1; i < N_x-1; i++){
        A[i] = AA;
        B[i] = BB;
        C[i] = CC;
    }
    A[N_x-1] = -betal;
    B[N_x-1] = alphal*h + betal;
    C[N_x-1] = 0.;
    
    vector<double> d(N_x);
    
    for(int n = 1; n < N_t; n++){
        d[0] = Gamma0(t[n])*h;
        for(int i = 1; i < N_x-1; i++){
            d[i] = -F(t[n-1], x[i])*tau - u[n-1][i] - a*a*sigma*tau/h/h*(u[n-1][i+1] - 2.*u[n-1][i] + u[n-1][i-1]);
        }
        d[N_x-1] = Gammal(t[n])*h;
        u[n]=progon(x, A, B, C, d, N_x);
    }
    
    return u;
}

vector<vector<double>> Sol2(vector<double> x, vector<double> t, vector<vector<double>> f, double tau, double h, const double a, const double alpha0, const double alphal, const double beta0, const double betal, const double sigma){
    const int N_x = x.size();
    const int N_t = t.size();
    vector<vector<double>> u(N_t, vector<double>(N_x));
    for (int i = 0; i < N_x; i++){
        u[0][i] = Phi(x[i]);    
    }
    
    vector<double> A(N_x);
    vector<double> B(N_x);
    vector<double> C(N_x);
    const double AA = a*a*(1.-sigma)*tau/h/h;
    const double BB = -2.*a*a*(1.-sigma)*tau/h/h-1.;
    const double CC = a*a*(1.-sigma)*tau/h/h;
    A[0] = 0.;
    B[0] = 2.*alpha0*h - 2.*beta0 - beta0*h*h/tau/a/a;
    C[0] = 2.*beta0;
    for(int i = 1; i < N_x-1; i++){
        A[i] = AA;
        B[i] = BB;
        C[i] = CC;
    }
    A[N_x-1] = -2.*betal;
    B[N_x-1] = 2.*alphal*h + 2.*betal + betal*h*h/tau/a/a;
    C[N_x-1] = 0.;
    
    vector<double> d(N_x);
    
    for(int n = 1; n < N_t; n++){
        d[0] = Gamma0(t[n])*h*2. - beta0*F(t[n], x[0])*h*h/a/a - u[n-1][0]*h*h*beta0/tau/a/a;
        for(int i = 1; i < N_x-1; i++){
            d[i] = -F(t[n-1], x[i])*tau - u[n-1][i] - a*a*sigma*tau/h/h*(u[n-1][i+1] - 2.*u[n-1][i] + u[n-1][i-1]);
        }
        d[N_x-1] = 2.*Gammal(t[n])*h + h*h*betal*u[n-1][N_x-1]/a/a/tau + betal*h*h*F(t[n], x[N_x-1])/a/a;
        u[n]=progon(x, A, B, C, d, N_x);
    }
    
    return u;
}


int main()
{
    double start_x = 0.;
    double end = 1.;
    double h = 0.2;
    int N_x = int ((end-start_x)/h);

    vector <double> x;
    for (int i = 0; i <= N_x; i++) {
        x.push_back(start_x + i * h);
    }


    double start_t = 0.;
    double T = 1.;
    double tau = 0.2;
    int N_t = int ((T-start_t)/tau);

    vector <double> t;
    for (int i = 0; i <= N_t; i++) {
        t.push_back(start_t + i * tau);
    }

    const double a = 1.;
    const double alpha0 = 1.;
    const double alphal = 1.;
    const double beta0 = -1.;
    const double betal = 1.;
    const double sigma = 0.5;

    vector<vector<double>> f(N_t+1, vector<double>(N_x+1));

    for(int i = 0; i <= N_t; i++){
        for(int j = 0; j <= N_x; j++){
            f[i][j] = F(t[i], x[j]);
        }
    }

    cout << "first precision" << endl;
    vector<vector<double>> sol1 = Sol1(x, t, f, tau, h, a, alpha0, alphal, beta0, betal, sigma);
    for(int i = 0; i < N_t; i++){
        for(int j = 0; j < N_x; j++){
            cout << sol1[i][j] << "|" << U0(t[i], x[j]) << "|" << abs(sol1[i][j] - U0(t[i], x[j])) << ", ";
        }
        cout << endl;
    }
    
    cout << "second precision" << endl;
    vector<vector<double>> sol2 = Sol2(x, t, f, tau, h, a, alpha0, alphal, beta0, betal, sigma);
    for(int i = 0; i <= N_t; i++){
        for(int j = 0; j <= N_x; j++){
            cout << sol2[i][j] << "|" << U0(t[i], x[j]) << "|" << abs(sol2[i][j] - U0(t[i], x[j])) << ", ";
        }
        cout << endl;
    }
    

    return 0;
}