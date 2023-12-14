#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double F(double t, double x){
    return (-2 + (2*x*x - t*t)*sinh(x*t))/2.;
}

double Phi(double x){
    return x*x;
}

double Psi(double x){
    return x;
}

double Gamma0(double t){
    return t;
}

double Gammal(double t){
    return 3. + sinh(t) + t*cosh(t);
}

double U0(double t, double x){
    return x*x + sinh(x*t);
}

vector<vector<double>> Sol1(vector<double> x, vector<double> t, vector<vector<double>> f, double tau, double h, const double a, const double alpha0, const double alphal, const double beta0, const double betal){
    int Nx = x.size();
    int Nt = t.size();
    vector<vector<double>> u(Nt+1, vector<double>(Nx+1));

    for(int i = 0; i < Nx; i++){
        u[0][i] = Phi(x[i]);
        u[1][i] = u[0][i] + tau*Psi(x[i]);
    }

    for( int i = 1; i < Nt-1; i++){
        for(int j = 1; j < Nx-1; j++){
            u[i+1][j] = 2.*u[i][j] - u[i-1][j] + pow((a*tau/h), 2)*(u[i][j+1] - 2.*u[i][j] + u[i][j-1]) + tau*tau*f[i][j];
        }
        u[i+1][0] = 1./(alpha0*h - beta0)*(h*Gamma0(t[i+1]) - beta0*u[i+1][1]);
        u[i+1][Nx-1] = 1./(alphal*h + betal)*(h*Gammal(t[i+1]) + betal*u[i+1][Nx-2]);
    }

    return u;
}

vector<vector<double>> Sol2(vector<double> x, vector<double> t, vector<vector<double>> f, double tau, double h, const double a, const double alpha0, const double alphal, const double beta0, const double betal){
    int Nx = x.size();
    int Nt = t.size();
    vector<vector<double>> u(Nt, vector<double>(Nx));

    for(int i = 1; i < Nx-1; i++){
        u[0][i] = Phi(x[i]);
        u[1][i] = 2.*u[0][i] - u[1][i] + pow((a*tau/h), 2)*(u[0][i+1] - 2.*u[0][i] + u[0][i-1]) + tau*tau*f[0][i] + 2.*tau*Psi(x[i]);
    }

    u[0][0] = Phi(x[0]);
    u[0][Nx-1] = Phi(x[Nx-1]);
    u[1][0] = 1./(2.*alpha0*h - 3.*beta0)*(2.*h*Gamma0(t[1]) - 4.*beta0*u[1][1] + beta0*u[1][2]);
    u[1][Nx-1] = 1./(2.*alphal*h + 3.*betal)*(2.*h*Gammal(t[1]) + 4.*betal*u[1][Nx-2] - betal*u[1][Nx-3]);

    for( int i = 1; i < Nt-1; i++){
        for(int j = 1; j < Nx-1; j++){
            u[i+1][j] = 2.*u[i][j] - u[i-1][j] + pow((a*tau/h), 2)*(u[i][j+1] - 2.*u[i][j] + u[i][j-1]) + tau*tau*f[i][j];
        }
        u[i+1][0] = 1./(2.*alpha0*h - 3.*beta0)*(2.*h*Gamma0(t[i+1]) - 4.*beta0*u[i+1][1] + beta0*u[i+1][2]);
        u[i+1][Nx-1] =  1./(2.*alphal*h + 3.*betal)*(2.*h*Gammal(t[i+1]) + 4.*betal*u[i+1][Nx-2] - betal*u[i+1][Nx-3]);
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

    const double a = 1/sqrt(2);
    const double alpha0 = 0.;
    const double alphal = 1.;
    const double beta0 = 1.;
    const double betal = 1.;

    vector<vector<double>> f(N_t, vector<double>(N_x));

    for(int i = 0; i < N_t; i++){
        for(int j = 0; j < N_x; j++){
            f[i][j] = F(t[i], x[j]);
        }
    }

    cout << "first precision" << endl;
    vector<vector<double>> sol1 = Sol1(x, t, f, tau, h, a, alpha0, alphal, beta0, betal);
    for(int i = 0; i < N_t; i++){
        for(int j = 0; j < N_x; j++){
            cout << sol1[i][j] << "|";
            cout << abs(sol1[i][j] - U0(t[i], x[j])) << ", ";
        }
        cout << endl;
    }
    
    cout << "second precision" << endl;
    vector<vector<double>> sol2 = Sol2(x, t, f, tau, h, a, alpha0, alphal, beta0, betal);
    for(int i = 0; i <= N_t; i++){
        for(int j = 0; j <= N_x; j++){
            //cout << sol2[i][j] << ", ";
            cout << abs(sol2[i][j] - U0(t[i], x[j])) << ", ";
        }
        cout << endl;
    }
    

    return 0;
}