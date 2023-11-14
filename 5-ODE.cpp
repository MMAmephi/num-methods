#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

double solution(double x){
    return (cos(x) + 2*log(x + 1));
}

double dsolution(double x){
    return(-sin(x) + 2/(x+1));
}

double U(double x, double u, double v){
    return (2*tan(x)/(1+x)*log(1+x) - cos(x) - tan(x)/(1+x)*u - v/(1+x));
}

double V(double v){
    return v;
}

void Runge(vector<double> x, vector<double> &u, vector<double> &v, double u0, double v0, double h){
    u.push_back(u0);
    v.push_back(v0);
    double q0, q1, q2, q3, k0, k1, k2, k3;
    for(int i=0; i < x.size(); i++){
        q0 = U(x[i], u[i], v[i]);
        k0 = V(v[i]);
        q1 = U(x[i] + h/2, u[i] + k0*h/2, v[i] + q0*h/2);
        k1 = V(v[i] + q0*h/2);
        q2 = U(x[i] + h/2, u[i] + k1*h/2, v[i] + q1*h/2);
        k2 = V(v[i] + q1*h/2);
        q3 = U(x[i] + h, u[i] + k2*h, v[i] + q2*h);
        k3 = V(v[i] + q2*h);
        u.push_back(u[i] + h/6*(k0 + 2*k1 + 2*k2 + k3));
        v.push_back(v[i] + h/6*(q0 + 2*q1 + 2*q2 + q3));
    }
}

void Adams3(vector<double> x, vector<double> &u, vector<double> &v, double u0, double v0, double h){
    double u1, u2, u3, v1, v2, v3;
    
    u.push_back(u0);
    v.push_back(v0);
    v1 = v0 + h * U(x[0], u0, v0);
    u1 = u0 + h * V(v0);
    u.push_back(u1);
    v.push_back(v1);
    v2 = v1 + h * (3/2*U(x[1], u1, v1) - 1/2*U(x[0], u0, v0));
    u2 = u1 + h * (3/2*V(v1) - 1/2*V(v0));
    u.push_back(u2);
    v.push_back(v2);

    for(int i=2; i < x.size()-1; i++){
        v3 = (v2 + h * (23/12*U(x[i], u2, v2) - 16/12*U(x[i-1], u1, v1) + 5/12*U(x[i-2], u0, v0)));
        u3 = (u2 + h * (23/12*V(v2) - 16/12*V(v1) + 5/12*V(v0)));
        u.push_back(u3);
        v.push_back(v3);
        u0 = u1;
        u1 = u2;
        u2 = u3;
        v0 = v1;
        v1 = v2;
        v2 = v3;
    }
}

void Adams2(vector<double> x, vector<double> &u, vector<double> &v, double u0, double v0, double h){
    double u1, u2, v1, v2;
    
    u.push_back(u0);
    v.push_back(v0);
    v1 = v0 + h * U(x[0], u0, v0);
    u1 = u0 + h * V(v0);
    u.push_back(u1);
    v.push_back(v1);

    for(int i=1; i < x.size(); i++){
        v2 = v1 + h * (3/2*U(x[i], u1, v1) - 1/2*U(x[i-1], u0, v0));
        u2 = u1 + h * (3/2*V(v1) - 1/2*V(v0));
        u.push_back(u2);
        v.push_back(v2);
        u0 = u1;
        u1 = u2;
        v0 = v1;
        v1 = v2;
    }
}

void Adams1(vector<double> x, vector<double> &u, vector<double> &v, double u0, double v0, double h){
    double u1, v1;
    
    u.push_back(u0);
    v.push_back(v0);

    for(int i=0; i < x.size(); i++){
        v1 = v0 + h * U(x[i], u0, v0);
        u1 = u0 + h * V(v0);
        u.push_back(u1);
        v.push_back(v1);
        u0 = u1;
        v0 = v1;
    }
}


int main()
{
    double a = 0.;
    double b = 1.;
    double h = 0.05;
    double u0 = 1.;
    double du0 = 2.;

    vector<double> x_grid;
    vector<double> u1;
    vector<double> v1;

    vector<double> u2;
    vector<double> v2;

    for(int i = 0; i <= (b-a)/h; i++){
        x_grid.push_back(h * i);
    }

    Runge(x_grid, u1, v1, u0, du0, h);

    for(int i = 0; i <= (b-a)/h; i++){
        cout << setprecision(10) << "Runge: " << u1[i] << ", analytical: " << solution(x_grid[i]) << ", error: " << abs(u1[i] - solution(x_grid[i])) << "\n";
    }

    vector<double> u2h;
    vector<double> v2h;

    Runge(x_grid, u2h, v2h, u0, du0, 2*h);

    for(int i = 0; i <= (b-a)/h; i++){
        cout << setprecision(10) << "error evaluation: " << abs(u2h[i] - u1[i])/15 << "\n";
    }

    Adams3(x_grid, u2, v2, u0, du0, h);

    for(int i = 0; i <= (b-a)/h; i++){
        cout << setprecision(10) << "Adams: " << u2[i] << ", analytical: " << solution(x_grid[i]) << ", error: " << abs(u2[i] - solution(x_grid[i])) << "\n";
    }
    for(int i = 0; i <= (b-a)/h; i++){
        cout << setprecision(10) << "Adams: " << v2[i] << ", analytical: " << dsolution(x_grid[i]) << ", error: " << abs(v2[i] - dsolution(x_grid[i])) << "\n";
    }

    return 0;
}