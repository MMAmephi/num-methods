#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double P(double x){
    return cos(x);
}

double Q(double x){
    return sin(x);
}

double F(double x){
    return (1-sin(x)-cos(x));
}

double U0(double x){
    return cos(x)+sin(x);
}

vector <double> progon (vector<double> X, vector<double> A, vector<double> B, vector<double> C, vector<double> d, int N){
    vector<double> u(N, 0);
    vector<double> a;
    vector<double> b;

    a.push_back(-C[0]/B[0]);
    b.push_back(d[0]/B[0]);

    for(int i = 1; i <= N; i++){
        a.push_back(-C[i]/(A[i]*a[i-1]+B[i]));
        b.push_back((d[i]-A[i]*b[i-1])/(A[i]*a[i-1]+B[i]));
    }
    /*
    for(int i = 0; i < a.size(); i++){
        cout << i << " " << a[i] << " " << b[i] << " " << d[i]  << endl;   
    }
    */
    u.push_back((d[N]-A[N]*b[N-1])/(B[N]+A[N]*a[N-1]));

    for(int i = N-1; i >= 0; i--){
        u[i] = a[i] * u[i+1] + b[i];
    }

    return u;
}

vector<double> sol1(vector<double> X, double gammaa, double gammab, double nua, double nub, double etaa, double etab, int N, double h){
    vector<double> p;
    for(int i = 0; i <= N; i++){
        p.push_back(P(X[i]));
    }

    vector<double> q;
    for(int i = 0; i <= N; i++){
        q.push_back(Q(X[i]));
    }

    vector<double> f;
    f.push_back(gammaa);
    for(int i = 1; i < N; i++){
        f.push_back(F(X[i]));
    }
    f.push_back(gammab);


    vector<double> A;
    vector<double> B;
    vector<double> C;

    A.push_back(0);
    B.push_back(nua - etaa/h);
    C.push_back(etaa/h);

    for(int i = 0; i < N-1; i++){
        A.push_back(1/h/h - p[i]/2/h);
        B.push_back(q[i]-2/h/h);
        C.push_back(1/h/h+p[i]/2/h);
    }

    A.push_back(-etab/h);
    B.push_back(nub + etab/h);
    C.push_back(0);

    vector<double> sol = progon(X, A, B, C, f, N);
    return sol;
}

vector<double> sol2(vector<double> X, double gammaa, double gammab, double nua, double nub, double etaa, double etab, int N, double h){
    vector<double> p;
    for(int i = 0; i <= N; i++){
        p.push_back(P(X[i]));
    }

    vector<double> q;
    for(int i = 0; i <= N; i++){
        q.push_back(Q(X[i]));
    }

    vector<double> f;
    f.push_back(gammaa + F(X[0])*etaa*h/(2.-p[0]*h));
    for(int i = 1; i < N; i++){
        f.push_back(F(X[i]));
    }
    f.push_back(gammab - etab*F(X[N])*h/(2.+p[N]*h));


    vector<double> A;
    vector<double> B;
    vector<double> C;

    A.push_back(0);
    B.push_back(nua - etaa/2./h*((4.-2*q[0]*h*h)/(2.-p[0]*h)));
    C.push_back(etaa/2./h*(1.+(2.+p[0]*h)/(2.-p[0]*h)));

    for(int i = 0; i < N-1; i++){
        A.push_back(1/h/h - p[i]/2./h);
        B.push_back(q[i]-2./h/h);
        C.push_back(1/h/h+p[i]/2./h);
    }

    A.push_back(-etab/2./h*(1.-(p[N]*h-2.)/(p[N]*h+2.)));
    B.push_back(nub + etab/2./h*((4.-2.*q[N]*h*h)/(2.+p[N]*h)));
    C.push_back(0);

    vector<double> sol = progon(X, A, B, C, f, N);
    return sol;
}

int main()
{
    double start = 0.; //левая граница
    double fin = 1.; //правая граница
    double h = 0.1; // шаг
    int N = int ((fin-start)/h);

    vector <double> X;
    for (int i = 0; i <= N; i++) {
        X.push_back(start + i * h);
    } 

    const double nua = 1.;
    const double etaa = -1.;
    const double nub = 1.;
    const double etab = 0.;
    const double gammaa = 0;
    const double gammab = 1.3818;


    /*
    ofstream csv_file;

    cout << "X, U" << endl;   

    const char csv_file_name1[64] = "U1.csv";
    csv_file.open(csv_file_name1);
    csv_file << "X,U\n";
    for (int i = 0; i < N; i++) {
        csv_file  << i << " " << X[i] << ", " << U[i] <<"\n";
    }
    */

    vector<double> ans1 = sol1(X, gammaa, gammab, nua, nub, etaa, etab, N, h);
    vector<double> ans2 = sol2(X, gammaa, gammab, nua, nub, etaa, etab, N, h);

    for(int i = 0; i < ans1.size(); i++){
        cout << U0(X[i]) << " " << ans1[i] << " " << abs(U0(X[i])-ans1[i]) << endl;   
    }
    cout << "second precision" << endl;
    for(int i = 0; i < ans2.size(); i++){
        cout << U0(X[i]) << " " << ans2[i] << " " << abs(U0(X[i])-ans2[i]) << endl;   
    }
    

    return 0;
}