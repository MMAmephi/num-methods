#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

double p (double x) {return (cos (x));} //коэфы перед U'

double q (double x) {return (sin (x));} //коэфы перед U

double U0 (double x) {return sin(x) + cos(x);} //решение

double W (double x, double u, double v) {return (1 - cos(x) - sin(x) - p(x)*v - q(x)*u);}  //Функция для U"

void RungeKutta(const vector<double> x, double h, vector<double> &u, vector<double> &v, double u0, double v0) // решение задачи Коши
{
    u.push_back(u0);
    v.push_back(v0);

    double k0, k1, k2, k3;
    double m0, m1, m2, m3;

    for(int i = 0; i < x.size(); ++i)
    {
        m0 = W(x.at(i), u.at(i), v.at(i));
        k0 = v.at(i);

        m1 = W(x.at(i) + h/2., u.at(i) + k0 * h/2., v.at(i) + m0 * h/2.);
        k1 = v.at(i) + m0 * h/2.;

        m2 = W(x.at(i) + h/2., u.at(i) + k1 * h/2., v.at(i) + m1 * h/2.);
        k2 = v.at(i) + m1 * h/2.;

        m3 = W(x.at(i) + h, u.at(i) + k2 * h, v.at(i) + m2 * h);
        k3 = v.at(i) + m2 * h;

        u.push_back(u.at(i) + h/6 * (k0 + 2*k1 + 2*k2 + k3));
        v.push_back(v.at(i) + h/6 * (m0 + 2*m1 + 2*m2 + m3));
    }
}

double podgon(const double (&A)[2], const double (&B)[2], double realB) { // поправка на u(0)
    return A[1] - ((B[1] - realB) * (A[1] - A[0]))/((B[1] - realB) - (B[0] - realB)); }

vector <double> strelba (vector <double> X, double h, double gammaa, double nua, double etaa, double right, int n, double eps) {  // собственно сама стрельба
    double A[2]{21635456.2465, 3241654465.168498},  B[2]{};  //это Ak, Ak+1, Bk и Bk+1
    vector<double> U;
    vector<double> V;
    while(true)
    {
        double u0 =  A[1]; //поправляем U(0)
        double du0 = (gammaa - nua*A[1])/etaa; //поправляем U'(0)

        RungeKnutta(X, h, U, V, u0, du0); //Решаем задачу Коши с поправленными нач данными

        B[1] = U.at(n); 

        if (abs(B[1] -  right) < eps) break; // проверка точности

        double tmp = podgon(A, B, right); //Ak+1

        A[0] = A[1];
        A[1] = tmp;
        B[0] = B[1];

        U.clear();
        V.clear();
    }
    return U;
}

int main()
{
    double a = 0.;
    double b = 1.;
    double h = 0.05;
    int n = (b - a) / h + 1;

    const double nua = 1.;
    const double etaa = -1.;
    const double gammaa = 0;
    const double nub = 1.;
    const double etab = 0.;
    const double gammab = 1.3818;

    double right = gammab/nub;

    vector<double> X;
    vector<double> Y;

    for(int i = 0; i <= n; i++) { 
        X.push_back(i * h);
        Y.push_back(U0 (X[i]));
    }
    double eps = h;

    vector <double> U = strelba (X, h, gammaa, nua, etaa, right, n, eps);
    
    ofstream csv_file;

    const char csv_file_name1[64] = "Us1.csv";
    csv_file.open(csv_file_name1);
    csv_file << "X,U,Y\n";
    for (int i = 0; i < n; i++) {
        csv_file << X[i] << ", " << U[i] << ", " << Y[i] <<endl;
    }
    csv_file.close();

    vector <long double> H; // массив шагов
    vector <long double> O; // массив ошибок
    vector <long double> Hl; // массив шагов отлагорифмированных
    vector <long double> Ol; // массив ошибок отлагорифмированнаых

    for (int i = 1; i <= 5; i++) { //сюда вписать до какого шага мы идём (10^-наше число)
        long double h1 = pow(10, -i);
        double eps1 = h1;
        H.push_back(h1);
        int N1 = int (double(b-a)/pow(10, -i));

        vector <double> X1;
        vector <double> Y1;
        for (int i = 0; i <= N1; i++) {
            X1.push_back(a + i * h1);
        }
        for (int i = 0; i <= N1; i++) {
            Y1.push_back(U0 (X1[i]));            
        }

        vector <double> U1 = strelba (X1, h1, gammaa, nua, etaa, right, N1, eps1);;

        long double maxDiff = 0;
        long double diff = 0;

        for (int i = 0; i <= N1; i++) {
            diff = abs(Y1[i] - U1[i]);
            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }
        O.push_back(maxDiff);
        Hl.push_back(log(h1));
        Ol.push_back(log(maxDiff));
    }

    const char csv_file_name2[64] = "Us2.csv";
    csv_file.open(csv_file_name2);
    csv_file << "H,O,Hl,Ol,\n";
    for (int i = 0; i < H.size(); i++) {
        csv_file << H[i] << ", " << O[i] << ", " << Hl[i] << ", " << Ol[i] << endl;
        cout << H[i] << ", " << O[i] << ", " << Hl[i] << ", " << Ol[i] << endl;
    }
    csv_file.close();
    

    return 0;
}