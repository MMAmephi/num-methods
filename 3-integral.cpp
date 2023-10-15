#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

double func(double x){
    return 1/pow(cosh(x), 2);
}

double true_integral(double x){
    return tanh(x);
}

double left_rectangle(double a, double b, int n){
    double I = 0.;
    double h = (b - a) / (n - 1);
    for (int i = 0; i < n; i++){
        I += func(a + h * i);
    }
    return I * h;
}

double right_rectangle(double a, double b, int n){
    double I = 0.;
    double h = (b - a) / (n - 1);
    for (int i = 1; i <= n; i++){
        I += func(a + h * i);
    }
    return I * h;
}

double rectangle(double a, double b, int n){
    double I = 0.;
    double h = (b - a) / (n - 1);
    for (int i = 0; i < n; i++){
        I += func(a + h * (i + 0.5));
    }
    return I * h;
}

double trapezoid(double a, double b, int n){
    double I = 0.;
    double h = (b - a) / (n - 1);
    for (int i = 0; i < n; i++){
        I += func(a + h * (i + 1)) + func(a + h * i);
    }
    return I * h / 2;
}

double Simpson(double a, double b, int n){
    double I = 0.;
    double h = (b - a) / (n - 1);
    for (int i = 1; i < n; i++){
        I += func(a + h * (i - 1)) + 4 * func(a + h * i) + func(a + h * (i + 1));
    }
    return I * h / 6;
}

int main()
{
    int n = 100000; // n - кол-во узлов
    double a = -1., b = 3.;
    double h = (b-a)/(n-1);

    cout << "true value: " << true_integral(b) - true_integral(a) << '\n';
    cout << "left rectangle: " << left_rectangle(a, b, n) << '\n';
    cout << "right rectangle: " << right_rectangle(a, b, n) << '\n';
    cout << "mid rectangle: " << rectangle(a, b, n) << '\n';
    cout << "trapezoid: " << trapezoid(a, b, n) << '\n';
    cout << "Simpson: " << Simpson(a, b, n) << '\n';

    cout << "error left rectangle = " << fabs(left_rectangle(a, b, n) - true_integral(b) + true_integral(a)) << '\n';
    cout << "error right rectangle = " << fabs(right_rectangle(a, b, n) - true_integral(b) + true_integral(a)) << '\n';
    cout << "error mid rectangle = " << fabs(rectangle(a, b, n) - true_integral(b) + true_integral(a)) << '\n';
    cout << "error trapezoid = " << fabs(trapezoid(a, b, n) - true_integral(b) + true_integral(a)) << '\n';
    cout << "error Simpson   = " << fabs(Simpson(a, b, n) - true_integral(b) + true_integral(a)) << '\n';

    return 0;
}