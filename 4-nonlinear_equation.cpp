#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
 
double func(double x){
    return sin(x+1);
}

double dfunc(double x){
    return cos(x+1);
}

double ddfunc(double x){
    return -sin(x+1);
}

double bisection(double a, double b, double eps){
    double x0 = (a + b)/2;
    int k = 0;
    while((b - a >= eps) && (abs(func(x0)) >= eps)){
        k++;
        x0 = (a + b) / 2;
        if(func(b)*func(x0) < 0)
            a = x0;
        else
            b = x0;
    }
    cout << "bisection: " << setprecision(10) << x0 << "\nnumber of iterrations - " << k <<"\n";
    return x0;
}

double Newton(double a, double b, double eps){
    double x0 = b;
    double x1 = x0 - func(x0)/dfunc(x0);
    int k = 1;
    while((abs(x1 - x0) >= eps) && (abs(func(x1)) >= eps)){
        k++;
        x0 = x1;
        x1 = x0 - func(x0) / dfunc(x0);
    }
    cout << "Newton: " << setprecision(10) << x1 << "\nnumber of iterrations - " << k <<"\n";
    return x1;    
}

double mod_Newton(double a, double b, double eps){
    int k = 1;
    double c, x0;
    if(func(a)*ddfunc(a) > 0){
        x0 = a;
    }
    else{
        x0 = b;
    }
    double temp = dfunc(x0);
    double x1 = x0 - func(x0)/temp;
    while((abs(x1 - x0) >= eps) && (abs(func(x1)) >= eps)){
        k++;
        x0 = x1;
        x1 = x0 - func(x0) / temp;
    }
    cout << "modified Newton: " << setprecision(10) << x1 << "\nnumber of iterrations - " << k <<"\n";
    return x1;    
}

double chord(double a, double b, double eps){
    int k = 1;
    double c, x0;
    if(func(a)*ddfunc(a) > 0){
        c = a;
        x0 = b;
    }
    else{
        c = b;
        x0 = a;
    }
    double x1 = x0 - func(x0)*(x0 - c) / (func(x0) - func(c));
    while((abs(x1 - x0) >= eps) && (abs(func(x1)) >= eps)){
        k++;
        x0 = x1;
        x1 = x0 - func(x0)*(x0 - c) / (func(x0) - func(c));
    }
    cout << setprecision(10) << "chord: " << x1 << "\nnumber of iterrations - " << k <<"\n";
    return x1;
}

int main()
{
    double a = 0;
    double b = 3;
    const double eps = pow(10, -3);
    bisection(a, b, eps);
    Newton(a, b, eps);
    mod_Newton(a, b, eps);
    chord(a, b, eps);
    return 0;
}