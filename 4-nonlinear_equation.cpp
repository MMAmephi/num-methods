#include <iostream>
#include <cmath>

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
    double ans = (a + b)/2;
    while((b - a >= eps) && (abs(func(ans)) >= eps)){
        ans = (a + b) / 2;
        if(func(b)*func(ans) < 0)
            a = ans;
        else
            b = ans;
    }
    return ans;
}

double Newton(double a, double b, double eps){
    double x0 = b;
    double x1 = x0 - func(x0)/dfunc(x0);
    while((abs(x1 - x0) >= eps) && (abs(func(x1)) >= eps)){
        x0 = x1;
        x1 = x0 - func(x0) / dfunc(x0);
    }
    return x1;    
}

double mod_Newton(double a, double b, double eps){
    double x0 = b;
    double temp = dfunc(x0);
    double x1 = x0 - func(x0)/temp;
    while((abs(x1 - x0) >= eps) && (abs(func(x1)) >= eps)){
        x0 = x1;
        x1 = x0 - func(x0) / temp;
    }
    return x1;    
}

double chord(double a, double b, double eps){
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
        x0 = x1;
        x1 = x0 - func(x0)*(x0 - c) / (func(x0) - func(c));
    }
    return x1;
}

int main()
{
    double a = 0;
    double b = 3;
    const double eps = pow(10, -10);
    cout << bisection(a, b, eps) << "\n";
    cout << Newton(a, b, eps) << "\n";
    cout << mod_Newton(a, b, eps) << "\n";
    cout << chord(a, b, eps) << "\n";
    return 0;
}