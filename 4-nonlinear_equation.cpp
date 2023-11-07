#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;
 
double func(double x){
    return exp(sin(x+2)) - 1;
}

double dfunc(double x){
    return exp(sin(x+2))*cos(x+2);
}

double ddfunc(double x){
    return -sin(x+1);
}

bool sign(double x){
    if(x > 0) return 1;
    else return 0;
}

void boundary(double a, double b, int n, vector<double> &boundaries){
    double temp;
    double prev = a;
    bool sign_f;
    bool sign_df;
    for(int i = 1; i <= n; i++){
        temp = a + i*(b - a)/n;
        if((sign(dfunc(temp)) == sign(dfunc(prev))) && (sign(func(temp)) != sign(func(prev)))){
            boundaries.push_back(prev);
            boundaries.push_back(temp);
        }
        prev = temp;
    }
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
    double b = 10;
    const double eps = pow(10, -3);
    int n = 10;
    vector<double> boundaries;

    boundary(a, b, n, boundaries);
    for(int i = 0; i < boundaries.size(); i=i+2){
        cout << i/2+1 << " solution" << "\n";
        bisection(boundaries[i], boundaries[i+1], eps); 
        Newton(boundaries[i], boundaries[i+1], eps);
        mod_Newton(boundaries[i], boundaries[i+1], eps);
        chord(boundaries[i], boundaries[i+1], eps);   
    }
    /*bisection(a, b, eps);
    Newton(a, b, eps);
    mod_Newton(a, b, eps);
    chord(a, b, eps);*/
    /*for(int i=0; i < boundaries.size(); i++){
        cout << boundaries[i] <<"\n";
    };*/
    return 0;
}