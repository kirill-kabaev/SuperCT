#include <stdio.h>
#include <math.h>
#include <iostream>
#include <ctime>

using namespace std;

#define eps 0.000001
#define dx 0.00001
double fx(double x) {
    return sqrt(x+1)+x;
}
double dfx(double x) {
    return (fx(x + dx) - fx(x)) / dx;
}

typedef double(*function)(double x);

double solve(function fx, double x0) {
    int start = clock();
    double x1 = x0 - fx(x0) / dfx(x0);
    while (fabs(x1 - x0) > eps) {
        x0 = x1;
        x1 = x0 - fx(x0) / dfx(x0);
    }
    return x1;
}

/*
int main() {
    double start = clock();

    double x0 = 0.;
    double solving;

    solving = solve(fx, x0);


    double end = clock();
    double t = end - start; // / CLOCKS_PER_SEC;

    cout << solving << endl<< t <<endl;
    return 0;
}
*/
