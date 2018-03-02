#include <iostream>
#include "Forsythe.h"
#include <conio.h>
#include <math.h>
#include <iomanip>
#include <fstream>

void f(double x, double y[], double dydx[]) {

    dydx[0] = -45*y[0] - 60*y[1] + sin(x+1);
    dydx[1] = 70*y[0] -110*y[1] + cos(1-x) + x + 1;
}

void rk2(void(*f)(double x, double y[], double dydx[]),
         double t, double tout, double Zn[],  double h);

int main() {

    const int x1 = 5;
    const int x2 = -1;

    double x[2] = {x1, x2};

    int left_border = 0;
    int right_border = 1;


    //rkf45
    const int neqn = 2;
    unsigned char work[6*(neqn*sizeof(Float)) + sizeof(struct rkf_inside)];

    rkf ARG;
    ARG.f = f;
    ARG.neqn = neqn;
    ARG.re = 0.0001;
    ARG.ae = 0.0001;
    ARG.work = work;
    ARG.flag = 1;
    ARG.Y = x;
    ARG.t = 0;

    double step=0.05;

    std::ofstream out("out.txt");

    out <<"RKF45\n\n";
    out <<"   t |     x[0]  |     x[1]  |  Flag\n";
    out <<"_______________________________________\n";
    out.setf(std::ios::fixed);

    for(double h=(left_border + step); h<(right_border + step); h+=0.05) {

        ARG.tout = h;
        rkf45(&ARG);
        out << std::setw(5)  << std::setprecision(2) << ARG.t
            << std::setw(11) << std::setprecision(6) << x[0]
            << std::setw(11) << std::setprecision(6) << x[1]
            << std::setw(4)  << std::setprecision(0) << ARG.flag
            << std::endl;
    }
    out <<"_______________________________________\n\n";
}