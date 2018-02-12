#include <iostream>
#include <cmath>
#include "Forsythe.h"


int main() {

    double alpha = 2;

    double H[4][4] = {{alpha,         1,             0,             0},
                      {pow(alpha, 2), alpha,         1,             0},
                      {pow(alpha, 3), pow(alpha, 2), alpha,         1},
                      {pow(alpha, 4), pow(alpha, 3), pow(alpha, 2), alpha}};

    double E[4][4] = {{1, 0, 0, 0},
                      {0, 1, 0, 0},
                      {0, 0, 1, 0},
                      {0, 0, 0, 1}};


    double beta[3] = {0.1, 0.01, 0.001};
    double A1[4][4];



// Вычисления для betta1
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            A1[i][j] = E[i][j] * beta[0] + H[i][j];
        }
    }

// DECOMP и SOLVE
    const int size = 4;
    double cond;
    double newA[size * size];
    int ipvt[size];
    int k = 0;
    std::cout.precision(2);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {

            newA[k] = A1[i][j];
            k++;
        }
    }


    Decomp(size, newA, &cond, ipvt);
    std::cout << "После DECOMP";
    std::cout << "\n";
    for (int i = 0; i < size * size; i++) {
        if (i % 4 == 0) {

            std::cout << "\n";
        }
        std::cout << newA[i] << " |";
    }

    std::cout << "После SOLVE";
    std::cout << "\n";

    double B[size];

    for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++) {
            B[k] = 0;
        }
        B[i] = 1;
        Solve(size, newA, B, ipvt);

    }

    std::cout << "\n";
    std::cout << "После A^-1";
 // ВЫВОД РЕЗУЛЬТАТА A^-1
    for (int i = 0; i < size * size; i++) {
        if (i % 4 == 0) {

            std::cout << "\n";
        }
        std::cout << B[i] << " |";
    }
    std::cout << "\n";
    // ВЫВОД R
    std::cout << "\n";
    std::cout << "Вычисление R";
    double R[4][4];
    int t;
    for (int i=0; i<4 ;i++) {
        std::cout << "\n";
        for (int j = 0; j < 4; j++) {
            t++;
            R[i][j]=B[t]*A1[i][j]-E[i][j];
            std::cout << "|" << R[i][j];
        }
    }


    std::cout << "\n";
    std::cout << "Норма матрицы R";

    double norm=0;
    double tmp=0;

    for(int i = 0; i<size; i++){
        tmp=0;
        for(int j = 0; j<size; j++){
            tmp=tmp+abs(R[j][i]);
        }
        if(tmp>norm) {norm=tmp;}
    }

    std::cout << "\n";
    std::cout << norm;




}
