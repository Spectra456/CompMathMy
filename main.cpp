#include "stdio.h"
#include "math.h"
#include "rkf45.h"

void f(float t, float *x, float *dx) {
    dx[0] = -45 * x[0] + 60 * x[1] + sin(1+t);
    dx[1] = 70 * x[0] - 110 * x[1] +cos(1-t) + t + 1;
    return;
}
int main() {
    int n = 2;
    float x0[2] = { 5, -1 };
    float t = 0, tout = 0;
    float re = 0.0001, ae = 0.0001;
    int iflag = 1;
    float work[15];
    int iwork[30];
    float h = 0.05;
    while (tout < 1.01) {
        RKF45(f, n, x0, &t, &tout, &re, &ae, &iflag, work, iwork);
        printf("t = %.2f, x = { %f; %f }, iflag = %d\n", t, x0[0], x0[1], iflag);
        tout += h;
    }

    float zn[2], k1[2], k2[2], k3[2], k4[2], ktmp[2];
    printf("\nRK 4 degree of accuracy, h = 0.025\n");
    x0[0] = 0; x0[1] = 1;
    t = 0;
    h = 0.025;
    for (int i;t < 1.01;i++) {
        f(t, x0, k1);
        k1[0] *= h;
        k1[1] *= h;
        ktmp[0] = x0[0] + k1[0]/3;
        ktmp[1] = x0[1] + k1[1]/3;
        f(t + h/3, ktmp, k2);
        k2[0] *= h;
        k2[1] *= h;
        ktmp[0] = x0[0] - k1[0]/3 + k2[0];
        ktmp[1] = x0[1] - k1[1]/3 + k2[1];
        f(t + 2*h/3, ktmp, k3);
        k3[0] *= h;
        k3[1] *= h;
        ktmp[0] = x0[0] + k1[0] - k2[0] + k3[0];
        ktmp[1] = x0[1] + k1[1] - k2[1] + k3[1];
        f(t+h, ktmp, k4);
        k4[0] *= h;
        k4[1] *= h;
        zn[0] = x0[0] + (k1[0] + 3*k2[0] + 3*k3[0] + k4[0])/8;
        zn[1] = x0[1] + (k1[1] + 3*k2[1] + 3*k3[1] + k4[1])/8;
        x0[0] = zn[0];
        x0[1] = zn[1];
        if (!(i%2)) {
            printf("t = %.2f, x = { %f; %f }\n", t, zn[0], zn[1]);
        }
        t += h;
    }
    printf("\nRK 4 degree of accuracy, h = 0.001\n");
    x0[0] = 0; x0[1] = 1;
    t = 0;
    h = 0.001;
    for (int i; t < 1.005; i++) {
        f(t, x0, k1);
        k1[0] *= h;
        k1[1] *= h;
        ktmp[0] = x0[0] + k1[0]/3;
        ktmp[1] = x0[1] + k1[1]/3;
        f(t + h/3, ktmp, k2);
        k2[0] *= h;
        k2[1] *= h;
        ktmp[0] = x0[0] - k1[0]/3 + k2[0];
        ktmp[1] = x0[1] - k1[1]/3 + k2[1];
        f(t + 2*h/3, ktmp, k3);
        k3[0] *= h;
        k3[1] *= h;
        ktmp[0] = x0[0] + k1[0] - k2[0] + k3[0];
        ktmp[1] = x0[1] + k1[1] - k2[1] + k3[1];
        f(t+h, ktmp, k4);
        k4[0] *= h;
        k4[1] *= h;
        zn[0] = x0[0] + (k1[0] + 3*k2[0] + 3*k3[0] + k4[0])/8;
        zn[1] = x0[1] + (k1[1] + 3*k2[1] + 3*k3[1] + k4[1])/8;
        x0[0] = zn[0];
        x0[1] = zn[1];
        if (!(i%20)) {
            printf("t = %.2f, x = { %f; %f }\n", t, zn[0], zn[1]);
        }
        t += h;
    }
}