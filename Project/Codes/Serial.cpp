#include <iostream>
#include <math.h>
#include <string.h>
#include <omp.h>

void getIndVel(double g, double x, double y, double x0, double y0, double uv[2]) {
    double den = pow((x - x0), 2) + pow((y - y0), 2);
    uv[0] = (g*(y - y0))/(2*M_PI*den);
    uv[1] = -(g*(x - x0))/(2*M_PI*den);
}

int main () {
    int i, j;           // Indices

    int N = 10;         // No. of collocation points
    double u   = 20;    // Freestream velocity
    double c   = 5;     // Chord length of the flat plate
    double alp = 0;     // Angle of attack of the plate
    double rho = 1.225; // Density
    double dx  = c/N;   // Spacing between two points

    double sn = sin(alp*M_PI/180);
    double cn = cos(alp*M_PI/180);

    double vorloc[N], colcloc[N], plate[N+1];

    for (i = 0; i < N; i++) {
        vorloc[i]  = (i + 0.25)*dx;
        colcloc[i] = (i + 0.75)*dx;
        plate[i]   = i*dx;
    }
    plate[N] = N*dx;

    double A[N][N], B[N], gIVRes[2];

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            getIndVel(1, colcloc[i], 0, vorloc[j], 0, gIVRes);
            A[i][j] = gIVRes[1];
        }
        getIndVel(1, colcloc[i], 0, (c + 0.1*dx), 0, gIVRes);
        B[i] = gIVRes[1];
    }

    double C[N+1][N+1] {};

    for (i = 0; i <= N; i++) {
        for (j = 0; j <= N; j++) {
            
        }
    }

    return 0;
}