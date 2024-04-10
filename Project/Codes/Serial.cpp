#include <iostream>
#include <math.h>
#include <omp.h>

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

    double A[N][N], B[N];

    return 0;
}