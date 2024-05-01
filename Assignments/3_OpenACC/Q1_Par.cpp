#include <iostream>
#include <fstream>
#include <math.h>

#define N 10 // Problem size
#define Ng 1 // No. of gangs

// Functions
double myfunc (double x) {
    double f = sin(5*x);
    return f;
}

double myderv (double x) {
    double fdot = 5*cos(5*x);
    return fdot;
}

void LUfunc (double A[N][N], double L[N][N], double U[N][N]) {

}

void printVector(double a[], int n) {
   for (int i = 0; i < n; i++) {
        printf("%.4f\n",a[i]);
   }
}

void printMatrix (double** a, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.4f\t",a[i][j]);
        }
        printf("\n");
    }
}



int main (int argc, char* argv[]) {
    int i {}, j {};
    printf("\n");

    double A[N][N] {}, xgrid[N] {};
    double b[N] {}, x[N] {}, y[N] {};
    double fvec[N] {}, fdvec[N] {};
    double L[N][N] {}, U[N][N] {};
    double h = 3/(double(N-1));

    for (i = 0; i < N; i++) {
        xgrid[i] = i*h;
        fvec[i] = myfunc(xgrid[i]);
        fdvec[i] = myderv(xgrid[i]);
    }

    // Initializing A and b
    for (int i = 0; i < N; i++)  {
        if (i == 0) {
            A[i][0] = 1;
            A[i][1] = 2;

            b[i] = (-2.5*fvec[0] + 2*fvec[1] + 0.5*fvec[2])/h;
        }
        else if (i == N - 1) {
            A[i][N-2] = 2;
            A[i][N-1] = 1;

            b[i] = (2.5*fvec[N-1] - 2*fvec[N-2] - 0.5*fvec[N-3])/h;
        }
        else {
            A[i][i-1] = 1;
            A[i][i] = 4;
            A[i][i+1] = 1;

            b[i] = 3*(fvec[i+1] - fvec[i-1])/h;
        }
    }

    // Ax = b => LUx = b => Ly = b & Ux = y

    printf("\nx = \n");
    printVector(x, N);
    printf("\ny = \n");
    printVector(y, N);

    // The following code is to output to a .txt file to plot the solution
    
    // std::ofstream oFile("./Res/LU.txt");

    // if (oFile.is_open()) {
    //     for (i = 0; i < N; i++) {
    //         oFile << x[i] << "," << fdvec[i] << "\n";
    //     }
    //     oFile.close();
    //     printf("Saved in file ./Res/LU.txt\n");
    // }
    // else {
    //     printf("Error opening file\n");
    // }

    return 0;
}