#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

// Problem size
#define N 1000
// No. of gangs
#define Ng 10
// Tolerance
#define tol 1e-6

// Given function
double myfunc (double x) {
    double f = sin(5*x);
    return f;
}

// Actual derivative function
double myderv (double x) {
    double fdot = 5*cos(5*x);
    return fdot;
}

// Function to make all elements of a matrix zero
void nullFunc (double A[N][N]) {
    #pragma acc parallel loop num_gangs(Ng) collapse(2) default(present)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 0;
            }
        }
}

// Initialization Function
void initFunc (double A[N][N], double b[N], double xgrid[N], double fvec[N], double fdvec[N], double h) {
    #pragma acc parallel num_gangs(Ng) default(present)  present(fvec[0:N])
    {
        #pragma acc loop
            for (int i = 0; i < N; i++) {
                xgrid[i] = i*h;
                fvec[i] = myfunc(xgrid[i]);
                fdvec[i] = myderv(xgrid[i]);
            }

        // Initializing A and b. Compiler thought it needed fvec[-1:12] here.
        #pragma acc loop
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
    }
}

// LU Decomposition for Tridiagonal systems (Check notes)
void LUfunc (double A[N][N], double L[N], double U[N]) {
    #pragma acc serial default(present)
    {    
        U[0] = A[0][0];
        L[0] = 0;
        
        for (int i = 1; i < N; i++) {
            L[i] = A[i][i-1] / U[i-1];
            U[i] = A[i][i] - L[i] * A[i-1][i];
        }
    }
}

// Forward and backward substitution to get final solution
void subsFunc (double A[N][N], double b[N], double L[N], double U[N], double x[N], double y[N]) {
    #pragma acc serial default(present)
    {
        // Ly = b
        y[0] = b[0];
        
        for (int i = 1; i < N; i++) {
            y[i] = b[i] - L[i] * y[i-1];
        }

        // Ux = y
        x[N-1] = y[N-1] / U[N-1];

        for (int i = N-2; i >= 0; i--) {
            x[i] = (y[i] - A[i][i+1] * x[i+1]) / U[i];
        }
    }
}

// Vector printing
void printVector(double a[N]) {
    if (N > 10) {
        printf("Avoiding printing of large vector...\n");
    }
    else {
        for (int i = 0; i < N; i++) {
                if (abs(a[i]) < tol) {
                    printf("0      ");
                }
                else {
                    printf("%.4f ", a[i]);
                }
        }
        printf("\n");
    }
}

// Matrix printing
void printMat (double a[N][N]) {
    if (N > 10) {
        printf("Avoiding printing of large matrix...\n");
    }
    else {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (abs(a[i][j]) < tol) {
                    printf("0      ");
                }
                else {
                    printf("%.4f ", a[i][j]);
                }
            }
            printf("\n");
        }
    }
}



int main (int argc, char* argv[]) {
    clock_t t;

    printf("N = %d\n\n",N);

    double A[N][N], xgrid[N];
    double b[N], x[N], y[N];
    double fvec[N], fdvec[N];
    double L[N] {}, U[N] {};
    double h = 3/(double(N-1));
    int fs = 0;

    t = clock();

    #pragma acc data create(A[0:N][0:N], b[0:N], xgrid[0:N], y[0:N], fvec[0:N]) \
    copyin(h) copyout(L[0:N], U[0:N], x[0:N], fdvec[0:N])
    {
        nullFunc(A);
        initFunc(A, b, xgrid, fvec, fdvec, h);

        // Ax = b => LUx = b => Ly = b & Ux = y
        LUfunc(A, L, U);
        subsFunc(A, b, L, U, x, y);
    }

    t = clock() - t;

    // printf("A = \n");
    // printMat(A);
    // printf("\n");

    printf("L = \n");
    printVector(L);
    printf("\n");

    printf("U = \n");
    printVector(U);
    printf("\n");

    printf("\nx = \n");
    printVector(x);
    printf("\n");

    // printf("\ny = \n");
    // printVector(y);
    // printf("\n");

    double val;
    if (N <= 10) {
        printf("res = \n");

        for (int i = 0; i < N; i++) {
            val = x[i] = fdvec[i];
            
            if (abs(val) < tol) {
                printf("0\n");
            }
            else {
                printf("%.4f\n",val);
            }
        }
    }
    else {
        val = 0;

        for (int i = 0; i < N; i++) {
            val += pow(x[i] - fdvec[i], 2);
        }

        printf("Norm of res = %.4f (log = %.4f)\n",sqrt(val),0.5*log10(val));
    }

    // The following code is to output to a .txt file to plot the solution
    if (fs == 1) {
        printf("\n");
        std::ofstream oFile("./ParLU.txt");

        if (oFile.is_open()) {
            for (int i = 0; i < N; i++) {
                oFile << x[i] << "," << fdvec[i] << "\n";
            }
            oFile.close();
            printf("Saved in file ./ParLU.txt\n");
        }
        else {
            printf("Error opening file\n");
        }
    }
    else {printf("Not saving in file...\n");}


    double time_taken = double(t) / double(CLOCKS_PER_SEC);
    printf("\nTime taken = %.6f secs\n", time_taken);
    return 0;
}