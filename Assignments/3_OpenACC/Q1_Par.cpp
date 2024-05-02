#include <iostream>
#include <fstream>
#include <math.h>

// Problem size
#define N 750
// No. of gangs
#define Ng 10
// Tolerance
#define tol 1e-6

// Functions
double myfunc (double x) {
    double f = sin(5*x);
    return f;
}

double myderv (double x) {
    double fdot = 5*cos(5*x);
    return fdot;
}

void initFunc (double A[N][N], double b[N], double xgrid[N], double fvec[N], double fdvec[N], double h) {
    #pragma acc parallel loop num_gangs(Ng) default(present)
        for (int i = 0; i < N; i++) {
            xgrid[i] = i*h;
            fvec[i] = myfunc(xgrid[i]);
            fdvec[i] = myderv(xgrid[i]);
        }

    
    // Initializing A and b. Compiler thought it needed fvec[-1:12] here.
    #pragma acc parallel loop num_gangs(Ng) default(present) present(fvec[0:N])
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

void LUfunc (double A[N][N], double L[N][N], double U[N][N]) {
    #pragma acc serial default(present)
    {
        U[0][0] = A[0][0];
        U[0][1] = A[0][1];
        L[0][0] = 1;
    
        for (int i = 1; i < N; i++) {
            L[i][i] = 1;
            L[i][i-1] = A[i][i-1] / U[i-1][i-1];
            U[i][i] = A[i][i] - L[i][i-1] * A[i-1][i];
            if (i < N-1) {
                U[i][i+1] = A[i][i+1];
            }
        }
    }
}

void subsFunc (double b[N], double L[N][N], double U[N][N], double x[N], double y[N]) {
    #pragma acc serial default(present)
    {
        // Ly = b
        y[0] = b[0] / L[0][0];
        
        for (int i = 1; i < N; i++) {
            y[i] = b[i] - L[i][i-1] * y[i-1];
        }

        // Ux = y
        x[N-1] = y[N-1] / U[N-1][N-1];

        for (int i = N-2; i >= 0; i--) {
            x[i] = (y[i] - U[i][i+1] * x[i+1]) / U[i][i];
        }
    }
}

void multiplyMatrices (double A[N][N], double B[N][N], double C[N][N]) {
    #pragma acc parallel loop num_gangs(Ng) collapse(2) default(present)
        for (int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j ++) {
                double sum = 0;
                #pragma acc loop vector reduction(+: sum)
                    for (int k = 0; k < N; k ++) {
                        sum += A[i][k]*B[k][j];
                    }
                
                C[i][j] = sum;
            }
        }
}

void printVector(double a[N]) {
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
    double A[N][N] {}, xgrid[N];
    double b[N], x[N], y[N];
    double fvec[N], fdvec[N];
    double L[N][N] {}, U[N][N] {};
    double C[N][N];
    double h = 3/(double(N-1));

    #pragma acc data copy(A[0:N][0:N], L[0:N][0:N], U[0:N][0:N]) create(b[0:N], xgrid[0:N], y[0:N], fvec[0:N]) \
    copyin(h) copyout(x[0:N], fdvec[0:N], C[0:N][0:N])
    {
        initFunc(A, b, xgrid, fvec, fdvec, h);

        // Ax = b => LUx = b => Ly = b & Ux = y
        LUfunc(A, L, U);
        subsFunc(b, L, U, x, y);

        // C = L*U
        multiplyMatrices(L, U, C);
    }

    printf("L = \n");
    printMat(L);
    printf("\n");

    printf("U = \n");
    printMat(U);
    printf("\n");

    // printf("C = \n");
    // printMat(C);
    // printf("\n");

    printf("A = \n");
    printMat(A);
    printf("\n");

    double val;
    if (N <= 10) {
        printf("res = \n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                val = A[i][j] - C[i][j];
                if (abs(val) < tol) {
                    printf("0      ");
                }
                else {
                    printf("%.4f ",val);
                }
            }
            printf("\n");
        }
    }
    else {
        val = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                val += pow(A[i][j] - C[i][j], 2);
            }
        }
        printf("Norm of res = %.4f",sqrt(val));
    }
    

    

    // printf("\nx = \n");
    // printVector(x);
    // printf("\ny = \n");
    // printVector(y);

    // The following code is to output to a .txt file to plot the solution
    // printf("\n");
    // std::ofstream oFile("./ParLU.txt");

    // if (oFile.is_open()) {
    //     for (int i = 0; i < N; i++) {
    //         oFile << x[i] << "," << fdvec[i] << "\n";
    //     }
    //     oFile.close();
    //     printf("Saved in file ./ParLU.txt\n");
    // }
    // else {
    //     printf("Error opening file\n");
    // }

    return 0;
}