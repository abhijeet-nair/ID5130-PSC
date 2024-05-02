#include <iostream>
#include <fstream>
#include <math.h>

#define N 10     // Problem size
#define Ng 1     // No. of gangs
#define tol 1e-6 // Tolerance

// Functions
double myfunc (double x) {
    double f = sin(5*x);
    return f;
}

double myderv (double x) {
    double fdot = 5*cos(5*x);
    return fdot;
}

// void LUfunc1 (double A[N][N], double L[N][N], double U[N][N]) {
//     for (int k = 0; k < N; k++) {
//         L[k][k] = 1;
//         for (int i = k + 1; i < N; i++) {
//             L[i][k] = A[i][k]/A[k][k]; // Multiplier (L part)
//         }

//         for (int i = k + 1; i < N; i++) {
//             for (int j = k + 1; j < N; j++) {
//                 U[i][j] += -L[i][k]*A[k][j]; // U part
//             }
//         }
//     }
// }

void LUfunc2 (double A[N][N], double L[N][N], double U[N][N]) {
    int i = 0, j = 0, k = 0;
    // for (i = 0; i < N; i++)
    // {
    //     for (j = 0; j < N; j++)
    //     {
    //         if (j < i)
    //             L[j][i] = 0;
    //         else
    //         {
    //             L[j][i] = A[j][i];

    //             for (k = 0; k < i; k++)
    //             {
    //                 L[j][i] = L[j][i] - L[j][k] * U[k][i];
    //             }
    //         }
    //     }
    //     for (j = 0; j < N; j++)
    //     {
    //         if (j < i)
    //             U[i][j] = 0;
    //         else if (j == i)
    //             U[i][j] = 1;
    //         else
    //         {
    //             U[i][j] = A[i][j] / L[i][i];
    //             for (k = 0; k < i; k++)
    //             {
    //                 U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
    //             }
    //         }
    //     }
    // }
    for (i = 0; i < N; i++) {
        for (j = 0; j < i; )
    }
}

void multiplyMatrices (double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i ++) {
		for (int j = 0; j < N; j ++) {
			C[i][j] = 0;
			for (int k = 0; k < N; k ++) {
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void printVector(double a[N]) {
   for (int i = 0; i < N; i++) {
        if (abs(a[i]) < tol) {
            printf("0\n");
        }
        else {
            printf("%.4f\n", a[i]);
        }
   }
}

void printMat (double a[N][N]) {
    if (N > 10) {
        printf("Avoiding printing of large matrix...\n");
    }
    else {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (abs(a[i][j]) < tol) {
                    printf("0\t");
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
    for (i = 0; i < N; i++)  {
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
    LUfunc2(A, L, U);

    double C[N][N];
    multiplyMatrices(L, U, C);

    printf("L = \n");
    printMat(L);
    printf("\n");

    printf("U = \n");
    printMat(U);
    printf("\n");

    printf("C = \n");
    printMat(C);
    printf("\n");

    printf("A = \n");
    printMat(A);
    printf("\n");

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%.4f ", A[i][j] - C[i][j]);
        }
        printf("\n");
    }


    printf("\nx = \n");
    printVector(x);
    printf("\ny = \n");
    printVector(y);

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