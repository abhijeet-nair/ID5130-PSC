#include <iostream>
#include <math.h>
#include <time.h>

// Datatype
#define TYPE float
// Problem size
#define N 10
// A small value
#define sval 0.001
// Tolerance
#define tol 1e-6

// Initialization Function
void initmult (TYPE mat[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            mat[i][j] = (i + j) / pow(N,2);
            mat[j][i] = mat[i][j];
        }

        mat[i][i] = 1;
    }
}

// Matrix printing
void printMat (TYPE a[N][N]) {
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

// Cholesky Decomposition. Also, checks for Positive Definiteness
void cholesky (TYPE a[N][N], TYPE l[N][N]) {
    for (int i = 0; i < N; i++) {        
        for (int j = 0; j < i; j++) {
            if (abs(a[i][j] - a[j][i]) > tol) {
                printf("Matrix is not symmetric...\n");
                return;
            }

            for (int k = 0; k < j; k++) {
                l[i][j] += l[i][k] * l[j][k];
            }

            l[i][j] = a[i][j] - l[i][j];
            
            l[i][j] /= (l[j][j] > sval ? l[j][j] : 1);
        }

        for (int k = 0; k < i; k++) {
            l[i][i] += pow(l[i][k],2);
        }

        if (a[i][i] - l[i][i] < 0) {
            printf("Matrix is not positive definite...\n");
            return;
        }

        l[i][i] = sqrt(a[i][i] - l[i][i]);
    }
}


int main () {
    clock_t t;
    
    TYPE a[N][N], l[N][N] {};
    t = clock();

    initmult(a);
    cholesky(a, l);

    t = clock() - t;

    printf("A = \n");
    printMat(a);
    printf("\n");
    
    printf("L = \n");
    printMat(l);

    double time_taken = double(t) / double(CLOCKS_PER_SEC);
    printf("\nTime taken = %.6f secs\n", time_taken);
    return 0;
}