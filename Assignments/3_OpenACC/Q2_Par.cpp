#include <iostream>
#include <math.h>
#include <time.h>

#define TYPE float
#define N 10
#define sval 0.001
#define tol 1e-6

// OG FUNCTION
void initmult (TYPE mat[][N]) {
    #pragma acc parallel loop present(mat[:][:])
        for (int i = 0; i < N; i++) {
            #pragma acc loop vector
                for (int j = 0; j < i; j++) {
                    // mat[i][j] = (i + j) / pow(N,2);
                    mat[i][j] = (i + j) / (N * N);
                    mat[j][i] = mat[i][j];
                }

            mat[i][i] = 1;
        }
}

// // TEST FUNCTION ==> TEST PASSED
// void initmult (TYPE mat[][N]) {
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < i; j++) {
//             mat[i][j] = abs(i - j) / pow(N,2);
//             mat[j][i] = mat[i][j];
//         }

//         // PD
//         // mat[i][i] = i + 1;

//         // PSD
//         // if (i < N-1) {
//         //     mat[i][i] = i + 1;
//         // }
//         // else {
//         //     mat[i][i] = 0;
//         // }

//         // ND
//         // mat[i][i] = -(i + 1);
//     }
// }

void printMat (TYPE a[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.4f ", a[i][j]);
        }
        printf("\n");
    }
}

void cholesky (TYPE a[N][N], TYPE l[N][N], int stat[2]) {
    stat[0] = 0;
    stat[1] = 0;
    // #pragma acc parallel present(a[:][:], l[:][:], stat[:])
        for (int i = 0; i < N; i++) {        
            // #pragma acc loop worker
            for (int j = 0; j < i; j++) {
                if (abs(a[i][j] - a[j][i]) > tol) {
                    // printf("Matrix is not symmetric...\n");
                    // return;
                    stat[0] = 1;
                }

                // can be parallelized
                for (int k = 0; k < j; k++) {
                    l[i][j] += l[i][k] * l[j][k];
                }

                l[i][j] = a[i][j] - l[i][j];
                
                l[i][j] /= (l[j][j] > sval ? l[j][j] : 1);
            }

            // can be parallelized
            for (int k = 0; k < i; k++) {
                l[i][i] += pow(l[i][k],2);
            }

            if (a[i][i] - l[i][i] < 0) {
                // printf("Matrix is not positive definite...\n");
                // return;
                stat[1] = 1;
            }

            l[i][i] = sqrt(a[i][i] - l[i][i]);
        }
}


int main () {
    TYPE a[N][N];
    TYPE l[N][N] {};
    int stat[2] {};

    #pragma acc data create(a[:][:]) copyin(l[:][:], stat[:]) copyout(a[:][:], l[:][:], stat[:])
    {
        initmult(a);
        // cholesky(a, l, stat);
    }
    
    printMat(a);
    printf("\n");
    if (stat[0] == 1) {
        printf("Matrix is not symmetric...\n");
    }
    if (stat[1] == 1) {
        printf("Matrix is not positive definite...\n");
    }
    if (stat[0] == 0 && stat[1] == 0) {
        printMat(l);
    }

    return 0;
}