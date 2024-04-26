#include <iostream>
#include <math.h>
#include <time.h>

#define TYPE float
#define N 10
#define sval 0.001

void initmult (TYPE mat[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            mat[i][j] = (i + j) / pow(N,2);
            mat[j][i] = mat[i][j];
        }

        mat[i][i] = 1;
    }
}

void printMat (TYPE a[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.4f ", a[i][j]);
        }
        printf("\n");
    }
}

void cholesky (TYPE a[N][N], TYPE l[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < j; k++) {
                l[i][j] += l[i][k] * l[j][k];
            }
            l[i][j] = a[i][j] - l[i][j];
            
            l[i][j] /= (l[j][j] > sval ? l[j][j] : 1);
        }

        for (int k = 0; k < i; k++) {
            l[i][i] += pow(l[i][k],2);
        }

        l[i][i] = sqrt(a[i][i] - l[i][i]);
    }
}


int main () {
    TYPE a[N][N];

    initmult(a);
    printMat(a);
    printf("\n");

    TYPE l[N][N] {};
    cholesky(a, l);
    printMat(l);

    return 0;
}