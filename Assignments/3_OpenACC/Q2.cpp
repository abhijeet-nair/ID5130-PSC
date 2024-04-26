#include <iostream>
#include <math.h>
#include <time.h>

#define TYPE float
#define N 10
#define sval 0.001
#define tol 1e-6

void initmult (TYPE mat[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            mat[i][j] = (i + j) / pow(N,1);
            mat[j][i] = mat[i][j];
        }

        // if (i != N) {
        //     mat[i][i] = 10*(i+1);
        // }
        // else {
        //     mat[i][i] = 0;
        // }
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
        for (int k = 0; k < i; k++) {
            l[i][i] += pow(l[i][k],2);
        }

        if (a[i][i] - l[i][i] <=0) {
            printf("Matrix is possibly not positive definite...\n");
            return;
        }

        l[i][i] = sqrt(a[i][i] - l[i][i]);

        for (int j = 0; j < i; j++) {
            if (abs(a[i][j] - a[j][i]) > tol) {
                // printf("(i,j) = (%d, %d)\tval = %.4f - %.4f = %.4f (log = %.4f)\n",i,j,a[i][j],a[j][i],abs(a[i][j] - a[j][i]),log10(abs(a[i][j] - a[j][i])));
                printf("Matrix is not symmetric...\n");
                return;
            }

            for (int k = 0; k < j; k++) {
                l[i][j] += l[i][k] * l[j][k];
            }

            l[i][j] = a[i][j] - l[i][j];
            
            l[i][j] /= (l[j][j] > sval ? l[j][j] : 1);
        }
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