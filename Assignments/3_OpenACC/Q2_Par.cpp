#include <iostream>
#include <math.h>
#include <time.h>

#define TYPE float
#define N 10
#define sval 0.001
#define tol 1e-6

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

TYPE myabs(TYPE a) {
    if (a >= 0) {
        return a;
    }
    else {
        return -a;
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

void cholesky (TYPE a[N][N], TYPE l[N][N], int stat[2]) {
    // Ensure that stat has both terms as 0
    
    #pragma acc parallel loop gang default(present)
        for (int i = 0; i < N; i++) {        
            #pragma acc loop vector
                for (int j = 0; j < i; j++) {
                    if (myabs(a[i][j] - a[j][i]) > tol) {
                        stat[0] = 1;
                    }
                }
        }
    
    if (stat[0] == 1) {
        printf("Matrix is not symmetric...\n");
        return;
    }
    double sum = 0;

    #pragma acc data present(a[:][:], l[:][:], stat[:]) copyin(stat)
    {
        #pragma acc parallel num_gangs(1)
        {
            for (int i = 0; i < N; i++) {        
                for (int j = 0; j < i; j++) {
                    sum = 0;
                    // can be parallelized
                    #pragma acc loop vector reduction(+: sum)
                        for (int k = 0; k < j; k++) {
                            sum += l[i][k] * l[j][k];
                        }

                    l[i][j] = a[i][j] - sum;
                    
                    l[i][j] /= (l[j][j] > sval ? l[j][j] : 1);
                }

                sum = 0;
                // can be parallelized
                #pragma acc loop vector reduction(+: sum)
                    for (int k = 0; k < i; k++) {
                        sum += l[i][k] * l[i][k];
                    }

                if (a[i][i] - sum < 0) {
                    stat[1] = 1;
                }

                l[i][i] = sqrt(a[i][i] - sum);
            }
        }
    }
}


int main () {
    TYPE a[N][N];
    TYPE l[N][N] {};
    int stat[2] {};

    #pragma acc data copyout(a[:][:], l[:][:], stat[:])
    {
        initmult(a);
        cholesky(a, l, stat);
    }

    printMat(a);
    printf("\n");
    printMat(l);

    return 0;
}