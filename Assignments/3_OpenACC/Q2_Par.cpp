#include <iostream>
#include <math.h>
#include <time.h>

#define TYPE float // Datatype
#define N 500      // Problem size
#define sval 0.001 // A small value
#define tol 1e-6   // Tolerance

void initmult (TYPE mat[N][N]) {
    #pragma acc parallel loop present(mat[:][:])
        for (int i = 0; i < N; i++) {
            #pragma acc loop vector
                for (int j = 0; j < i; j++) {
                    // mat[i][j] = (i + j) / pow(N,2);
                    mat[i][j] = TYPE(i + j) / (N * N);
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

void printMat (TYPE a[N][N]) {
    if (N > 10) {
        printf("Avoiding printing of large matrix...\n");
    }
    else {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (a[i][j] < tol) {
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

void cholesky (TYPE a[N][N], TYPE l[N][N], int stat) {
    // Ensure that stat is 0
    
    // Symmetric Matrix Check
    #pragma acc parallel loop gang default(present) copyout(stat)
        for (int i = 0; i < N; i++) {        
            #pragma acc loop vector
                for (int j = 0; j < i; j++) {
                    if (myabs(a[i][j] - a[j][i]) > tol) {
                        stat = 1;
                    }
                }
        }
    
    if (stat == 1) {
        printf("Matrix is not symmetric...\n");
        return;
    }

    stat = 0; // Re-initializing, if modified
    double sum = 0;

    // Main Cholesky Decomposition Code
    #pragma acc data present(a[:][:], l[:][:]) copyin(sum) copy(stat)
    {
        #pragma acc parallel num_gangs(1)
        {
            for (int i = 0; i < N; i++) {        
                for (int j = 0; j < i; j++) {
                    sum = 0;
                    // Can be parallelized
                    #pragma acc loop vector reduction(+: sum)
                        for (int k = 0; k < j; k++) {
                            sum += l[i][k] * l[j][k];
                        }

                    l[i][j] = a[i][j] - sum;
                    
                    l[i][j] /= (l[j][j] > sval ? l[j][j] : 1);
                }

                sum = 0;
                // Can be parallelized
                #pragma acc loop vector reduction(+: sum)
                    for (int k = 0; k < i; k++) {
                        sum += l[i][k] * l[i][k];
                    }

                // Positive Definite Matrix Check
                if (a[i][i] - sum < 0) {
                    stat = 1;
                }

                l[i][i] = sqrt(a[i][i] - sum);
            }
        }
    }
}


int main () {
    TYPE a[N][N];
    TYPE l[N][N] {};
    int stat {};

    #pragma acc data copyout(a[:][:], l[:][:]) copy(stat)
    {
        initmult(a);
        cholesky(a, l, stat);
    }
    
    if (stat == 1) {
        printf("Matrix is not positive definite...\n");
    }
    else {
        printMat(a);
        printf("\n");
        printMat(l);
        printf("\n");
    }

    int r1, r2;
    for (int i = 0; i < 3; i++) {
        r1 = rand()%(N-1);
        r2 = rand()%(N-1);
        printf("(%d, %d) ==>\t", r1, r2);
        printf("a = %.4f\tl = %.4f\n",a[r1][r2],l[r1][r2]);
    }
    
    return 0;
}