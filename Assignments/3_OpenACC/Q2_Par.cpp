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

// Function to make all elements of a matrix zero
void nullFunc (TYPE A[N][N]) {
    #pragma acc parallel loop num_gangs(10) collapse(2) default(present)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 0;
            }
        }
}

// Custom absolute value function (to avoid header issue in GPU)
TYPE myabs(TYPE a) {
    if (a >= 0) {
        return a;
    }
    else {
        return -a;
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
void cholesky (TYPE a[N][N], TYPE l[N][N], int stat[2], TYPE sum) {
    // Ensure that stat is 0
    
    // Symmetric Matrix Check
    #pragma acc parallel loop num_gangs(10) default(present)
        for (int i = 0; i < N; i++) {        
            #pragma acc loop vector
                for (int j = 0; j < i; j++) {
                    if (myabs(a[i][j] - a[j][i]) > tol) {
                        stat[0] = 1;
                    }
                }
        }
    
    // Main Cholesky Decomposition Code
    #pragma acc parallel num_gangs(1) default(present)
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
                stat[1] = 1;
            }

            l[i][i] = sqrt(a[i][i] - sum);
        }
    }
}


int main () {
    clock_t t;

    TYPE a[N][N];
    TYPE l[N][N];
    int stat[2] {};
    TYPE sum;

    t = clock();

    #pragma acc data copy(stat) create(sum, a[:][:], l[:][:])
    {
        initmult(a);
        nullFunc(l);
        cholesky(a, l, stat, sum);
    }

    t = clock() - t;
    
    if (stat[0] == 1) {
        printf("Matrix is not symmetric...\n");
    }
    else if (stat[1] == 1) {
        printf("Matrix is not positive definite...\n");
    }
    else {
        printf("L = \n");
        printMat(l);
    }

    if (N > 10) {
        printf("\nPrinting random locations for check...\n");
        int r1, r2;
        for (int i = 0; i < 3; i++) {
            r1 = rand()%(N-1);
            r2 = rand()%(N-1);
            printf("(%d, %d) ==>\t", r1, r2);
            // printf("a = %.4f\tl = %.4f\n",a[r1][r2],l[r1][r2]);
            printf("l = %.4f\n",l[r1][r2]);
        }
    }
    
    double time_taken = double(t) / double(CLOCKS_PER_SEC);
    printf("\nTime taken = %.6f secs\n", time_taken);
    return 0;
}