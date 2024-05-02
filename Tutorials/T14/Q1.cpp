#include <iostream>
#include <math.h>
#include <random>

#define M 5
#define N 4
#define tol 1e-6

double myabs (double x) {
    if (x > 0) {
        return x;
    }
    else {
        return (-x);
    }
}

void multiplyMatrices (double A[M][N], double B[N][N], double C[M][N]) {
    for (int i = 0; i < M; i ++) {
		for (int j = 0; j < N; j ++) {
			C[i][j] = 0;
			for (int k = 0; k < N; k ++) {
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

int main () {
    double A[M][N], Q[M][N], R[N][N] {};
    std::uniform_real_distribution<double> unif(0,M-1);
    std::default_random_engine re;
    
    for (int i = 0; i < M; i++) {
        for (int j =0; j < N; j++) {
            A[i][j] = unif(re);
            Q[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            R[i][i] += Q[j][i] * Q[j][i];
        }
        R[i][i] = sqrt(R[i][i]);

        for (int j = 0; j < M; j++) {
            Q[j][i] = Q[j][i] / R[i][i];
        }

        for (int j = i+1; j < N; j++) {
            for (int k = 0; k < M; k++) {
                R[i][j] += Q[j][k] * Q[k][i];
            }

            for (int k = 0; k < M; k++) {
                Q[k][j] = Q[k][j] - R[i][j] * Q[k][i];
            }
        }
    }

    double res[M][N], val;
    multiplyMatrices(Q, R, res);

    printf("A = \n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            val = A[i][j];
            if (myabs(val) < tol) {
                printf("0      ");
            }
            else {
                printf("%.4f ",val);
            }
        }
        printf("\n");
    }
    printf("\n");

    printf("Q = \n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            val = Q[i][j];
            if (myabs(val) < tol) {
                printf("0      ");
            }
            else {
                printf("%.4f ",val);
            }
        }
        printf("\n");
    }
    printf("\n");

    printf("R = \n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            val = R[i][j];
            if (myabs(val) < tol) {
                printf("0      ");
            }
            else {
                printf("%.4f ",val);
            }
        }
        printf("\n");
    }
    printf("\n");

    printf("diff = \n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            val = A[i][j] - res[i][j];
            if (myabs(val) < tol) {
                printf("0      ");
            }
            else {
                printf("%.4f ",val);
            }
        }
        printf("\n");
    }
    printf("\n");
}