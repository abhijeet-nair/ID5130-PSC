#include "functions.h"
#include <iostream>
#include <iomanip>

void printMatrix (double** a, int m, int n) {
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << a[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void printVector (double* a, int m) {
    // std::cout << std::setprecision(4);
    for (int i = 0; i < m; i++) {
        std::cout << a[i] << std::endl;
    }
}

double** addTwoMatrices (double** a, double** b, int m, int n) {
    double** sum = new double*[m];
    for (int i = 0; i < m; i++) {
        sum[i] = new double[n];
    }

    for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
            sum[i][j] = a[i][j] + b[i][j];
        }
    }
    return sum;
}

double** multiplyMatrices (double** a, double** b, int m, int n, int p) {
    double** prod = new double*[m];
    for (int i = 0; i < m; i++) {
        prod[i] = new double[p];
    }

    for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
			prod[i][j] = 0;
			for (int k = 0; k < p; k ++) {
				prod[i][j] += a[i][k]*b[k][j];
			}
		}
	}
    return prod;
}

double** transposeMatrix (double** a, int m, int n) {
    double** aT = new double*[m];
    for (int i = 0; i < m; i++) {
        aT[i] = new double[n];
    }
    for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
            aT[i][j] = a[j][i];
        }
    }
    return aT;
}

int compareTwoMatrices (double** a, double** b, int m, int n) {
    // int cnt {0};
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (a[i][j] != b[i][j]) {                
                return 0;
            }
        }
    }
    return 1;
}

double* multiplyMatVec (double** a, double* b, int m, int n) {
    double* prod = new double[m] {};

    for (int i = 0; i < m; i ++) {
        for (int j = 0; j < n; j ++) {
            prod[i] += a[i][j]*b[j];
        }
	}
    return prod;
}