#include <iostream>
#include <math.h>
#include <cstdio>
using namespace std;

void multiplyTwoMatrices (double a[][100], double b[][100], double c[][100], int m, int n, int p) {
	for (int i = 1; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
			c[i][j] = 0;
			for (int k = 0; k < p; k ++) {
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return;
}

void printMatrix (double** arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << arr[i][j] << "\t";
        }
        cout << "\n";
    }
    return;
}

int main() {
    int N;
    cout << "Enter number of elements (< 10) of the matrices: ";
    cin >> N;

    if (N >= 10)
    {
        cout << "N Not Valid!!!";
    }
    else
    {
        double const PI = 3.14;
        double** mat1 = new double*[N];
        double** mat2 = new double*[N];

        for (int i = 0; i<N; i++) {
		mat1[i] = new double[N];
		mat2[i] = new double[N];
	}

        for (int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++) {
                mat1[i][j] = pow(0.5,0.5*i) * sin(i*j*PI/(N+1));
                mat2[j][i] = pow(0.5,0.5*i) * cos(i*j*PI/(N+1));
            }
        }

        printMatrix(mat1, N);
        printMatrix(mat2, N);
    }
}