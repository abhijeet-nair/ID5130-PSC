#include <iostream>
#include <math.h>
using namespace std;

void printMatrix (double** arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << arr[i][j] << "\t";
        }
        cout << "\n";
    }
}

void addTwoMatrices (double** a, double** b, double** sum, int m, int n) {
    for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
            sum[i][j] = a[i][j] + b[i][j];
        }
    }
}

void multiplyTwoMatrices (double** a, double** b, double** c, int m, int n, int p) {
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
			c[i][j] = 0;
			for (int k = 0; k < p; k ++) {
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
}

void transposeMatrix (double** a, double** aT, int m, int n) {
    for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
            aT[i][j] = a[j][i];
        }
    }
}

void compareTwoMatrices (double** a, double** b, int m, int n) {
    // int cnt {0};
    cout << "The given matrices are ";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (a[i][j] != b[i][j]) {                
                cout << "not same..." << endl;
                return;
            }
        }
    }
    cout << "same..." << endl;
}

int main() {
    int N {};
    cout << "Enter number of elements (< 10) of the matrices: ";
    cin >> N;

    if (N >= 10) {
        cout << "Invalid N!!!";
    }
    else {
        cout << "Valid N..." << endl;
        double const PI = 3.14;
        double** A = new double*[N];
        double** B = new double*[N];
        double** AB = new double*[N];

        double** A_T = new double*[N];
        double** B_T = new double*[N];
        double** AB_T = new double*[N];
        double** B_TA_T = new double*[N];
        double** ApA_T = new double*[N];
        double** ApA_T_T = new double*[N];

        for (int i = 0; i < N; i++) {
		    A[i] = new double[N];
		    B[i] = new double[N];
            AB[i] = new double[N];

            A_T[i] = new double[N];
		    B_T[i] = new double[N];
            AB_T[i] = new double[N];
            B_TA_T[i] = new double[N];
            ApA_T[i] = new double[N];
            ApA_T_T[i] = new double[N];
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = pow(0.5,0.5*i) * sin(i*j*PI/(N+1));
                B[j][i] = pow(0.5,0.5*i) * cos(i*j*PI/(N+1));
                // A[i][j] = 1;
                // B[j][i] = 1;
            }
	    }

        cout << endl << "Question 1:" << endl;
        cout << "------------------" << endl;
        cout << "Matrix 1 is:" << endl;
        printMatrix(A, N);

        cout << endl;
        cout << "Matrix 2 is:" << endl;
        printMatrix(B, N);

        // cout << endl << "Multiplying A*B..." << endl;
        multiplyTwoMatrices(A, B, AB, N, N, N);
        transposeMatrix(AB, AB_T, N, N);
        // cout << "Product is..." << endl;
        // printMatrix(AB, N);

        transposeMatrix(A, A_T, N, N);
        transposeMatrix(B, B_T, N, N);
        multiplyTwoMatrices(B_T, A_T, B_TA_T, N, N, N);

        cout << endl;
        compareTwoMatrices(AB_T, B_TA_T, N, N);

        cout << endl << "Question 2:" << endl;
        cout << "------------------" << endl;

        addTwoMatrices(A, A_T, ApA_T, N, N);
        transposeMatrix(ApA_T, ApA_T_T, N, N);
        compareTwoMatrices(ApA_T, ApA_T_T, N, N);

        cout << endl << "Question 3:" << endl;
        cout << "------------------" << endl;

    }

    cout << endl;
    return 0;
}