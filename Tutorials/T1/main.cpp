#include <bits/stdc++.h>
#include "functions.h"

using namespace std;

int main () {
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

        for (int i = 0; i < N; i++) {
		    A[i] = new double[N];
		    B[i] = new double[N];
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
        printMatrix(A, N, N);

        cout << endl;
        cout << "Matrix 2 is:" << endl;
        printMatrix(B, N, N);

        /*
        double** sum = addTwoMatrices(A, B, N, N);
        cout << "Sum is:" << endl;
        printMatrix(sum, N, N);
        // printMatrix(addTwoMatrices(A, B, N, N), N, N);

        double** prod = multiplyMatrices(A, B, N, N, N);
        cout << "Product is:" << endl;
        printMatrix(prod, N, N);
        */

        double** res1 = transposeMatrix(multiplyMatrices(A, B, N, N, N), N, N);
        double** res2 = multiplyMatrices(transposeMatrix(B, N, N), transposeMatrix(A, N, N), N, N, N);

        cout << endl;
        int stat1 = compareTwoMatrices(res1, res2, N, N);
        
        if (stat1 == 1){cout << "Matrices are same!!!" << endl;}
        else {cout << "Matrices are not same!!!" << endl;}

        cout << endl << "Question 2:" << endl;
        cout << "------------------" << endl;

        double** sum = addTwoMatrices(A, transposeMatrix(A, N, N), N, N);
        int stat2 = compareTwoMatrices(sum, transposeMatrix(sum, N, N), N, N);
        
        if (stat2 == 1){cout << "Matrices are same!!!" << endl;}
        else {cout << "Matrices are not same!!!" << endl;}

        cout << endl << "Question 3:" << endl;
        cout << "------------------" << endl;

        
    }
    return 0;
}