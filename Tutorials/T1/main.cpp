// #include <bits/stdc++.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "functions.h"

int main () {
    int N {};
    std::cout << "Enter number of elements (> 0) of the matrices: ";
    std::cin >> N;

    // if (N >= 10) {
    //     std::cout << "Invalid N!!!";
    // }
    // else {
        std::cout << "Valid N..." << std::endl;
        double const PI = 3.14;
        double** A = new double*[N];
        double** B = new double*[N];
        double* C = new double[N];

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
            C[i] = 0.5*i;
	    }

        /*
        std::cout << std::endl << "Question 1:" << std::endl;
        std::cout << "------------------" << std::endl;
        // std::cout << "Matrix 1 is:" << std::endl;
        // printMatrix(A, N, N);

        // std::cout << std::endl;
        // std::cout << "Matrix 2 is:" << std::endl;
        // printMatrix(B, N, N);
        */

        /*
        double** sum = addTwoMatrices(A, B, N, N);
        std::cout << "Sum is:" << std::endl;
        printMatrix(sum, N, N);
        // printMatrix(addTwoMatrices(A, B, N, N), N, N);

        double** prod = multiplyMatrices(A, B, N, N, N);
        std::cout << "Product is:" << std::endl;
        printMatrix(prod, N, N);
        */

        /*
        double** res1 = transposeMatrix(multiplyMatrices(A, B, N, N, N), N, N);
        double** res2 = multiplyMatrices(transposeMatrix(B, N, N), transposeMatrix(A, N, N), N, N, N);

        std::cout << std::endl;
        int stat1 = compareTwoMatrices(res1, res2, N, N);
        
        if (stat1 == 1){std::cout << "Matrices are same!!!" << std::endl;}
        else {std::cout << "Matrices are not same!!!" << std::endl;}

        std::cout << std::endl << "Question 2:" << std::endl;
        std::cout << "------------------" << std::endl;

        double** sum = addTwoMatrices(A, transposeMatrix(A, N, N), N, N);
        int stat2 = compareTwoMatrices(sum, transposeMatrix(sum, N, N), N, N);
        
        if (stat2 == 1){std::cout << "Matrices are same!!!" << std::endl;}
        else {std::cout << "Matrices are not same!!!" << std::endl;}
        */

        std::cout << std::endl << "Question 3:" << std::endl;
        std::cout << "------------------" << std::endl;

        double T {};

        for(int i = 0; i < 10; i++) {
            clock_t t;
            t = clock();
            double* prod1 = multiplyMatVec(A, C, N, N);
            t = clock()-t;
            T += ((float)t)/CLOCKS_PER_SEC;
            std::cout << "i = " << i << ", T = " << ((float)t)/CLOCKS_PER_SEC << std::endl;
        }
        T = T/10;
        std::cout << "Time taken for N = " << N << " is: " << T << " sec..." << std::endl;
        // t = ((float)t)/CLOCKS_PER_SEC;
        // std::cout << "Result is:" << std::endl;
        // printVector(prod1, N);
        // std::cout << "Time taken for N = " << N << " is: " << t << " sec..." << std::endl;
        // printf ("Multiplication took %d clicks (%f seconds).\n\n",t,((float)t)/CLOCKS_PER_SEC);
    // }
    return 0;
}