#include <iostream>
#include <fstream>
#include <math.h>

double q (double x, double y) {
    double val = 2*(2 - pow(x, 2) - pow(y, 2));
    return val;
}

double phiSol (double x, double y) {
    double phi = (pow(x, 2) - 1) * (pow(y, 2) - 1);
    return phi;
}

void printVector(double a[], int n) {
   for (int i = 0; i < n; i++) {
        printf("%.4f\n",a[i]);
   }
}

void printMatrix (double** a, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.4f\t",a[i][j]);
        }
        printf("\n");
    }
}

double norm (double A[], int n) {
    double res {};
    for (int i =0; i < n; i++) {
        res += pow(A[i],2);
    }
    res = sqrt(res);
    return res;
}


int main (int argc, char* argv[]) {
    int N {}, i {}, j {};
    double del {};

    std::cout << "Enter grid spacing: ";
    std::cin >> del;

    double del2 = pow(del, 2);

    // std::cout << "Enter y-grid spacing: ";
    // std::cin >> dely;

    N = int(2/del) + 1;
    // N = del;

    double xi[N] {};
    double yi[N] {};

    double** phik = new double*[N];
    double** phik1 = new double*[N];
    double** solMat = new double*[N];
    double** qij = new double*[N];

    for (i = 0; i < N; i++) {
        phik[i] = new double[N] {};
        phik1[i] = new double[N] {};
        qij[i] = new double[N] {};
        solMat[i] = new double[N] {};

        xi[i] = -1 + i*del;
        yi[i] = -1 + i*del;
    }

    // printf("xi = \n");
    // printVector(xi,N);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            qij[i][j] = q(xi[i], yi[j]);
            solMat[i][j] = phiSol(xi[i], yi[j]);
        }
    }

    // printf("solMat = \n");
    // printMatrix(solMat,N,N);

    double err = 1, eps = 1e-2;
    double errvec[N*N];
    int cnt = 1;

    int lim = 1e5;

    while ((err > eps) && (cnt < lim)) {      
        for (i = 1; i < N-1; i++) {
            for (j = 1; j < N-1; j++) {
                phik1[i][j] = 0.25*(phik[i+1][j] + phik1[i-1][j] + phik[i][j+1] + phik1[i][j-1] + del2*qij[i][j]);
            }
        }

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                errvec[N*i + j] = phik1[i][j] - solMat[i][j];
                phik[i][j] = phik1[i][j];
            }
        }
        
        err = norm(errvec, N*N);
        if (cnt % 10 == 0) {
            printf("cnt = %d  err = %.2f\n",cnt,err);
        }
        cnt += 1;
    }

    double numSolVec[N] {};
    double actSolVec[N] {};

    if (err > eps) {
        printf("Crossed iteration limit of 1e%d\n",log10(lim));
    }
    else {
        printf("Converged to required tolerance\nNo. of iterations = %d\n",cnt);
        int yInd = int(0.5*N);

        for (i = 0; i < N; i++) {
            numSolVec[i] = phik[i][yInd];
            actSolVec[i] = solMat[i][yInd];
        }
    }

    std::ofstream oFile("./Res/GS_Ser.txt");

    if (oFile.is_open()) {
        for (i = 0; i < N; i++) {
            oFile << numSolVec[i] << "," << actSolVec[i] << "\n";
        }
        oFile.close();
        printf("Saved in file ./Res/GS_Ser.txt\n");
    }
    else {
        printf("Error opening file\n");
    }

    return 0;
}