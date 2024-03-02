#include <iostream>
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
    int thrd_cnt = 1;

    if (argc == 2) {
        thrd_cnt = strtol(argv[1], NULL, 10);
    }
    else {
        printf("\n A command line argument other than name of the executable is required... Exiting the program...\n");
        return 1;
    }

    int N {}, Nd {};
    int i {}, j {}, l {};
    int ibeg {}, iend {};
    double del {};

    std::cout << "Enter grid spacing: ";
    std::cin >> del;

    double del2 = pow(del, 2);

    // std::cout << "Enter y-grid spacing: ";
    // std::cin >> dely;

    N = int(2/del) + 1;
    Nd = 2*N - 1;
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

    // Diagonal numbered from 1 to 2N - 1
    while ((err > eps) && (cnt < lim)) {
        for (l = 2; l <= Nd - 1; l++) {
            if (l < N) {
                ibeg = 1;
                iend = l - 1;
            }
            else {
                ibeg = l - N + 2;
                iend = N - 2;
            }

            for (i = ibeg; i <= iend; i++) {
                j = l - i;

                // printf("l, i, j, ib, ie, Nd, N = %d, %d, %d, %d, %d, %d, %d\n",l,i,j,ibeg,iend,Nd,N);
                
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

    if (err > eps) {
        printf("Crossed iteration limit of 1e%d\n",log10(lim));
    }
    else {
        printf("Converged to required tolerance\nNo. of iterations = %d\n",cnt);
        int yInd = int(0.5*N);

        double solVec[N] {};

        for (i = 0; i < N; i++) {
            solVec[i] = phik[i][yInd];
        }
    }

    return 0;
}