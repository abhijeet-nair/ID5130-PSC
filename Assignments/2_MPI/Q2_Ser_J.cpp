#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

// RHS function
double q (double x, double y) {
    double val = pow(x, 2) + pow(y, 2);
    return val;
}

// Norm function
double norm (double A[], int n) {
    double res {};
    for (int i =0; i < n; i++) {
        res += pow(A[i],2);
    }
    res = sqrt(res);
    return res;
}


int main (int argc, char* argv[]) {
    int i, j;
    double del = 0.1, del2 = pow(del, 2);

    int N = int(2/del) + 1;

    double xi[N] {};
    double yi[N] {};

    double** phik = new double*[N];
    double** phik1 = new double*[N];
    double** qij = new double*[N];

    for (i = 0; i < N; i++) {
        phik[i] = new double[N] {};
        phik1[i] = new double[N] {};
        qij[i] = new double[N] {};

        xi[i] = -1 + i*del;
        yi[i] = -1 + i*del;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            qij[i][j] = q(xi[i], yi[j]);
        }
    }

    double err = 1, eps = 1e-4;
    double errvec[N*N];
    int cnt = 1;
    int lim = 1e7;

    clock_t t;
    t = clock();

    while ((err > eps) && (cnt < lim)) {
        for (i = 1; i < N-1; i++) {
            phik1[0][i] = sin(2*M_PI*yi[i]);
            phik1[i][0] = 0;
            phik1[i][N-1] = 0;
            for (j = 1; j < N-1; j++) {
                phik1[i][j] = 0.25*(phik[i+1][j] + phik[i-1][j] + phik[i][j+1] + phik[i][j-1] + del2*qij[i][j]);
            }
        }

        for (i = 1; i < N-1; i++) {
            phik1[N-1][i] = (4*phik1[N-2][i] - phik1[N-3][i])/3;
        }

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                errvec[N*i + j] = phik1[i][j] - phik[i][j];
                phik[i][j] = phik1[i][j];
            }
        }
        
        err = norm(errvec, N*N);
        if (cnt % 500 == 0) {printf("cnt = %d  err = %.6f\n",cnt,err);}
        cnt += 1;
    }
    t = clock() - t;
    double totT = double(t) / double(CLOCKS_PER_SEC);

    printf("\ncnt = %d\terr = %.6f\tTime taken = %.4f s\n\n",cnt-1,err,totT);

    double phivsx0[N] {}, phivsy0[N] {};
    int rind = int(0.5*N); // Index for (N/2)+1-th element

    for (i = 0; i < N; i++) {
        phivsx0[i] = phik1[i][rind];
        phivsy0[i] = phik1[rind][i];
    }

    // // Saving in a file for plotting
    // char fname[25];
    // sprintf(fname, "./Res/Q2_Ser_J_e%d.txt",-int(log10(del)));
    
    // std::ofstream oFile(fname);

    // if (oFile.is_open()) {
    //     for (i = 0; i < N; i++) {
    //         oFile << phivsx0[i] << "," << phivsy0[i] << "\n";
    //     }

    //     oFile.close();
    //     printf("Saved in file %s\n",fname);
    // }
    // else {
    //     printf("Error opening file\n");
    // }

    delete phik1;
    delete phik;
    delete qij;

    return 0;
}