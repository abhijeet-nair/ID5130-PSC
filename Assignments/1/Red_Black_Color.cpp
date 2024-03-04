#include <iostream>
#include <fstream>
#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

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

    int N {};
    int i {}, j {}, it {};
    double del {};

    std::cout << "Enter grid spacing: ";
    std::cin >> del;
    printf("\n");

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

    #pragma omp parallel for num_threads(thrd_cnt) collapse(2) default(none) shared(qij, solMat, N, xi, yi) private(i, j)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                qij[i][j] = q(xi[i], yi[j]);
                solMat[i][j] = phiSol(xi[i], yi[j]);
            }
        }

    // printf("solMat = \n");
    // printMatrix(solMat,N,N);

    double myerr = 1, eps = 1e-2;
    int cnt = 1;
    int N2 = N*N;

    int lim = 1e7;
    double errvec[N*N] {};
    double numSolVec[N] {};
    double actSolVec[N] {};

    double errtmp {};
    double t = omp_get_wtime();
    #pragma omp parallel num_threads(thrd_cnt) default(none) shared(phik, phik1, qij, del2, N, N2, lim, solMat, myerr, errtmp, errvec, eps, cnt) private(i, j)
    {
        while ((myerr > eps) && (cnt < lim)) {
            #pragma omp for collapse(2)
                for (i = 1; i < N - 1; i++) {
                    for (j = 1; j < N - 1; j++) {
                        if ((i + j) % 2 == 1) {
                            phik1[i][j] = 0.25*(phik[i+1][j] + phik[i-1][j] + phik[i][j+1] + phik[i][j-1] + del2*qij[i][j]);
                            // phik1[i][j] = 0.25*(phik[i+1][j] + phik1[i-1][j] + phik[i][j+1] + phik1[i][j-1] + del2*qij[i][j]);
                        }
                    }
                }

            #pragma omp for collapse(2)
                for (i = 1; i < N - 1; i++) {
                    for (j = 1; j < N - 1; j++) {
                        if ((i + j) % 2 == 0) {
                            phik1[i][j] = 0.25*(phik1[i+1][j] + phik1[i-1][j] + phik1[i][j+1] + phik1[i][j-1] + del2*qij[i][j]);
                            // phik1[i][j] = 0.25*(phik[i+1][j] + phik1[i-1][j] + phik[i][j+1] + phik1[i][j-1] + del2*qij[i][j]);
                        }
                    }
                }
            
            // errtmp = 0.0;
            #pragma omp for collapse(2) // reduction(+:errtmp)
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        errvec[N*i + j] = phik1[i][j] - solMat[i][j];
                        // errtmp += pow(phik1[i][j] - solMat[i][j],2);
                        phik[i][j] = phik1[i][j];
                    }
                }
            
            errtmp = 0.0;
            // int N2 = N*N;
            #pragma omp barrier
            #pragma omp for reduction(+:errtmp)
                for (i = 0; i < N2; i++) {
                    errtmp += pow(errvec[i],2);
                }

            #pragma omp barrier

            if (omp_get_thread_num() == 0) {
                myerr = sqrt(errtmp);
                if (cnt % 500 == 0) {
                    printf("cnt = %d  err = %.4f\n",cnt,myerr);
                }
                cnt += 1;
            }

            #pragma omp barrier

            // printf("cnt = %d  err = %.2f\n",cnt,err);
            // #pragma omp single
            // {
            //     // err = norm(errvec, N*N);
            //     myerr = sqrt(errtmp);
            //     if (cnt % 5 == 0) {
            //         printf("cnt = %d  err = %.4f\n",cnt,myerr);
            //     }
            //     cnt += 1;
            // }
        }
    }
    t = omp_get_wtime() - t;

    if (myerr > eps) {
        printf("\nCrossed iteration limit of 1e%2.0f\n",log10(lim));
    }
    else {
        printf("\nConverged to required tolerance\nError = %.3f\nNo. of iterations = %d\n",myerr,cnt);
        printf("Time Taken = %.6f s\n",t);

        double errvec[N*N];
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                errvec[N*i + j] = phik1[i][j] - solMat[i][j];
            }
        }
        printf("ERR = %.4f\n",norm(errvec, N*N));

        int yInd = int(0.5*N);

        #pragma omp parallel for num_threads(thrd_cnt) default(none) \
        shared(numSolVec, actSolVec, phik, solMat, N, yInd) private(i)
            for (i = 0; i < N; i++) {
                numSolVec[i] = phik[i][yInd];
                actSolVec[i] = solMat[i][yInd];
            }
    }

    // std::ofstream oFile("./Res/GS_RBC.txt");

    // if (oFile.is_open()) {
    //     for (i = 0; i < N; i++) {
    //         oFile << numSolVec[i] << "," << actSolVec[i] << "\n";
    //     }
    //     oFile.close();
    //     printf("Saved in file ./Res/GS_RBC.txt\n");
    // }
    // else {
    //     printf("Error opening file\n");
    // }

    return 0;
}