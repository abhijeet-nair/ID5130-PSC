#include <iostream>
#include <fstream>
#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

double myfunc (double x) {
    double f = sin(5*x);
    return f;
}

double myderv (double x) {
    double fdot = 5*cos(5*x);
    return fdot;
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

int main (int argc, char* argv[]) {
    int thrd_cnt = 1;

    if (argc == 2) {
        thrd_cnt = strtol(argv[1], NULL, 10);
    }
    else {
        printf("\n A command line argument other than name of the executable is required... Exiting the program...\n");
        return 1;
    }
    
    int N {}, i {}, j {};
    std::cout << "Enter number of elements (> 0) of the matrices: ";
    std::cin >> N;
    printf("\n");

    double** A = new double*[N];
    double yi_k[N] {};
    double fvec[N] {};
    double fdvec[N] {};
    double h = 3/(double(N-1));
    double xi {};

    for (i = 0; i < N; i++) {
        A[i] = new double[N] {};
    }

    #pragma omp parallel num_threads(thrd_cnt) default(none) shared(A, yi_k, fvec, fdvec, N, h) private(i, xi)
    {
        #pragma omp for
        for (i = 0; i < N; i++) {
            xi = i*h;
            fvec[i] = myfunc(xi);
            fdvec[i] = myderv(xi);
        }


        // Initializing A
        #pragma omp for
        for (i = 0; i < N; i++)  {
            if (i == 0) {
                A[i][0] = 1;
                A[i][1] = 2;

                yi_k[i] = (-2.5*fvec[0] + 2*fvec[1] + 0.5*fvec[2])/h;
            }
            else if (i == N - 1) {
                A[i][N-2] = 2;
                A[i][N-1] = 1;

                yi_k[i] = (2.5*fvec[N-1] - 2*fvec[N-2] - 0.5*fvec[N-3])/h;
            }
            else {
                A[i][i-1] = 1;
                A[i][i] = 4;
                A[i][i+1] = 1;

                yi_k[i] = 3*(fvec[i+1] - fvec[i-1])/h;
            }
        }
    }

    double ai_k[N] {};
    double bi_k[N] {};
    double ci_k[N] {};

    #pragma omp parallel for num_threads(thrd_cnt) default(none) shared(A, ai_k, bi_k, ci_k, N) private(i)
        for (i = 0; i < N - 1; i++) {
            ai_k[i+1] = A[i+1][i];
            bi_k[i] = A[i][i];
            ci_k[i] = A[i][i+1];
        }

    bi_k[N-1] = A[N-1][N-1];

    // Initializing variables for RD Alg
    int n_RD = ceil(log2(N));
    double ai_k1[N] {};
    double bi_k1[N] {};
    double ci_k1[N] {};
    double yi_k1[N] {};
    double alpi_k[N] {};
    double beti_k[N] {};
    int l1 {}, l2 {};
    double t {};

    // RD Algorithm
    t = omp_get_wtime();
    #pragma omp parallel num_threads(thrd_cnt) \
    shared(A, ai_k, bi_k, ci_k, yi_k, ai_k1, bi_k1, ci_k1, yi_k1, alpi_k, beti_k, N, n_RD) \
    private(i, l1, l2)
    {
        for (int k = 1; k <= n_RD; k++) {
            l1 = pow(2, k-1);
            l2 = pow(2, k);

            #pragma omp for
                for (i = 0; i < N; i++) {

                    bi_k1[i] = bi_k[i];
                    yi_k1[i] = yi_k[i];

                    if (i <= N - l1 - 1) {
                        beti_k[i] = -ci_k[i]/bi_k[i+l1];

                        bi_k1[i] += beti_k[i]*ai_k[i+l1];
                        yi_k1[i] += beti_k[i]*yi_k[i+l1];
                    }
                    else {
                        beti_k[i] = 0;
                    }

                    if (i >= l1) {
                        alpi_k[i] = -ai_k[i]/bi_k[i-l1];

                        bi_k1[i] += alpi_k[i]*ci_k[i-l1];
                        yi_k1[i] += alpi_k[i]*yi_k[i-l1];
                    }
                    else {
                        alpi_k[i] = 0;
                    }

                    if (i <= N - l2 - 1) {
                        ci_k1[i] = beti_k[i]*ci_k[i+l1];
                    }
                    else {
                        ci_k1[i] = 0;
                    }

                    if (i >= l2) {
                        ai_k1[i] = alpi_k[i]*ai_k[i-l1];
                    }
                    else {
                        ai_k1[i] = 0;
                    }
                }

            if (k <= n_RD - 1) {
                #pragma omp for
                    for (i = 0; i < N; i++) {
                        ai_k[i] = ai_k1[i];
                        bi_k[i] = bi_k1[i];
                        ci_k[i] = ci_k1[i];
                        yi_k[i] = yi_k1[i];
                    }
            }
        }
    }

    double xsol[N] {};
    
    #pragma omp parallel for num_threads(thrd_cnt) default(none) shared(xsol, yi_k1, bi_k1, N) private(i)
    for (i = 0; i < N; i++) {
        xsol[i] = yi_k1[i]/bi_k1[i];
    }
    t = omp_get_wtime() - t;

    printf("Completed...\n");
    printf("Time Taken = %.6f s\n",t);
    // printf("xi = \n");
    // printVector(xsol, N);

    // std::ofstream oFile("./Res/RDA.txt");

    // if (oFile.is_open()) {
    //     for (i = 0; i < N; i++) {
    //         oFile << xsol[i] << "," << fdvec[i] << "\n";
    //     }
    //     oFile.close();
    //     printf("Saved in file ./Res/RDA.txt\n");
    // }
    // else {
    //     printf("Error opening file\n");
    // }
}