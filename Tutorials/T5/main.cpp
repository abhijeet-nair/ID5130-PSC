#include <iostream>
#include <math.h>
#include <time.h>
#include <iomanip>
#ifdef _OPENMP
    #include <omp.h>
#endif


double norm (double A[], int n) {
    double res {};
    for (int i =0; i < n; i++) {
        res += pow(A[i],2);
    }
    res = sqrt(res);
    return res;
}

void printMatrix (double** a, int m, int n) {
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << a[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void printVector (double a[], int n) {
    std::cout << std::fixed;
    for (int i = 0; i < n; i++) {
        std::cout << a[i] << std::endl;
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

    double** A = new double*[N];

    for (i = 0; i < N; i++) {
        A[i] = new double[N];
    }

    // printf("Defining A...\n");

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++) {
            // printf("(i,j) = (%d,%d)\n",i,j);
            if (i == j) {
                A[i][j] = i + j + 2;
            }
            else if (i == 0 && j == N - 1) {
                A[i][j] = 1;
            }
            else if (i == N - 1 && j == 0) {
                A[i][j] = 2*N - 1;
            }
            else {
                A[i][j] = 1/(double(N));
            }
        }
    }
    double b[N] {};
    b[0] = 1;

    // printf("Matrix A is:\n");
    // printMatrix(A, N, N);
    // printf("Vector b is:\n");
    // printVector(b, N);
    // printf("\n");

    double xk[N] {}, xk_1[N] {};

    double err = 1.0;
    double eps = 1e-6;
    double sum {};
    double errvec[N] {};
    int cnt {};

    // printf("All matrices defined...\n");
    for (i = 0; i < N; i++) {
        xk_1[i] = 0.5*i;
    }
    
    printf("Question 1: Serial code for Jacobi Iteration:\n\n");
    double t_ser = omp_get_wtime();
    while (err > eps)
    {
       cnt += 1;
       if (cnt % 10 == 0) {printf("Count = %d\n", cnt);}

       // Transferring the old xk+1 to new xk
        for (i = 0; i < N; i++) {
            xk[i] = xk_1[i];
        }
        // printf("Transfer complete...\n");

        // Finding new xk+1
        for (i = 0; i < N; i++) {
            sum = 0.0;
            for (j = 0; j < N; j++) {
                if (j != i) {
                    sum += A[i][j]*xk[j];
                }
            }
            xk_1[i] = (b[i] - sum)/A[i][i];
        }
        // printf("Found new xk+1...\n");

        // Evaluating new error
        for (i = 0; i < N; i++) {
            errvec[i] = xk_1[i] - xk[i];
        }
        err = norm(errvec, N)/norm(xk, N);
        // printf("Evaluated error...\n");
    }
    t_ser = omp_get_wtime() - t_ser;
    
    printf("\nTotal count = %d\n", cnt);
    // printf("Solution x is:\n");

    double solSer[N] {};
    for (i = 0; i < N; i++){
        // printf("%.4f\n", xk_1[i]);
        solSer[i] = xk_1[i];
    }


    printf("\nQuestion 2: Parallel code for Jacobi Iteration:\n\n");
    err = 1.0;
    cnt = 0;

    for (i = 0; i < N; i++) {
        xk_1[i] = 0.5*i;
    }

    double t_par = omp_get_wtime();
    #pragma omp parallel num_threads(thrd_cnt) default(none) shared(A, b, xk, xk_1, N, errvec, cnt, eps, err) private(i, j, sum)
    {
        while (err > eps)
        {
            #pragma omp single
            {
                cnt += 1;
                if (cnt % 10 == 0) {printf("Count = %d\n", cnt);}
            }

            // Transferring the old xk+1 to new xk
            #pragma omp for
                for (i = 0; i < N; i++) {
                    xk[i] = xk_1[i];
                }
            // printf("Transfer complete...\n");

            // Finding new xk+1
            #pragma omp for
                for (i = 0; i < N; i++) {
                    sum = 0.0;
                    for (j = 0; j < N; j++) {
                        if (j != i) {
                            sum += A[i][j]*xk[j];
                        }
                    }
                    xk_1[i] = (b[i] - sum)/A[i][i];
                }
            // printf("Found new xk+1...\n");

            // Evaluating new error
            #pragma omp for
                for (i = 0; i < N; i++) {
                    errvec[i] = xk_1[i] - xk[i];
                }

            #pragma omp single
                err = norm(errvec, N)/norm(xk, N);
            // printf("Evaluated error...\n");
        }
    }
    t_par = omp_get_wtime() - t_par;

    printf("\nTotal count = %d\n", cnt);
    // printf("Solution x is:\n");

    double solPar[N] {};
    for (i = 0; i < N; i++){
        // printf("%.4f\n", xk_1[i]);
        solPar[i] = xk_1[i];
        errvec[i] = solPar[i] - solSer[i];
    }

    err = norm(errvec, N);

    printf("\nResults:\n");
    printf("Serial Time = %.6f s\n", t_ser);
    printf("Parallel Time = %.6f s\n", t_par);
    printf("Time Saved: %.6f s (%3.2f\%)\n", t_ser - t_par, (t_ser - t_par)*100/t_ser);
    printf("Norm(error) among two = %.6f\n", err);

    return 0;
}