#include <iostream>
#include <math.h>
#include <time.h>
#include <iomanip>
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
   std::cout << std::fixed << std::setprecision(4);
   for (int i = 0; i < n; i++) {
        std::cout << a[i] << std::endl;
   }
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


int mai_kn (int argc, char* argv[]) {
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

    double** A = new double*[N+1];
    double b[N+1] {};
    double fvec[N+1] {};
    double fdvec[N+1] {};
    double h = 3/(double(N));
    double xi {};

    for (i = 0; i <= N; i++) {
        A[i] = new double[N+1] {};

        xi = i*h;
        fvec[i] = myfunc(xi);
        fdvec[i] = myderv(xi);
    }

    // Initializing A
    for (i = 0; i <= N; i++)  {
        // printf("(%d,%d)\n",i,N);

        if (i == 0) {
            A[i][0] = 1;
            A[i][1] = 2;

            b[i] = (-2.5*fvec[0] + 2*fvec[1] + 0.5*fvec[2])/h;
        }
        else if (i == N) {
            A[i][N-1] = 2;
            A[i][N] = 1;

            b[i] = (2.5*fvec[N] - 2*fvec[N-1] - 0.5*fvec[N-2])/h;
        }
        else {
            A[i][i-1] = 1;
            A[i][i] = 4;
            A[i][i+1] = 1;

            b[i] = 3*(fvec[i+1] - fvec[i-1])/h;
        }
    }

    double ai_k[N] {};
    double bi_k[N+1] {};
    double ci_k[N] {};

    for (i = 0; i < N; i++) {
        ai_k[i] = A[i+1][i];
        bi_k[i] = A[i][i];
        ci_k[i] = A[i][i+1];
    }
    bi_k[N] = A[N][N];

    // Serial LU Decomposition
    // Finding LU
    double u[N+1] {};
    // double l[N+1] {}; // Will not use l[0]. Do something later.
    double l[N] {};

    u[0] = b[0];

    for (i = 0; i < N; i++) {
        l[i] = ai_k[i]/u[i];  // l_{i+1} = a_{i+1}/u_i = A_{i+1,i}/u_i
        u[i+1] = bi_k[i+1] - l[i]*ci_k[i];
        // printf("l(%d), u(%d) = %.4f, %.4f\n",i,i+1,l[i],u[i+1]);
    }

    // Forward substitution for z
    double z[N+1] {};
    z[0] = b[0];

    for (i = 0; i < N; i++) {
        z[i+1] = b[i+1] - l[i]*z[i];
    }

    // printf("\nf(x) = \n");
    // printVector(fvec,N+1);
    // printf("\nb = \n");
    // printVector(b,N+1);
    // printf("\nz = \n");
    // printVector(z,N+1);

    double x[N+1] {};
    x[N] = z[N]/u[N];

    for (i = N - 1; i >= 0; i--) {
        x[i] = (z[i] - ci_k[i]*x[i+1])/u[i];
    }
    // printf("\nx = \n");
    // printVector(x,N+1);

    // std::cout << std::fixed << std::setprecision(4);
    // for (int i = 0; i <= N; i++) {
    //     std::cout << x[i] << "\t" << fdvec[i] << std::endl;
    // }

    // printMatrix(A, N+1, N+1);


    // Parallel Recursive Doubling
    int n_RD = ceil(log2(N));
    // printf("n_RD = %d\n",n_RD);

    
    double ai_kp1[N] {};
    double bi_kp1[N+1] {};
    double ci_kp1[N] {};
    double alp_i_k {};
    double bet_i_k {};

    for (int k = 1; k <= n_RD; k++) {
        for (int i = 0; i <= N; i++) {
            
        }
    }


    return 0;
}