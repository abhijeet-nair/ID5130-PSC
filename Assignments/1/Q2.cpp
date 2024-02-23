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

void printMatrix (double** a, int m, int n) {
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << a[i][j] << "\t";
        }
        std::cout << std::endl;
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

    double** A = new double*[N+1];
    double b[N+1] {};
    double fvec[N+1] {};
    double fdvec[N+1] {};
    double h = 3/(double(N));
    double xi {};

    for (i = 0; i <= N; i++) {
        A[i] = new double[N+1];

        xi = i*h;
        fvec[i] = myfunc(xi);
        fdvec[i] = myderv(xi);
    }

    for (i = 0; i <= N; i++)  {
        // printf("(%d,%d)\n",i,N);

        if (i == 0) {
            A[i][0] = 1;
            A[i][1] = 2;

            b[i] = (-2.5*fvec[0] + 2*fvec[2] + 0.5*fvec[2])/h;
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

    

    // printMatrix(A, N+1, N+1);
}