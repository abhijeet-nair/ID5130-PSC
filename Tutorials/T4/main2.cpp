#include <iostream>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

double** addTwoMatrices (double** a, double** b, int m, int n) {
    double** sum = new double*[m];
    for (int i = 0; i < m; i++) {
        sum[i] = new double[n];
    }

    for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
            sum[i][j] = a[i][j] + b[i][j];
        }
    }
    return sum;
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

    printf("Question 1:\n------------------------------------------\n");
    int N {};
    std::cout << "Enter number of elements (> 0) of the matrices: ";
    std::cin >> N;

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
        }
    }

    double** sum = addTwoMatrices(A, B, N, N);

    return 0;
}