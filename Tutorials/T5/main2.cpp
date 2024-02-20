#include <iostream>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

// void printMatrix (double** a, int m, int n) {
//     std::cout << std::fixed; // << std::setprecision(4);
//     for (int i = 0; i < m; i++) {
//         for (int j = 0; j < n; j++) {
//             std::cout << a[i][j] << "\t";
//         }
//         std::cout << std::endl;
//     }
// }

double** addTwoMatrices (double** a, double** b, int m, int n, int thrd_cnt, double& t) {
    double** sum = new double*[m];
    for (int i = 0; i < m; i++) {
        sum[i] = new double[n];
    }

    int i {}, j {};

    t = omp_get_wtime();
    #pragma omp parallel for collapse(2) default(none) shared(a, b, sum, m, n) private(i, j) num_threads(thrd_cnt)
        for (i = 0; i < m; i ++) {
            for (j = 0; j < n; j ++) {
                sum[i][j] = a[i][j] + b[i][j];
            }
        }
    
    t = omp_get_wtime() - t;
    return sum;
}


double** multTwoMatrices (double** a, double** b, int m, int n, int p, int thrd_cnt, double& t) {
    double** prod = new double*[m];
    for (int i = 0; i < m; i++) {
        prod[i] = new double[p] {};
    }

    // printMatrix(prod, m, p);
    int i {}, j {}, k {};

    t = omp_get_wtime();
    #pragma omp parallel for collapse(2) default(none) shared(a, b, prod, m, n, p) private(i,j, k) num_threads(thrd_cnt)
        for (i = 0; i < m; i ++) {
            for (j = 0; j < n; j ++) {
                prod[i][j] = 0;
                for (k = 0; k < p; k ++) {
                    prod[i][j] += a[i][k]*b[k][j];
                }
            }
        }

    t = omp_get_wtime() - t;
    return prod;
}

// Question 3 Functions
double myfunc (double x) {return 7.0 - x*tan(x);}

double myderv (double x) {return -tan(x) - x/pow(cos(x),2);}

double forward (double fi1, double fi, double delx) {return (fi1 - fi)/delx;}

double backward (double fi, double fi_1, double delx) {return (fi - fi_1)/delx;}

double central2nd (double fi1, double fi_1, double delx) {return 0.5*(fi1 - fi_1)/delx;}

double central4th (double fi2, double fi1, double fi_1, double fi_2, double delx) {return (-fi2 + 8*fi1 - 8*fi_1 + fi_2)/(12*delx);}


int main (int argc, char* argv[]) {
    int thrd_cnt = 1;

    if (argc == 2) {
        thrd_cnt = strtol(argv[1], NULL, 10);
    }
    else {
        printf("\n A command line argument other than name of the executable is required... Exiting the program...\n");
        return 1;
    }

    // printf("Question 1:\n------------------------------------------\n");
    // int N {};
    // std::cout << "Enter number of elements (> 0) of the matrices: ";
    // std::cin >> N;

    // double const PI = 3.14;
    // double** A = new double*[N];
    // double** B = new double*[N];
    // int i {};

    // for (i = 0; i < N; i++) {
    //     A[i] = new double[N];
    //     B[i] = new double[N];
    // }

    // for (i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         A[i][j] = pow(0.5,0.5*i) * sin(i*j*PI/(N+1));
    //         B[j][i] = pow(0.5,0.5*i) * cos(i*j*PI/(N+1));
    //     }
    // }

    // // // double t1 = omp_get_wtime();
    // double tadd {}, tmp {};
    // double** sum;
    // for (i = 0; i < 10; i++) {
    //     sum = addTwoMatrices(A, B, N, N, thrd_cnt, tmp);
    //     tadd += tmp;
    //     // std::cout << tmp << std::endl;
    // }
    // tadd *= 0.1;
    // // double t2 = omp_get_wtime();
    // printf("\nAddition Done!!!\nTime taken: %.8f sec\n",tadd);

    
    // printf("\n\nQuestion 2:\n------------------------------------------\n");

    // double tmult {};
    // tmp = 0.0;
    // double** prod;
    // for (i = 0; i < 10; i++) {
    //     prod = multTwoMatrices(A, B, N, N, N, thrd_cnt, tmp);
    //     tmult += tmp;
    // }
    // tmult *= 0.1;
    // printf("Multiplication Done!!!\nTime taken: %.8f sec\n",tmult);

    
    printf("\n\nQuestion 3:\n------------------------------------------\n");
    double delx = 0.01;
    int Ngrid = int(2/delx) + 1;
    
    // LESSON LEARNT: ARRAY WAS DEFINED WITH ONE LESS ELEMENT.
    // DEFINED WITH Ngrid INSTEAD OF Ngrid + 1.
    // CHANGED Ngrid DEFINITION NOW.

    double grid[Ngrid] {};
    // double grid {};
    double funcvec[Ngrid] {};
    double dervvec[Ngrid] {};
    int i {};

    // Initializing the grid and function values
    #pragma omp parallel for num_threads(thrd_cnt) default(none) shared(grid, funcvec, dervvec, delx, Ngrid) private(i)
        for (i = 0; i < Ngrid; i++) {
            grid[i] = -1 + delx*i;
            funcvec[i] = myfunc(grid[i]);
            // printf("%d\t%.4f\n",i,funcvec[i]);
            dervvec[i] = myderv(grid[i]);
        }
    
    // printf("%.4f\t%.4f\n",funcvec[Ngrid-1],funcvec[Ngrid]);

    // 2nd Order Central Difference Scheme
    double dervvec2FD[Ngrid] {};
    double t2FD = omp_get_wtime();
    #pragma omp parallel for num_threads(thrd_cnt) default(none) shared(dervvec2FD, funcvec, delx, Ngrid) private(i)
        for (i = 0; i < Ngrid; i++) {
            if (i == 0) { // x = -1
                dervvec2FD[i] = forward (funcvec[i+1], funcvec[i], delx);
                // printf("2-1: %d, %.4f, %.4f, %.4f, %.4f\n", i, funcvec[i+1], funcvec[i], delx, (funcvec[i+1] - funcvec[i])/delx);
            }
            else if (i == Ngrid - 1) { // x = 1
                dervvec2FD[i] = backward (funcvec[i], funcvec[i-1], delx);
            }
            else {
                dervvec2FD[i] = central2nd (funcvec[i+1], funcvec[i-1], delx);
                // dervvec[i] = central4th (funcvec[i+2], funcvec[i+1], funcvec[i-1], funcvec[i-2], delx);
            }
        }
    t2FD = omp_get_wtime() - t2FD;

    // 4th Order Central Difference Scheme
    double dervvec4FD[Ngrid] {};
    double t4FD = omp_get_wtime();
    #pragma omp parallel for num_threads(thrd_cnt) default(none) shared(dervvec4FD, funcvec, delx, Ngrid) private(i)
        for (i = 0; i < Ngrid; i++) {
            if (i == 0) { // x = -1
                dervvec4FD[i] = forward (funcvec[i+1], funcvec[i], delx);
            }
            else if (i == 1 || i == Ngrid - 2) { // x = -1 + delx, 1 - delx
                dervvec4FD[i] = central2nd (funcvec[i+1], funcvec[i-1], delx);
            }
            else if (i == Ngrid - 1) { // x = 1
                dervvec4FD[i] = backward (funcvec[i], funcvec[i-1], delx);
            }
            else {
                dervvec4FD[i] = central4th (funcvec[i+2], funcvec[i+1], funcvec[i-1], funcvec[i-2], delx);
            }
        }
    t4FD = omp_get_wtime() - t4FD;

    printf("Index\tActual\t\t2nd Order\t4th Order\n");
    for (i = 0; i < Ngrid; i++) {
        printf("%d\t%.4f\t\t%.4f\t\t%.4f\n",i,dervvec[i],dervvec2FD[i],dervvec4FD[i]);
    }
    printf("\n");
    printf("Time taken for 2nd order scheme = %.8f s\n",t2FD);
    printf("Time taken for 4th order scheme = %.8f s\n",t4FD);

    return 0;
}