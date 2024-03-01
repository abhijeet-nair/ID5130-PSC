#include <iostream>
#include <math.h>

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
    int N {}, i {}, j {};
    std::cout << "Enter number of elements (> 0) of the matrices: ";
    std::cin >> N;

    double** A = new double*[N];
    double yi_k[N] {};
    double fvec[N] {};
    double fdvec[N] {};
    double h = 3/(double(N-1));
    double xi {};

    for (i = 0; i < N; i++) {
        A[i] = new double[N] {};

        xi = i*h;
        fvec[i] = myfunc(xi);
        fdvec[i] = myderv(xi);
    }

    // printf("h = %.4f\nf = ",h);
    // printVector(fvec, N);

    // Initializing A
    for (i = 0; i < N; i++)  {
        // printf("(%d,%d)\n",i,N);

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

    // A[0][0] = 3;
    // A[0][1] = -1;
    // A[1][0] = -1;
    // A[1][1] = 3;
    // A[1][2] = -1;
    // A[2][1] = -1;
    // A[2][2] = 3;
    // A[2][3] = -1;
    // A[3][2] = -1;
    // A[3][3] = 3;

    // yi_k[0] = 2;
    // yi_k[1] = 1;
    // yi_k[2] = 1;
    // yi_k[3] = 2;

    printMatrix(A, N, N);
    printf("\n");
    printVector(yi_k, N);
    printf("\n");

    double ai_k[N] {};
    double bi_k[N] {};
    double ci_k[N] {};

    for (i = 0; i < N-1; i++) {
        ai_k[i+1] = A[i+1][i];
        bi_k[i] = A[i][i];
        ci_k[i] = A[i][i+1];
    }
    bi_k[N-1] = A[N-1][N-1];

    int n_RD = ceil(log2(N));
    double ai_k1[N] {};
    double bi_k1[N] {};
    double ci_k1[N] {};
    double yi_k1[N] {};
    double alpi_k[N] {};
    double beti_k[N] {};
    int l1 {}, l2 {};

    for (int k = 1; k <= n_RD; k++) {
        l1 = pow(2, k-1);
        l2 = pow(2, k);

        for (i = 0; i < N; i++) {
            printf("%d, %d, %d, %d\n",i,k,l1,l2);

            bi_k1[i] = bi_k[i];
            yi_k1[i] = yi_k[i];

            if (i <= N - l1 - 1) {
                // printf("betIf\n");
                beti_k[i] = -ci_k[i]/bi_k[i+l1];

                bi_k1[i] += beti_k[i]*ai_k[i+l1];
                yi_k1[i] += beti_k[i]*yi_k[i+l1];
            }
            else {
                beti_k[i] = 0;
            }

            if (i >= l1) {
                // printf("alpIf\n");
                alpi_k[i] = -ai_k[i]/bi_k[i-l1];

                bi_k1[i] += alpi_k[i]*ci_k[i-l1];
                yi_k1[i] += alpi_k[i]*yi_k[i-l1];
            }
            else {
                alpi_k[i] = 0;
            }

            if (i <= N - l2 - 1) {
                // printf("cIf\n");
                ci_k1[i] = beti_k[i]*ci_k[i+l1];
            }
            else {
                ci_k1[i] = 0;
            }

            if (i >= l2) {
                // printf("aIf\n");
                ai_k1[i] = alpi_k[i]*ai_k[i-l1];
            }
            else {
                ai_k1[i] = 0;
            }

            // printf("alpi = \n");
            // printVector(alpi_k, N);
            // printf("ai = \n");
            // printVector(ai_k1, N);
            // printf("beti = \n");
            // printVector(beti_k, N);
            // printf("ci = \n");
            // printVector(ci_k1, N);
            // printf("bi = \n");
            // printVector(bi_k1, N);
            // printf("yi = \n");
            // printVector(yi_k1, N);
            // printf("\n");
        }

        if (k <= n_RD - 1) {
            for (i = 0; i < N; i++) {
                ai_k[i] = ai_k1[i];
                bi_k[i] = bi_k1[i];
                ci_k[i] = ci_k1[i];
                yi_k[i] = yi_k1[i];
            }
        }
    }

    double xsol[N] {};

    for (i = 0; i < N; i++) {
        xsol[i] = yi_k1[i]/bi_k1[i];
    }

    printf("xi = \n");
    printVector(xsol, N);

    // printf("ai = \n");
    // printVector(ai_k1, N);
    // printf("bi = \n");
    // printVector(bi_k1, N);
    // printf("ci = \n");
    // printVector(ci_k1, N);
    // printf("yi = \n");
    // printVector(yi_k1, N);
    // printf("A = \n");
    // printMatrix(A, N, N);
}