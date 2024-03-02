#include <iostream>
#include <fstream>
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
    double y[N] {};
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

            y[i] = (-2.5*fvec[0] + 2*fvec[1] + 0.5*fvec[2])/h;
        }
        else if (i == N - 1) {
            A[i][N-2] = 2;
            A[i][N-1] = 1;

            y[i] = (2.5*fvec[N-1] - 2*fvec[N-2] - 0.5*fvec[N-3])/h;
        }
        else {
            A[i][i-1] = 1;
            A[i][i] = 4;
            A[i][i+1] = 1;

            y[i] = 3*(fvec[i+1] - fvec[i-1])/h;
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

    // y[0] = 2;
    // y[1] = 1;
    // y[2] = 1;
    // y[3] = 2;

    printMatrix(A, N, N);
    printf("\n");
    printVector(y, N);
    printf("\n");

    double ai[N] {};
    double bi[N] {};
    double ci[N] {};

    for (i = 0; i < N-1; i++) {
        ai[i+1] = A[i+1][i];
        bi[i] = A[i][i];
        ci[i] = A[i][i+1];
    }
    bi[N-1] = A[N-1][N-1];

    // printf("ai = \n");
    // printVector(ai, N);
    // printf("bi = \n");
    // printVector(bi, N);
    // printf("ci = \n");
    // printVector(ci, N);

    double u[N] {};
    double l[N] {};

    u[0] = bi[0];

    for (i = 1; i < N; i++) {
        l[i] = ai[i]/u[i-1];
        u[i] = bi[i] - l[i]*ci[i-1];
    }

    // printf("l = \n");
    // printVector(l, N);
    // printf("u = \n");
    // printVector(u, N);

    // Forward substitution for z
    double z[N] {};
    z[0] = y[0];

    for (i = 1; i < N; i++) {
        z[i] = y[i] - l[i]*z[i-1];
    }

    printf("z = \n");
    printVector(z, N);

    // Backward substitution for x
    double x[N] {};
    x[N-1] = z[N-1]/u[N-1];

    for (i = N - 2; i >= 0; i--) {
        x[i] = (z[i] - ci[i]*x[i+1])/u[i];
    }

    printf("\nx = \n");
    printVector(x, N);

    std::ofstream oFile("./Res/LU.txt");

    if (oFile.is_open()) {
        for (i = 0; i < N; i++) {
            oFile << x[i] << "," << fdvec[i] << "\n";
        }
        oFile.close();
        printf("Saved in file ./Res/LU.txt\n");
    }
    else {
        printf("Error opening file\n");
    }

    return 0;
}