#include <iostream>
#include <math.h>
#include <mpi.h>
#include <fstream>
#include <string.h>

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
    int i, j, myid, np;
    double del = 0.1, del2 = pow(del, 2);    

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Row-wise block decomposition
    int nx = int(2/del) + 1, ny = nx, lny;

    if (myid == 0) {lny = ny - int(ny/np)*(np - 1);}
    else {lny = int(ny/np);};

    double** phik = new double*[nx];
    double** phik1 = new double*[nx];
    double** qij = new double*[nx];

    for (i = 0; i < nx; i++) {
        phik[i] = new double[lny+2] {};
        phik1[i] = new double[lny+2] {};
        qij[i] = new double[lny] {};
    }

    int cnts[np] {}, dsplc[np] {};
    for (i = 1; i < np; i++) {
        cnts[i] = int(ny/np);
        dsplc[i] = ny - (np - i)*int(ny/np);
    }
    cnts[0] = dsplc[1];

    for (i = 0; i < nx; i++) {
        for (j = 0; j < lny; j++) {
            qij[i][j] = q(-1 + i*del, -1 + (dsplc[myid] + j)*del);
        }
    }

    // For Upwind:
    // Usable range of values --> [0,lny+1]
    // Owned range of values  --> [1,lny]

    int tagd1 = myid*11 - 1;
    int tagd2 = tagd1 + 1;
    int tagu1 = myid*11 + 10;
    int tagu2 = tagu1 + 1;
    int prtnr[2] {myid - 1, myid + 1};
    int is = 1, ie = lny;

    if (myid == 0) {
        prtnr[0] = MPI_PROC_NULL;
        is = 2;
    }
    else if (myid == (np - 1)) {
        prtnr[1] = MPI_PROC_NULL;
        ie = lny - 1;
    }

    double err = 1, eps = 1e-4;
    double errvec[nx*lny];
    int cnt = 1;
    int lim = 1e7;

    while ((err > eps) && (cnt < lim)) {
        if (myid % 2 == 0) {
            // Upsend
            MPI_Sendrecv(&phik[nx-1][],,MPI_DOUBLE,prtnr[1],tagu1,&ui_t[lnx+2],1,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Downsend
            MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,prtnr[0],tagd2,&ui_t[1],1,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else {
            // Downsend
            MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,prtnr[0],tagd1,&ui_t[1],1,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Upsend
            MPI_Sendrecv(&ui_t[lnx+1],1,MPI_DOUBLE,prtnr[1],tagu2,&ui_t[lnx+2],1,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }


    MPI_Finalize();
    return 0;
}