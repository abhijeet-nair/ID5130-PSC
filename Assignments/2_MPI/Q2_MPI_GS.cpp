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
double norm2 (double A[], int n) {
    double res {};
    for (int i =0; i < n; i++) {
        res += pow(A[i],2);
    }
    return res;
}


int main (int argc, char* argv[]) {
    int i, j, myid, np;
    double del = 0.1, del2 = pow(del, 2);    

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Row-wise block decomposition
    int nx = int(2/del) + 1, ny = nx, lnx;

    if (myid == 0) {lnx = nx - int(nx/np)*(np - 1);}
    else {lnx = int(nx/np);};

    double** phik = new double*[lnx+2];
    double phik1[lnx+2][ny];
    double** qij = new double*[lnx];

    for (i = 0; i < lnx; i++) {
        phik[i] = new double[ny] {};
        qij[i] = new double[ny] {};
    }
    phik[lnx] = new double[ny] {};
    phik[lnx+1] = new double[ny] {};

    int cnts[np] {}, dsplc[np] {};
    for (i = 1; i < np; i++) {
        cnts[i] = int(nx/np);
        dsplc[i] = nx - (np - i)*int(nx/np);
    }
    cnts[0] = dsplc[1];

    for (i = 0; i < lnx; i++) {
        for (j = 0; j < ny; j++) {
            qij[i][j] = q(-1 + (dsplc[myid] + i)*del, -1 + j*del);
        }
    }

    // Usable range of values --> [0,lnx+1]
    // Owned range of values  --> [1,lnx]

    int tagd1 = myid*11 - 1;
    int tagd2 = tagd1 + 1;
    int tagu1 = myid*11 + 10;
    int tagu2 = tagu1 + 1;
    int prtnr[2] {myid - 1, myid + 1};
    int is = 1, ie = lnx;

    if (myid == 0) {
        prtnr[0] = MPI_PROC_NULL;
        is = 2;
    }
    else if (myid == (np - 1)) {
        prtnr[1] = MPI_PROC_NULL;
        ie = lnx - 1;
    }

    double err = 1, eps = 1e-4;
    double lerr {};
    double errvec[lnx*ny];
    int cnt = 1;
    int lim = 1e7;

    MPI_Datatype mtype_1, mtype_2;
    int CT = int(0.5*ny) + (ny%2);
    int BL = 1;
    int ST = 2;
    MPI_Type_vector(CT,BL,ST,MPI_DOUBLE,&mtype_1);
    MPI_Type_commit(&mtype_1);
    CT = int(0.5*ny);
    MPI_Type_vector(CT,BL,ST,MPI_DOUBLE,&mtype_2);
    MPI_Type_commit(&mtype_2);

    while ((err > eps) && (cnt < lim)) {

    }

    MPI_Finalize();
    return 0;
}