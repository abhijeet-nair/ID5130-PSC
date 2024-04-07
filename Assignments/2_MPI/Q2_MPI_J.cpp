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
    double** phik1 = new double*[lnx+2];
    double** qij = new double*[lnx];

    for (i = 0; i < lnx; i++) {
        phik[i] = new double[ny] {};
        phik1[i] = new double[ny] {};
        qij[i] = new double[ny] {};
    }
    phik[lnx] = new double[ny] {};
    phik[lnx+1] = new double[ny] {};
    phik1[lnx] = new double[ny] {};
    phik1[lnx+1] = new double[ny] {};

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

    // For Upwind:
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
    
    while ((err > eps) && (cnt < lim)) {
        if (myid % 2 == 0) {
            // Upsend
            MPI_Sendrecv(&phik[lnx][0],ny,MPI_DOUBLE,prtnr[1],tagu1,&phik[lnx+1][0],ny,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Downsend
            MPI_Sendrecv(&phik[1][0],ny,MPI_DOUBLE,prtnr[0],tagd2,&phik[0][0],ny,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else {
            // Downsend
            MPI_Sendrecv(&phik[1][0],ny,MPI_DOUBLE,prtnr[0],tagd1,&phik[0][0],ny,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Upsend
            MPI_Sendrecv(&phik[lnx][0],ny,MPI_DOUBLE,prtnr[1],tagu2,&phik[lnx+1][0],ny,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        for (i = 1; i < lnx+1; i++) {
            phik1[i][0] = 0;
            phik1[i][ny-1] = 0;
        }
        
        if (myid == 0) {
            if (cnt % 100 == 0) {printf("cnt = %d\terr = %.5f\n",cnt,err);}
            for (j = 0; j < ny; j++) {
                phik1[1][j] = sin(2*M_PI*(-1 + j*del));
            }
        }
        
        for (i = is; i <= ie; i++) {
            for (j = 1; j < (ny - 1); j++) {
                phik1[i][j] = 0.25*(phik[i+1][j] + phik[i-1][j] + phik[i][j+1] + phik[i][j-1] + del2*qij[i-1][j]);
            }
        }

        if (myid == (np - 1)) {
            for (j = 0; j < ny; j++) {
                phik1[lnx][j] = (4*phik1[lnx-1][j] - phik1[lnx-2][j])/3;
            }
        }

        for (i = 1; i <= lnx; i++) {
            for (j = 0; j < ny; j++) {
                errvec[(i-1)*ny+j] = phik1[i][j] - phik[i][j];
                phik[i][j] = phik1[i][j];
            }
        }
        lerr = norm2(errvec, lnx*ny);

        err = 0;
        MPI_Allreduce(&lerr, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        err = sqrt(err);
        cnt += 1;
    }

    int rind = int(0.5*nx) + 1;
    int val = rind - dsplc[myid];

    double* phivsx0;
    if (myid == 0) {
        phivsx0 = (double *)malloc(nx*sizeof(double));
        double phivsy0[ny] {};
        int src {};

        for (i = 0; i < np; i++) {
            // printf("val = %d\n",rind - dsplc[i]);
            if (rind - dsplc[i] < cnts[i]) {
                src = i;
                break;
            }
        }
        // printf("src = %d\n",src);
        if (src != 0) {
            MPI_Recv(&phivsy0[0],ny,MPI_DOUBLE,src,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }

    if ((val >= 0) && (val < lnx) && (myid != 0)) {
        // printf("myid = %d\n",myid);
        MPI_Send(&phik1[val][0],ny,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
    }

    cnt = lnx;
    int BL = 1;
    int ST = ny;

    MPI_Datatype my_type;
    MPI_Type_vector(cnt,BL,ST,MPI_DOUBLE,&my_type);
    MPI_Type_commit(&my_type);

    MPI_Gatherv(&phik1[1][rind],1,my_type,&phivsx0[0],cnts,dsplc,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (myid == 0) {
        for (i = 0; i < nx; i++) {
            printf("val[%d] = %.4f\n",i,phivsx0[i]);
        }
    }

    MPI_Finalize();
    return 0;
}