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
    for (int i = 0; i < n; i++) {
        res += pow(A[i],2);
    }
    return res;
}


int main (int argc, char* argv[]) {
    int i, j, myid, np;
    double del = 0.01, del2 = pow(del, 2);    

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // i along x and j along y.
    // x, y in the usual sense as real world.
    // Row-wise block decomposition
    int nx = int(2/del) + 1, ny = nx, lnx;

    if (myid == 0) {lnx = nx - int(nx/np)*(np - 1);}
    else {lnx = int(nx/np);};

    double phik[lnx+2][ny] {};
    double phik1[lnx+2][ny] {};
    double** qij = new double*[lnx];

    for (i = 0; i < lnx; i++) {
        qij[i] = new double[ny] {};
    }

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

    // Tags and partners for each processor
    int tagd1 = myid*11 - 1;  // -1, 10, 21, 32, 43, 54, 65, 76...
    int tagd2 = tagd1 + 1;    //  0, 11, 22, 33, 44, 55, 66, 77...
    int tagu1 = myid*11 + 10; // 10, 21, 32, 43, 54, 65, 76, 87...
    int tagu2 = tagu1 + 1;    // 11, 22, 33, 44, 55, 66, 77, 88...
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

    // Making custom data-type to send alternate data in a row.
    // mtype_1 starts at the 1st point (j = 0)
    // mtype_2 starts at the 2nd point (j = 1)
    MPI_Datatype mtype_1; // End-type
    MPI_Datatype mtype_2; // Mid-type
    int CT = int(0.5*ny) + (ny%2);
    int BL = 1;
    int ST = 2;
    MPI_Type_vector(CT,BL,ST,MPI_DOUBLE,&mtype_1);
    MPI_Type_commit(&mtype_1);
    CT = int(0.5*ny);
    MPI_Type_vector(CT,BL,ST,MPI_DOUBLE,&mtype_2);
    MPI_Type_commit(&mtype_2);

    // Starting index of P0 (0) is red.
    // So, if (dsplc[myid] + i + j) is even, it is red (i-row & j-col index) and so on.

    double runT = MPI_Wtime();
    // Main loop
    while ((err > eps) && (cnt < lim)) {
        // Dirichlet Boundary Conditions on y = -1, 1 faces
        for (i = 1; i < lnx+1; i++) {
            phik1[i][0] = 0;
            phik1[i][ny-1] = 0;
        }

        // Dirichlet Boundary Conditions on x = -1 face (Only on proc 0)
        if (myid == 0) {
            for (j = 0; j < ny; j++) {
                phik1[1][j] = sin(2*M_PI*(-1 + j*del));
            }
        }
        
        // Blue-send for red points
        if (myid % 2 == 0) {
            // Upsend
            if (dsplc[myid+1] % 2 == 1) { // Next proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik[lnx][1],1,mtype_2,prtnr[1],tagu1,&phik[lnx+1][0],1,mtype_1,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                        // Next proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik[lnx][0],1,mtype_1,prtnr[1],tagu1,&phik[lnx+1][1],1,mtype_2,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            // Downsend
            if (dsplc[myid] % 2 == 1) { // This proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik[1][0],1,mtype_1,prtnr[0],tagd2,&phik[0][1],1,mtype_2,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                      // This proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik[1][1],1,mtype_2,prtnr[0],tagd2,&phik[0][0],1,mtype_1,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else {
            // Downsend
            if (dsplc[myid] % 2 == 1) { // This proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik[1][0],1,mtype_1,prtnr[0],tagd1,&phik[0][1],1,mtype_2,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                      // This proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik[1][1],1,mtype_2,prtnr[0],tagd1,&phik[0][0],1,mtype_1,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            // Upsend
            if (dsplc[myid+1] % 2 == 1) { // Next proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik[lnx][1],1,mtype_2,prtnr[1],tagu2,&phik[lnx+1][0],1,mtype_1,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                        // Next proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik[lnx][0],1,mtype_1,prtnr[1],tagu2,&phik[lnx+1][1],1,mtype_2,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        
        // Red-calculations
        for (i = is; i <= ie; i++) {
            for (j = 1; j < (ny - 1); j++) {
                if ((dsplc[myid] + i - 1 + j) % 2 == 0) {
                    phik1[i][j] = 0.25*(phik[i+1][j] + phik[i-1][j] + phik[i][j+1] + phik[i][j-1] + del2*qij[i-1][j]);
                }
            }
        }

        // Red-send for blue points
        if (myid % 2 == 0) {
            // Upsend
            if (dsplc[myid+1] % 2 == 1) { // Next proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[lnx][0],1,mtype_1,prtnr[1],tagu1,&phik1[lnx+1][1],1,mtype_2,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                        // Next proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[lnx][1],1,mtype_2,prtnr[1],tagu1,&phik1[lnx+1][0],1,mtype_1,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            // Downsend
            if (dsplc[myid] % 2 == 1) { // This proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[1][1],1,mtype_2,prtnr[0],tagd2,&phik1[0][0],1,mtype_1,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                      // This proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[1][0],1,mtype_1,prtnr[0],tagd2,&phik1[0][1],1,mtype_2,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else {
            // Downsend
            if (dsplc[myid] % 2 == 1) { // This proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[1][1],1,mtype_2,prtnr[0],tagd1,&phik1[0][0],1,mtype_1,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                      // This proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik[1][0],1,mtype_1,prtnr[0],tagd1,&phik1[0][1],1,mtype_2,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            // Upsend
            if (dsplc[myid+1] % 2 == 1) { // Next proc has a blue at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[lnx][0],1,mtype_1,prtnr[1],tagu2,&phik1[lnx+1][1],1,mtype_2,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {                        // Next proc has a red at (i,j) = (0,0)
                MPI_Sendrecv(&phik1[lnx][1],1,mtype_2,prtnr[1],tagu2,&phik1[lnx+1][0],1,mtype_1,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }

        // Blue calculations
        for (i = is; i <= ie; i++) {
            for (j = 1; j < (ny - 1); j++) {
                if ((dsplc[myid] + i - 1 + j) % 2 == 1) {
                    phik1[i][j] = 0.25*(phik1[i+1][j] + phik1[i-1][j] + phik1[i][j+1] + phik1[i][j-1] + del2*qij[i-1][j]);
                }
            }
        }

        // Neumann Boundary Conditions on x = 1 face (Only on proc np-1)
        if (myid == (np - 1)) {
            for (j = 0; j < ny; j++) {
                phik1[lnx][j] = (4*phik1[lnx-1][j] - phik1[lnx-2][j])/3;
            }
        }

        // Error calculations and array reallocations
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
        if (myid == 0) {if (cnt % 500 == 0) {printf("cnt = %5.0f\terr = %.6f\n",double(cnt),err);}}
        cnt += 1;
    }

    // Getting data for x = 0
    int rind = int(0.5*nx); // Index for (nx/2)+1-th element
    int val = rind - dsplc[myid];

    double* phivsx0,* phivsy0;
    if (myid == 0) {
        phivsx0 = (double *)malloc(nx*sizeof(double));
        phivsy0 = (double *)malloc(ny*sizeof(double));
        int src {};

        // Finding owner of x = 0 data
        for (i = 0; i < np; i++) {
            if (rind - dsplc[i] < cnts[i]) {
                src = i;
                break;
            }
        }

        // Collecting x = 0 data
        if (src != 0) {
            MPI_Recv(&phivsy0[0],ny,MPI_DOUBLE,src,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else {
            for (i = 0; i < ny; i++) {
                phivsy0[i] = phik1[rind+1][i]; // +1 because of the halo row 0
            }
        }
    }

    if ((val >= 0) && (val < lnx) && (myid != 0)) {
        val += 1; // +1 because of the halo row 0. Essentially, MATLAB indexing now.
        MPI_Send(&phik1[val][0],ny,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
    }

    // Getting data for y = 0
    rind = int(0.5*ny); // Index for (ny/2)+1-th element
    CT = lnx;
    BL = 1;
    ST = ny;

    // Making custom data-type to send data from same column, all rows.
    MPI_Datatype mtype_3;
    MPI_Type_vector(CT,BL,ST,MPI_DOUBLE,&mtype_3);
    MPI_Type_commit(&mtype_3);

    MPI_Gatherv(&phik1[1][rind],1,mtype_3,&phivsx0[0],cnts,dsplc,MPI_DOUBLE,0,MPI_COMM_WORLD);
    runT = MPI_Wtime()- runT;

    if (myid == 0) {
        printf("Iterations = %d\terr = %.8f\tTime = %.4f s\n",cnt-1,err,runT);

        // // Saving in a file for plotting
        // char fname[25];
        // sprintf(fname, "./Res/Q2_MPI_GS%d.txt", np);

        // std::ofstream oFile(fname);

        // if (oFile.is_open()) {
        //     for (i = 0; i < nx; i++) {
        //         oFile << phivsx0[i] << "," << phivsy0[i] << "\n";
        //     }

        //     oFile.close();
        //     printf("Saved in file %s\n",fname);
        // }
        // else {
        //     printf("Error opening file\n");
        // }
    }

    delete qij;
    MPI_Type_free(&mtype_1);
    MPI_Type_free(&mtype_2);
    MPI_Type_free(&mtype_3);

    MPI_Finalize();
    return 0;
}