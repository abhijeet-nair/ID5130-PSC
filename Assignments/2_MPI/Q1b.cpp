#include<iostream>
#include<math.h>
#include<mpi.h>
#include <fstream>
#include <string.h>

double u0 (double x) {
    double res;
    if ((x <= 0.5) && (x >= 0)) {
        res = sin(4*M_PI*x);
    }
    else {
        res = 0;
    }

    return res;
}

int main (int argc, char* argv[]) {
    double c = 1.0, L = 2.0, xi {}, ti {};
    double delx = 0.002, delt = 0.001;
    double ct_x = c*delt/delx;
    double ct_x_1 = 1 - ct_x;
    double ct_x8;
    double t0 = 0, tf = 1;
    double tr[3] = {0, 0.5, 1};

    int nx = L/delx + 1;         // No. of grid points
    int nt = (tf - t0)/delt + 1; // No. of time steps
    int i, j, myid, np, lnx, lnxe;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    if (myid == 0) {lnx = nx - int(nx/np)*(np - 1);}
    else {lnx = int(nx/np);}
    
    lnxe = lnx + 3;

    // Usable range of values --> [1,lnxe+2]
    // Owned range of values  --> [2,lnxe+1]
    double ui_t[lnxe] {}, ui_t1[lnxe] {};
    double l_uMat[3*lnx] {};
    double l_uMat_UP[2*lnx] {};
    double l_uMat_QK[2*lnx] {};

    double *uMat, *uMat_UP, * uMat_QK;
    int cnts[np] {}, dsplc[np] {};
    if (myid == 0) {
        uMat = (double *)malloc(3*nx*sizeof(double));
        uMat_UP = (double *)malloc(2*nx*sizeof(double));
        uMat_QK = (double *)malloc(2*nx*sizeof(double));
    }
    
    for (i = 1; i < np; i++) {
        cnts[i] = int(nx/np);
        dsplc[i] = nx - (np - i)*int(nx/np);
    }
    cnts[0] = dsplc[1];

    // Initial conditions
    double res;
    for (i = 0; i < lnx; i++) {
        xi = (dsplc[myid] + i)*delx;
        res = u0(xi);
        ui_t[i+2] = res;
        l_uMat[i] = res;
    }

    int cnt = 1;
    int tagd1 = myid*11 - 1;
    int tagd2 = tagd1 + 1;
    int tagu1 = myid*11 + 10;
    int tagu2 = tagu1 + 1;
    int prtnr[2] {};
    int is, ie;

    if (myid == 0) {
        prtnr[0] = MPI_PROC_NULL;
        prtnr[1] = 1;

        is = 3;
        ie = lnx + 1;
    }
    else if (myid == (np - 1)) {
        prtnr[0] = np - 2;
        prtnr[1] = MPI_PROC_NULL;

        is = 2;
        ie = lnx;
    }
    else {
        prtnr[0] = myid - 1;
        prtnr[1] = myid + 1;

        is = 2;
        ie = lnx + 1;
    }

    // First-order upwind scheme with Euler explicit time discretization
    for (j = 1; j < nt; j++) {
        if (myid == 0) {
            if(j % 100 == 0) {printf("UP - j = %d\n",j);}
            ui_t1[2] = 0;
        }
        else if (myid == np - 1) {
            ui_t1[lnxe-2] = 0;
        }
        
        if (myid % 2 == 0) {
            // Upsend
            MPI_Sendrecv(&ui_t[lnx+1],1,MPI_DOUBLE,prtnr[1],tagu1,&ui_t[lnx+2],1,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Downsend
            MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,prtnr[0],tagd2,&ui_t[1],1,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else {
            // Downsend
            MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,prtnr[0],tagd1,&ui_t[1],1,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Upsend
            MPI_Sendrecv(&ui_t[lnx+1],1,MPI_DOUBLE,prtnr[1],tagu2,&ui_t[lnx+2],1,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        for (i = is; i <= ie; i++) {
            ui_t1[i] = ct_x_1*ui_t[i] + ct_x*ui_t[i-1];
        }

        ti = t0 + j*delt;
        if (ti == tr[cnt]) {
            for (i = 0; i < lnx; i++) {
                xi = (dsplc[myid] + i)*delx;
                l_uMat[lnx*cnt+i] = u0(xi - c*ti);
                l_uMat_UP[lnx*(cnt-1)+i] = ui_t1[i+2];
            }
            cnt += 1;
        }
        memcpy(ui_t, ui_t1, lnxe*sizeof(double));
    }

    // Reinitializing the array...
    for (i = 0; i < lnx; i++) { ui_t[i+2] = u0((dsplc[myid] + i)*delx); }

    delt  = 1e-4;
    nt    = (tf - t0)/delt + 1; // No. of time steps
    ct_x  = c*delt/delx;
    ct_x8 = ct_x/8;
    cnt   = 1;

    // QUICK scheme with Euler explicit time discretization
    for (j = 1; j < nt; j++) {
        if (myid == 0) {
            if(j % 100 == 0) {printf("QK - j = %d\n",j);}
            ui_t1[2] = 0;
            ui_t1[3] = ct_x_1*ui_t[3] + ct_x*ui_t[2];
            is = 4;
        }
        else if (myid == np - 1) {
            ui_t1[lnxe-2] = 0;
        }
        
        if (myid % 2 == 0) {
            // Upsend
            MPI_Sendrecv(&ui_t[lnx],2,MPI_DOUBLE,prtnr[1],tagu1,&ui_t[lnx+2],1,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Downsend
            MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,prtnr[0],tagd2,&ui_t[0],2,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else {
            // Downsend
            MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,prtnr[0],tagd1,&ui_t[0],2,MPI_DOUBLE,prtnr[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // Upsend
            MPI_Sendrecv(&ui_t[lnx],2,MPI_DOUBLE,prtnr[1],tagu2,&ui_t[lnx+2],1,MPI_DOUBLE,prtnr[1],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        for (i = is; i <= ie; i++) {
            ui_t1[i] = ui_t[i] - ct_x8*(3*ui_t[i] - 7*ui_t[i-1] + ui_t[i-2] + 3*ui_t[i+1]);
        }

        ti = t0 + j*delt;
        if (ti == tr[cnt]) {
            for (i = 0; i < lnx; i++) {
                l_uMat_QK[lnx*(cnt-1)+i] = ui_t1[i+2];
            }
            cnt += 1;
        }
        memcpy(ui_t, ui_t1, lnxe*sizeof(double));
    }

    int sts = 0;
    if (myid == 0) {
        for (i = 0; i < lnx; i++) {
            printf("uMatQK[%d] = %.4f\n",i,l_uMat_QK[i]);
        }
        MPI_Send(&sts, 1, MPI_INT, myid+1,10,MPI_COMM_WORLD);
    }
    else if (myid == np - 1) {
        MPI_Recv(&sts, 1, MPI_INT,myid-1,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for (i = 0; i < lnx; i++) {
            printf("uMatQK[%d] = %.4f\n",dsplc[myid]+i,l_uMat_QK[i]);
        }
    }
    else {
        MPI_Recv(&sts, 1, MPI_INT,myid-1,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for (i = 0; i < lnx; i++) {
            printf("uMatQK[%d] = %.4f\n",dsplc[myid]+i,l_uMat_QK[i]);
        }
        MPI_Send(&sts, 1, MPI_INT, myid+1,10,MPI_COMM_WORLD);
    }

    MPI_Gatherv(&l_uMat[0], lnx, MPI_DOUBLE, &uMat[0], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&l_uMat[lnx], lnx, MPI_DOUBLE, &uMat[nx], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&l_uMat[2*lnx], lnx, MPI_DOUBLE, &uMat[2*nx], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&l_uMat_UP[0], lnx, MPI_DOUBLE, &uMat_UP[0], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&l_uMat_UP[lnx], lnx, MPI_DOUBLE, &uMat_UP[nx], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&l_uMat_QK[0], lnx, MPI_DOUBLE, &uMat_QK[0], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&l_uMat_QK[lnx], lnx, MPI_DOUBLE, &uMat_QK[nx], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        char fname[20] = "./Res/Q1_MPI.txt";
        std::ofstream oFile(fname);

        if (oFile.is_open()) {
            for (i = 0; i < nx; i++) {
                oFile << i << "," << uMat[i] << "," << uMat[i] << "," << uMat[i] << "\n";
            }
            oFile << "\n";

            for (i = 0; i < nx; i++) {
                oFile << i << "," << uMat[nx+i] << "," << uMat_UP[i] << "," << uMat_QK[i] << "\n";
            }
            oFile << "\n";

            for (i = 0; i < nx; i++) {
                oFile << i << "," << uMat[2*nx+i] << "," << uMat_UP[nx+i] << "," << uMat_QK[nx+i] << "\n";
            }
            oFile << "\n";

            oFile.close();
            printf("Saved in file %s\n",fname);
        }
        else {
            printf("Error opening file\n");
        }

    }

    MPI_Finalize();
    return 0;
}