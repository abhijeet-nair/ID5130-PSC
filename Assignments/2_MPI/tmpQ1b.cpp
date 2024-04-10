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

    if (myid == 0) {
        lnx = nx - int(nx/np)*(np - 1);
        // lnxe = lnx + 1;
    }
    else if (myid == np - 1) {
        lnx = int(nx/np);
        // lnxe = lnx + 2;
    }
    else {
        lnx = int(nx/np);
        // lnxe = lnx + 3;
    }
    lnxe = lnx + 3;

    double ui_t[lnxe] {}, ui_t1[lnxe] {};
    double l_uMat[3][lnx] {};
    double l_uMat_UP[2][lnx] {};
    double l_uMat_QK[2][lnx] {};

    // // Initial conditions
    // double res;
    // for (i = 0; i < lnx; i++) {
    //     xi = i*delx;
    //     res = u0(xi);
    //     // if (myid == 0) {ui_t[i] = res;}
    //     // else {ui_t[i+2] = res;}
    //     ui_t[i+2] = res;
        
    //     l_uMat[0][i] = res;
    // }

    // int cnt   = 1;
    // int tag11 = myid*11 - 1;  // -1, 10, 21, 32, 43, 54, 65, 76...
    // int tag12 = tag11 + 1;    //  0, 11, 22, 33, 44, 55, 66, 77...
    // int tag21 = myid*11 + 10; // 10, 21, 32, 43, 54, 65, 76, 87...
    // int tag22 = tag21 + 1;    // 11, 22, 33, 44, 55, 66, 77, 88...
    // int upprt = myid + 1;
    // int dnprt = myid - 1;
    // int is, ie;

    // // First-order upwind scheme with Euler explicit time discretization
    // for (j = 1; j < nt; j++) {
    //     if (myid == 0) {
    //         ui_t1[2] = 0;
    //         is = 3;
    //         ie = lnx + 1;

    //         // Upsend
    //         MPI_Sendrecv(&ui_t[lnx+1],1,MPI_DOUBLE,upprt,tag21,&ui_t[lnx+2],1,MPI_DOUBLE,upprt,tag22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    //     }
    //     else if (myid == np - 1) {
    //         ui_t1[lnxe-2] = 0;
    //         is = 2;
    //         ie = lnx;

    //         // Downsend
    //         MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,dnprt,tag11,&ui_t[1],1,MPI_DOUBLE,dnprt,tag12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     }
    //     else {
    //         is = 2;
    //         ie = lnx + 1;

    //         // Downsend
    //         MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,dnprt,tag11,&ui_t[1],1,MPI_DOUBLE,dnprt,tag12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //         // Upsend
    //         MPI_Sendrecv(&ui_t[lnx+1],1,MPI_DOUBLE,upprt,tag21,&ui_t[lnx+2],1,MPI_DOUBLE,upprt,tag22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     }

    //     // Usable range of values --> [1,lnxe+2]
    //     // Owned range of values  --> [2,lnxe+1]

    //     for (i = is; i <= ie; i++) {
    //         ui_t1[i] = ct_x_1*ui_t[i] + ct_x*ui_t[i-1];
    //     }

    //     ti = t0 + j*delt;
    //     if (ti == tr[cnt]) {
    //         for (i = 2; i <= lnx + 1; i++) {
    //             xi = i*delx;
    //             l_uMat[cnt][i-2] = u0(xi - c*ti);
    //             l_uMat_UP[cnt-1][i-2] = ui_t1[i];
    //             // if (cnt == 1) { printf("l_uMat[%d] = %.4f\tl_uMatUP[%d] = %.4f\tval = %.4f\n",i,l_uMat[cnt][i],i,l_uMat_UP[cnt-1][i],xi - c*ti); }
    //         }
    //         cnt += 1;
    //     }
    //     memcpy(ui_t, ui_t1, lnxe*sizeof(double));
    // }


    // // Reinitializing the array...
    // for (i = 0; i < lnx; i++) { ui_t[i+2] = u0((i-2)*delx); }

    // delt  = 1e-4;
    // nt    = (tf - t0)/delt + 1; // No. of time steps
    // ct_x  = c*delt/delx;
    // ct_x8 = ct_x/8;
    // cnt   = 1;

    // // QUICK scheme with Euler explicit time discretization
    // for (j = 1; j < nt; j++) {
    //     if (myid == 0) {
    //         ui_t1[2] = 0;
    //         ui_t1[3] = ct_x_1*ui_t[3] + ct_x*ui_t[2];
    //         is = 3;
    //         ie = lnx + 1;

    //         // Upsend
    //         MPI_Sendrecv(&ui_t[lnx],2,MPI_DOUBLE,upprt,tag21,&ui_t[lnx+2],1,MPI_DOUBLE,upprt,tag22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    //     }
    //     else if (myid == np - 1) {
    //         ui_t1[lnxe-2] = 0;
    //         is = 2;
    //         ie = lnx;

    //         // Downsend
    //         MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,dnprt,tag11,&ui_t[0],2,MPI_DOUBLE,dnprt,tag12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     }
    //     else {
    //         is = 2;
    //         ie = lnx + 1;

    //         // Upsend
    //         MPI_Sendrecv(&ui_t[2],1,MPI_DOUBLE,dnprt,tag11,&ui_t[0],2,MPI_DOUBLE,dnprt,tag12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //         // Downsend
    //         MPI_Sendrecv(&ui_t[lnx],2,MPI_DOUBLE,upprt,tag21,&ui_t[lnx+2],1,MPI_DOUBLE,upprt,tag22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     }
    //     // CAN REWRITE ABOVE THING WITH MPI_PROC_NULL

    //     // Usable range of values --> [1,lnxe+2]
    //     // Owned range of values  --> [2,lnxe+1]

    //     for (i = is; i <= ie; i++) {
    //         ui_t1[i] = ui_t[i] - ct_x8*(3*ui_t[i] - 7*ui_t[i-1] + ui_t[i-2] + 3*ui_t[i+1]);
    //     }

    //     ti = t0 + j*delt;
    //     if (ti == tr[cnt]) {
    //         for (i = 2; i <= lnx + 1; i++) {
    //             l_uMat_QK[cnt-1][i-2] = ui_t1[i];
    //             // if (cnt == 1) { printf("l_uMat[%d] = %.4f\tl_uMatUP[%d] = %.4f\tval = %.4f\n",i,l_uMat[cnt][i],i,l_uMat_UP[cnt-1][i],xi - c*ti); }
    //         }
    //         cnt += 1;
    //     }
    //     memcpy(ui_t, ui_t1, lnxe*sizeof(double));
    // }

    // double *uMat, *uMat_UP, * uMat_QK;
    // int cnts[np] {}, dsplc[np] {};
    // if (myid == 0) {
    //     uMat = (double *)malloc(3*nx*sizeof(double));
    //     uMat_UP = (double *)malloc(3*nx*sizeof(double));
    //     uMat_QK = (double *)malloc(3*nx*sizeof(double));

    //     cnts[0] = lnx;
    //     for (i = 1; i < np; i++) {
    //         cnts[i] = int(nx/np);
    //         dsplc[i] = lnx + (i-1)*int(nx/np);
    //     }

    //     double uMat[3][nx] {}, uMat_UP[3][nx] {}, uMat_QK[3][nx] {};
    // }

    // MPI_Gatherv(&l_uMat[2], lnx, MPI_DOUBLE, &uMat[0], cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // if (myid == 0) {
    //     for (i = 0; i < nx; i++) {
    //         printf("uMat[%d] = %.4f\n",i,uMat[i]);
    //     }
    // }
    


    // if (myid == -1) {
    //     char fname[20] = "./Res/Q1_MPI.txt";

    //     std::ofstream oFile(fname);

    //     if (oFile.is_open()) {
    //         for (i = 0; i < nx; i++) {
    //             oFile << i << "," << uMat[0][i] << "," << uMat[0][i] << "," << uMat[0][i] << "\n";
    //         }
    //         oFile << "\n";
            
    //         for (j = 1; j < 3; j++) {
    //             for (i = 0; i < nx; i++) {
    //                 oFile << i << "," << uMat[j][i] << "," << uMat_UP[j-1][i] << "," << uMat_QK[j-1][i] << "\n";
    //             }
    //             oFile << "\n";
    //         }
    //         oFile.close();
    //         printf("Saved in file %s\n",fname);
    //     }
    //     else {
    //         printf("Error opening file\n");
    //     }
    // }

    MPI_Finalize();
    return 0;
}