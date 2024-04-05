#include <iostream>
#include <math.h>
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
    int i, j;
    double c = 1.0, L = 2.0, xi {}, ti {};
    double delx = 0.002, delt = 0.001;
    const double ct_x = c*delt/delx, ct_x_1 = 1 - ct_x;
    double t0 = 0, tf = 1;

    int nx = L/delx + 1;     // No. of grid points
    int nt = (tf - t0)/delt + 1; // No. of time steps
    printf("nt = %d\tnx = %d\n",nt,nx);

    double ui_t[nx] {}, ui_t1[nx] {};
    // double u_xt_UP[nt][nx] {}, u_xt[nt][nx] {};
    printf("Bruh...\n");
    double uMat[3][nx] {};
    double uMat_UP[2][nx] {};
    double uMat_QK[2][nx] {};

    
    printf("Starting Initialization...\n");
    // Initial conditions
    double res;
    for (i = 0; i < nx; i++) {
        xi = i*delx;
        res = u0(xi);
        ui_t[i] = res;
        uMat[0][i] = res;
    }

    printf("Completed Initialization...\n");

    double tr[3] = {0, 0.5, 1};
    int cnt = 1;

    // First-order upwind scheme with Euler explicit time discretization
    for (j = 1; j < nt; j++) {
        ui_t1[0] = 0;
        ui_t1[nx-1] = 0;

        // Boundary Conditions are not put, because they are also zero.

        for (i = 1; i < nx - 1; i++) {
            ui_t1[i] = ct_x_1*ui_t[i] + ct_x*ui_t[i-1];
            // u_xt_UP[j][i] = ct_x_1*u_xt_UP[j-1][i] + ct_x*u_xt_UP[j-1][i-1];
        }

        ti = t0 + j*delt;
        // printf("ti = %.4f\n",ti);

        if (ti == tr[cnt]) {
            for (i = 0; i < nx; i++) {
                xi = i*delx;
                uMat[cnt][i] = u0(xi - c*ti);
                uMat_UP[cnt-1][i] = ui_t1[i];
                if (cnt == 1) { printf("uMat[%d] = %.4f\tuMatUP[%d] = %.4f\tval = %.4f\n",i,uMat[cnt][i],i,uMat_UP[cnt-1][i],xi - c*ti); }
            }
            cnt += 1;
        }
        memcpy(ui_t, ui_t1, nx*sizeof(double));
    }

    // double u_xt_QK[nt][nx] {};
    for (i = 0; i < nx; i++) { ui_t[i] = u0(i*delx); }

    double ct_x8 = ct_x/8;
    delt = 1e-4;
    nt = (tf - t0)/delt + 1; // No. of time steps
    cnt = 0;
    // QUICK scheme with Euler explicit time discretization
    for (j = 1; j < nt; j++) {
        // ui_t1[0] = 0;
        // ui_t1[nx-1] = 0;

        // Boundary Conditions are not put, because they are also zero.

        ui_t1[1] = ct_x_1*ui_t[1] + ct_x*ui_t[0];
        for (i = 2; i < nx - 1; i++) {
            ui_t1[i] = ui_t[i] - ct_x8*(3*ui_t[i] - 7*ui_t[i-1] + ui_t[i-2] + 3*ui_t[i+1]);
            // u_xt_QK[j][i] = u_xt_QK[j-1][i] - ct_x8*(3*u_xt_QK[j-1][i] - 7*u_xt_QK[j-1][i-1] + u_xt_QK[j-1][i-2] + 3*u_xt_QK[j-1][i+1]);
        }

        ti = t0 + j*delt;
        // printf("ti = %.4f\n",ti);

        if (ti == tr[cnt]) {
            for (i = 0; i < nx; i++) {
                uMat_QK[cnt][i] = ui_t1[i];
            }
            cnt += 1;
        }
        memcpy(ui_t, ui_t1, nx*sizeof(double));
    }

    char fname[20] = "./Res/Q1_Ser.txt";
    // char fname2[20] = "./Res/Q1_Ser_QK.txt";

    std::ofstream oFile(fname);

    if (oFile.is_open()) {
        for (i = 0; i < nx; i++) {
            oFile << uMat[0][i] << "," << uMat[0][i] << "," << uMat[0][i] << "\n";
        }
        oFile << "\n";
        
        for (j = 1; j < 3; j++) {
            for (i = 0; i < nx; i++) {
                oFile << uMat[j][i] << "," << uMat_UP[j][i] << "," << uMat_QK[j][i] << "\n";
            }
            oFile << "\n";
        }
        oFile.close();
        printf("Saved in file %s\n",fname);
    }
    else {
        printf("Error opening file\n");
    }

    return 0;
}