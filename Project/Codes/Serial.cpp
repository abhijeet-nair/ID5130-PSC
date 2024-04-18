#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

void getIndVel(double g, double x, double y, double x0, double y0, double uv[2]) {
    double den = pow((x - x0), 2) + pow((y - y0), 2);
    uv[0] = (g*(y - y0))/(2*M_PI*den);
    uv[1] = -(g*(x - x0))/(2*M_PI*den);
}

double h (double t, int f1, int f2, int a) {
    double T1 = 1/double(f1);
    double T  = 0.5*(T1 + 1/double(f2));
    double tr = remainder(t, T);

    double res;
    if (tr < 0.5*T1) {
        res = a*cos(2*M_PI*f1*tr);
    }
    else {
        res = a*cos(2*M_PI*f2*(tr - 0.5*T1) + M_PI);
    }

    return res;
}

double hdot (double t, int f1, int f2, int a) {
    double T1 = 1/double(f1);
    double T  = 0.5*(T1 + 1/double(f2));
    double tr = remainder(t, T);

    double res;
    if (tr < 0.5*T1) {
        res = -a*2*M_PI*f1*sin(2*M_PI*f1*tr);
    }
    else {
        res = -a*2*M_PI*f2*sin(2*M_PI*f2*(tr - 0.5*T1) + M_PI);
    }

    return res;
}

double deg2rad (double x) {
    return x*M_PI/180;
}



int main () {
    int i, j, k, m, p;  // Indices

    int Nl     = 10;    // No. of collocation points
    double u   = 20;    // Freestream velocity
    double c   = 5;     // Chord length of the flat plate
    double alp = 0;     // Angle of attack of the plate
    double rho = 1.225; // Density
    double dx  = c/Nl;  // Spacing between two points

    double sn = sin(deg2rad(alp));
    double cs = cos(deg2rad(alp));

    double dt = 0.01;        // Time step
    double tf = 5;           // Final time
    int Nt = int(tf/dt) + 1; // No. of time steps

    int n  = 1;
    int f1 = 1;
    int f2 = n*f1;
    int a  = 1;

    printf("Nl = %d\nNt = %d\n",Nl,Nt);

    double vorloc[Nl], colcloc[Nl], plate[Nl+1];

    for (i = 0; i < Nl; i++) {
        vorloc[i]  = (i + 0.25)*dx;
        colcloc[i] = (i + 0.75)*dx;
        plate[i]   = i*dx;
    }
    plate[Nl] = Nl*dx;

    double A[Nl][Nl], B[Nl], gIVRes[2];

    for (i = 0; i < Nl; i++) {
        for (j = 0; j < Nl; j++) {
            getIndVel(1, colcloc[i], 0, vorloc[j], 0, gIVRes);
            A[i][j] = gIVRes[1];
        }
        getIndVel(1, colcloc[i], 0, (c + 0.1*dx), 0, gIVRes);
        B[i] = gIVRes[1];
    }

    double C[Nl+1][Nl+1] {}, U[Nl+1] {}, U1[Nl+1] {};

    for (i = 0; i <= Nl; i++) {
        for (j = 0; j <= Nl; j++) {
            if (i == Nl) {
                C[i][j] = 1;
            }
            else if (j == Nl) {
                C[i][j] = B[i];
            }
            else {
                C[i][j] = A[i][j];
            }
        }
    }

    double ydot, t_cur, extflow;
    // double x0[Nt], y0[Nt], L[Nt], D[Nt];
    double* x0 = new double[Nt];
    double* y0 = new double[Nt];
    double* L = new double[Nt];
    double* D = new double[Nt];
    // double xw[Nt], yw[Nt];
    double** xw = new double* [Nt];
    double** yw = new double* [Nt];
    double clx, cly;
    double* R = new double[Nl+1];
    double* gw = new double[Nt];
    double** gbm = new double* [2];
    gbm[0] = new double[Nl];
    gbm[1] = new double[Nl];
    // double R[Nl+1], gw[Nt], gbm[2][Nl];
    double Bjs, uw, vw;
    double dgbm_dt, vindw[Nl];

    double err, eps = 1e-6, lim = 1e7, sum;
    int cnt;

    for (m = 0; m < Nt; m++) {
        xw[m] = new double[Nt] {};
        yw[m] = new double[Nt] {};
    }

    // For loop for time marching
    for (m = 0; m < 2; m++) {
        t_cur = m*dt;
        ydot = hdot(t_cur, f1, f2, a);

        extflow = u*sn - ydot*cs;

        // Origin location at time t
        x0[m] = -u*t_cur;
        y0[m] = h(t_cur, f1, f2, a);

        // New wake location
        xw[m][m] = x0[m] + (c + 0.1*dx)*cs;
        yw[m][m] = y0[m] - (c + 0.1*dx)*sn;

        for (i = 0; i < Nl; i++) {
            // Collocation point location
            clx = colcloc[i]*cs + x0[m];
            cly = colcloc[i]*sn + y0[m];

            // B vectors and R (RHS) vector calculations
            Bjs = 0;
            for (j = 0; j < m; j++) {
                getIndVel(1, clx, cly, xw[m-1][j], yw[m-1][j], gIVRes);
                // Bj[i][j] = gIVRes[0]*sn + gIVRes[1]*cs;
                Bjs += (gIVRes[0]*sn + gIVRes[1]*cs)*gw[j];
            }
            R[i] = -extflow - Bjs;
        }

        for (j = 0; j < m; j++) {
            R[Nl] += gw[j];
        }
        R[Nl] = -R[Nl];

        // Computation of unknown gbm and gwm
        // Gauss-Seidel
        err = 1;
        cnt = 1;
        while (err > eps and cnt < lim) {
            for (i = 0; i <= Nl; i++) {
                sum = 0;
                for (j = 0; j <= Nl; j++) {
                    if (j < i) {
                        sum += C[i][j]*U1[j];
                    }
                    else if (j > i) {
                        sum += C[i][j]*U[j];
                    }
                }
                U1[i] = (R[i] - sum)/C[i][i];
            }

            err = 0;
            for (i = 0; i <= Nl; i++) {
                err += pow(U1[i] - U[i], 2);
                U[i] = U1[i];
            }
            

            err = sqrt(err);
            cnt += 1;
            // if (cnt%10 == 0) {printf("cnt = %d\n",cnt);}
        }
        printf("m = %d\titer = %d\n",m,cnt-1);
        for (i = 0; i <= Nl; i++) {
                printf("%.4f\n",U[i]);
            }
            printf("\n");

        // Solved body and wake circulations
        for (i = 0; i < Nl; i++) {
            gbm[1][i] = U[i];
        }
        gw[m] = U[Nl];

        // Wake Position Update
        for (k = 0; k < m; k++) {
            uw = 0;
            vw = 0;

            for (i = 0; i < Nl; i++) {
                // Due to body vortices on wakes
                getIndVel(gbm[1][i], xw[m-1][k], yw[m-1][k], vorloc[i]*cs + x0[m], -vorloc[i]*sn + y0[m], gIVRes);
                uw += gIVRes[0];
                vw += gIVRes[1];
            }
            
            for (p = 0; p < m; p++) {
                // Other wakes on TE wake
                if (p != k) {
                    getIndVel(gw[p], xw[m-1][k], yw[m-1][k], xw[m-1][p], yw[m-1][p], gIVRes);
                    uw += gIVRes[0];
                    vw += gIVRes[1];
                }
            }
            xw[m][k] += uw*dt;
            yw[m][k] += vw*dt;
        }

        // Aerodynamic Load Calculations
        dgbm_dt = 0;
        if (m == 1) {
            for (i = 0; i < Nl; i++) {
                dgbm_dt += gbm[1][i];
            }
            dgbm_dt /= dt;
        }
        else {
            for (i = 0; i < Nl; i++) {
                dgbm_dt += gbm[1][i] - gbm[0][i];
            }
            dgbm_dt /= dt;
        }

        L[m] = 0;
        for (i = 0; i < Nl; i++) {
            vindw[i] = 0;
            for (p = 0; p < m; p++) {
                getIndVel(gw[p], xw[m-1][k], yw[m-1][k], xw[m][p], yw[m][p], gIVRes);
                vindw[i] += gIVRes[0]*sn + gIVRes[1]*cs;
            }

            L[m] += gbm[1][i];
            D[m] += vindw[i]*gbm[1][i];
        }
        L[m] = rho*(u*L[m] + dgbm_dt*c);
        D[m] = rho*(D[m] + dgbm_dt*c*deg2rad(alp));
    }
    // Data needed finally:
    // At each instant:
    // x0, y0, platex, platey, xw, yw, L, D
    
    // Common:
    // times, rho, u, c

    // platex and platey will be calculated in python
    // So, instead send plate, cs and sn also as common
    // times is not there. So, send Nt and dt

    // Final List:
    // At each instant:
    // x0[Nt], y0[Nt], xw[Nt][Nt], yw[Nt][Nt], L[Nt], D[Nt]

    // Common:
    // Nt, Nl, dt, dx, rho, u, c, alp

    // // File writing
    // char fname[25];
    // sprintf(fname, "./Res/Ser_1.txt");

    // std::ofstream oFile(fname);

    // if (oFile.is_open()) {
    //     oFile << Nl << "\n" << Nt << "\n";
    //     oFile << dt << "\n" << dx << "\n";
    //     oFile << rho << "\n" << u << "\n";
    //     oFile << c << "\n" << alp << "\n";

    //     oFile << "\n";

    //     for (i = 0; i < Nt; i++) {
    //         oFile << x0[i] << "," << y0[i] << "\n";
    //     }

    //     oFile << "\n";

    //     for (i = 0; i < Nt; i++) {
    //         oFile << L[i] << "," << D[i] << "\n";
    //     }

    //     oFile << "\n\n";

    //     for (i = 0; i < Nt; i++) {
    //         for (j = 0; j < Nt; j++) {
    //             oFile << xw[i][j] << "," << yw[i][j] << "\n";
    //         }
    //         oFile << "\n";
    //     }

    //     oFile.close();
    //     printf("Saved in file %s\n",fname);
    // }
    // else {
    //     printf("Error opening file\n");
    // }

    delete x0, y0, L, D, xw, yw, R, gw, gbm;
    return 0;
}