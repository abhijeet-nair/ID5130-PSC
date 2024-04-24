#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

void getIndVel(double g, double x, double y, double x0, double y0, double uv[2]) {
    double den = pow((x - x0), 2) + pow((y - y0), 2);
    // printf("den = %.4f\n",den);
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



int main (int argc, char* argv[]) {
    int np = 1, cnt;

    if (argc == 2) {
        np = strtol(argv[1], NULL, 10);
    }
    printf("Running in code with %d processors...\n",np);
    
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
    plate[Nl] = Nl*dx;

    // double A[Nl][Nl], B[Nl], gIVRes[2];
    // double C[Nl+1][Nl+1] {}, U[Nl+1] {}, U1[Nl+1] {};
    double gIVRes[2], dgbm_dt, vindw[Nl], xw[Nt] {}, yw[Nt] {};
    double ydot, t_cur, extflow, sum, Bjs, uw, vw, clx, cly;
    double* B = new double[Nl];
    double* U = new double[Nl+1];
    double* U1 = new double[Nl+1];
    double** A = new double* [Nl];
    double** C = new double* [Nl+1];

    for (i = 0; i <= Nl; i++) {
        if (i < Nl) {
            A[i] = new double[Nl] {};
            C[i] = new double[Nl+1] {};
        }
        C[i] = new double[Nl+1] {};
    }

    // double x0[Nt], y0[Nt], L[Nt], D[Nt];
    double* x0 = new double[Nt];
    double* y0 = new double[Nt];
    double* L = new double[Nt];
    double* D = new double[Nt];
    // double xw[Nt], yw[Nt];
    double** xwM = new double* [Nt];
    double** ywM = new double* [Nt];
    double* R = new double[Nl+1];
    double* gw = new double[Nt];
    double** gbm = new double* [2];
    gbm[0] = new double[Nl];
    gbm[1] = new double[Nl];
    // double R[Nl+1], gw[Nt], gbm[2][Nl];

    // double err, eps = 1e-6, lim = 1e7;

    for (m = 0; m < Nt; m++) {
        xwM[m] = new double[Nt] {};
        ywM[m] = new double[Nt] {};
    }


    #pragma omp parallel num_threads(np) default(shared) private(i, j, gIVRes)
    {
        #pragma omp for
        for (i = 0; i < Nl; i++) {
            vorloc[i]  = (i + 0.25)*dx;
            colcloc[i] = (i + 0.75)*dx;
            plate[i]   = i*dx;
        }

        #pragma omp barrier

        #pragma omp for
        for (i = 0; i < Nl; i++) {
            for (j = 0; j < Nl; j++) {
                getIndVel(1, colcloc[i], 0, vorloc[j], 0, gIVRes);
                A[i][j] = gIVRes[1];
            }
            getIndVel(1, colcloc[i], 0, (c + 0.1*dx), 0, gIVRes);
            B[i] = gIVRes[1];
        }

        #pragma omp barrier

        #pragma omp for collapse(2)
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
    }

    for (k = 0; k <= Nl; k++) {
        for (i = k+1; i <= Nl; i++) {
            C[i][k] = C[i][k]/C[k][k];

            for (j = k+1; j <= Nl; j++) {
                C[i][j] += -C[i][k]*C[k][j];
            }
        }
    }

    // For loop for time marching
    for (m = 0; m < Nt; m++) {
        printf("m = %3.0f\t",double(m));
        t_cur = m*dt;
        ydot = hdot(t_cur, f1, f2, a);

        extflow = u*sn - ydot*cs;

        // Origin location at time t
        x0[m] = -u*t_cur;
        y0[m] = h(t_cur, f1, f2, a);

        // New wake location
        xw[m] = x0[m] + (c + 0.1*dx)*cs;
        yw[m] = y0[m] - (c + 0.1*dx)*sn;

        for (i = 0; i < Nl; i++) {
            // Collocation point location
            clx = colcloc[i]*cs + x0[m];
            cly = colcloc[i]*sn + y0[m];

            // B vectors and R (RHS) vector calculations
            Bjs = 0;
            for (j = 0; j < m; j++) {
                getIndVel(1, clx, cly, xw[j], yw[j], gIVRes);
                // Bj[i][j] = gIVRes[0]*sn + gIVRes[1]*cs;
                Bjs += (gIVRes[0]*sn + gIVRes[1]*cs)*gw[j];
            }
            R[i] = -extflow - Bjs;
        }

        R[Nl] = 0;
        for (j = 0; j < m; j++) {
            R[Nl] += gw[j];
            // printf("gw[%d] = %.4f\n",j,gw[j]);
        }
        R[Nl] = -R[Nl];
    }

    return 0;
}