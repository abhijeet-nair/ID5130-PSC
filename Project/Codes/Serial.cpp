#include <iostream>
#include <math.h>
#include <string.h>
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



int main () {
    int i, j, m;        // Indices

    int Nl = 10;        // No. of collocation points
    double u   = 20;    // Freestream velocity
    double c   = 5;     // Chord length of the flat plate
    double alp = 0;     // Angle of attack of the plate
    double rho = 1.225; // Density
    double dx  = c/Nl;   // Spacing between two points

    double sn = sin(alp*M_PI/180);
    double cs = cos(alp*M_PI/180);

    double dt = 0.01;    // Time step
    double tf = 5;       // Final time
    int Nt = int(tf/Nt); // No. of time steps

    int n  = 1;
    int f1 = 1;
    int f2 = n*f1;
    int a  = 1;

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

    double C[Nl+1][Nl+1] {};

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

    double ydot, t_cur, extflow, x0, y0;
    double xw[Nt], yw[Nt];
    double clx, cly;
    double R[Nl+1], gw[Nt];
    double Bjs;

    for (m = 0; m < Nt; m++) {
        t_cur = m*dt;
        ydot = hdot(t_cur, f1, f2, a);

        extflow = u*sn - ydot*cs;

        x0 = -u*t_cur;
        y0 = h(t_cur, f1, f2, a);

        xw[m] = x0 + (c + 0.1*dx)*cs;
        yw[m] = y0 - (c + 0.1*dx)*sn;

        for (i = 0; i < Nl; i++) {
            clx = colcloc[i]*cs + x0;
            cly = colcloc[i]*sn + y0;

            Bjs = 0;
            for (j = 0; j < m; j++) {
                getIndVel(1, clx, cly, xw[j], yw[j], gIVRes);
                // Bj[i][j] = gIVRes[0]*sn + gIVRes[1]*cs;
                Bjs += (gIVRes[0]*sn + gIVRes[1]*cs)*gw[j];
            }
            R[i] = -extflow - Bjs;
        }
        for (j = 0; j < m; j++) {
            R[Nl] += gw[j];
        }
        R[Nl] = -R[Nl];

    }

    return 0;
}