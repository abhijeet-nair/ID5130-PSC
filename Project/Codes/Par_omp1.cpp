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

void printVector(double a[], int n) {
   for (int i = 0; i < n; i++) {
        printf("%.4f\n",a[i]);
   }
}

void printMatrix (double** a, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.4f\t",a[i][j]);
        }
        printf("\n");
    }
}



int main (int argc, char* argv[]) {
    int np = 1, cnt, fs = 1;

    if (argc == 2) {
        np = strtol(argv[1], NULL, 10);
    }
    printf("Running in code with %d processors...\n",np);
    
    int i, j, k, m, p;  // Indices

    int Nl     = 1000;    // No. of collocation points
    double u   = 20;    // Freestream velocity
    double c   = 10;     // Chord length of the flat plate
    double alp = 0;     // Angle of attack of the plate
    double rho = 1.225; // Density
    double dx  = c/Nl;  // Spacing between two points

    double sn = sin(deg2rad(alp));
    double cs = cos(deg2rad(alp));

    double dt = 0.005;        // Time step
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
    double gIVRes[2], dgbm_dt, vindw, xw[Nt] {}, yw[Nt] {};
    double t_cur, sum, Bjs, uw, vw, clx, cly;
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
    double* extflow = new double[Nt];
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

    // LU Decomposition and save in the same matrix
    for (k = 0; k <= Nl; k++) {
        for (i = k+1; i <= Nl; i++) {
            C[i][k] = C[i][k]/C[k][k];

            for (j = k+1; j <= Nl; j++) {
                C[i][j] += -C[i][k]*C[k][j];
            }
        }
    }

    // For loop for time marching
    #pragma omp parallel num_threads(np) default(shared) private(i, j, k, m, p, t_cur, gIVRes, clx, cly, Bjs)
    {
        int myid = omp_get_thread_num();

        // Origin location at time t
        #pragma omp for
            for (i = 0; i < Nt; i++) {
                t_cur = i*dt;
                x0[i] = -u*t_cur;
                y0[i] = h(t_cur, f1, f2, a);
                extflow[i] = u*sn - hdot(t_cur, f1, f2, a)*cs;
            }
        

        for (m = 0; m < Nt; m++) {
            // if (myid == 0) {
            //     if (m % 50 == 0) {printf("m = %3.0f\n",double(m));}

            //     // New wake location
            //     xw[m] = x0[m] + (c + 0.1*dx)*cs;
            //     yw[m] = y0[m] - (c + 0.1*dx)*sn;

            //     R[Nl] = 0;
            //     printf("xw = %.4f\tyw = %.4f\n",xw[m],yw[m]);
            // }
            #pragma omp single
            {
                if (m % 50 == 0) {printf("m = %3.0f\n",double(m));}

                // New wake location
                xw[m] = x0[m] + (c + 0.1*dx)*cs;
                yw[m] = y0[m] - (c + 0.1*dx)*sn;

                R[Nl] = 0;
                // printf("xw = %.4f\tyw = %.4f\n",xw[m],yw[m]);
            }

            // #pragma omp barrier

            #pragma omp for
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
                    // printf("i = %d\tBjs = %.4f\n",i,Bjs);
                    R[i] = -extflow[m] - Bjs;
                }
            
            #pragma omp barrier

            #pragma omp for reduction(+: R[Nl])
                for (j = 0; j < m; j++) {
                    R[Nl] += gw[j];
                    // printf("gw[%d] = %.4f\n",j,gw[j]);
                }
            
            // if (myid == 0) {
            //     R[Nl] = -R[Nl];
            //     // printVector(R, Nl+1);

            //     // Computation of unknown gbm and gwm
            //     // Forward substitution
            //     U1[0] = R[0];
            //     for (i = 1; i <= Nl; i++) {
            //         sum = 0;
            //         for (j = 0; j < i; j++){
            //             sum += C[i][j]*U1[j];
            //         }
            //         U1[i] = R[i] - sum;
            //     }

            //     // Backward substitution
            //     U[Nl] = U1[Nl]/C[Nl][Nl];
            //     for (i = Nl-1; i >= 0; i--) {
            //         sum = 0;
            //         for (j = i+1; j <= Nl; j++) {
            //             sum += C[i][j]*U[j];
            //         }
            //         U[i] = (U1[i] - sum)/C[i][i];
            //     }

            //     // Solved wake circulation
            //     gw[m] = U[Nl];
            //     printVector(U, Nl+1);
            // }

            #pragma omp single
            {
                R[Nl] = -R[Nl];
                // printf("R = \n");
                // printVector(R, Nl+1);

                // Computation of unknown gbm and gwm
                // Forward substitution
                U1[0] = R[0];
                for (i = 1; i <= Nl; i++) {
                    sum = 0;
                    for (j = 0; j < i; j++){
                        sum += C[i][j]*U1[j];
                    }
                    U1[i] = R[i] - sum;
                }

                // Backward substitution
                U[Nl] = U1[Nl]/C[Nl][Nl];
                for (i = Nl-1; i >= 0; i--) {
                    sum = 0;
                    for (j = i+1; j <= Nl; j++) {
                        sum += C[i][j]*U[j];
                    }
                    U[i] = (U1[i] - sum)/C[i][i];
                }

                // Solved wake circulation
                gw[m] = U[Nl];
                // printf("U = \n");
                // printVector(U, Nl+1);
            }

            // #pragma omp barrier

            // Solved body circulations
            #pragma omp for
                for (i = 0; i < Nl; i++) {
                    gbm[1][i] = U[i];
                }
            
            for (k = 0; k <= m; k++) {
                #pragma omp single
                {
                    uw = 0;
                    vw = 0;
                }

                // #pragma omp barrier

                #pragma omp for reduction(+: uw, vw)
                    for (i = 0; i < Nl; i++) {
                        getIndVel(gbm[1][i], xw[k], yw[k], vorloc[i]*cs + x0[m], -vorloc[i]*sn + y0[m], gIVRes);
                        uw += gIVRes[0];
                        vw += gIVRes[1];
                    }
                
                #pragma omp for reduction(+: uw, vw)
                    for (p = 0; p <= m; p++) {
                        if (p != k) {
                            getIndVel(gw[p], xw[k], yw[k], xw[p], yw[p], gIVRes);
                            uw += gIVRes[0];
                            vw += gIVRes[1];
                        }
                    }
                    
                #pragma omp single
                {
                    xw[k] += uw*dt;
                    yw[k] += vw*dt;
                    xwM[m][k] = xw[k];
                    ywM[m][k] = yw[k];
                    // printf("uw = %.4f\tvw = %.4f\n",uw,vw);
                    // printVector(xw, 4);
                    // printf("\n");
                    // printVector(yw, 4);
                    // printf("\n\n");
                }
            }

            // Aerodynamic Load Calculations
            if (myid == 0) {dgbm_dt = 0;}

            #pragma omp barrier

            if (m == 1) {
                #pragma omp for reduction(+: dgbm_dt)
                    for (i = 0; i < Nl; i++) {
                        dgbm_dt += gbm[1][i];
                    }
            }
            else {
                #pragma omp for reduction(+: dgbm_dt)
                    for (i = 0; i < Nl; i++) {
                        dgbm_dt += gbm[1][i] - gbm[0][i];
                    }
            }
            
            if (myid == 0) {
                dgbm_dt /= dt;
                L[m] = 0;
                D[m] = 0;
            }

            #pragma omp barrier

            #pragma omp for reduction(+: L[m], D[m])
                for (i = 0; i < Nl; i++) {
                    gbm[0][i] = gbm[1][i];
                    vindw = 0;

                    // Collocation point location
                    clx = colcloc[i]*cs + x0[m];
                    cly = colcloc[i]*sn + y0[m];

                    for (p = 0; p <= m; p++) {
                        getIndVel(gw[p], clx, cly, xw[p], yw[p], gIVRes);
                        vindw += gIVRes[0]*sn + gIVRes[1]*cs;
                    }


                    L[m] += gbm[1][i];
                    D[m] += vindw*gbm[1][i];
                }
            
            #pragma omp single
            {
                L[m] = rho*(u*L[m] + dgbm_dt*c);
                D[m] = rho*(D[m] + dgbm_dt*c*deg2rad(alp));
            }

            // L[m] = rho*(u*L[m] + dgbm_dt*c);
            // D[m] = rho*(D[m] + dgbm_dt*c*deg2rad(alp));
        } // End of time loop
    } // End of parallel

    // printf("extflow = \n");
    // printVector(extflow, 4);
    // printMatrix(xwM, 10, 10);

    if (fs == 1) {
        // File writing
        printf("Saving in file...\n");
        char fname[25];
        sprintf(fname, "./Res/Par_OMP.txt");

        std::ofstream oFile(fname);

        if (oFile.is_open()) {
            oFile << Nl << "\n" << Nt << "\n";
            oFile << dt << "\n" << dx << "\n";
            oFile << rho << "\n" << u << "\n";
            oFile << c << "\n" << alp << "\n";

            oFile << "\n";

            for (i = 0; i < Nt; i++) {
                oFile << x0[i] << "," << y0[i] << "\n";
            }

            oFile << "\n";

            for (i = 0; i < Nt; i++) {
                oFile << L[i] << "," << D[i] << "\n";
            }

            oFile << "\n\n";

            for (i = 0; i < Nt; i++) {
                for (j = 0; j < Nt; j++) {
                    oFile << xwM[i][j] << "," << ywM[i][j] << "\n";
                }
                oFile << "\n";
            }

            oFile.close();
            printf("Saved in file %s\n",fname);
        }
        else {
            printf("Error opening file\n");
        }
    }
    else {printf("Not saving in file...\n");}

    delete B, U, U1, A, C, x0, y0, extflow, L, D, xwM, ywM, R, gw, gbm;
    return 0;
}