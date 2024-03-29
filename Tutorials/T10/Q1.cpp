#include<iostream>
#include<math.h>
#include<mpi.h>

double myfunc (double x) {return x*tan(x);}

double myderv (double x) {return (tan(x) + x/pow(cos(x),2));}

double forward (double fi1, double fi, double delx) {return (fi1 - fi)/delx;}

double backward (double fi, double fi_1, double delx) {return (fi - fi_1)/delx;}

double central2nd (double fi1, double fi_1, double delx) {return 0.5*(fi1 - fi_1)/delx;}


int main (int argc, char* argv[]) {
    int myid, np, i, n {}, ln {};
    double dat[2] {}, la {}, lb {}, delx {}, xi {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int cnts[np] {}, dsplc[np] {};

    if (myid == 0) {
        double a {}, b{};
        printf("Enter start: ");
        std::cin >> a;

        printf("Enter stop : ");
        std::cin >> b;

        printf("Enter delx : ");
        std::cin >> delx;

        n = (b - a)/delx + 1;
        ln = n/np;
        dat[0] = delx;
        dat[1] = ln;

        for (i = 1; i < np; i++) {
            la = a + delx*(n - ln*(np - i));
            MPI_Send(&la, 1, MPI_DOUBLE, i, i*10, MPI_COMM_WORLD);
            cnts[i] = ln;
            dsplc[i] = n - ln*(np - i);
        }
        la  = a;
        ln = n - (np - 1)*ln;
        cnts[0] = ln;
    }
    else {
        MPI_Recv(&la, 1, MPI_DOUBLE, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    MPI_Bcast(&dat, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myid != 0) {
        delx = dat[0];
        ln   = int(dat[1]);
    }

    double fvals[ln] {};
    double fdvals[ln] {};
    double fdvact[ln] {};
    double prtnr[2] {};
    double fextra[2] {};

    for (i = 0; i < ln; i++) {
        xi = la + i*delx;
        fvals[i] = myfunc(xi);
        fdvact[i] = myderv(xi);
    }

    for (i = 1; i < ln - 1; i++) {
        fdvals[i] = central2nd(fvals[i+1], fvals[i-1], delx);
    }

    if (myid == 0) {
        prtnr[0] = MPI_PROC_NULL;
        prtnr[1] = 1;
    }
    else if (myid == (np - 1)) {
        prtnr[0] = np - 2;
        prtnr[1] = MPI_PROC_NULL;
    }
    else {
        prtnr[0] = myid - 1;
        prtnr[1] = myid + 1;
    }

    MPI_Sendrecv(&fvals[0], 1, MPI_DOUBLE, prtnr[0], 0, &fextra[0], 1, MPI_DOUBLE, prtnr[0], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&fvals[ln-1], 1, MPI_DOUBLE, prtnr[1], 0, &fextra[1], 1, MPI_DOUBLE, prtnr[1], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (myid == 0) {
        fdvals[0] = forward(fvals[1], fvals[0], delx);
        fdvals[ln-1] = central2nd(fextra[1], fvals[ln-2], delx);
    }
    else if (myid == (np - 1)) {
        fdvals[0] = central2nd(fvals[1], fextra[0], delx);
        fdvals[ln-1] = backward(fvals[ln-1], fvals[ln-2], delx);
    }
    else {
        fdvals[0] = central2nd(fvals[1], fextra[0], delx);
        fdvals[ln-1] = central2nd(fextra[1], fvals[ln-2], delx);
    }

    double *res, *act;
    if (myid == 0) {
        res = (double *)malloc(n*sizeof(double));
        act = (double *)malloc(n*sizeof(double));
    }
    
    MPI_Gatherv(&fdvals, ln, MPI_DOUBLE, res, cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&fdvact, ln, MPI_DOUBLE, act, cnts, dsplc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("cnt\tcalcv\t\tactv\n");

        for (i = 0; i < n; i++) {
            printf("%3.0f\t%.4f\t\t%.4f\n",double(i),res[i],act[i]);
        }    
    }

    MPI_Finalize();
    return 0;
}