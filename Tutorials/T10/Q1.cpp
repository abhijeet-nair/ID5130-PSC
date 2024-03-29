#include<iostream>
#include<math.h>
#include<mpi.h>

double myfunc (double x) {return x*tan(x);}

double myderv (double x) {return (tan(x) + x/pow(cos(x),2));}

double forward (double fi1, double fi, double delx) {return (fi1 - fi)/delx;}

double backward (double fi, double fi_1, double delx) {return (fi - fi_1)/delx;}

double central2nd (double fi1, double fi_1, double delx) {return 0.5*(fi1 - fi_1)/delx;}


int main (int argc, char* argv[]) {
    int myid, np, i, j, n {}, ln {};
    double dat[2] {}, la {}, lb {}, delx {}, xi {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (myid == 0) {
        double a {}, b{};
        printf("Enter start: ");
        std::cin >> a;

        printf("Enter stop: ");
        std::cin >> b;

        printf("Enter delx: ");
        std::cin >> delx;

        n = (b - a)/delx + 1;
        ln = n/np;
        dat[0] = delx;
        dat[1] = ln;

        // printf("i = 0\ta = 1\tln = %d\n",ln);
        for (i = 1; i < np; i++) {
            la = a + delx*(n - ln*(np - i));
            // dat[0] = la;
            // MPI_Send(&dat, 3, MPI_DOUBLE, i, i*10, MPI_COMM_WORLD);
            MPI_Send(&la, 1, MPI_DOUBLE, i, i*10, MPI_COMM_WORLD);
        }
        la  = a;
        ln = n - (np - 1)*ln;
        // lb  = a + delx*ln;
    }
    else {
        // MPI_Recv(&dat, 3, MPI_DOUBLE, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // la   = dat[0];
        // delx = dat[1];
        // ln   = int(dat[2]);
        // lb   = la + delx*ln;
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

    for (i = 0; i < ln; i++) {
        xi = la + i*delx;
        fvals[i] = myfunc(xi);
        fdvact[i] = myderv(xi);
    }

    if (myid == 0) {
        fdvals[0] = forward(fvals[1], fvals[0], delx);

        prtnr[0] = MPI_PROC_NULL;
        prtnr[1] = 1;
    }
    else if (myid == (np - 1)) {
        fdvals[ln] = backward(fvals[ln-1], fvals[ln-2], delx);

        prtnr[1] = MPI_PROC_NULL;
        prtnr[0] = np - 2;
    }

    MPI_Finalize();
    return 0;
}