#include<iostream>
#include<math.h>
#include<mpi.h>

double myfunc (double x) {
    return x*tan(x);
}

double myderv (double x) {
    return (tan(x) + x/pow(cos(x),2));
}

double forward (double fi1, double fi, double delx) {return (fi1 - fi)/delx;}

double backward (double fi, double fi_1, double delx) {return (fi - fi_1)/delx;}

double central2nd (double fi1, double fi_1, double delx) {return 0.5*(fi1 - fi_1)/delx;}

int main (int argc, char* argv[]) {
    int n;
    int myid, np, i, j;
    double a {}, b{}, delx {};
    int ln {}, ls {}, dat[2] {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (myid == 0) {
        printf("Enter start: ");
        std::cin >> a;

        printf("Enter stop: ");
        std::cin >> b;

        printf("Enter delx: ");
        std::cin >> delx;

        int n = (b - a)/delx + 1;
        ln = n/np;
        dat[0] = ln;

        // printf("i = 0\ta = 1\tln = %d\n",ln);
        for (i = 1; i < np; i++) {
            ls = n - ln*(np - i);
            dat[1] = ls;
            MPI_Send(&dat, 2, MPI_INT, i, i*10, MPI_COMM_WORLD);
        }

    }
    else {
        MPI_Recv(&dat, 2, MPI_INT, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Finalize();
    return 0;
}