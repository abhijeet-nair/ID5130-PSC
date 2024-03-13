#include<iostream>
#include<math.h>
#include<mpi.h>

double func (double x) {
    return (0.5*sin(x)/pow(x,3));
}

double mySimps (double a, double h, int is, int ie) {
    double sum {}, x {};

    for (int i = is; i <= ie; i++) {
        x = a + i*h;

        if (i % 2 == 0) {
            sum += 2*func(x);
        }
        else {
            sum += 4*func(x);
        }
    }

    return sum;
}

int main (int argc, char* argv[]) {
    double a, b, fSum {}, lsum, h;
    int myid, np, i, n, is, ie, ln;
    double dat[4] {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (myid == 0) {
        printf("Enter start location a: ");
        std::cin >> a;

        printf("Enter end location b  : ");
        std::cin >> b;

        printf("Enter no. of trapz n  : ");
        std::cin >> n;

        h = (b - a)/n;
        ln = int((n - 2)/np);
        dat[0] = a;
        dat[1] = h;

        for (i = 1; i < np; i++) {
            is = (i - 1)*ln + 1;
            ie = i*ln;

            dat[2] = is;
            dat[3] = ie;

            MPI_Send(&dat, 4, MPI_DOUBLE, i, i*10, MPI_COMM_WORLD);
        }
        
        is = (np - 1)*ln + 1;
        ie = n - 1;
    }
    else {
        MPI_Recv(&dat, 4, MPI_DOUBLE, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        a  = dat[0];
        h  = dat[1];
        is = int(dat[2]);
        ie = int(dat[3]);
    }

    lsum = mySimps(a, h, is, ie);

    MPI_Reduce(&lsum, &fSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    fSum += func(a) + func(b);
    fSum *= h/3;

    if (myid == 0) {
        printf("\nResult = %.6f\n",fSum);
    }

    MPI_Finalize();
    return 0;
}