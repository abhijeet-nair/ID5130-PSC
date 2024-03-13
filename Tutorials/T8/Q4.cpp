#include<iostream>
#include<math.h>
#include<mpi.h>

double func (double x) {
    return (0.5*sin(x)/pow(x,3));
}

double myTrapz (double la, double lb, double h, int ln) {
    double sum {}, x {};
    int i;

    sum = 0.5*(func(la) + func(lb));

    for (i = 1; i <= ln - 1; i++) {
        x = la + i*h;
        sum += func(x);
    }

    sum *= h;
    return sum;
}


int main (int argc, char* argv[]) {
    double a, b, fSum {}, la, lb, lsum, h;
    int myid, np, i, n, ln;
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

        ln = int(n/np);
        h = (b - a)/n;

        for (i = 1; i < np; i++) {
            dat[0] = a + (i - 1)*ln*h;
            dat[1] = a + i*ln*h;
            dat[2] = h;
            dat[3] = double(ln);

            MPI_Send(&dat, 4, MPI_DOUBLE, i, i*10, MPI_COMM_WORLD);
        }

        la = a + (np - a)*ln*h;
        lb = b;
        ln = n - ln*(np - 1);
    }
    else {
        MPI_Recv(&dat, 4, MPI_DOUBLE, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        la = dat[0];
        lb = dat[1];
        h  = dat[2];
        ln = int(dat[3]);
    }

    lsum = myTrapz(la, lb, h, ln);

    MPI_Reduce(&lsum, &fSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("\nResult = %.6f\n",fSum);
    }

    MPI_Finalize();
    return 0;
}