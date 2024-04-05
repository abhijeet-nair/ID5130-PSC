#include<iostream>
#include<math.h>
#include<mpi.h>


int main (int argc, char* argv[]) {
    int myid, np, i, l {}, m {}, n {};
    double dat[3] {};
    int lli {}, lmi {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (myid == 0) {
        printf("Enter no. of rows of A: ");
        std::cin >> l;

        printf("Enter no. of cols of A: ");
        std::cin >> m;

        printf("Enter no. of cols of B: ");
        std::cin >> n;

        dat[0] = l;
        dat[1] = m;
        dat[2] = n;
        printf("\n");
    }

    MPI_Bcast(&dat, 3, MPI_INT, 0, MPI_COMM_WORLD);

    l = dat[0];
    m = dat[1];
    n = dat[2];

    lli = int(l/np);
    lmi = int(m/np);

    MPI_Finalize();
    return 0;
}