#include<iostream>
#include<math.h>
#include<mpi.h>

int main (int argc, char* argv[]) {
    int n;
    int myid, np, i, j;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // MPI_Reduce already used in T8, Q4.cpp. So, not doing here.

    // Use && in the run scripts to avoid run when compilation throws errors.

    if (myid == 0) {
        printf("Enter size of the array: ");
        std::cin >> n;
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int *send, *recv;

    if (myid == 0) {
        send = new int[n];

        for (i = 0; i < n; i++) {
            send[i] = i + 1;
            printf("send[%d]: %d\n",i, send[i]);
        }        
    }
    int ln = n/np;
    MPI_Scatter(send, ln, MPI_INT, recv, ln, MPI_INT, 0, MPI_COMM_WORLD);

    for (i = 0; i < np; i++) {
        if (myid == i) {
            for (j = 0; j < ln; j++) {
                printf("arr[%d, %d]:%d\n", myid, j, recv[i]);
            }
        }
    }

    MPI_Finalize();
    return 0;
}