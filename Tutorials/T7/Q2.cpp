#include<iostream>
#include<string>
#include<mpi.h>


int main(int argc, char** argv)
{
    int i, myid, thrd_cnt;
    int val {}, final_sum {};

    MPI_Status status; 		/* data type that is defined in mpi.h... */

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &thrd_cnt); /* tells about the number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* will return the rank or id of the process that called it.  */

    val = myid*10;

    if (myid == 0) {
        final_sum = val;
        for (int i = 1; i < thrd_cnt; i++) {
            MPI_Recv(&val, 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            final_sum += val;
        }

        printf("Final sum = %d\n",final_sum);
    }
    else {
        MPI_Send(&val, 1, MPI_INT, 0, myid, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}