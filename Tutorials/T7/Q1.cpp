#include<iostream>
#include<string>
#include<mpi.h>


int main(int argc, char** argv)
{
    int i, myid, thrd_cnt, tag = 100, msg_size = 100;
    char msg_send[msg_size], msg_recv[msg_size];
    // char message_send[50], message_recv[50];
    MPI_Status status; 		/* data type that is defined in mpi.h... */

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &thrd_cnt); /* tells about the number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* will return the rank or id of the process that called it.  */

    if (myid == 0) {
        for (int i = 1; i < thrd_cnt; i++) {
            sprintf(msg_send, "Hello from process %d!!! Received this from process 0.\n", i);
            // sprintf(msg_send, )
            MPI_Send(msg_send, msg_size, MPI_CHAR, i, i*tag, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(msg_recv, msg_size, MPI_CHAR, 0, myid*tag, MPI_COMM_WORLD, &status);
        printf("\n%s", msg_recv);
        printf("Process id: %d... Confirmed!!!\n", myid);
    }

    MPI_Finalize();

    return 0;
}