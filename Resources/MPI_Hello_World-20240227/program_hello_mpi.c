#include<stdio.h>
#include<string.h>
#include<mpi.h> 		/* new and contains the function prototypes, macros, variables.. */

int main(int argc, char** argv)
{
  int i, myid, size, tag=100;
  char message_send[50], message_recv[50];
  MPI_Status status; 		/* data type that is defined in mpi.h... */

  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* tells about the number of processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* will return the rank or id of the process that called it.  */

  /* printf("%d", size); */

  if (myid != 0) 		/* not the designated process to receive..send the data.. */
    {
      sprintf(message_send, " Hello from process id: %d \n ", myid);
      MPI_Send(message_send, 50, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
    }
  else
    {				/* only by process id = 0, myrank = 0, myid = 0 */
      for (i = 1; i < size; i++)
	{
	  MPI_Recv(message_recv, 50, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);
	  printf("\n %s", message_recv); 
	}

      sprintf(message_send, " Hello from process id: %d \n ", myid);
      printf("\n %s", message_send); 

    }

  MPI_Finalize();

  return 0;
}
