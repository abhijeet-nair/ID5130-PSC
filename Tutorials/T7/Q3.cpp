#include<iostream>
#include<string>
#include<mpi.h>


int main(int argc, char** argv)
{
    int i, j, myid, p_cnt, ln;
    int sum {}, final_sum {}, N {};

    MPI_Status status; 		/* data type that is defined in mpi.h... */

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &p_cnt); /* tells about the number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* will return the rank or id of the process that called it.  */

    if (myid == 0) {
        printf("Enter size of array: ");
        std::cin >> N;
        double a[N] {};
        ln = N/p_cnt;
        double la[ln] {};

        for (i = 1; i < p_cnt; i++) {
            MPI_Send(&ln, 1, MPI_INT, i, i*10, MPI_COMM_WORLD);

            // for (j = 0; j < ln; j++) {
            //     la[j] = a[(i-1)*ln + j];
            // }
            memcpy(la, a + (i-1)*ln, ln*sizeof(double));
            MPI_Send(&la, ln, MPI_DOUBLE, i, i*11, MPI_COMM_WORLD);
        }

        for (i = 1; i < p_cnt; i++) {
            MPI_Recv(&la, ln, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // for (j = 0; j < ln; j++) {
            //     a[(i-1)*ln + j] = la[j];
            // }
            memcpy(a + (i-1)*ln, la, ln*sizeof(double));
        }

        for (i = (p_cnt - 1)*ln; i < N; i++) {
            a[i] = myid + 10;
        }

        printf("Updated array is:\n");
        for (i = 0; i < N; i++) {
            printf("%.4f\n",a[i]);
        }
    }
    else {
        MPI_Recv(&ln, 1,  MPI_INT, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        double la[ln] {};
        MPI_Recv(&la, ln, MPI_DOUBLE, 0, myid*11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (i = 0; i < ln; i++) {
            la[i] = myid + 10;
        }
        
        MPI_Send(&la, ln, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}