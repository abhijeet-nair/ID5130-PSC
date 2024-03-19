#include<iostream>
#include<math.h>
#include<mpi.h>

int main (int argc, char* argv[]) {
    int n;
    int myid, np, i;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // MPI_Reduce already used in T8, Q4.cpp. So, not doing here.

    // Use && in the run scripts to avoid run when compilation throws errors.

    if (myid == 0) {
        printf("Enter size of the array: ");
        std::cin >> n;

        int arr[n] {};
        // printf("Enter elements of the array:\n");
        for (i = 0; i < n; i++) {
            // printf("arr[%d]: ", i);
            // std::cin >> arr[i];
            arr[i] = i + 1;
        }

        MPI_Scatter(&arr,v )
    }

    MPI_Finalize();
    return 0;
}