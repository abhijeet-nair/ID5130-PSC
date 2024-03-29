#include<iostream>
#include<math.h>
#include<mpi.h>

int main (int argc, char* argv[]) {
    int myid, np, i, j, m{}, n {}, lm {}, lmi {};
    int dat[2] {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (myid == 0) {
        printf("Enter no. of rows: ");
        std::cin >> m;

        printf("Enter no. of cols: ");
        std::cin >> n;

        dat[0] = m;
        dat[1] = n;
    }
    
    MPI_Bcast(&dat, 2, MPI_INT, 0, MPI_COMM_WORLD);

    m = dat[0];
    n = dat[1];
    lmi = int(m/np);

    if (myid == 0) {lm = m - lmi*(np - 1);}
    else {lm = lmi;};

    int lA[lm][n] {};
    int lv[lm] {};

    // Considering each processor fills the data also...
    // Can be modified for the other case as well...
    for (i = 0; i < lm; i++) {
        lv[i] = myid*10;
        for (j = 0; j < n; j++){
            lA[i][j] = myid*(j+1)*10;
        }
    }

    int v[m] {}, cnts[np] {}, displs[np] {};

    // printf("lm = %d\n",lm);
    for (i = 0; i < np; i++) {
        if (i == 0) {
            cnts[i] = m - lmi*(np - 1);
        }
        else {
            cnts[i] = lmi;
            displs[i] = cnts[0] + (i - 1)*lmi;
        }
        // if (myid == 0) {
        //     printf("cnts[%d] = %d\tdispls[%d] = %d\n",i,cnts[i],i,displs[i]);
        // }
    }
    MPI_Allgatherv(&lv, lm, MPI_INT, &v, cnts, displs, MPI_INT, MPI_COMM_WORLD);

    if (myid == 0) {
        for (i = 0; i < m; i++) {
            printf("v[%2.0f] = %d\n",double(i),v[i]);
        }
    }
    MPI_Finalize();
    return 0;
}