#include<iostream>
#include<math.h>
#include<mpi.h>

int main (int argc, char* argv[]) {
    int myid, np, i, j, k, m{}, n {}, lm {}, lmi {}, ln {}, lni {};
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
        printf("\n");
    }
    
    MPI_Bcast(&dat, 2, MPI_INT, 0, MPI_COMM_WORLD);

    m = dat[0];
    n = dat[1];
    lmi = int(m/np);
    lni = int(n/np);

    if (myid == 0) {
        lm = m - lmi*(np - 1);
        ln = n - lni*(np - 1);
    }
    else {
        lm = lmi;
        ln = lni;
    }

    int lA[lm][n] {};
    int lv[ln] {};
    int lr[lm] {};

    // Considering each processor fills the data also...
    // Can be modified for the other case as well...
    for (i = 0; i < lm; i++) {
        for (j = 0; j < n; j++){
            lA[i][j] = (myid + 1)*(j+1)*10;
        }
    }
    
    for (i = 0; i < np; i++) {
        if (myid == i) {
            for (j = 0; j < lm; j++) {
                for (k = 0; k < n; k++) {
                    printf("%2.0f\t",double(lA[j][k]));
                }
                printf("\n");
            }
        }
    }

    for (i = 0; i < ln; i++) {
        lv[i] = myid*10;
    }

    int v[n] {}, cnts[np] {}, displs[np] {}, r[m] {};

    // printf("lm = %d\n",lm);
    for (i = 0; i < np; i++) {
        if (i == 0) {
            cnts[i] = n - lni*(np - 1);
        }
        else {
            cnts[i] = lni;
            displs[i] = cnts[0] + (i - 1)*lni;
        }
        // if (myid == 0) {
        //     printf("cnts[%d] = %d\tdispls[%d] = %d\n",i,cnts[i],i,displs[i]);
        // }
    }
    MPI_Allgatherv(&lv, ln, MPI_INT, &v, cnts, displs, MPI_INT, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("\n");
        for (i = 0; i < n; i++) {
            printf("v[%2.0f] = %d\n",double(i),v[i]);
        }
    }

    for (i = 0; i < lm; i++) {
        for (j = 0; j < n; j++) {
            lr[i] += lA[i][j]*v[j];
        }
    }

    // for (i = 0; i < np; i++) {
    //     if (myid == i) {
    //         for (j = 0; j < lm; j++) {
    //             printf("r[%2.0f] = %d\n",double(i),lr[j]);
    //         }
    //     }
    // }

    for (i = 0; i < np; i++) {
        if (i == 0) {
            cnts[i] = m - lmi*(np - 1);
        }
        else {
            cnts[i] = lmi;
            displs[i] = cnts[0] + (i - 1)*lmi;
        }
    }

    MPI_Allgatherv(&lr, lm, MPI_INT, &r, cnts, displs, MPI_INT, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("\n");
        for (i = 0; i < m; i++) {
            printf("r[%2.0f] = %d\n",double(i),r[i]);
        }
    }

    MPI_Finalize();
    return 0;
}