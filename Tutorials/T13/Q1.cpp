#include <iostream>
#include <time.h>

int main () {
    int N = 100;
    int probzero = 33;

    srand(time(NULL));

    int a[N][N]; // cities X localities
    int x[N];    // sum of infected people for all the localities of a city
    int y[N];    // number of localities in a city
    int rv[N] {}; // rand vector

    for (int ii = 0; ii < N; ii++) {
        rv[ii] = rand();
    }

    #pragma acc parallel loop
        for (int ii = 0; ii < N; ++ii) {
            x[ii] = 0;
            for (int jj = 0; jj < N; ++jj) {
                a[ii][jj] = rv[ii] % N;	//ii * N + jj;
                if (rv[ii] % 100 < probzero) a[ii][jj] = 0;	// add some zeros to make life interesting.
		    }
        }

    #pragma acc parallel loop collapse(2) reduction(+:x)
        for (int ii = 0; ii < N; ++ii) {
            for (int jj = 0; jj < N; ++jj) {
                x[ii] += a[ii][jj];
                if (a[ii][jj] > 0) y[ii]++;
            }
        }
    
    for (int ii = 0; ii < N; ++ii)
		printf("%.0f ", x[ii] * 100.0 / (y[ii] * N));
	printf("\n");

    return 0;
}