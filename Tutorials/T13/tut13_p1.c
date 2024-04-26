#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 		100
#define PROBZERO	33

int main() {
	int a[N][N];	// cities X localities
	int x[N], 	// sum of infected people for all the localities of a city
	    y[N];	// number of localities in a city

	srand(time(NULL));

	#pragma acc parallel loop
	for (int ii = 0; ii < N; ++ii) {
		x[ii] = 0;
		for (int jj = 0; jj < N; ++jj) {
			a[ii][jj] = rand() % N;	//ii * N + jj;
			if (rand() % 100 < PROBZERO) a[ii][jj] = 0;	// add some zeros to make life interesting.
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
