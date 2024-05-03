/* red black solver - openacc parallel version */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 1000
#define BILLION 1000000000L
#define MIN(a,b) ( ((a)<(b))?(a):(b) )

void initialize(double matrix[][N])
{
  
#pragma acc parallel present(matrix[0:N][0:N])
  {
#pragma acc loop collapse(2) 
    for(int i = 0; i < N; i++)
      for(int j = 0; j < N; j++)
	matrix[i][j] = 100.0;

    //#pragma acc parallel loop default(present)
#pragma acc loop
    for (int i = 0; i < N; i++)
      {
	matrix[i][0] = 0.0;
	matrix[i][N-1] = 0.1*i;
      }

    //#pragma acc parallel loop default(present)
#pragma acc loop
    for (int j = 0; j < N; j++)
      {
	matrix[0][j] = 0.0;
	matrix[N-1][j] = 0.1*j;
      }
  }
  
    return;
}

void nullify(double matrix[][N]){
#pragma acc parallel loop present(matrix[0:N][0:N]) collapse(2)
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      matrix[i][j] = 0.0;
  
  return;
}


double red_black_solver(double A[][N])
{
  double error;

  error = 0.0;
#pragma acc data copy(error)
  {
#pragma acc parallel loop collapse(2) present(A[0:N][0:N]) present(error) reduction(+:error)
    for (int i = 1; i < N-1; i++)
      for (int j = 1; j < N-1; j++)
	if ((i+j)%2 == 0)
	  {
	    double temp  = A[i][j];
	    A[i][j] = 0.25*(A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1]);
	    error += fabs(A[i][j] - temp);
	  }
  
#pragma acc parallel loop collapse(2) present(A[0:N][0:N]) present(error) reduction(+:error)
    for (int i = 1; i < N-1; i++)
      for (int j = 1; j < N-1; j++)
	if ((i+j)%2 == 1)
	  {
	    double temp  = A[i][j];
	    A[i][j] = 0.25*(A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1]);
	    error += fabs(A[i][j] - temp);
	  }
  }
  
  error = error/(N*N);
  
  
  //printf("\n Error = %lf", error);

  return error;
  
}

void printmatrix(double matrix[N][N]){

  if (N<=10)
#pragma acc serial loop present(matrix[0:N][0:N])
    for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++)
	printf("%lf ", matrix[i][j]);
      printf("\n");
    }
  else
    printf("Skipping printing of large matrix...\n");
  
    
  return;
}

int main()
{
  struct timespec start_t, end_t;
  uint64_t diff;

  double T[N][N]; 
  double error;
  int iter = 0; 

  clock_gettime(CLOCK_MONOTONIC, &start_t);	/* mark start time */  

#pragma acc data create(T[0:N][0:N]) copyout(T[0:N][0:N])
  {
    initialize(T);
    
    do
      {
	error = red_black_solver(T);
	iter = iter + 1 ;
	
	if (iter%100==0)
	  printf ("iter = %d, error = %15.5f \n", iter, error);
	
      } while(iter < 3000);

  }
    
  printf ("iter = %d, error = %15.5f \n", iter, error);

  clock_gettime(CLOCK_MONOTONIC, &end_t);	/* mark the end time */
  
  printf("\n sanity check value for T[1][1] = %lf \n", T[1][1]);
  printf("\n sanity check value for T[N-5][N-5] = %lf \n", T[N-5][N-5]);

  diff = BILLION * (end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
  printf("elapsed time = %lf seconds \n", (double) diff/1000000000.0);

  return 0;
}
 
  


