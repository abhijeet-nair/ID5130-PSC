/* matrix multiplication - openacc parallel version */
#include <stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include <time.h>

#define N 32
#define BILLION 1000000000L
#define MIN(a,b) ( ((a)<(b))?(a):(b) )

void initialize(double matrix[][N]){
#pragma acc parallel loop present(matrix[0:N][0:N]) collapse(2)
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      matrix[i][j] = (i+1)*N+(j+1);
  
  return;
}

void nullify(double matrix[][N]){
#pragma acc parallel loop present(matrix[0:N][0:N]) collapse(2)
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      matrix[i][j] = 0.0;
  
  return;
}


void multiplymatrices(double A[][N], double B[][N], double C[][N])
{
  int tile = 32;
  
#pragma acc parallel present(A[0:N][0:N],B[0:N][0:N],C[0:N][0:N]) num_workers(8) vector_length(32)
#pragma acc loop gang collapse(2)
    for (int i = 0; i < N; i+= tile)
      {
       for (int j = 0; j < N; j+= tile)
	 {
#pragma acc loop worker
	  for (int ii = i; ii < MIN(i+tile, N); ii++)
	    {
#pragma acc loop vector
	     for (int jj = j; jj < MIN(j+tile, N); jj++)
	       {
#pragma acc loop seq
		for (int k = 0; k < N; k++)
		  C[ii][jj] += A[ii][k] * B[k][jj];
	       }
	    }
	 }
      }
  
  return;
}

void printmatrix(double matrix[N][N]){

  if (N<=10)
#pragma acc serial loop present(matrix[0:N][0:N])
    for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++)
	printf("%5d ", matrix[i][j]);
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

  double A[N][N]={0.0}, B[N][N]={0.0}, C[N][N]={0.0};

  clock_gettime(CLOCK_MONOTONIC, &start_t);	/* mark start time */  

  
#pragma acc data create(A[0:N][0:N],B[0:N][0:N],C[0:N][0:N]) copyout(C[0:N][0:N])
  {
    initialize(A);
    initialize(B);
    nullify(C);
    
    printf("\n Matrix A is:\n");
    printmatrix(A);
    printf("\n Matrix B is:\n");
    printmatrix(B);
    
    multiplymatrices(A, B, C);
    
    printf("\n Matrix C is:\n");
    printmatrix(C);
  }

  clock_gettime(CLOCK_MONOTONIC, &end_t);	/* mark the end time */
  
  printf("\n sanity check value for C[0][0] = %lf \n", C[0][0]);

  diff = BILLION * (end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
  //  printf("elapsed time = %llu nanoseconds \n", (long long unsigned int) diff);
  printf("elapsed time = %lf seconds \n", (double) diff/1000000000.0);

  return 0;
}
 
  


