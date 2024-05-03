/* matrix multiplication - serial version */
#include <stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include <time.h>

#define N 32
#define BILLION 1000000000L

void initialize(double matrix[][N]){
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      matrix[i][j] = (i+1)*N+(j+1);
  
  return;
}

void multiplymatrices(double A[][N], double B[][N], double C[][N])
{

  double temp;
  
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      temp = 0;
      for (int k = 0; k < N; k++)
	temp += A[i][k] * B[k][j];
      C[i][j] = temp;
    }

  return;
}

void printmatrix(double matrix[N][N]){

  if (N<=10)
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

  double A[N][N]={0.0}, B[N][N]={0.0}, C[N][N]={0.0};
  
  clock_gettime(CLOCK_MONOTONIC, &start_t);	/* mark start time */  
  
  initialize(A);
  initialize(B);
  
  printf("\n Matrix A is:\n");
  printmatrix(A);
  printf("\n Matrix B is:\n");
  printmatrix(B);
  
  multiplymatrices(A, B, C);
    
  printf("\n Matrix C is:\n");
  printmatrix(C);
  
  clock_gettime(CLOCK_MONOTONIC, &end_t);	/* mark the end time */
  
  printf("\n sanity check value for C[0][0] = %lf \n", C[0][0]);

  diff = BILLION * (end_t.tv_sec - start_t.tv_sec) + end_t.tv_nsec - start_t.tv_nsec;
  //  printf("elapsed time = %llu nanoseconds \n", (long long unsigned int) diff);
  printf("elapsed time = %lf seconds \n", (double) diff/1000000000.0);

  return 0;
}
 
  


