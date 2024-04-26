/* serial code for Cholesky decomposition */
/* make sure that the init function setups a  */
/* symmetric and positive definite matrix  */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE		float
#define N		100
#define SMALLVALUE	0.001

void initmult(TYPE mat[][N])
{
  for (int ii = 0; ii < N; ++ii)
    for (int jj = 0; jj < N && jj < ii; ++jj)
      {	mat[ii][jj] = (ii + jj) / (float)N / N;
	mat[jj][ii] = (ii + jj) / (float)N / N;}

  for (int ii = 0; ii < N; ++ii)
    mat[ii][ii] = 1.0;
}
			
			
void printMat(TYPE a[][N])
{
  for (int ii = 0; ii < N; ++ii)
    {
      for (int jj = 0; jj < N; ++jj)
	printf("%.2f ", a[ii][jj]);
      printf("\n");
    }
}

void cholesky(TYPE a[][N])
{
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < ii; ++jj) {
      for (int kk = 0; kk < jj; ++kk)
	a[ii][jj] += a[ii][kk] * a[jj][kk];
      a[ii][jj] /= (a[jj][jj] > SMALLVALUE ? a[jj][jj] : 1);
      //a[ii][jj] /= a[jj][jj];	// divide by zero.
    }
    for (int kk = 0; kk < ii; ++kk)
      a[ii][ii] += -a[ii][kk] * a[ii][kk];
    a[ii][ii] = sqrt(a[ii][ii]);
  }
}

int main()
{
  TYPE a[N][N];

  init(a);
  cholesky(a);
  printMat(a);

  return 0;
}
