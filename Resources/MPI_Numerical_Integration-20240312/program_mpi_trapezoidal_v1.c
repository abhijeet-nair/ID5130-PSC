/* MPI parallel version of trapezoidal rule */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

#define PI 3.14159265358

double func(double x)
{
  return (1.0 + sin(x));
}

double trapezoidal_rule(double la, double lb, double ln, double h)
{
  double total;
  double x;
  int i;

  total = (func(la) + func(lb))/2.0;
  for(i = 1; i <= ln-1; i++) /* sharing the work, use only local_n */
    {
      x = la + i*h;
      total += func(x);
    }
  total = total * h;

  return total;			/* total for each thread, private */
}


int main(int argc, char* argv[])
{
  double a, b, final_result, la, lb, lsum, h;
  int myid, nprocs, proc;
  int n, ln;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);	  /* myrank of the process */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* size of the communicator */

  n = 1024;			/* number of trapezoids.. */
  a = 0.0;
  b = PI;			/* hard-coded.. */
  final_result = 0.0;

  h = (b-a)/n;
  ln = n/nprocs; 		/* nprocs evenly divides number of trapezoids */

  la = a + myid*ln*h;
  lb = la + ln*h;
  lsum = trapezoidal_rule(la, lb, ln, h); /* every process calls this function... */

  if (myid != 0)
    {
      MPI_Send(&lsum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  else				/* process 0 */
    {
      final_result = lsum;
      for (proc = 1; proc < nprocs; proc++)
	{
	  MPI_Recv(&lsum, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  final_result += lsum;
	}
    }

  if (myid == 0) 		/* output is only printed by process 0 */
    {
      printf("\n The area under the curve (1+sin(x)) between 0 to PI is equal to %lf \n\n", final_result);
    }
  
  MPI_Finalize();
  return 0;  
}
