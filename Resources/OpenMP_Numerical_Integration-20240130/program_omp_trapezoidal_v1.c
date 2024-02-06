/* OpenMP parallel version of trapezoidal rule */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#define PI 3.14159265358

double func(double x)
{
  return (1.0 + sin(x));
}

void trapezoidal_rule(int n, double a, double b, double *result)
{
  double h, x, total; 		/* private  */
  int i;
  int my_rank = omp_get_thread_num(); /* private  */
  int thread_count = omp_get_num_threads();
  int local_n;			
  double local_a, local_b;	/* private scope  */

  h = (b-a)/n;			

  local_n = n/thread_count;	
  local_a = a + my_rank*local_n*h;
  local_b = local_a + local_n*h;
  
  total = (func(local_a) + func(local_b))/2.0;
  for(i = 1; i <= local_n-1; i++) 
    {
      x = local_a + i*h;
      total += func(x);
    }
  total = total * h;

#pragma omp critical
  *result += total; 		/* race condition avoided... */

}


int main(int argc, char* argv[])
{
  double a, b, final_result;
  int n;
  int thread_count = 1;

  if (argc == 2)
    {
      thread_count = strtol(argv[1], NULL, 10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..");
      return 1;
    }

  n = 124;			/* number of trapezoids.. */
  a = 0.0;			/* shared  */
  b = PI;
  final_result = 0.0;		/* shared  */

#pragma omp parallel num_threads(thread_count)
  trapezoidal_rule(n, a, b, &final_result);

  printf("\n The area under the curve (1+sin(x)) between 0 to PI is equal to %lf \n\n", final_result);

  return 0;  
}
