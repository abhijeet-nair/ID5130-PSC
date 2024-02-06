/* OpenMP parallel version of Fibonacci series -- ? */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#define N 10

int main(int argc, char* argv[])
{
  int i, n, fib[N];
  int thread_count; 

  if (argc == 2)
    {
      thread_count = strtol(argv[1], NULL, 10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..");
      return 1;
    }

  fib[0] = 1;
  fib[1] = 1;

#pragma omp parallel for num_threads(thread_count)
  for (i = 2; i < N; i++)
    {
      fib[i] = fib[i-1] + fib[i-2];
    }

  for (i = 0; i < N; i++)
    printf("\t%d", fib[i]);

  return 0;  
}
