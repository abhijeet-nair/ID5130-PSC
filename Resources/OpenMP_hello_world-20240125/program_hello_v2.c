/* A little more sophisticated version of hello world */
#include<stdio.h>
#include<stdlib.h>
#ifdef _OPENMP
#include<omp.h>
#endif

void Hello(void)
{
#ifdef _OPENMP
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
#else
  int my_rank = 0;
  int thread_count = 1;
#endif

  printf("\n Hello from thread %d out of total threads of %d", my_rank, thread_count);
}

int main(int argc, char* argv[])
{
  int thread_count = 1;

  if (argc == 2)
    {
      thread_count = strtol(argv[1], NULL, 10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }
  
#pragma omp parallel num_threads(thread_count)
  Hello();

  return 0;
}
