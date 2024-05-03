/* OpenMP program to highlight for loop exceptions */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#define N 50

int main(int argc, char* argv[])
{
  int a[N];
  int key, i, thread_count = 1;

  if (argc == 2)
    {
      thread_count = strtol(argv[1], NULL, 10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..");
      return 1;
    }
  
  /* fill some data */
# pragma omp parallel for num_threads(thread_count)
  for(i = 0 ; i < N; i++ )	
    {
      a[i] = i+i^2;
    }
  
  key = 5 + 5^2;

  /* search for a key */
#pragma omp parallel for num_threads(thread_count)
  for(i = 0; i < N; i++)
    if (a[i] == key) exit;	/* not a structured block of code.. */
  
  printf("\n found the key at location, i = %d", i);

}
