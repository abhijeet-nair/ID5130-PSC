#include <iostream>
#ifdef _OPENMP
    #include <omp.h>
#endif

void hello(void) {
    #ifdef _OPENMP
        int my_rank = omp_get_thread_num();
        int thrd_cnt = omp_get_num_threads();
    #else
        int my_rank = 0;
        int thrd_cnt = 1;
    #endif

    // std::cout << "Hello from thread " << my_rank << " out of total threads of " << thrd_cnt << std::endl;
    printf("Hello from thread %d out of total threads of %d\n", my_rank, thrd_cnt);
}

int main (int argc, char* argv[]) {
    int thrd_cnt = 1;

    if (argc == 2) {
        thrd_cnt = strtol(argv[1], NULL, 10);
    }
    else {
        // std::cout << "A command line argument other than name of the executable is required... exiting the program...";
        // std::cout << std::endl;
        printf("\n A command line argument other than name of the executable is required... Exiting the program...\n");
        return 1;
    }
    printf("\n");
    #pragma omp parallel num_threads(thrd_cnt)
        hello();

    return 0;
}