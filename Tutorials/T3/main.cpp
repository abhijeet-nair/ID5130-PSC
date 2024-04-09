#include <iostream>
#include <math.h>
#include <time.h>
#include <omp.h>

void swap (int a[], int i, int j) {
    int tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
}

void printArray(int a[], int n) {
    int i = 0;
    std::cout << "[";
    for (i = 0; i < n; i++) {
        std::cout << a[i];
        if (i < n-1) {std::cout << ", ";}
    }
    std::cout << "]" << std::endl;
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
    
    printf("Question 1:\n------------------------------------------\n");
    // Subtraction
    printf("Subtraction:\n");
    double arr[10] = {1, -2.5, 3.4, -5, -6.9, 9.8, -0.45, 1, 3.14, -0.981};
    double total = 0.0;
    int i;

    double act = 0.0;
    for (i = 0; i < 10; i++) {
        act += arr[i];
    }
    act = -act;

    #pragma omp parallel for num_threads(thrd_cnt) reduction(-:total)
        for (i = 0; i < 10; i++) {
            total -= arr[i];
        }

    printf("Result with reduction operator: %.4f\n", total);
    printf("Correct answer:                 %.4f\n", act);

    // Multiplication


    // Division - Doesn't work
    // printf("\nDivision:\n");
    // total = 1.0;
    // act = 1.0;

    // for (i = 0; i < 10; i++) {
    //     act /= arr[i];
    // }

    // #pragma omp parallel for num_threads(thrd_cnt) reduction(/:total)
    //     for (i = 0; i < 10; i++) {
    //         total /= arr[i];
    //     }

    // printf("Result with reduction operator: %.4f\n", total);
    // printf("Correct answer:                 %.4f\n", act);


    printf("\n\nQuestion 2:\n------------------------------------------\n");
    int n = 0;
    std::cout << "Enter the size of the array: ";
    std::cin >> n;
    int arr1[n] {};

    printf("Array before sorting:\n");
    // std::cout << "[";
    for (i = 0; i < n; i++) {
        arr1[i] = rand() % 100;
        // std::cout << "i = " << i << ", arr[i] = " << arr1[i] << std::endl;
        // std::cout << arr1[i];
        // if (i < n-1) {std::cout << ", ";}
    }
    // std::cout << "]" << std::endl;
    // printArray(arr1, n);
   
    // Odd-Even Transposition Sort
    int pass = 0;

    clock_t t;
    t = clock();
    // time_t start, end;
    // time(&start);

    #pragma omp parallel num_threads(thrd_cnt) default(none) shared(arr1, n) private(i, pass)
    {
        for (pass = 0; pass < n; pass++) {
            if (pass%2 == 0) {
                #pragma omp for
                    for (i = 1; i < n; i += 2) {
                        if (arr1[i-1] > arr1[i]) {swap(arr1, i-1, i);}
                    }
            }
            else {
                #pragma omp for
                    for (i = 1; i < n-1; i += 2) {
                        if (arr1[i] > arr1[i+1]) {swap(arr1, i, i+1);}
                    }
            }
        }
    }
    t = clock() - t;
    // time(&end); 
 
    // double time_taken = double(end - start);
    double time_taken = double(t) / double(CLOCKS_PER_SEC);

    printf("\nArray after sorting:\n");
    // std::cout << "[";
    // for (i = 0; i < n; i++) {
    //     std::cout << arr1[i];
    //     if (i < n-1) {std::cout << ", ";}
    // }
    // std::cout << "]" << std::endl;
    // printArray(arr1, n);

    // printf("\nTime taken for sorting: %f seconds\n",((float)t)/CLOCKS_PER_SEC);
    printf("\nTime taken for sorting: %.6f seconds\n",time_taken);
    // printf("Done\n");

    /*
    printf("\n\nQuestion 3:\n------------------------------------------\n");
    int a = 10;

    // private variable
    printf("Example of private variable:\n");
    printf("Value before parallel block, a = %d. Address = ", a);
    std::cout << &a << "\n\n";

    #pragma omp parallel num_threads(thrd_cnt) private(a)
    {
        printf("Value from thread %d, a = %d\n", omp_get_thread_num(), a);
        a += 10;
        printf("Value from thread %d, a = %d\n", omp_get_thread_num(), a);
    }

    printf("\nValue after parallel block, a = %d. Address = ", a);
    std::cout << &a << "\n\n";
    
    // firstprivate variable
    printf("\nExample of firstprivate variable:\n");
    a = 10;

    printf("Value before parallel block, a = %d. Address = ", a);
    std::cout << &a << "\n\n";

    #pragma omp parallel num_threads(thrd_cnt) firstprivate(a)
    {
        printf("Value from thread %d, a = %d\n", omp_get_thread_num(), a);
        a += 10;
        printf("Value from thread %d, a = %d\n", omp_get_thread_num(), a);
    }

    printf("\nValue after parallel block, a = %d. Address = ", a);
    std::cout << &a << "\n\n";

    // threadprivate variable
    printf("\nExample of threadprivate variable:\n");
    static int b = 10;

    #pragma omp threadprivate(b)
    #pragma omp parallel num_threads(thrd_cnt) copyin(b)
    {
        printf("Value from thread %d, b = %d\n", omp_get_thread_num(), b);
        b += 10;
        printf("Value from thread %d, b = %d\n", omp_get_thread_num(), b);
    }

    printf("\nFirst parallel block done...\n\n");

    // #pragma omp threadprivate(b)
    #pragma omp parallel num_threads(thrd_cnt)
    {
        printf("Value from thread %d, b = %d\n", omp_get_thread_num(), b);
        b += 10;
        printf("Value from thread %d, b = %d\n", omp_get_thread_num(), b);
    }

    printf("\nSecond parallel block done...\n");
    */

    return 0;
}