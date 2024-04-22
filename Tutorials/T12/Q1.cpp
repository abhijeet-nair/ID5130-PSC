#include <iostream>

int main (int argc, char* argv[]) {
    int ngangs = 1;
    
    if (argc == 2) {
        ngangs = atoi(argv[1]);
    }
    else {
        printf("Run the program as ./a.out 10...");
        return 1;
    }

    printf("\nngangs = %d\n", ngangs);

    #pragma acc parallel num_gangs(ngangs)
    {
        printf("Hello world\n");
        printf("Bye world\n");
    }

    printf("Host\n");

    #pragma acc parallel num_gangs(ngangs/2)
        printf("Second pragma\n");

    return 0;
}