#include<iostream>

// THIS CODE IS FOR TESTING THE ROW-MAJOR ORDER OF C++.
// RUN THE FILE AND YOU WILL SEE HOW C++ STORES 2D ARRAYS.
int main() {
    printf("\nTHIS CODE IS FOR TESTING THE ROW-MAJOR ORDER OF C++...\n\n\n");

    int A[5][5] {}, i {}, j {};
    int* pntr = &A[0][0];
    double val;

    int cnt = 1;
    printf("   (i, j)\tdat\t    Address\t\tcnt\tdat\n");
    printf("  ----------------------------------------------------------\n");
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            A[i][j] = cnt;
            val = *(pntr + cnt - 1);
            printf("   (%d, %d)\t%2.0f\t", i, j, double(A[i][j]));
            std::cout << &A[i][j] << "\t\t";
            printf("%2.0f\t%2.0f\n", double(cnt), val);
            cnt += 1;
        }
    }
    printf("\t*Note that first dat is through direct access...\n");
    printf("\tSecond dat is through pointer... Just to verify...\U0001F600\n");
    printf("\nMatrix in usual order:\n");
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            printf("\t%2.0f", double(A[i][j]));
        }
        printf("\n");
    }
    printf("\nMatrix in Sir's order:\n");
    for (i = 4; i >= 0; i--) {
        for (j = 0; j < 5; j++) {
            printf("\t%2.0f", double(A[i][j]));
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}