#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void printMatrix (double** a, int m, int n);
void printVector (double* a, int m);
double** addTwoMatrices (double** a, double** b, int m, int n);
double** multiplyMatrices (double** a, double** b, int m, int n, int p);
double** transposeMatrix (double** a, int m, int n);
int compareTwoMatrices (double** a, double** b, int m, int n);
double* multiplyMatVec (double** a, double* b, int m, int n);

#endif