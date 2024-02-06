#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
#include <cstdio>
using namespace std;

void addTwoMatrices (double mat1[][100], double mat2[][100], double sum[][100], int m, int n) {
	// Add two matrices
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			sum[i][j] = mat1[i][j]+mat2[i][j];
 		}
	}

	return;
}

void multiplyTwoMatrices (double a[][100], double b[][100], double c[][100], int m, int n, int p) {
	for (int i = 1; i < m; i ++) {
		for (int j = 0; j < n; j ++) {
			c[i][j] = 0;
			for (int k = 0; k < p; k ++) {
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return;
}

void printMatrix (double mat[][100], int m, int n) {
	
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			printf("%f\t",mat[i][j]);
		}
		cout << "\n";
	}
	return;
}

int main () {

	const double PI = 3.14;
	double mat1[100][100], mat2[100][100];
	double sum[100][100], prod[100][100];
    
	// Get matrix size
	int m, n;
	cout << "Enter the number of rows (<100) and cols in matrix ";
    cin >> m >> n;

	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			mat1[i][j] = pow(0.5,0.5*i) * sin(i*j*PI/(m+1));
			mat2[j][i] = pow(0.5,0.5*i) * cos(i*j*PI/(m+1));
		}
	}
	
	clock_t t; 
	t = clock();
	addTwoMatrices(mat1,mat2,sum,m,n);
	t = clock()-t;
    printf ("Addition took %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	cout << "Sum of the matrices: \n\n";
	// printMatrix(sum,m,n);
	
	int p = m; // We are multiplying square matrices
	
	t = clock();
	multiplyTwoMatrices(mat1,mat2,prod,m,n,p);
	t = clock()-t;
	printf ("Multiplication took %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	cout << "Product of the matrices: \n\n";
	// printMatrix(prod,m,n);
	
	return 0;
}
