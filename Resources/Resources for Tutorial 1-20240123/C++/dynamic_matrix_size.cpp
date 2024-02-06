#include <iostream>
#include <math.h>       /* pow */
#include <time.h>       /* clock_t */
using namespace std;

void addTwoMatrices (double** mat1, double** mat2, double** sum, int m, int n) {
	// Add two matrices
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			sum[i][j] = mat1[i][j]+mat2[i][j];
 		}
	}
	return;
}

void multiplyTwoMatrices (double** a, double** b, double** c, int m, int n, int p) {
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

void printMatrix (double** mat, int m, int n) {
	
	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			printf("%f\t",mat[i][j]);
		}
		cout << "\n";
	}
	return;
}

double matrix_operations (int m, int n) {
	
	// Print out general info
	printf("Matrix size %dx%d\n",m,n);
	// Declare matrices mat1, mat2
	double** mat1 = new double*[m];
	double** mat2 = new double*[m];
	double** sum = new double*[m];
	double** prod = new double*[m];

	for (int i = 0; i<m; i++) {
		mat1[i] = new double[n];
		mat2[i] = new double[n];
		sum[i] = new double[n];
		prod[i] = new double[n];
	}

	for (int i = 0; i < m; i ++) {
		for (int j = 0; j < n; j++) {
			mat1[i][j] = 1.0; //pow(0.5,0.5*i) * sin(i*j*PI/(m+1));
			mat2[j][i] = 2.0; //pow(0.5,0.5*i) * cos(i*j*PI/(m+1));
		}
	}
	
	clock_t t;
	t = clock();
	addTwoMatrices(mat1,mat2,sum,m,n);
	t = clock()-t;
	printf ("Addition took %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	// cout << "Sum of the matrices: \n\n";
	// printMatrix(sum,m,n);
	
	int p = m; // We are multiplying square matrices

	t = clock();
	multiplyTwoMatrices(mat1,mat2,prod,m,n,p);
	t = clock()-t;
	printf ("Multiplication took %d clicks (%f seconds).\n\n",t,((float)t)/CLOCKS_PER_SEC);
	double time_taken = ((float)t)/CLOCKS_PER_SEC;
	
	// cout << "Product of the matrices: \n\n";
	// printMatrix(prod,m,n);
	
	for (int i = 0; i<m; i++) {
		delete mat1[i];
		delete mat2[i];
		delete sum[i];
		delete prod[i];		
	}
	
	return time_taken;
}

int main () {

	const double PI = 3.14;
    
	// Get matrix size
	int m, n;
	cout << "Enter the number of rows in square matrix ";
    cin >> m;
	n = m; // nRows = nCols
	matrix_operations(m,n);
	
	// Give a predefined m,n
	int mArray[5] = {256, 512, 1024, 2048, 4096};
	int nArray[5];
	int nTrials = 5;
	for (int i = 0; i < nTrials; i ++) {
		nArray[i] = mArray[i];
		matrix_operations(mArray[i],nArray[i]);
	}
	
		
	return 0;
}