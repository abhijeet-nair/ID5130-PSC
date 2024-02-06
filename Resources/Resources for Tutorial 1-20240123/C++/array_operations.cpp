#include <iostream>
#include <string>
#include <math.h> /* pow */
using namespace std;

// Function declarations
// Compute dot product of two arrays, whose length is given
double dotProduct (double arr1[], double arr2[], int n) {
	
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += arr1[i]*arr2[i];
	}
	
	return sum;	
}


int main () {

	// Create and retrieve elements of array
	// Create two arrays of size 5
	int n = 5;
	double arr1[5] = {1000.0, 2.0, 3.4, 17.0, 50.0};
	double arr2[5] = {100.0, 0.20, 0.34, 1.7, 5.0};
	string names[5] = {"first","second","third","fourth","fifth"};
	
	for (int i = 0; i< n; i++) {
		cout << "The " << names[i] << " element is: " << arr1[i] <<"\n";
	}
	
	// Dynamic array - read size from user
	int m;
	cout << "Enter the number of elements in array ";
    cin >> m;
	double* arr3 = new double[m];
	
	// Populate with entries
	cout << "Dynamic array entries are: \n";
	for (int i=0; i<m; i++) {
		arr3[i] = pow(i,m);
		cout << arr3[i] << "\t";
	}
	cout << "\n";
	delete arr3;
	
	// Compute sum of arrays
	double s;
	cout << "Sum of arrays: ";
	for (int i = 0; i < n; i++) {
		s = arr1[i]+arr2[i];
		cout << s << "\t";
	}
	cout << "\n";
		
	// Compute dot product of arrays
	cout << "Dot product: " <<  dotProduct(arr1,arr2,n) << "\n";
	return 0;
}
