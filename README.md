# ID5130-PSC

Codes written for the institute course ID5130: Parallel Scientific Computingy. A summary given below of the contents covered through coding practice.

## Tutorial 1
Basic C++ practice. Adapted functions provided by the professor for matrix operations. Generated a plot with Python.

## Tutorial 2
Numerical integration with *Trapezoidal* and *Simpson's Rule*. Used *OpenMP* to parallelize the calculations. Used *manual work-splitting* and *parallel for*.

## Tutorial 3
Experimentation with *reduction* clause in OpenMP with subtraction and division. Implemented *Odd-Even Transposition Sort* with parallelized inner loops. Tested the behaviour of *private*, *firstprivate* and *threadprivate* variable clauses.

## Tutorial 4
Starting Numerical Methods with OpenMP. Addition, multiplication and implicit & explicit finite differences.

## Tutorial 5
Implemented *Iterative Jacobi Method* for solving system of linear equations.

## Assignment 1
Completely covered numerical methods for OpenMP.

Wrote a *serial LU decomposition program*, which is used to find solution to a system of linear equation, generated using *implicit finite difference scheme* (Pade scheme). Implemented a parallel program for *Recursive-Doubling algorithm* for the same.

Implemented *serial Gauss-Seidel algorithm* for finding solution to the Poisson equation. Parallelized it using two methods, *diagonal approach* and *red-black coloring approach*. Compared performance with different grid sizes and number of threads.

## Tutorial 6
Installed *OpenMPI*. Run a sample Hello World program for MPI provided by the professor for testing.

## Tutorial 7
Wrote *MPI* programs to send and receive messages between the processes. Sent integers, strings and arrays.

## Tutorial 8
Numerical integration with *Trapezoidal* and *Simpson's Rule*. Used *OpenMPI* to parallelize the calculations. Made sure that the important values can be taken from the user and then passed on to other processes. Used *MPI_Reduce* function as well as an alternative to *MPI_Send*/*MPI_Recv*.

## Tutorial 9
Wrote an example program for the usage of *MPI* Collective Communication (CC) functions, like *MPI_Reduce*, *MPI_Allreduce*, *MPI_Scatter*, *MPI_Gather*, *MPI_Allgather*.

Developed an *MPI* program which performs matrix-vector multiplication, with the matrix and vector being block-decomposed along rows. Used only *MPI* CC functions. Verified the results.

## Tutorial 10
Developed an *MPI* program for computing derivatives of a function using explicit first and second-order accurate schemes. Verified the results.

Learnt about new functions such as *MPI_Gatherv*, *MPI_Scatterv*, *MPI_Dims_create*, *MPI_Cart_create*, *MPI_Cart_rank*, *MPI_Cart_coords*.