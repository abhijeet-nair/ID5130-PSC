#include<iostream>
#include<math.h>
#include<mpi.h>

double func(double x)
{
  return (0.5*sin(x)/pow(x,3));
}

double trapezoidal_rule(double la, double lb, double ln, double h)
{
  double total;
  double x;
  int i;

  total = (func(la) + func(lb))/2.0;
  for(i = 1; i <= ln-1; i++) /* sharing the work, use only local_n */
    {
      x = la + i*h;
      total += func(x);
    }
  total = total * h;

  return total;			/* total for each thread, private */
}

int main(int argc, char* argv[]) {
    double a, b, fSum, la, lb, lsum, h;
    int myid, np, i;
    int n, ln;
    double dat[3] {};

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);	  /* myrank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &np); /* size of the communicator */

    if (myid == 0) {
        printf("Enter start location a   : ");
        std::cin >> a;

        printf("Enter end location b     : ");
        std::cin >> b;

        printf("Enter no. of trapezoids n: ");
        std::cin >> n;

        ln = n/np;
        h = (b - a)/n;

        for (i = 1; i < np; i++){
            dat[0] = a + (i-1)*ln*h;
            dat[1] = a + i*ln*h;
            dat[2] = ln;

            MPI_Send(&dat, 3, MPI_DOUBLE, i, i*10, MPI_COMM_WORLD);
        }
        
        la = a + (np - 1)*ln*h;
        lb = b;
        ln = n - ln*(np - 1);
    }
    else {
        MPI_Recv(&dat, 3, MPI_DOUBLE, 0, myid*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        la = dat[0];
        lb = dat[1];
        ln = int(dat[2]);
    }

    // for (i = 0; i < np; i++) {
    //     if (myid == i) {
    //         printf("ID:%d\nla:%.4f\nlb:%.4f\nln:%d\n\n",i,la,lb,ln);
    //     }
    // }
    lsum = trapezoidal_rule(la, lb, ln, h);

    if (myid == 0) {
        fSum += lsum;
        
        for (i = 1; i < np; i++) {
            MPI_Recv(&lsum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            fSum += lsum;
        }

        printf("Result = %.4f\n",fSum);
    }
    else {
        MPI_Send(&lsum, 1, MPI_DOUBLE, 0, myid*11, MPI_COMM_WORLD);
    }
    
    // printf("ID:%d\na:%.0f\nb:%.2f\nn:%d\n",myid,a,b,n);

    MPI_Finalize();
    return 0;
}