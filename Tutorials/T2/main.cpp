#include <iostream>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358

double func (double x) {
    return sin(x)/(2.0*pow(x, 3.0));
}

void trapz_serial (int n, double a, double b, double* result) {
    double h, x;
    int i;
    h = (b - a)/n;

    *result += (func(a) + func(b)) / 2.0;
    for (i = 1; i <= n-1; i++) {
        x = a + i*h;
        *result += func(x);
    }
    *result = *result * h;
}

void trapz_para_crit(int n, double a, double b, double *result)
{
    double h, x, total; 		/* private  */
    int i;
    int my_rank = omp_get_thread_num(); /* private  */
    int thread_count = omp_get_num_threads();
    int local_n;			
    double local_a, local_b;	/* private scope  */

    h = (b-a)/n;			

    local_n = n/thread_count;	
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;

    total = (func(local_a) + func(local_b))/2.0;
    for(i = 1; i <= local_n-1; i++) 
        {
            x = local_a + i*h;
            total += func(x);
        }
    total = total * h;

    #pragma omp critical
        *result += total; 		/* race condition avoided... */

}

void trapz_para_for (int n, double a, double b, double *result, int thrd_cnt) {
    double h, x, total; 		/* private  */
    int i;

    h = (b - a)/n;

    total = (func(a) + func(b)) / 2.0;
    #pragma omp parallel for num_threads(thrd_cnt) reduction (+:total)
        for (i = 1; i <= n-1; i++) {
            x = a + i*h;
            total += func(x);
        }
    
    total *= h;

    *result = total;
}

void simps_para_for (int n, double a, double b, double* result, int thrd_cnt) {
    double h, x, total; 		/* private  */
    int i;
    h = (b - a)/n;
    total = func(a) + func(b);

    #pragma omp parallel for num_threads(thrd_cnt) reduction(+:total)
        for (i = 1; i <= n-1; i++) {
            x = a + i*h;
            if (i%2 == 0) {
                total += 2*func(x);
            }
            else {
                total += 4*func(x);
            }
        }
    
    total *= h/3;
    *result = total;
}

int main (int argc, char* argv[]) {
    // double ans = func(1.0);
    // std::cout << ans << std::endl;
    double a, b;
    int n;
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
    // printf("\n");

    // n = 32;			/* number of trapezoids.. */
    // a = 1;			/* shared  */
    std::cout << "Enter start point: ";
    std::cin >> a;

    std::cout << "Enter no. of trapezoids: ";
    std::cin >> n;
    b = PI;
    double res_serial =  0.0;
    double res_para_crit = 0.0;		/* shared  */
    double res_para_for = 0.0;
    double res_simps_for = 0.0;

    printf("Question 1,2:\n------------------------------------");
    // Serial
    trapz_serial(n, a, b, &res_serial);

    // Parallel Critical
    #pragma omp parallel num_threads(thrd_cnt)
    {
        trapz_para_crit(n, a, b, &res_para_crit);
    }

    // Parallel for
    trapz_para_for(n, a, b, &res_para_for, thrd_cnt);

    // Parallel for Simpsons
    simps_para_for(n, a, b, &res_simps_for, thrd_cnt);


    double act_ans = 0.198557;
    printf("\n                    Answer\tDifference");
    printf("\nActual           : %lf \t %lf", act_ans, act_ans-act_ans);
    printf("\nSerial           : %lf \t %lf", res_serial, res_serial-act_ans);
    printf("\nParallel Critical: %lf \t %lf", res_para_crit, res_para_crit-act_ans);
    printf("\nParallel For Loop: %lf \t %lf", res_para_for, res_para_for-act_ans);
    printf("\nSimpson's parfor : %lf \t %lf \n\n", res_simps_for, res_simps_for-act_ans);

    return 0;
}