! a little more sophisticated version of hello world
module data
  use OMP_LIB
  implicit none

  double precision, parameter :: PI = 3.14159265358

contains

  double precision function func(x)
    implicit none
    double precision, intent(in) :: x

    func = (1.0 + dsin(x))
    
  end function func

  subroutine trapezoidal_rule(n, a, b, final_result)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: a, b
    double precision, intent(out) :: final_result
    double precision :: h, x, total
    integer :: local_n, i
    double precision :: local_a, local_b
    integer :: my_rank, thread_count

    my_rank = OMP_GET_THREAD_NUM()
    thread_count = OMP_GET_NUM_THREADS()

    h = (b-a)/n
    local_n = n/thread_count;
    local_a = a + my_rank*local_n*h
    local_b = local_a + local_n*h

    total = (func(local_a) + func(local_b))/2.0
    do i = 1, local_n-1
       x = local_a + i*h
       total = total + func(x)
    end do
    total = total*h
    
    !$OMP CRITICAL
    final_result = final_result + total
    !$OMP END CRITICAL

  end subroutine trapezoidal_rule
  
end module data

program main
  use data
  implicit none
  character(100) :: numchar
  character(100) :: name
  double precision :: a, b, final_result
  integer :: n
  integer :: thread_count = 1

  if (COMMAND_ARGUMENT_COUNT() .ne. 1) then
     write(*, *) 'Error, one command line arguments is required, stopping the program'
     STOP
  end if

  call GET_COMMAND_ARGUMENT(0,name)
  write(*, *) name
  
  call GET_COMMAND_ARGUMENT(1, numchar)
  READ(numchar, *) thread_count

  n = 100;
  a = 0.0;
  b = PI
  final_result = 0.0;

  !$OMP PARALLEL NUM_THREADS(thread_count)
  call trapezoidal_rule(n, a, b, final_result)
  !$OMP END PARALLEL

  write(*, *) 'The final result is ', final_result

end program main
