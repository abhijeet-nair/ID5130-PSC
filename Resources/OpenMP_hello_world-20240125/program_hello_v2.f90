! a little more sophisticated version of hello world
module data
#ifdef _OPENMP
  use OMP_LIB
#endif
  implicit none

contains

  subroutine hello()
    implicit none
    integer :: my_rank, thread_count

#ifdef _OPENMP
    my_rank = OMP_GET_THREAD_NUM()
    thread_count = OMP_GET_NUM_THREADS()
#else
    my_rank = 0
    thread_count = 1 
#endif

    write(*, *) 'Hello from thread', my_rank, 'out of total threads of', thread_count

    !$OMP PARALLEL NUM_THREADS(4)
    write(*, *) 'Hi'
    !$OMP END PARALLEL

  end subroutine hello
  
end module data

program main
  use data
  implicit none
  integer :: thread_count = 4 
  character(100) :: numchar
  character(100) :: name

  if (COMMAND_ARGUMENT_COUNT() .ne. 1) then
     write(*, *) 'Error, one command line arguments is required, stopping the program'
     STOP
  end if

  

  call GET_COMMAND_ARGUMENT(0,name)
  write(*, *) name
  
  call GET_COMMAND_ARGUMENT(1, numchar)
  READ(numchar, *) thread_count

  !  call OMP_SET_MAX_ACTIVE_LEVELS(3)
  call OMP_SET_NESTED(.false.)
  
  !$OMP PARALLEL NUM_THREADS(thread_count)
  call hello()
  !$OMP END PARALLEL
    

end program main
