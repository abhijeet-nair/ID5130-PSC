program main
  implicit none
  include "mpif.h"
  integer :: i, myid, size, mpierror, tag, status(MPI_STATUS_SIZE)
  character(len=50) :: message_send, message_recv

  call MPI_INIT(mpierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, mpierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

  tag = 100

  if (myid /= 0) then 
     write(message_send, *), 'Hello from process id:', myid
     call MPI_SEND(message_send, 50, MPI_CHARACTER, 0, tag, MPI_COMM_WORLD, mpierror)
  else
     do i = 1, size-1
        call MPI_RECV(message_recv, 50, MPI_CHARACTER, i, tag, MPI_COMM_WORLD, status, mpierror)
        write(*, *), message_recv
     end do
     write(message_send, *), 'Hello from process id:', myid     
     write(*, *), message_send
  end if

  call MPI_FINALIZE(mpierror)

end program main
