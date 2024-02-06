!My first program 
program main 
  ! tell fortran compiler to assume nothing 
  implicit none
  ! declare a character array to hold my name
  character(len=128) :: myname

  write(*, *), 'Please enter your name: '
  read (*, *), myname
  write(*, *), 'Hello ',trim(myname), '!, Welcome to ID5130 !'

end program main
