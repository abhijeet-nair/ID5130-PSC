!program to perform some basic array operations

! subroutine to compute dot product of 
! two vectors with given length
subroutine dotProduct(n, a, b, aDotb)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(1:n), intent(in) :: a, b
  double precision, intent(out) :: aDotb
  integer :: i 
  
  ! initialize the dot product value 
  aDotb = 0.0

  ! compute the dot product
  do i = 1, n
     aDotb = aDotb + a(i)*b(i)
  end do

end subroutine dotProduct

program main
  implicit none
  ! declare and define two vectors of length 5 
  integer, parameter :: n = 5 
  double precision, dimension(1:n) :: arr1 = (/ 1000.0, 2.0, 3.4, 17.0,  50.0 /)
  double precision, dimension(1:n)  :: arr2 = (/ 100.0, 0.20, 0.34, 1.7, 5.0 /)
  character(len=128), dimension(1:n) :: names = (/ 'first ', 'second', 'third ', 'fourth', 'fifth ' /)
  integer :: i 
  double precision :: product

  ! retrieve and print the elements of the array 
  write(*, *), 'Elements of the first array are: '
  do i = 1, n
     write(*, *), 'The ', trim(names(i)), ' element is: ', arr1(i)
  end do

  write(*, *), 'Elements of the second array are: '
  do i = 1, n
     write(*, *), 'The ', trim(names(i)), 'elements is: ', arr2(i)
  end do

  ! compute sum of the two arrays
  write(*, *), 'Sum of arrays: '
  do i = 1, n 
     write(*, *), arr1(i)+arr2(i)
  end do

  ! compute dot product of arrays by calling 
  ! the respective subroutine 
  call dotProduct(n, arr1, arr2, product)
  
  write(*, *), 'Dot product is: ', product

end program main
