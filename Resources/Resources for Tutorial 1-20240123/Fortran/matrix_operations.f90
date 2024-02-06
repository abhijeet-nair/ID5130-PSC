! program to compute addition and multiplication of 
! two matrices 

subroutine addTwoMatrices(m, n, mat1, mat2, matSum)
  implicit none
  integer, intent(in) :: m, n
  double precision, dimension(1:m, 1:n), intent(in) :: mat1, mat2
  double precision, dimension(1:m, 1:n), intent(out) :: matSum
  integer :: i, j
  double precision :: t1, t2 

  ! time at the beginning
  call cpu_time(t1)

  ! add the two matrices
  do j = 1, n
     do i = 1, m
        matSum(i, j) = mat1(i, j) + mat2(i, j)
     end do
  end do

  ! time at the end 
  call cpu_time(t2)

  write(*, *), 'Addition of matrices took ', t2-t1, ' seconds'

end subroutine addTwoMatrices

subroutine multiplyTwoMatrices(m, n, mat1, mat2, matProduct)
  implicit none
  integer, intent(in) :: m, n
  double precision, dimension(1:m, 1:n), intent(in) :: mat1, mat2
  double precision, dimension(1:m, 1:n), intent(out) :: matProduct
  integer :: i, j, k
  double precision :: t1, t2 

  ! check the dimensions
  if (m /= n) then 
     write(*, *), 'Column size of mat1 is not equal to Row size of mat2'
     write(*, *), 'Cannot multiply the matrices, please choose m = n'
     stop
  end if

  ! time at the beginning
  call cpu_time(t1)

  ! initialize elements of matProduct
  do j = 1, n
     do i = 1, m 
        matProduct(i, j) = 0.0
     end do
  end do

  ! multiply the two matrices
  do j = 1, n
     do i = 1, m
        do k = 1, m
           matProduct(i, j) = matProduct(i, j) + mat1(i, k) * mat2(k, j)
        end do
     end do
  end do
  
  ! time at the end 
  call cpu_time(t2)

  write(*, *), 'Multiplication of matrices took ', t2-t1, ' seconds'

end subroutine multiplyTwoMatrices

subroutine printMatrix(m, n, mat)
  implicit none
  integer, intent(in) :: m, n 
  double precision, dimension(1:m, 1:n), intent(in) :: mat
  integer :: i, j 
  
  write(*, *) 'Printing matrix:'
  
  do i = 1, m
     write(*, '(100f10.4)'), (mat(i, j), j = 1, n)
  end do
  
end subroutine printMatrix

program main 
  double precision :: PI = 4.0*atan(1.0)
  double precision, dimension(:, :), allocatable :: mat1, mat2, matSum, matProduct
  integer :: m, n, i, j, ii, jj

  ! read the dimension of the matrix from the user
  write(*, *), 'Please enter the number of rows and number of columns:'
  read(*, *), m, n

  ! check 
  if (m <= 0 .or. n <= 0) then 
     write(*, *), 'Matrix dimension can only be positive'
     stop
  end if

  ! allocate memory
  allocate (mat1(1:m, 1:n), mat2(1:m, 1:n))
  allocate (matSum(1:m, 1:n), matProduct(1:m, 1:n))

  ! populate the matrices
  do j = 1, n
     do i = 1, m 

        ii = i-1
        jj = j-1

        mat1(i, j) = 0.5**(0.5*ii) * sin(ii*jj*PI/(m+1))
        mat2(j, i) = 0.5**(0.5*ii) * cos(ii*jj*PI/(m+1))

     end do
  end do

  call printMatrix(m, n, mat1)
  call printMatrix(m, n, mat2)

  ! add the two matrices together
  call addTwoMatrices(m, n, mat1, mat2, matSum)
  call printMatrix(m, n, matSum)

  ! compute product of the two matrices
  call multiplyTwoMatrices(m, n, mat1, mat2, matProduct)
  call printMatrix(m, n, matProduct)

  ! dellocate memory
  deallocate (mat1, mat2, matSum, matProduct)

end program main
