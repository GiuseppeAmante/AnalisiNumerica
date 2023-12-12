subroutine generation_matrix(matrix_a,n,response)
    implicit none
    real(kind=8), dimension(n,n) :: matrix_a
    real(kind=8), allocatable :: x(:)
    integer :: n, i ,j
    real(kind=8), parameter :: zero=0.0d0, one = 1.0d0
    real(kind=8) :: conv, val, row_sum
    character(len=30), intent(out) :: response
!
    write(*,*) 'What type of square matrix do you want to generate?'
    write(*,*) '[hilbert, diagonally dominant, pascal, vandermode]'
    read(*,*) response
    if (trim(response) == 'hilbert') then
        ! Hilbert matrix h_ij = 1/i+j-1
        do i = 1, n
          do j = 1, n
            conv = real(i+j-1)
            matrix_a(i,j) = one / conv
          enddo
        enddo
    else if (trim(response) == 'diagonallydominant') then
        do i = 1, n
          matrix_a(i, i) = 2.0d0
          do j = 1, n
            if (i /= j) then
                val = matrix_a(i, i)+matrix_a(i,i)
                matrix_a(i, j) = sqrt(val)/real(i+i)
            end if
          enddo
        enddo
    ! Test if the matrix is diagonally dominant
        do i = 1, n
          row_sum = zero
          do j = 1, n
            if (i /= j) then
              row_sum = row_sum + abs(matrix_a(i, j))
            end if
          end do
          if ( row_sum > abs( matrix_a(i, i) ) ) then
            write(*, *) 'The matrix is not diagonally dominant!'
          end if
        end do
    else if (trim(response) == 'pascal') then
        do i = 1, n
          do j = 1, i
            if (j == 1 .or. j == i) then
                matrix_a(i, j) = 1
            else
                matrix_a(i, j) = matrix_a(i-1, j-1) + matrix_a(i-1, j)
            end if
          end do
        end do
    else if (trim(response) == 'vandermode') then
        allocate(x(n))
        do i = 1, n
          x(i) = real(i-1)
        enddo
        matrix_a = one
        do i = 1, n
          do j = 2, n
            matrix_a(i, j) = x(i)**(j-1)
          end do
        end do
    else 
        write(*,*) 'Error!' 
        stop  
    endif
    end subroutine
