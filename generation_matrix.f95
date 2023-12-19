subroutine generation_matrix(matrix_a,n,response)
    implicit none
    real(kind=8), dimension(n,n) :: matrix_a
    real(kind=8), allocatable :: x(:)
    integer :: n, i ,j
    real(kind=8), parameter :: zero=0.0d0, one = 1.0d0
    character(len=10), intent(out) :: response
!
    write(*,*) 'What type of square matrix do you want to generate?'
    write(*,*) '[hilbert, lehmer, vandermode, toeplitz]'
    read(*,*) response
    if (trim(response) == 'hilbert') then
        ! Hilbert matrix h_ij = 1/i+j-1
        do i = 1, n
          do j = 1, n
            matrix_a(i,j) = one / real(i+j-1)
          enddo
        enddo
    else if (trim(response) == 'lehmer') then
        do i = 1, n
          do j = 1, n
              matrix_a(i,j) = real(min(i,j)) / real(max(i,j))
          enddo
        enddo
    else if (trim(response) == 'vandermode') then
        allocate(x(n))
        do i = 1, n
          x(i) = real(i)
        enddo
        matrix_a = one
        do i = 1, n
          do j = 2, n
            matrix_a(i, j) = x(i)**(j-1)
          end do
        end do
        deallocate(x)
    else if (trim(response) == 'toeplitz') then
        allocate(x(n))
        do i = 1, n
          x(i) = real(i)
        enddo
        do i = 1, n
          do j = 1, n
            matrix_a(i, j) = x(abs(i-j) + 1)
          end do
        end do
        deallocate(x)
    else 
        write(*,*) 'No matrix generated'
        stop  
    endif
    end subroutine
