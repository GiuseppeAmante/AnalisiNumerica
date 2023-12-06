subroutine lufact(A,L,U,g,n,response)
    implicit none
    real(kind=8), dimension(n,n) :: A
    real(kind=8), dimension(n,n), intent(out) :: L, U
    reaL(kind=8), allocatable :: copy(:,:), matrix_G(:,:), matrix_error(:,:)
    real(kind=8) :: g, val_1, val_2, mach_precision, norm_matrix_1, norm_matrix_2, nrb_error
    integer :: n, i, j, k, ierr
    real(kind=8), parameter :: zero = 0.0d0, one=1.0d0
    character(len=30) :: response 
!    
open(unit=1, file='factLU.txt', status='replace', action='write', iostat=ierr)
  if (ierr.ne.0) then
    write (*,'(A)') 'Error! on open the file factLU.txt'
    stop
  endif
!
    mach_precision = epsilon(1.0d0)
    g = zero
    L = zero
    U = zero
    allocate(copy(n,n))
    copy = A
! 
    call print_matrix(1,A,n,n,response)

    do k = 1, n
      L(k, k) = one
      do i = k + 1, n
        L(i, k) = A(i, k) / A(k, k)
        if (L(i,k) < mach_precision) then
          write(1,'(A)') 'Error! Factorizzatio LU terminated'
          write(1,'(A)') 'due to division by a quantity smaller than machine precision'
          stop
        endif
        do j = k, n
          A(i, j) = A(i, j) - L(i, k)*A(k, j)
        end do
      end do
      do j = k, n
        U(k, j) = A(k, j)
      end do
    end do

    response = 'Lower Triangular'
    call print_matrix(1,L,n,n,response)
    response = 'Upper Triangular'
    call print_matrix(1,U,n,n,response)
!  G = |L| x |U|
    allocate (matrix_G(n,n), matrix_error(n,n))
    matrix_G = zero
    matrix_error = zero
    do i = 1, n
      do j = 1, n
        matrix_G(i,j) = zero
        do k =1, n
          matrix_G(i,j) = matrix_G(i,j) + L(i,k)*U(k,j)
        enddo
        matrix_error(i,j) = copy(i,j) - matrix_G(i,j)
      enddo
    enddo
!  growth factor g = max( |L| x |U| ) / max( |A| )
    val_1 = maxval(abs(matrix_G))
    val_2 = maxval(abs(copy))
    g = val_1/val_2
    if ( g < mach_precision ) then
       write(1,'(A)') 'Error: Division by a quantity smaller than machine precision '
       stop
    else
       write(1,*) ' Growth factor is :', g
    endif
! normwise relative backward error
    call norm_inf_matrix(matrix_error,n,n,norm_matrix_1)
    write(1,*) 'Norm inf A - LU :', norm_matrix_1
    call norm_inf_matrix(copy,n,n,norm_matrix_2) 
    write(1,*) 'Norm inf A :', norm_matrix_2
    nrb_error = norm_matrix_1/norm_matrix_2 
    write(1,*) 'Norm inf A - LU/ Norm inf A :', nrb_error
    write(1,*) 'Machine precision :', mach_precision
    if ( nrb_error < mach_precision ) then
       write(1,'(A)') 'Error: Division by a quantity smaller than machine precision '
       stop
    else
       write(1,*) ' Normwise relative backward error is :',  nrb_error
    endif
    deallocate(copy,matrix_G,matrix_error)
close(1)
    end subroutine
