subroutine wilkin(n,matrix)
    implicit none
    real(kind=8), dimension(n,n) :: matrix
    integer :: n
!    
    integer :: i , j, ierr, info
    real(kind=8) :: rmse, mach_precision, error_relative, norm_1, norm_2,&
                    backward_error, backward_error_ge
    character(len=30) :: type_matrix
    real(kind=8), allocatable :: e(:), b(:), y(:), x(:), L(:,:), U(:,:),&
        x_ge(:), work_1(:), work_2(:), copy(:,:)
    integer, allocatable :: pivots(:)
    real(kind=8), parameter :: zero = 0.0d0, one=1.0d0
!
open(unit=2, file='wilkin.txt', status='replace', action='write', iostat=ierr)
  if (ierr .ne. 0) then
    write (*,'(A)') 'Error:  open the file wilkin.txt !'
    stop
  endif
   mach_precision= epsilon(1.0d0)
!  Generate the nxn Wilkinson matrix
    matrix= zero
    do i = 1, n
      matrix(i,i) = one 
      matrix(i,n) = one
    enddo
    do i = 2, n
      do j = 1, i-1
        matrix(i,j) = -one
      enddo
    enddo
!
    type_matrix = 'Wilkinson'
    call print_matrix(2,matrix,n,n,type_matrix)
!
! Numerical experiment
!
    allocate(e(n),b(n),y(n),x(n),L(n,n),U(n,n))
    e = one
    b = zero
    y = zero
    x = zero
    L = zero
    U = zero
!A=LU
    call lufact(matrix,n,type_matrix,L,U)
    type_matrix = 'Lower Triangular'
    call print_matrix(2,L,n,n,type_matrix)
    type_matrix = 'Upper Triangular'
    call print_matrix(2,U,n,n,type_matrix)
!b = Ae
    call dgemv('n',n,n,one,matrix,n,e,1,zero,b,1)
!Ly = b
    y = b
    do i = 1, n-1
        y(i) = y(i) / L(i,i)
        y(i+1:n) = y(i+1:n) - y(i) * L(i+1:n,i)
    enddo
    y(n) = y(n) / L(n,n)
!Ux=y  
    x=y
    do i = n,2,-1
      x(i) = x(i) / U(i,i)
      x(1:i-1) = x(1:i-1) - x(i) * U(1:i-1,i)
    enddo 
    x(1) = x(1) / U(1,1)
deallocate(y,L,U)
!
!!
!    
    allocate(work_1(n),work_2(n))
    work_1 = zero
    work_2 = zero
! work_1 = Ax
    call dgemv('n',n,n,one,matrix,n,x,1,zero,work_1,1)
    work_2 = work_1 - b 
    call norm_inf_vector(work_2,n,norm_1)
    call norm_inf_vector(b,n,norm_2)
    backward_error_ge = norm_1 / norm_2
    deallocate(b,work_1,work_2)
    rmse = sqrt( sum( (x-e)**2 ) / real(n, kind=8) )
    write(2,*) 'GE Root Mean Square                    is :', rmse
    write(2,*) 'GE Backward error                      is :', backward_error_ge
    allocate(x_ge(n))
    x_ge = x
    deallocate(e,x)
    !
    allocate(e(n),b(n),x(n))
    e = one
    x = zero
    b = zero
    call dgemv('n',n,n,one,matrix,n,e,1,zero,x,1)
    allocate(pivots(n))
    allocate(copy(n,n))
    copy = matrix
    b = x
    call dgesv(n,1,copy,n,pivots,x,n,info)
    if (info .eq. 0) then
        write(2,'(A)') ' The Linear Equation System: successful !'
    else
        write(2,'(A)') 'Error : The System of equations has not been solved !'
    endif
    !write(2,'(A)') 'Solution W_n x = b '
    !do i = 1, n 
    !  write(2,*) x(i)
    !enddo 
    deallocate(copy)
    allocate(work_1(n),work_2(n))
    work_1 = zero
    work_2 = zero
! work_1 = Ax
    call dgemv('n',n,n,one,matrix,n,x,1,zero,work_1,1)
! work_2 = Ax - b    
    work_2 = work_1 - b 
    call norm_inf_vector(work_2,n,norm_1)
    call norm_inf_vector(b,n,norm_2)
    deallocate(b,work_1,work_2)
    backward_error = norm_1 / norm_2
    rmse = sqrt( sum( (x-e)**2 ) / real(n, kind=8) )
    error_relative  = abs(backward_error_ge - backward_error)/abs(backward_error_ge)
    write(2,*) 'GEPP Backward error                    is :', backward_error
    write(2,*) 'GEPP Root Mean Square                  is :', rmse
    write(2,*) 'Relative Backward error                is :', error_relative
    deallocate(e,x,x_ge)
!    
close(2)
end subroutine

