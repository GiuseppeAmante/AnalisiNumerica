subroutine wilkin(n,matrix)
    implicit none
    real(kind=8), dimension(n,n) :: matrix
    integer :: n, i , j, ierr, info
    real(kind=8) :: error_norm, rmse
    character(len=30) :: type_matrix
    real(kind=8), allocatable :: e(:), x(:) 
    integer, allocatable :: pivots(:)
    real(kind=8), parameter :: zero = 0.0d0, one=1.0d0
!
open(unit=2, file='wilkin.txt', status='replace', action='write', iostat=ierr)
  if (ierr .ne. 0) then
    write (*,'(A)') 'Error:  open the file wilkin.txt !'
    stop
  endif
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
    allocate(e(n),x(n))
    e = one
    x = zero
    call dgemv('n',n,n,one,matrix,n,e,1,zero,x,1)
    allocate(pivots(n))
    call dgesv(n,1,matrix,n,pivots,x,n,info)
    if (info .eq. 0) then
        write(2,'(A)') ' The Linear Equation System: successful !'
    else
        write(2,'(A)') 'Error : The System of equations has not been solved !'
    endif
    error_norm = sqrt(sum((e - x)**2))
    rmse = sqrt(sum((e - x)**2) / real(n))
    write(2,*) 'Euclidean Norm of the Error (e - x) is : ', error_norm
    write(2,*) 'Root Mean Square (x - e)            is : ', rmse
    deallocate(e,x)
!    
close(2)
end subroutine

