subroutine intpol
    implicit none
    real(kind=8), dimension(9)   :: x, y
    real(kind=8), dimension(3)   :: sol_ne
    real(kind=8), dimension(9,3) :: V
    real(kind=8), allocatable    :: array(:), array_2(:),  matrix(:,:), L(:,:)
    integer :: i, j, info, lwork, ierr
    real(kind=8), parameter :: one=1.0d0, zero=0.0d0   
!
open(unit=3, file='intpol.txt', status='replace', action='write', iostat=ierr)
  if (ierr .ne. 0) then
    write (*,'(A)') 'Error! on open the file intpol.txt'
    stop
endif
!
    x(1) = 8d0
    x(2) = 10d0
    x(3) = 12d0
    x(4) = 16d0
    x(5) = 20d0
    x(6) = 30d0
    x(7) = 40d0
    x(8) = 60d0
    x(9) = 100d0
!
    y(1) = 0.88d0
    y(2) = 1.22d0
    y(3) = 1.64d0
    y(4) = 2.72d0
    y(5) = 3.96d0
    y(6) = 7.66d0
    y(7) = 11.96d0
    y(8) = 21.56d0
    y(9) = 43.16d0
!
! First way: normal equation  
!
    V(:,1) = one
    V(:,2) = x(:)
    allocate(array(9))
    array = x*x
    V(:,3) = array(:)
    deallocate(array)
    allocate(matrix(3,3),array(3))
    matrix = zero
write(3,*) 'Matrice A'
do i = 1, 9
    write(3,*) (V(i,j),j=1,3)
enddo
! matrix = A^T A    
    call dgemm('T','N',3,3,9,one,V,9,V,9,zero,matrix,3)
write(3,*) 'A^T A'
do i = 1, 3
    write(3,*) (matrix(i,j),j=1,3)
enddo
! array = A^T y    
    call dgemv('T',9,3,one,V,9,y,1,zero,array,1)
write(3,*) 'b = A^T y'
do i = 1, 3
    write(3,*) array(i)
enddo
! cholesky A^T A = L L^T where L is lower triangular    
    allocate(L(3,3))
    L = matrix
    call dpotrf('L',3,L,3,info)
write(3,*) 'Cholesky Uscita dpotrf'
do i = 1, 3
    write(3,*) (L(i,j),j=1,3)
enddo
    if (info == 0) write(*,*) 'Cholesky factorization successful exit'
    do i = 1, 3
      do j = i+1,3
        L(i,j)= zero
      enddo
    enddo
write(3,*) 'Cholesky sistemata'
do i = 1, 3
    write(3,*) (L(i,j),j=1,3)
enddo
! L array_2 = array( A^T y )    
    allocate(array_2(3))
    array_2= array
    do i = 1, 3-1
      array_2(i) = array_2(i)/L(i,i)
      array_2(i+1:3) = array_2(i+1:3) -array_2(i)*L(i+1:3,i)
    enddo
    array_2(3) = array_2(3)/L(3,3)
write(3,*) 'L y* = b'
do i = 1, 3
    write(3,*) array_2(i)
enddo
! L^T sol_ne = array_2    
    sol_ne= array_2
    do i = 1, 3-1
      sol_ne(i) = sol_ne(i)/L(i,i)
      sol_ne(i+1:3) = sol_ne(i+1:3) -sol_ne(i)*L(i,i+1:3)
    enddo
    sol_ne(3) = sol_ne(3)/L(3,3)
write(3,*) 'L^T x = y*'
do i = 1, 3
    write(3,*) sol_ne(i)
enddo
    deallocate(matrix,array,array_2,L)    
    write(3,'(A)') '---Solution Normal Equation---'
    do i = 1, 3
      write(3,*) sol_ne(i)
    enddo
! Second way : pivoted QR factorization    
    lwork = 3+3
    allocate(array(lwork))    
    call dgels('N',9,3,1,V,9,y,9,array,lwork,info) 
    if (info == 0) write(*,*) 'QR factorization successful exit'
    deallocate(array)
    write(3,'(A)') '---Solution QR FACTORIZATION---'
    do i = 1, 3
      write(3,*) y(i)
    enddo
close(3)
end subroutine intpol
    
