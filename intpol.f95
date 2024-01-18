subroutine intpol
    implicit none
    real(kind=8), dimension(9)   :: x, y
    real(kind=8), dimension(3)   :: sol_ne, alpha, x_propto, x_exact
    real(kind=8), dimension(9,3) :: V
    real(kind=8), allocatable    :: array(:), array_2(:), array_3(:), array_4(:), array_5(:), work(:)
    real(kind=8), allocatable    :: matrix(:,:), matrix_2(:,:), L(:,:)
    integer                      :: i, j, info, lwork, ierr
    integer, allocatable         :: ipiv(:)
    real(kind=8), parameter      :: one=1.0d0, zero=0.0d0   
    real(kind=8)  :: norm_r, norm_d, norm_x, norm_x_diff,&
                     norm_C, norm_C_inv, k, k_inv, cond, rcond
!
open(unit=3, file='intpol.txt', status='replace', action='write', iostat=ierr)
  if (ierr .ne. 0) then
    write (*,'(A)') 'Error! on open the file intpol.txt'
    stop
endif
!
! Start : Sub-problem3
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
!write(3,*) 'Matrice A'
!do i = 1, 9
!    write(3,*) (V(i,j),j=1,3)
!enddo

! matrix = A^T A    
    call dgemm('T','N',3,3,9,one,V,9,V,9,zero,matrix,3)
!write(3,*) 'A^T A'
!do i = 1, 3
!    write(3,*) (matrix(i,j),j=1,3)
!enddo

! array = A^T y    
    call dgemv('T',9,3,one,V,9,y,1,zero,array,1)
!write(3,*) 'b = A^T y'
!do i = 1, 3
!    write(3,*) array(i)
!enddo

! cholesky A^T A = L L^T where L is lower triangular    
    allocate(L(3,3))
    L = matrix
    call dpotrf('L',3,L,3,info)
!write(3,*) 'Cholesky Uscita dpotrf'
!do i = 1, 3
!    write(3,*) (L(i,j),j=1,3)
!enddo
    if (info == 0) write(*,*) 'Cholesky factorization successful exit'
    do i = 1, 3
      do j = i+1,3
        L(i,j)= zero
      enddo
    enddo
!write(3,*) 'Cholesky sistemata'
!do i = 1, 3
!    write(3,*) (L(i,j),j=1,3)
!enddo

! L array_2 = array( A^T y )    
    allocate(array_2(3))
    array_2= array
    do i = 1, 3-1
      array_2(i) = array_2(i)/L(i,i)
      array_2(i+1:3) = array_2(i+1:3) -array_2(i)*L(i+1:3,i)
    enddo
    array_2(3) = array_2(3)/L(3,3)
!write(3,*) 'L y* = b'
!do i = 1, 3
!    write(3,*) array_2(i)
!enddo

! L^T sol_ne = array_2    
    sol_ne= array_2
    do i = 3, 1, -1
      sol_ne(i) = sol_ne(i)/L(i,i)
      sol_ne(1:i-1) = sol_ne(1:i-1) -sol_ne(i)*L(i,1:i-1)
    enddo
!write(3,*) 'L^T x = y*'
!do i = 1, 3
!    write(3,*) sol_ne(i)
!enddo
    deallocate(matrix,array,array_2,L)    
    write(3,'(A)') '---Solution Normal Equation---'
    do i = 1, 3
      write(3,*) sol_ne(i)
    enddo
! Second way : pivoted QR factorization    
    lwork = 3+3
    allocate(work(lwork),matrix(9,3))
    matrix = V    
    call dgels('N',9,3,1,matrix,9,y,9,work,lwork,info) 
    if (info == 0) write(*,*) 'QR factorization successful exit'
    deallocate(work,matrix)
    write(3,'(A)') '---Solution QR FACTORIZATION---'
    do i = 1, 3
      write(3,*) y(i)
    enddo
! End: Sub-problem 3
  
! Start: Sub-problem 4
write(3,*) 'Data set X'
do i = 1, 9
  write(3,*) i, x(i)
enddo
write(3,*) 'Matrix Vamdermode Type (A)'
do i = 1, 9
  write(3,*) (V(i,j), j=1,3)
enddo
    alpha(1) = -1.919d0
    alpha(2) = 0.2782d0
    alpha(3) = 0.001739d0
! f(x) = alpha(1) + alpha(2)*x + alpha(3)*x^2 where f(x_) = b
    allocate(array(9))
    array=zero
write(3,*) 'b_i = f(x_i) = a_1 + a_2*x + a_3*x**2'
    do i = 1, 9
      array(i) = alpha(1) + alpha(2)*x(i) + alpha(3)*(x(i)*x(i))
write(3,*) i, array(i)
    enddo        
! array_2[d] = A^T array[b]    
    allocate(array_2(3))
    array_2=zero
    call dgemv('T',9,3,one,V,9,array,1,zero,array_2,1)
write(3,*) 'd = A^T b'
do i = 1, 3
  write(3,*) array_2(i)
enddo
! C[matrix_2] = A^T A
    allocate(matrix_2(3,3))
    matrix_2=zero
    call dgemm('T','N',3,3,9,one,V,9,V,9,zero,matrix_2,3)
write(3,*) 'C = A^T * A'
do i = 1, 3
  write(3,*) (matrix_2(i,j), j=1,3)
enddo
! cholesky A^T A = L L^T where L is lower triangular    
    allocate(L(3,3))
    L = matrix_2
    call dpotrf('L',3,L,3,info)
    if (info == 0) write(*,*) 'Cholesky factorization successful exit'
    do i = 1, 3
      do j = i+1,3
        L(i,j)= zero
      enddo
    enddo
write(3,*) 'A --Cholesky--> L'
do i = 1, 3
  write(3,*) (L(i,j), j=1,3)
enddo
!  L y[array_3] = array_2[d]
    allocate(array_3(3))
    array_3 = array_2
    do i = 1, 3-1
      array_3(i) = array_3(i)/L(i,i)
      array_3(i+1:3) = array_3(i+1:3) -array_3(i)*L(i+1:3,i)
    enddo
    array_3(3) = array_3(3)/L(3,3)
write(3,*) 'L y = d'
do i = 1, 3
  write(3,*) array_3(i)
enddo 
!  L^T x[x_propto] = y[array_3]
    x_propto= array_3
    deallocate(array_3)
    do i = 3, 1, -1
      x_propto(i) = x_propto(i)/L(i,i)
      x_propto(1:i-1) = x_propto(1:i-1) - x_propto(i)*L(i,1:i-1)
    enddo
    deallocate(L)
    write(3,'(A)') '---x[propto]--'
    do i = 1, 3
      write(3,*) x_propto(i)
    enddo
!  x_exact[array]
    lwork = 3+3
    allocate(work(lwork),matrix(9,3))
    matrix = V    
    call dgels('N',9,3,1,matrix,9,array,9,work,lwork,info) 
    if (info == 0) write(*,*) 'QR factorization successful exit'
    x_exact = array(1:3)
    deallocate(work,matrix,array)
    write(3,'(A)') '---x[exact]---'
    do i = 1, 3
      write(3,*) x_exact(i)
    enddo
! r[array_3] = d[array_2] - C[matrix_2] x_propto
    allocate(array_3(3),array_4(3))
    call dgemv('N',3,3,one,matrix_2,3,x_propto,1,zero,array_4,1)
    array_3 = array_2 - array_4
write(3,*) 'Cx_propto', (array_4(i), i=1,3)
    deallocate(array_4)
write(3,*) 'C = A^T *A '
do i = 1, 3
write(3,*) (matrix_2(i,j), j=1,3)
enddo
write(3,*) 'd', (array_2(i), i=1,3)
    write(3,'(A)') 'Calculate Residual r = d - Cx_propto'
    do i = 1, 3
      write(3,*) array_3(i)
    enddo
    call norm_vector_2(3,array_3,norm_r)
write(3,*) 'NORM RESIDUAL', norm_r
    call norm_vector_2(3,array_2,norm_d)
write(3,*) 'NORM Vector d', norm_d
    call norm_vector_2(3,x_exact,norm_x)
write(3,*) 'NORM Exact x', norm_x
    allocate(array(3))
    array = x_exact-x_propto 
    call norm_vector_2(3,array,norm_x_diff)
write(3,*) 'NORM X_exact-X_propto', norm_x_diff
    deallocate(array)
    call dlange('E',3,3,matrix_2,3,norm_C)
write(3,*) 'NORM MATRIX C', norm_C
    lwork = 3
    allocate(ipiv(3),work(lwork))
    call dgetrf(3,3,matrix_2,3,ipiv,info)
    if(info == 0) write(*,*) 'LU decompision successful exit'
    call dgetri(3,matrix_2,3,ipiv,work,lwork,info)
    deallocate(ipiv,work)
    if(info == 0) write(*,*) 'Calculate Inverse Matrix successful exit'
    call dlange('E',3,3,matrix_2,3,norm_C_inv)
write(3,*) 'NORM MATRIX INVERSE C', norm_C_inv
    k = norm_C * norm_C_inv
write(3,*) 'k', k
    k_inv = one/k
write(3,*) 'k^-1', k_inv
    allocate(work(5*3),array(3),array_4(3), array_5(3))
    call dgesvd('N','N',3,3,matrix_2,3,array,array_4,3,array_5,3,work,5*3,info)
    if (info == 0) write(*,*) 'successful exit'
    write(3,*) 'Valori singolari'
    write(3,*) (array(i), i= 1, 3)
    cond = array(1)/array(3)
    rcond = one/cond
    write(3,*) 'Numero di condizionamento calcolato con SVD', cond
    deallocate(work,array,array_4,array_5)
    deallocate(array_2,array_3,matrix_2)
! stima dell'errore relativo tra la soluzione esatta e quella approssimata
    if ( ( rcond*norm_r/norm_d .le. norm_x_diff/norm_x ) .and.& 
         ( norm_x_diff/norm_x .le. cond*norm_r/norm_d ) ) then
      write(3,*) 'OK'
! La soluzione approssimata è accurata entro un centro regime
    else 
      write(3,*) 'A caga o cazz ei'
! La soluzione approssimata NON è accurata 
    endif
! End: Sub-problem 4
close(3)
end subroutine intpol
    
