subroutine lufact(A,L,U,g,n,response)
    implicit none
    real(kind=8), dimension(n,n) :: A
    real(kind=8), dimension(n,n), intent(out) :: L, U
    reaL(kind=8), allocatable :: copy(:,:), matrix_G(:,:), matrix_error(:,:), copy_copy(:,:), test(:,:), A_k(:,:), g_k(:)
    real(kind=8) :: g, val_1, val_2, mach_precision, norm_matrix_1, norm_matrix_2, nrb_error
    integer :: n, i, j, k, ierr, info
    integer, dimension(n) :: ipiv
    real(kind=8), parameter :: zero = 0.0d0, one=1.0d0
    character(len=30) :: response 
    logical :: pivot_null
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
    call print_matrix(1,A,n,n,response)
    allocate (copy(n,n))
    copy = A
!
    do k = 1, n
      L(k, k) = one
! control: pivot
      if (A(k,k).eq.zero) then
        write(*,'(A)') 'Null Pivot Found'
        exit
      endif 
      do i = k+1, n
        L(i,k) = A(i,k) / A(k,k)
        if (abs(L(i,k)) < mach_precision) then
          write(*,'(A)') 'Error! Factorizzation LU without pivoting terminated'
          write(*,'(A)') 'due to division by a quantity smaller than machine precision'
          stop
        endif
      enddo
      do j = k, n
        U(k,j) = A(k,j)
      enddo
! control: machine precision
      do i = k+1, n
        do j = k+1, n
          A(i,j) = A(i,j) - L(i,k)*U(k,j)
        enddo
      enddo
    end do
    write(*,'(A)') 'Factorizzation LU without pivoting successful'
    response = 'Lower Triangular'
    call print_matrix(1,L,n,n,response)
    response = 'Upper Triangular'
    call print_matrix(1,U,n,n,response)
! Test
    allocate(test(n,n)) 
    test=zero
    call dgemm('n','n',n,n,n,one,L,n,U,n,zero,test,n)
    response = 'Is it equal to origin matrix?'
    call print_matrix(1,test,n,n,response)
    do i = 1, n
      do j = 1, n
        test(i,j) = copy(i,j) - test(i,j)
        if (test(i,j).ne.zero) then 
          write(*,'(A)') 'Test factorizzation LU without pivoting FAILED!'
          stop
        endif
      enddo
    enddo
    write(1,'(A)') 'YES!'
    write(*,'(A)') 'Test factorizzation LU without pivoting successful'
    deallocate(test)
 
    write(1,'(A)') '====================================================================='
    write(1,'(A)') '                 Analysis step by step factorizzation LU'
    write(1,'(A)') '====================================================================='
!
    allocate(g_k(n-1),A_k(n,n))
    pivot_null = .FALSE.
    do k = 1, n-1
      L = zero
      U = zero
      call lu_step(k,n,copy,pivot_null,A_k)
      if (pivot_null) then
        write(*,*) 'Null Pivot Found at Step', k
        exit
      else 
        val_1 = maxval(abs(A_k))
        val_2 = maxval(abs(copy))
        g = val_1/val_2
        write(1,*) 'at Step', k, '  growth factor  :', g
        g_k(k) = g
      endif
    enddo
    deallocate(g_k,A_k)
!
!  LAPACK: DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
    L = zero
    U = zero
!
    response = 'Matrix before dgetrf[LAPACK]'
    call print_matrix(1,copy,n,n,response)
!
    allocate(copy_copy(n,n))
    copy_copy = copy
    CALL DGETRF(n, n, copy, n, ipiv, info)
    if (info == 0)  write(1,*) 'DGETRF: successful exit!'
    DO i = 1, n
      L(i, i) = one
      DO j = 1, i - 1
        L(i, j) = copy(i, j)
      END DO
      DO j = i, n
        U(i, j) = copy(i, j)
      END DO
    END DO
!
    response = 'Matrix after dgetrf[LAPACK]'
    call print_matrix(1,copy,n,n,response)
!
    response = 'Lower Triangular - LAPACK'
    call print_matrix(1,L,n,n,response)
    response = 'Upper Triangular - LAPACK'
    call print_matrix(1,U,n,n,response)

!  growth factor g = max( |U| ) / max( |A| )
    val_1 = maxval(abs(U))
    val_2 = maxval(abs(copy_copy))
    write(1,*) 'Maxval matrix |U|           :', val_1
    write(1,*) 'Maxval matrix |A|           :', val_2
    g = val_1/val_2
    if ( g < mach_precision ) then
       write(1,'(A)') 'Error: Division by a quantity smaller than machine precision for GROWTH FACTOR'
    else
       write(1,*) 'Growth factor (with partial pivoting) is :', g
    endif
!  A - LU
    allocate (matrix_G(n,n), matrix_error(n,n))
    matrix_G = zero
    matrix_error = zero
    do i = 1, n
      do j = 1, n
        matrix_G(i,j) = zero
        do k =1, n
          matrix_G(i,j) = matrix_G(i,j) + L(i,k)*U(k,j)
        enddo
        matrix_error(i,j) = copy_copy(i,j) - matrix_G(i,j)
      enddo
    enddo
    deallocate(matrix_G)
! normwise relative backward error
    call norm_inf_matrix(matrix_error,n,n,norm_matrix_1)
    write(1,*) 'Norm inf A - LU             :', norm_matrix_1
    call norm_inf_matrix(copy_copy,n,n,norm_matrix_2) 
    write(1,*) 'Norm inf A                  :', norm_matrix_2
    nrb_error = norm_matrix_1/norm_matrix_2 
    write(1,*) 'Norm inf A - LU/ Norm inf A :', nrb_error
    !write(1,*) 'Machine precision           :', mach_precision
    if ( nrb_error < mach_precision ) then
       write(1,'(A)') 'Error: Division by a quantity smaller than machine precision for NORWISE RELATIVE BACKWARD ERROR'
    else
       write(1,*) 'Normwise relative backward error is :',  nrb_error
    endif
    deallocate(matrix_error)
!
    deallocate(copy_copy)
close(1)
    end subroutine
