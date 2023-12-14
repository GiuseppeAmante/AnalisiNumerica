subroutine lufact(A,n,response_A)
    implicit none
    real(kind=8), dimension(n,n) :: A
    reaL(kind=8), allocatable :: L(:,:), U(:,:), copy(:,:), copy_copy(:,:),&
        test(:,:), ALU(:,:), A_k(:,:), g_k(:), work(:)
    integer, allocatable :: iwork(:)
    real(kind=8) :: g, val_1, val_2, val_3, val_4, mach_precision,&
       norm_matrix_1, norm_matrix_2, nrb_error, value_min, diff, max_diff,&
       rcond
    integer :: n, i, j, k, ierr, info
    integer, dimension(n) :: ipiv
    real(kind=8), parameter :: zero = 0.0d0, one=1.0d0
    character(len=10) :: response_A 
    character(len=30) :: response
    logical :: pivot_null, machine_prec, no_lufact
    real(kind=8), external :: dlange

    no_lufact = .false.
!    
open(unit=1, file='lufact.txt', status='replace', action='write', iostat=ierr)
  if (ierr.ne.0) then
    write (*,'(A)') 'Error! on open the file factLU.txt'
    stop
  endif
!
!START: FACTORIZATION LU
!
    write(*,'(A)') 'Welcome to LU factorization without pivoting.'
    write(1,'(A)') 'Welcome to LU factorization without pivoting.'
    allocate(L(n,n),U(n,n))
    mach_precision = epsilon(1.0d0)
    g = zero
    L = zero
    U = zero
    !call print_matrix(1,A,n,n,response_A)
    allocate (copy(n,n))
    copy = A
!
    outer_loop: do k = 1, n
      L(k, k) = one
! control: pivot
      if (A(k,k).eq.zero) then
        write(*,'(A)') 'Error : found null pivot'
        exit outer_loop 
      endif 
      if (k .eq. n) write(1,'(A)') 'Pivot control passed.' 
      do i = k+1, n
        L(i,k) = A(i,k) / A(k,k)
! control: machine precision
        !if (abs(L(i,k)) < mach_precision) then  ERRORE
        if (abs(A(k,k)) < mach_precision) then
          write(*,'(A)') 'Error : Division by a quantity smaller than machine precision'
          no_lufact= .true.
          exit outer_loop
        endif
    enddo
      if (k .eq. n) write(1,'(A)') 'Control on machine precision value passed.' 
      do j = k, n
        U(k,j) = A(k,j)
      enddo
      do i = k+1, n
        do j = k+1, n
          A(i,j) = A(i,j) - L(i,k)*U(k,j)
        enddo
      enddo
    end do outer_loop
!
!! 
!   
!!
!    
    if (no_lufact) then
      write(*,'(A)') 'LU factorization without pivoting: termined !'
    else    
      write(*,'(A)') 'LU factorization without pivoting: successful !'
      response = 'Lower Triangular'
      !call print_matrix(1,L,n,n,response)
      !response = 'Upper Triangular'
      !call print_matrix(1,U,n,n,response)
!  growth factor g = max( G ) / max( |A| ) where G = |L| x |U|
      val_1 = maxval(abs(L)) 
      val_2 = maxval(abs(U)) 
      val_3 = val_1*val_2
      val_4 = maxval(abs(copy))
      write(1,*) 'Maxval matrix G = |L| x |U|               :', val_3
      write(1,*) 'Maxval matrix |A|                         :', val_4
      if ( val_4 < mach_precision ) then
        write(1,'(A)') 'Error : Division by a quantity smaller than machine precision for GROWTH FACTOR'
      else
        g = val_3/val_4
        write(1,*) 'Growth factor (without pivoting) is       :', g
      endif
! A - LU
      allocate (test(n,n), ALU(n,n))
      test = zero
      ALU = zero
      call dgemm('n','n',n,n,n,one,L,n,U,n,zero,test,n)
      response = 'A'
      call print_matrix(1,copy,n,n,response)
      response = 'LU'
      call print_matrix(1,test,n,n,response)
      ALU = copy - test
      response = 'A - LU'
      call print_matrix(1,ALU,n,n,response)
      deallocate(test)
! normwise relative backward error
      call norm_inf_matrix(ALU,n,n,norm_matrix_1)
      write(1,*) 'Norm inf A - LU                           :', norm_matrix_1
      !
      allocate(work(n))
      work = zero
      norm_matrix_1 = zero
      norm_matrix_1 = dlange('i',n,n,ALU,n,work)
      write(1,*) 'Norm inf A - LU  [LAPACK]                 :', norm_matrix_1
      !do i =1, n
      !    write(1,*) work(i)
      !enddo  
      deallocate(work)
      !
      call norm_inf_matrix(copy,n,n,norm_matrix_2) 
      write(1,*) 'Norm inf A                                :', norm_matrix_2
      if ( norm_matrix_2 < mach_precision ) then
        write(1,'(A)') 'Error : Division by a quantity smaller than machine precision for NORWISE RELATIVE BACKWARD ERROR'
      else
        nrb_error = norm_matrix_1/norm_matrix_2 
        write(1,*) 'Normwise relative backward error (without pivoting) is :',  nrb_error
      endif
! TEST   
      max_diff = zero
      do i = 1, n
        do j = 1, n
          diff = abs(ALU(i,j))
          if (diff > max_diff) max_diff = diff
        end do
      end do
      deallocate(ALU)
  ! Soglia (ad esempio, 1e-10)
      value_min = 1e-10
  ! Valuta se la differenza fra il valore assoluto di A-LU e la soglia è accettabile
      if (max_diff < value_min) then
        write(*, '(A)') "LU factorization without pivoting is acceptable."
      else
        write(*, '(A)') "LU factorization without pivoting is not acceptable."
      end if
    endif
    write(1,'(A)') 'Goodbye from factorization LU without pivoting' 
    write(*,'(A)') 'Goodbye from factorization LU without pivoting' 
!
!END: FACTORIZZATION LU
!
    write(*,'(A)') 'Welcome from Analysis step by step LU factotization' 
    write(1,'(A)') '====================================================================='
    write(1,'(A)') '                 ANALYSIS STEP BY STEP LU FACTORIZATION'
    write(1,'(A)') '====================================================================='
!
    allocate(g_k(n-1),A_k(n,n))
    A_k = zero
    g_k = zero
    pivot_null = .FALSE.
    machine_prec = .FALSE.
    outer_loop_2: do k = 1, n-1
      call lu_step(k,n,copy,pivot_null,machine_prec,A_k)
      if (pivot_null) then
        write(*,*) 'Null Pivot Found at Step', k
        exit outer_loop_2
      elseif (machine_prec) then
        write(*,*) 'Error at Step', k, 'Division by a quantity smaller than machine precision'
        exit outer_loop_2
      else
        val_1 = maxval(abs(A_k))
        val_2 = maxval(abs(copy))
        g = val_1 / val_2
        write(1,*) 'at Step', k, '  growth factor  :', g
        g_k(k) = g
      endif
    enddo outer_loop_2
    deallocate(g_k,A_k)
    write(*,'(A)') 'Goodbye from Analysis step by step LU factotization' 
!
!  LAPACK: DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.

!
    write(*,'(A)') 'Welcome from LAPACK: LU factorization with partial pivoting' 
    write(1,'(A)') '====================================================================='
    write(1,'(A)') '           LAPACK: LU factorization with partial pivoting'
    write(1,'(A)') '====================================================================='
    L = zero
    U = zero
!
    !call print_matrix(1,copy,n,n,response_A)
!
    allocate(copy_copy(n,n))
    copy_copy = copy
    call dgetrf(n, n, copy, n, ipiv, info)
    if (info .ne. 0) then
      write(1,'(A)') 'LU factorization with partial pivoting [LAPACK] : failed !'
      stop
    endif
    do i = 1, n
      L(i, i) = one
      do j = 1, i - 1
        L(i, j) = copy(i, j)
      enddo
      do j = i, n
        U(i, j) = copy(i, j)
      enddo
    enddo
    deallocate(copy)
!
    !response = 'Lower Triangular [LAPACK]'
    !call print_matrix(1,L,n,n,response)
    !response = 'Upper Triangular [LAPACK]'
    !call print_matrix(1,U,n,n,response)
!    
!growth factor g = max( |U| ) / max( |A| )
    val_1 = maxval(abs(U))
    val_2 = maxval(abs(copy_copy))
    write(1,*) 'Maxval matrix |U|                         :', val_1
    write(1,*) 'Maxval matrix |A|                         :', val_2
    g = val_1/val_2
    if ( g < mach_precision ) then
       write(1,'(A)') 'Error : Division by a quantity smaller than machine precision for GROWTH FACTOR'
       write(*,'(A)') 'Error : Division by a quantity smaller than machine precision for GROWTH FACTOR'
    else
       write(1,*) 'Growth factor (with partial pivoting)  is :', g
    endif
!A - LU
    allocate (test(n,n), ALU(n,n))
    test = zero
    ALU = zero
    call dgemm('n','n',n,n,n,one,L,n,U,n,zero,test,n)
    ALU = copy_copy - test
    deallocate(test)
!normwise relative backward error
    call norm_inf_matrix(ALU,n,n,norm_matrix_1)
    write(1,*) 'Norm inf A - LU                           :', norm_matrix_1
    !
    allocate(work(n))
    work = zero
    norm_matrix_1 = zero
    norm_matrix_1 = dlange('i',n,n,ALU,n,work)
    write(1,*) 'Norm inf A - LU  [LAPACK]                 :', norm_matrix_1
    !do i =1, n
    !    write(1,*) work(i)
    !enddo  
    deallocate(work)
    deallocate(ALU)
    !
    call norm_inf_matrix(copy_copy,n,n,norm_matrix_2) 
    write(1,*) 'Norm inf A                                :', norm_matrix_2
    !
    allocate(work(n))
    work = zero
    norm_matrix_2 = zero
    norm_matrix_2 = dlange('i',n,n,copy_copy,n,work)
    write(1,*) 'Norm inf A  [LAPACK]                      :', norm_matrix_2
    !do i =1, n
    !    write(1,*) work(i)
    !enddo  
    deallocate(work)
    !
!CONDITION NUMBER    
! The reciprocal of the condition number of the matrix A, computed as RCOND = 1/(norm(A) * norm(inv(A)))
    allocate(work(4*n),iwork(n))
    work = zero
    iwork = 0
    call dgecon('I',n,copy_copy,n,norm_matrix_2,rcond,work,iwork,info)
    if (info .ne. 0) then
      write(1,'(A)') 'Condition number computation : failed !'
    endif
    write(1,*) 'Reciprocal of the condition number of the matrix ', response_A, 'is :', rcond
    deallocate(work,iwork)
!NORMWISE RELATIVE BACKWARD ERROR
    nrb_error = norm_matrix_1/norm_matrix_2 
    if ( nrb_error < mach_precision ) then
       write(1,'(A)') 'Error : Division by a quantity smaller than machine precision for NORWISE RELATIVE BACKWARD ERROR'
    else
       write(1,*) 'Normwise relative backward error (with partial pivoting)   is : ',  nrb_error
    endif
    deallocate(L, U, copy_copy)
    write(*,'(A)') 'Goodbye from LAPACK: LU factorization with partial pivoting' 
close(1)
    end subroutine
