subroutine lufact(A,n,response_A,L_np,U_np)
    implicit none
    real(kind=8), dimension(n,n) :: A
    real(kind=8), dimension(n,n), intent(out) :: L_np, U_np
    integer :: n
!    
    reaL(kind=8), allocatable :: L(:,:), U(:,:), copy_1(:,:), copy_2(:,:),&
                                 work(:)
    integer, allocatable :: iwork(:)
    real(kind=8) :: g, val_1, val_2, val_3, val_4, mach_precision,&
                    norm_matrix_1, norm_matrix_2, nrb_error,&
                    value_min, diff, max_diff, rcond, soglia_pivot
    integer :: i, j, k, ierr, info, ipivot
    real(kind=8), parameter :: zero = 0.0d0, one=1.0d0
    character(len=10) :: response_A 
    character(len=30) :: response
    logical :: no_lufact
!
    no_lufact = .false.
    ipivot= 1
!    
open(unit=1, file='lufact.txt', status='replace', action='write', iostat=ierr)
  if (ierr .ne. 0) then
    write (*,'(A)') 'Error! on open the file factLU.txt'
    stop
  endif
!
!START: FACTORIZATION LU without pivoting
!
    mach_precision = epsilon(1.0d0)
    !soglia_pivot   = epsilon(1.0d0)*10.0d0
    soglia_pivot   = zero
    write(*,'(A)') 'Welcome to LU factorization without pivoting.'
    write(1,'(A)') '====================================================================='
    write(1,'(A)') '                 LU factorization without pivoting'
    write(1,'(A)') '====================================================================='
    allocate(L(n,n),U(n,n))
    !call print_matrix(1,A,n,n,response_A)
    L = zero
    U = A
!
    outer_loop: do k = 1, n-1
      L(k, k) = one
! control: pivot
      if (abs(U(k,k)) < soglia_pivot) then
        write(*,'(A)') 'Error : found null pivot !'
        write(1,'(A)') 'Error : found null pivot !'
        no_lufact= .true.
        ipivot = 0
        exit outer_loop 
      endif 
      if (k .eq. n-1) write(1,'(A)') 'Pivot Control: passed !' 
      do i = k+1, n
        L(i,k) = U(i,k) / U(k,k)
        write(1,*) 'L', L(i,k)
        write(1,*) 'A_kk^k', U(k,k)
! control: machine precision
        if (abs(U(k,k)) < mach_precision) then 
          write(*,'(A)') 'Error : Division by a quantity smaller than machine precision'
          write(1,'(A)') 'Error : Division by a quantity smaller than machine precision'
          no_lufact= .true.
          exit outer_loop
        endif
        if (k .eq. n-1) write(1,'(A)') 'Control on machine precision value: passed !' 
        U(i,k+1:n) = U(i,k+1:n) - L(i,k) * U(k,k+1:n)
      enddo
    enddo outer_loop
    L(n,n) = one
    do i = 1, n
       do j = 1, i-1
           U(i,j) = zero
       enddo
   enddo
!
!! 
!   
!!
!    
    if (no_lufact) then
      write(*,'(A)') 'LU factorization without pivoting: failed !'
      write(1,'(A)') 'LU factorization without pivoting: failed !'
      if (ipivot == 0) write(1,'(A)') '                      due to a null pivot.' 
    else    
      write(*,'(A)') 'LU factorization without pivoting: successful !'
      write(1,'(A)') 'LU factorization without pivoting: successful !'
      !response = 'Lower Triangular'
      !call print_matrix(1,L,n,n,response)
      L_np = L
      !response = 'Upper Triangular'
      !call print_matrix(1,U,n,n,response)
      U_np = U
!  growth factor g = max( G ) / max( |A| ) where G = |L| x |U|
      val_1 = maxval(abs(L)) 
      val_2 = maxval(abs(U)) 
      val_3 = val_1*val_2
      val_4 = maxval(abs(A))
      write(1,*) 'Maxval matrix G = |L| x |U|               :', val_3
      write(1,*) 'Maxval matrix |A|                         :', val_4
      g = val_3/val_4
      write(1,*) 'Growth factor (without pivoting) is       :', g
! A - LU
      allocate (copy_1(n,n), copy_2(n,n))
      copy_1 = zero
      copy_2 = zero
      !copy_1 = L x U
      call dgemm('n','n',n,n,n,one,L,n,U,n,zero,copy_1,n)
      !copy_2 = A - LU
      copy_2 = A - copy_1
      deallocate(copy_1)
! normwise relative backward error
      call norm_inf_matrix(copy_2,n,n,norm_matrix_1)
      write(1,*) 'Norm inf A - LU                           :', norm_matrix_1
      call norm_inf_matrix(A,n,n,norm_matrix_2) 
      write(1,*) 'Norm inf A                                :', norm_matrix_2
      nrb_error = norm_matrix_1/norm_matrix_2 
      write(1,*) 'Normwise relative backward error (without pivoting) is :',  nrb_error
! TEST   
      max_diff = zero
      do i = 1, n
        do j = 1, n
          diff = abs(copy_2(i,j))
          if (diff > max_diff) max_diff = diff
        end do
      end do
      deallocate(copy_2)
      value_min = 1e-10
      if (max_diff < value_min) then
        write(*, '(A)') "LU factorization without pivoting is acceptable."
      else
        write(*, '(A)') "LU factorization without pivoting is not acceptable."
      endif
    endif
    deallocate(L,U)
    write(1,'(A)') 'Goodbye from factorization LU without pivoting' 
    write(*,'(A)') 'Goodbye from factorization LU without pivoting' 
!
!END: FACTORIZZATION LU without pivoting
!
!  LAPACK: DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
!
    write(*,'(A)') 'Welcome from LAPACK: LU factorization with partial pivoting' 
    write(1,'(A)') '====================================================================='
    write(1,'(A)') '           LAPACK: LU factorization with partial pivoting'
    write(1,'(A)') '====================================================================='
    allocate(L(n,n),U(n,n))
    L = zero
    U = zero
!
    !call print_matrix(1,copy,n,n,response_A)
!
    allocate(copy_1(n,n),iwork(n))
    copy_1 = A
    call dgetrf(n,n,copy_1,n,iwork,info)
    if (info .ne. 0) then
      write(1,'(A)') 'LU factorization with partial pivoting [LAPACK] : failed !'
      stop
    endif
    do i = 1, n
      L(i, i) = one
      do j = 1, i - 1
        L(i, j) = copy_1(i, j)
      enddo
      do j = i, n
        U(i, j) = copy_1(i, j)
      enddo
    enddo
    deallocate(copy_1,iwork)
!
    !response = 'Lower Triangular [LAPACK]'
    !call print_matrix(1,L,n,n,response)
    !response = 'Upper Triangular [LAPACK]'
    !call print_matrix(1,U,n,n,response)
!    
!growth factor g = max( |U| ) / max( |A| )
    val_1 = maxval(abs(U))
    val_2 = maxval(abs(A))
    write(1,*) 'Maxval matrix |U|                         :', val_1
    write(1,*) 'Maxval matrix |A|                         :', val_2
    g = val_1/val_2
    write(1,*) 'Growth factor (with partial pivoting)  is :', g
!A - LU
    allocate (copy_1(n,n), copy_2(n,n))
    copy_1 = zero
    copy_2 = zero
    call dgemm('n','n',n,n,n,one,L,n,U,n,zero,copy_1,n)
    copy_2 = A - copy_1
    deallocate(copy_1,L,U)
!normwise relative backward error
    call norm_inf_matrix(copy_2,n,n,norm_matrix_1)
    write(1,*) 'Norm inf A - LU                           :', norm_matrix_1
    deallocate(copy_2)
    call norm_inf_matrix(A,n,n,norm_matrix_2) 
    write(1,*) 'Norm inf A                                :', norm_matrix_2
    nrb_error = norm_matrix_1 / norm_matrix_2 
    write(1,*) 'Normwise relative backward error (with partial pivoting)   is : ', nrb_error
!CONDITION NUMBER    
! The reciprocal of the condition number of the matrix A, computed as RCOND = 1/(norm(A) * norm(inv(A)))
    allocate(work(4*n),iwork(n))
    work = zero
    iwork = 0
    call dgecon('I',n,A,n,norm_matrix_2,rcond,work,iwork,info)
    if (info .ne. 0) write(1,'(A)') 'Condition number computation : failed !'
    write(1,*) 'Reciprocal of the condition number of the matrix ', response_A, 'is :', rcond
    deallocate(work, iwork)
    write(*,'(A)') 'Goodbye from LAPACK: LU factorization with partial pivoting' 
    write(1,'(A)') 'Goodbye from LAPACK: LU factorization with partial pivoting' 
close(1)
    end subroutine
