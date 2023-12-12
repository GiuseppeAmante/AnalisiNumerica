subroutine lu_step(step,n,A,pivot_null,A_k)
    implicit none
    integer :: n, step
    real(kind=8), dimension(n,n), intent(in) :: A
    real(kind=8), dimension(n,n), intent(out) :: A_k
    logical :: pivot_null
!
    real(kind=8), parameter :: one=1.0d0, zero=0.0d0
    real(kind=8) :: mach_precision
    real(kind=8), allocatable :: M(:,:,:), temp(:,:), MA(:,:)
    integer :: i, j, k  
    character(len=30) :: response
!
    mach_precision = epsilon(1.0d0)
    allocate(M(n,n,step))
    M = zero
    do i = 1, n
      do j = 1, step
        M(i,i,j) = one
      enddo
    enddo
    do i = 1, step
      if (A(i,i).eq.zero) then
          pivot_null = .TRUE.
        exit
      endif
      do j = i+1, n
! PROBLEMA è QUI
        M(j,i,i) = - A(j,i) / A(i,i) 
      enddo
    enddo
    write(1,*)  '---Matrix that composes A---'
    do i = 1, step
      do j = 1, n
        write(1,*) (M(j,k,i), k =1,n)
      enddo
      write(1,*) '----------------------------'
    enddo
!
!    allocate (temp(n,n)) 
!    temp=zero
!    do k = step, 2, -1 
!      temp(:, :) = M(:, :, k)
!      M(:, :, k) = M(:, :, k - 1)
!      M(:, :, k - 1) = temp(:, :)
!    end do
!    deallocate(temp)
!    write(1,*)  '---Matrix that composes A (ORDER)---'
!    do i = 1, step
!      do j = 1, n
!        write(1,*) (M(j,k,i), k =1,n)
!      enddo
!      write(1,*) '----------------------------'
!    enddo
!
    allocate(MA(n,n))
    call dgemm('n','n',n,n,n,one,M(:,:,1),n,A,n,zero,MA,n)
        response='matrix uscita dgemm M1A'
        call print_matrix(1,MA,n,n,response)
    M(:,:,1) = zero
    M(:,:,1) = MA
    deallocate(MA)
    if (step >= 2) then
      allocate(temp(n,n))
     temp = zero
!----------canc-----------
        response='M2'
        temp = M(:,:,2)
        call print_matrix(1,temp,n,n,response)
!--------canc----------
      temp=zero
      do i = 1, step-1
        call dgemm('n','n',n,n,n,one,M(:,:,i+1),n,M(:,:,i),n,zero,temp,n)
        response='matrix uscita dgemm M2M1A'
        call print_matrix(1,temp,n,n,response)
        M(:,:,i+1) = temp
      enddo
      deallocate(temp)
    endif
    write(*,*) 'k-step   =', step
    write(1,*) 'k-step   =', step
    response = 'Matrix A at k-step'
    A_k = M(:,:,step)
    call print_matrix(1,A_k,n,n,response)
    deallocate(M)
end subroutine