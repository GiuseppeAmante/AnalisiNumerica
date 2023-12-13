subroutine lu_step(step,n,A,pivot_null,A_k)
    implicit none
    integer :: n, step
    real(kind=8), dimension(n,n), intent(in) :: A
    real(kind=8), dimension(n,n), intent(out) :: A_k
    logical :: pivot_null
!
    real(kind=8), parameter :: one=1.0d0, zero=0.0d0
    real(kind=8) :: mach_precision
    real(kind=8), allocatable :: M(:,:), A_step(:,:,:)
    integer :: i, j, k   
!
    mach_precision = epsilon(1.0d0)
    allocate(M(n,n))
    M = zero
    do i =1 , n
      M(i,i) =one
    enddo
!    
    A_k=A
    allocate(A_step(n,n,step))
    do k = 1, step
!control: Pivoting Null  
      if (A_k(k,k).eq.zero) then
          pivot_null = .TRUE.
        exit
      endif
!      
      do i = k + 1, n
        M(i, k) = A_k(i, k) / A_k(k, k)
        if ( abs(M(i,k)) < mach_precision ) then 
          write(1,'(A)') 'Error : Division by a quantity smaller than machine precision'  
        endif
        do j = k, n
          A_k(i, j) = A_k(i, j) - M(i, k) * A_k(k, j)
        end do
      end do
      A_step(:, :, k) = A_k
    end do

    !write(1,*)  '---Matrix that composes A---'
    !do k = 1, step
    !  do i = 1, n
    !    write(1,*) (A_step(i,j,k), j =1,n)
    !  enddo
    !  write(1,*) '----------------------------'
    !enddo
end subroutine
