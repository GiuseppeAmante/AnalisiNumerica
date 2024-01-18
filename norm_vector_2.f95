subroutine norm_vector_2(n,vector,norm)
    implicit none
    real(kind=8), dimension(n)            :: vector 
    real(kind=8)                          :: norm
    integer                               :: n
!
    integer :: i
    real(kind=8) :: sum_
    
    sum_ = 0.0d0
    do i = 1, n
      sum_ = sum_ + vector(i)*vector(i)
    enddo
    norm = sqrt(sum_)
end subroutine norm_vector_2
