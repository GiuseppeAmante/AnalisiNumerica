subroutine norm_inf_vector(x,n,norm)
    implicit none
    real(kind=8), dimension(n), intent(in) :: x
    integer :: n
    real(kind=8), intent(out) :: norm
!
    integer :: i
    norm = 0.0d0
    do i = 1, n
      norm = max( norm, abs( x(i) ) )
    end do
end subroutine norm_inf_vector
