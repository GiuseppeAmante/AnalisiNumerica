subroutine norm_inf_matrix(A,m,n,norm)
    implicit none
    real(kind=8), intent(in) :: A(m,n)
    integer, intent(in) :: m, n
    real(kind=8), intent(out) :: norm
    real(kind=8) :: max_sum, row_sum
    integer :: i, j

    max_sum = 0.0d0

    do i = 1, m
        row_sum = 0.0d0
        do j = 1, n
            row_sum = row_sum + abs(A(i, j))
        end do
        max_sum = max(max_sum, row_sum)
    end do

    norm = max_sum
end subroutine
