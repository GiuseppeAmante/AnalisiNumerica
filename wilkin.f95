subroutine wilkin(n,matrix)
    implicit none
    real(kind=8), dimension(n,n), intent(out) :: matrix
    integer :: n, i , j, ierr
    character(len=30) :: type_matrix
!
open(unit=2, file='wilkin.txt', status='replace', action='write', iostat=ierr)
  if (ierr.ne.0) then
    write (*,'(A)') 'Error! on open the file wilkin.txt'
    stop
  endif

    matrix= 0.0d0
    do i = 1, n
      matrix(i,i) = 1.0d0
      matrix(i,n) = 1.0d0
    enddo
    do i = 2, n
      do j = 1, i-1
        matrix(i,j) = -1.0d0
      enddo
    enddo

    type_matrix = 'Wilkinson'
    call print_matrix(2,matrix,n,n,type_matrix)
close(2)
end subroutine

