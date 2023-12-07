program main
      implicit none
      real(kind= 8), allocatable :: matrix_a(:,:), matrix_l(:,:), matrix_u(:,:)
      integer :: n, n_wilkin
      real(kind=8) :: g
      character(len=3) :: response
      character(len=30) :: type_matrix
!   
      write(*,*) 'Do you have a matrix available?'
      read(*,*) response
      if (trim(response) == 'yes') then
           !read 
      else if (trim(response) == 'no') then
           write(*,*) 'What is the size of the matrix you want to factorize?'
           read (*,*) n
           allocate(matrix_a(n,n)) 
           call generation_matrix(matrix_a,n,type_matrix)
      endif
! Subroutine lufact() computes the LU fattorization of a non-singular nxn matrix A without pivoting        
!    lufact(A,L,U,g)
! where, takes A as input and 
!        returns the unit lower triangular factor L,
!        the upper triangular factor U, and
!        the growth factor g (defined here as the larger entry in the matrix G= |L||U|, divided by the largest entry in |A|)
      allocate (matrix_l(n,n),matrix_u(n,n))
      call lufact(matrix_a,matrix_l,matrix_u,g,n,type_matrix)
      deallocate(matrix_a,matrix_l,matrix_u)
!
      write(*,*) 'What is the size of the matrix Wilkinson?'
      read(*,*) n_wilkin
      allocate(matrix_a(n_wilkin,n_wilkin))
      call wilkin(n_wilkin,matrix_a)
      deallocate(matrix_a)
      end program 
