program main
      implicit none
      real(kind= 8), allocatable :: matrix_a(:,:), L(:,:), U(:,:)
      integer :: n, n_wilkin, success, i, j, n_project
      character(len=3) :: response
      character(len=10) :: type_matrix, problem
!   
write(*,*) 'Welcome to my program, I am Giuseppe Amante. '
  write(*,*) 'To which Numerical Algebra Home Project are you referring?[1 or 2]'
  read (*,*) n_project
  if (n_project .eq. 1) then
      write(*,*) 'This program solves exercises on linear system solutions'
      write(*,*) 'Starting !'
!START:1      
      write(*,*) 'Which issue would you like me to address?[Problem1 or Problem2]'
      read(*,*) problem
      if (trim(problem) == 'Problem1' ) then
        write(*,*) 'Do you have a matrix available?[yes, no]'
        read(*,*) response
        if (trim(response) == 'yes') then
          write(*,*) "Enter the name of the file containing the matrix (make sure it is square):"
          read(*,*) type_matrix
          open(unit=10, file=trim(type_matrix), status='old', action='read', iostat=success)
          if (success .ne. 0) then
            write(*,*) "Error in opening the file."
            write(*,*) "Remember: Create a file and insert the elements of a matrix inside"  
            stop
          end if
          n=0 
          do
            read(10, *, iostat=success) 
            if (success .ne. 0) then
              exit
            else
              n = n + 1
            end if
          enddo
          rewind(unit=10)
          allocate(matrix_a(n,n))
          do i = 1, n
            read(10,*) (matrix_a(i,j), j =1 , n)
        enddo
          close(10)
        else if (trim(response) == 'no') then
          write(*,*) 'What is the size of the matrix you want to factorize?'
          read (*,*) n
          allocate(matrix_a(n,n)) 
          call generation_matrix(matrix_a,n,type_matrix)
        endif
        allocate(L(n,n),U(n,n))     
        call lufact(matrix_a,n,type_matrix,L,U)
        deallocate(matrix_a,L,U)
        write(*,*) 'Problem 1: Resolved !'
!END:1
      elseif (trim(problem) == 'Problem2') then  
!START 2
        write(*,*) 'What is the size of the matrix Wilkinson?'
        read(*,*) n_wilkin
        allocate(matrix_a(n_wilkin,n_wilkin))
        call wilkin(n_wilkin,matrix_a)
        deallocate(matrix_a)
        write(*,*) 'Problem 2: Resolved !'
!END:2      
      else 
        write(*,*) 'Ok, Goodbye !'
      endif
  else if (n_project .eq. 2) then 
    call intpol
  else 
    write(*,*) 'I apologize, but this project is not present in my code.'
  endif
      end program 
