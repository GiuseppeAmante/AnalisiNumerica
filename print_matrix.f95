  subroutine print_matrix(val_unit,matrix,m,n,type_matrix)
    real(8), intent(in) :: matrix(m,n)
    integer, intent(in) :: val_unit, m, n
    integer :: i, j
    character(len=10) :: type_matrix
!    
100 format('***********************************************************',a)
101 format('                       M A T R I X                         ',a)
102 format('  number of rows        :',I5)
103 format('  number of colons      :',I5)
104 format('  type of matrix        :',a)
105 format('==============================',a)
106 format('                                                          ',a)
!
write(val_unit,100)
write(val_unit,101)
write(val_unit,106)  
write(val_unit,105)   
write(val_unit,102) m
write(val_unit,103) n
write(val_unit,104) type_matrix
write(val_unit,105)   
write(val_unit,106)  
! 
  do i = 1, m
     if (val_unit.eq.1)  write(val_unit, '(*(F14.7))') (matrix(i, j), j=1, n)
     if (val_unit.eq.2)  write(val_unit, '(*(F7.1))') (matrix(i, j), j=1, n)
  end do
write(val_unit,100)
!
write(val_unit,106)  
write(val_unit,106)  
  end subroutine print_matrix
