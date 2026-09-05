!    The module checks the allocation status of arrays. 
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module error_status_module
  use, intrinsic :: iso_fortran_env
  use parameters_module
  implicit none
  integer, parameter :: ALLOCATION_ERR = 1
  
  contains
  
  subroutine check_allocation_status (istat)
    integer(sp),intent(in) :: istat
    print*, " "
    print*, "      PROGRAM STARTS!"    
    print*, "================================"
    if (istat/=0) then
      write ( unit = error_unit, fmt = '( A, I12 )' ) &
      '### Allocation Error: can not allocate array of size ', Nx,Ny
    stop 1
    else 
      print*, "arrays are allocated!"
    end if
    print*, "================================"
  end subroutine check_allocation_status

end module error_status_module