!    This module defines the precision and range for the program.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module precision_module
    use, intrinsic :: iso_fortran_env
    implicit none
    
    integer,parameter :: sp = int32
    integer,parameter :: dp = real64

end module precision_module