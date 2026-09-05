!    The module deals with the output routines.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

Module output_module
  use parameters_module
  use timestamp_module
  use, intrinsic :: iso_fortran_env
  implicit none
  
  contains
  
  subroutine print_counter_frequency ()
    if (mod(counter,frequency)==0) print*,"step ",counter
  end subroutine print_counter_frequency
  
  subroutine output ()
    call Report_time
    call Time_stamp ()
    call Report_compiler_properties ()
    !call Plot ()
    call Write_field (p,T)
  end subroutine output
  
  !Only relevant if dislin is installed!
  !subroutine plot ()
    ! to plot phi
    !call metafl ('pdf')
    !call disini
    !call qplclr (p,Nx,Ny)
    ! to plot temperature
    !call metafl ('pdf')
    !call disini        
    !call qplclr (T,Nx,Ny)
  !end subroutine plot
  
  !====
  
  subroutine Report_time ()
    time_computed = real(max(time_end-time_start,1_int64))/real(time_rate)
    print*, "================================"
    print*, " The code took", time_computed, "seconds"
  end subroutine Report_time
  
  !====
  
  subroutine Report_compiler_properties ()
    print*, "================================"
    print '(4a)', 'The code was compiled by ', &
    compiler_version(), ' with the options ', &
    compiler_options()
  end subroutine Report_compiler_properties
  
  !===
  
  subroutine write_field (phi,Tempr)
    real(dp),dimension(Nx,Ny),intent(in) :: phi,Tempr
    integer(sp)                           :: i,j 
    open(1, file = "phi.dat" )
    open(2, file = 'T.dat')
    do i=1,Nx
      write(1,10) ( phi(i,j),  j = 1,Ny )
      write(2,20) ( Tempr(i,j),j = 1,Ny )
    end do
    10 FORMAT(1000000F10.6)
    20 FORMAT(1000000F10.4)
    close (1)
    close (2)
  end subroutine write_field

end module output_module