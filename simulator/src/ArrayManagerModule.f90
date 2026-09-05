!    The module manages the memory.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module array_manage_module
  use Precision_Module
  use Parameters_Module
  use Error_status_Module
  implicit none
  
  contains
  subroutine allocate_arrays()
    implicit none
    ! allocate dynamic arrays
    allocate (p(Nx,Ny),    stat = istat)
    allocate (pnew(Nx,Ny), stat = istat)
    allocate (T(Nx,Ny),    stat = istat)
    allocate (Tnew(Nx,Ny), stat = istat)
    allocate (pdx(Nx,Ny),  stat = istat)
    allocate (pdy(Nx,Ny),  stat = istat)
    allocate (a(Nx,Ny),    stat = istat)
    allocate (da(Nx,Ny),   stat = istat)
    
    call Check_allocation_status (istat)
    
  end subroutine allocate_arrays
  
  !======
  
  subroutine free_memory ()
    deallocate (p,pnew,T,Tnew,pdx,pdy,a,da)
  end subroutine free_memory

end module array_manage_module