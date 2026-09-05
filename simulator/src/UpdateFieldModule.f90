!    The module updates the phi and temperature fields.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module update_module
  use parameters_module
  implicit none
  
  contains
  
  subroutine Update_field(phi,phi_new,Tempr,Tempr_new)
    real(dp),dimension(Nx,Ny),intent(in)  :: phi_new,Tempr_new
    real(dp),dimension(Nx,Ny),intent(out) :: phi,Tempr
    
    Phi   = Phi_new
    Tempr = Tempr_new
    
  end subroutine Update_field

end module update_module