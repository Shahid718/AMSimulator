!    The module solves the driving force.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module drivingForce_module
  use parameters_module
  implicit none
  
  contains
  
    pure function drivingForce (Tempr,i,j)
      
      ! list of arguments
      real(dp),dimension(Nx,Ny),intent(in) :: Tempr
      real(dp)                             :: drivingForce
      integer(sp),intent(in)               :: i,j
      
      drivingForce = -H*(Tempr(i,j)-Tm)/Tm
      
    end function drivingForce
  
  !===
  pure function DH (phi,i,j)
    
    ! list of arguments
    real(dp),dimension(Nx,Ny),intent(in) :: phi
    real(dp)                             :: DH
    integer(sp),intent(in)               :: i,j
    
    DH = 30.0_dp*phi(i,j)**2*(phi(i,j) - 1.0_dp)**2
    
  end function DH

end module drivingForce_module