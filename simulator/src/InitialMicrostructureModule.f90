!    The module introduces the initial nuclei.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module Initial_Module
  use Precision_Module
  use Parameters_Module
  implicit none
  
  contains
  
  subroutine Generate_initial_profiles (phi,Tempr)
    real(dp),dimension(Nx,Ny),intent(inout) :: phi,Tempr
    integer(sp)                             :: i,j
    real(dp)                                :: seed=3.0
    real(dp)                                :: d1,d2,d3,d4,d5,d6
    
    Tempr=T0
    phi=0.0_dp
    
    ! ====
    !  Single grain: To make profile at the center
    ! ====
    
    do j=1,Ny
      do i=1,Nx
        if ((i-0.5_dp*Nx)*(i-0.5_dp*Nx)+(j-0.5_dp*Ny)*(j-0.5_dp*Ny)&
        -seed <= 0.0_dp) then
          phi(i,j) = 1.0_dp       
        end if
      end do
    end do
    
  end subroutine Generate_initial_profiles
  
end module initial_module