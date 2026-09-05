!    The module introduces the secondary arms.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module noise_module
  use parameters_module
  implicit none
  
  contains
  
  subroutine Side_branching (phi,Tempr)
    real(dp),dimension(Nx,Ny),intent(in)  :: phi
    real(dp),dimension(Nx,Ny),intent(out) :: Tempr
    real(dp) :: random_no,noise,Tempr_noise = 1.0e-3 ! magnitude  
    integer(sp) :: i,j
    
    do j=1,Ny
      do i=1,Nx
        if (phi(i,j) .ge. 0.01 .and. phi(i,j) .le. 0.99) then 
          call random_number(random_no)
          noise = (2.0_dp*random_no - 1.0_dp) ! random no. btw -1 and +1
          Tempr(i,j) = Tempr(i,j)*(1.0_dp + Tempr_noise*noise) 
        end if
      end do
    end do
    
  end subroutine Side_branching

end module noise_module