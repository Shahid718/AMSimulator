!    The module solves the stencil operation.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module stencil_module
  use parameters_module
  implicit none
  
  contains
  
  pure function Laplace (field,i,j)
    
    ! list of arugments    
    real(dp),dimension(Nx,Ny),intent(in) :: field
    real(dp)                             :: Laplace
    integer(sp),intent(in)               :: i,j
    integer(sp)                          :: jp,jm,ip,im
    
    ! boundary conditions
    jp = j + 1
    jm = j - 1
    ip = i + 1
    im = i - 1
    
    if ( im==0 )        im = Nx
    if ( ip==(Nx+1) )   ip = 1
    if ( jm==0 )        jm = Ny
    if ( jp==(Ny+1) )   jp = 1 
    
    Laplace =  ( field(ip,j) + field(im,j) + field(i,jm) + &
    field(i,jp) - 4.0_dp*field(i,j)) / (dx*dy)
    
  end function LAPLACE
  
end module stencil_module