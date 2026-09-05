!    The module solves the anisotropy.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module anisotropy_module
  use parameters_module
  implicit none

  contains

  pure subroutine Compute_gradients(phi,dphidx,dphidy,a,da,i,j)
    real(dp),dimension(Nx,Ny),intent(in)    :: phi
    real(dp),dimension(Nx,Ny),intent(inout) :: dphidx,dphidy,a,da
    integer(sp),intent(in) :: i,j
    real(dp) :: theta
    integer(sp) :: jp,jm,ip,im

    jp = j+1;  jm = j-1;  ip = i+1;  im = i-1
    if (im==0)      im = Nx
    if (ip==(Nx+1)) ip = 1
    if (jm==0)      jm = Ny
    if (jp==(Ny+1)) jp = 1

    dphidx(i,j) = (phi(ip,j)-phi(im,j))/(2.0_dp*dx)
    dphidy(i,j) = (phi(i,jp)-phi(i,jm))/(2.0_dp*dy)
    theta = atan2(dphidy(i,j),dphidx(i,j))

    a(i,j)  =  a0*(1.0_dp+alpha*cos(aniso*(theta-theta0)))
    da(i,j) = -a0*aniso*alpha*sin(aniso*(theta-theta0))

  end subroutine Compute_gradients


  pure subroutine Solve_anisotropy(dphidx,dphidy,a,da,i,j,anisotropy)
    real(dp),dimension(Nx,Ny),intent(in) :: dphidx,dphidy,a,da
    integer(sp),intent(in)  :: i,j
    real(dp),intent(out)    :: anisotropy
    real(dp) :: term1,term2
    integer(sp) :: jp,jm,ip,im

    jp = j+1;  jm = j-1;  ip = i+1;  im = i-1
    if (im==0)      im = Nx
    if (ip==(Nx+1)) ip = 1
    if (jm==0)      jm = Ny
    if (jp==(Ny+1)) jp = 1

    term1 =  (a(i,jp)*da(i,jp)*dphidx(i,jp) - a(i,jm)*da(i,jm)*dphidx(i,jm))/dy
    term2 = -(a(ip,j)*da(ip,j)*dphidy(ip,j) - a(im,j)*da(im,j)*dphidy(im,j))/dx

    anisotropy = term1 + term2

  end subroutine Solve_anisotropy

end module anisotropy_module