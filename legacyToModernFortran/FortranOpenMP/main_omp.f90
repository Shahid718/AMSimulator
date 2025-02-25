program solidification_pure_Ni_omp
  use, intrinsic :: iso_fortran_env  ! Fortran 2008 module for precision
  use omp_lib                        ! OpenMP library 
  implicit none
  
  ! declare data
  integer, parameter::dp=real64,sp=int32 ! Fortran 2008 for precision
  integer(sp),parameter::Nx=100,Ny=100
  integer(sp)::nsteps,istep,i,j,ip,im,jp,jm,t1,t2,rate
  real(dp),dimension(Nx,Ny)::p,pnew,T,Tnew,pdx,pdy,a,da,T_lap,p_lap 
  real(dp)::dx,dy,delta,sigma,Tm,alpha,aniso,theta0,Tk,Cp,H,mu,T0,   &
  alambd,b,a0,w,pmob,term1,term2,f,theta,dt,x,y,dist,t3
  
  ! ====
  ! initialize data
  dx      = 20e-9
  dy      = 20e-9
  delta   = 4.0*dx
  sigma   = 0.37
  alpha   = 0.05
  aniso   = 4.0
  theta0  = 0.0
  Tm      = 1728.0
  Tk      = 84.01
  Cp      = 5.42e+06
  H       = 2.35e+09
  mu      = 2.0
  T0      = Tm-H/Cp*0.3
  alambd  = 0.1
  b       = 2.0*atanh(1.0-2.0*alambd)
  a0      = sqrt(3.0*delta*sigma/b)
  w       = 6.0*sigma*b/delta
  pmob    = b*Tm*mu/(3.0*delta*H)
  dt      = dx*dx/(5.0*Tk/Cp)
  nsteps  = 6000
  
  call system_clock (count=t1, count_rate=rate)  ! for multithread timing
  !====
  !initial microstructure
  T = T0
  p = 0
  do i = 1,Nx
    do j = 1,Ny
      x = (i-Nx/2.0)*(i-Nx/2.0)
      y = (j-Ny/2.0)*(j-Ny/2.0)
      dist = x+y
      if (dist-3 <= 0) then
        p(i,j) = 1.0
      end if
    end do
  end do
  !====
  !evolution      
  do istep = 1,nsteps
    !!$omp parallel do default (none)  & 
    !!$omp& private (...)              &
    !!$omp& shared (...)
    do i=1,Nx
      do j=1,Ny
        ! boundary conditions
        jp = j+1
        jm = j-1
        ip = i+1
        im = i-1
        if (im==0)        im=Nx
        if (ip==(Nx+1))   ip=1
        if (jm==0)        jm=Ny
        if (jp==(Ny+1))   jp=1        
        ! laplacian
        p_lap(i,j)=(p(ip,j)+p(im,j)+p(i,jm)+p(i,jp)-4.0*p(i,j))/(dx*dy)
        T_lap(i,j)=(T(ip,j)+T(im,j)+T(i,jm)+T(i,jp)-4.0*T(i,j))/(dx*dy)        
        ! gradients and angle
        pdx(i,j) = (p(ip,j)-p(im,j))/(2.0*dx)
        pdy(i,j) = (p(i,jp)-p(i,jm))/(2.0*dy)
        theta =  atan2(pdy(i,j),pdx(i,j))        
        ! anisotropy and its derivative
        a(i,j) = a0*(1.0+alpha*cos(aniso*(theta-theta0)))
        da(i,j)= -a0*aniso*alpha*sin(aniso*(theta-theta0))
        term1=(a(i,jp)*da(i,jp)*pdx(i,jp)-a(i,jm)*da(i,jm)*pdx(i,jm))/dy
        term2=-(a(ip,j)*da(ip,j)*pdy(ip,j)-a(im,j)*da(im,j)*pdy(im,j))/dx        
        ! driving force
        f = -H*(T(i,j)-Tm)/Tm       
        ! phi and temperature
        pnew(i,j) = p(i,j)+(a(i,j)**2*p_lap(i,j)+term1+term2         &
        & + 4.0*w*p(i,j)*(1.0-p(i,j))*(p(i,j)-0.5+15.0/(2.0*w)*f*    &
        & p(i,j)*(1.0-p(i,j))))*dt*pmob
        Tnew(i,j) = T(i,j)+(Tk/Cp)*T_lap(i,j)*dt+30.0*p(i,j)*p(i,j)* &
        & (1.0-p(i,j))*(1.0-p(i,j))*(H/Cp)*(pnew(i,j)-p(i,j) )
      end do
    end do
    !update fields using whole array
    T = Tnew
    p = pnew
  end do
  !get the end time
  call system_clock (count=t2)
  t3 = real(max(t2 - t1, 1_sp )) /real(rate)
  print*, " time", t3, " seconds"
  !call qplclr (p,Nx,Ny)
end program solidification_pure_Ni_omp