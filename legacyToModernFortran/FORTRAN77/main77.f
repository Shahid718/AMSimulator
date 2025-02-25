C     declare data      
!      implicit double precision (A-H,O-Z)
      parameter ( Nx = 400, Ny = 400 )
      dimension T(Nx,Ny),Tnew(Nx,Ny),T_lap(Nx,Ny),p_lap(Nx,Ny),a(Nx,Ny)
      dimension p(Nx,Ny),pnew(Nx,Ny),pdx(Nx,Ny),pdy(Nx,Ny),da(Nx,Ny)      
C     ====
C     initialize data    
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
      mu      = 2
      T0      = Tm-H/Cp*0.3
      alambd  = 0.1
      b       = 2.0*atanh(1.0-2.0*alambd)
      a0      = sqrt(3.0*delta*sigma/b)
      w       = 6.0*sigma*b/delta
      pmob    = b*Tm*mu/(3.0*delta*H)
      dt      = dx*dx/(5.0*Tk/Cp)
      nsteps  = 6000      
      call cpu_time (t1)      
C     ====
C     initial microstructure    
      T = T0
      p = 0      
      do i = 1, Nx
        do j = 1, Ny
          x = (i - Nx/2.0)*(i - Nx/2.0)
          y = (j - Ny/2.0)*(j - Ny/2.0)
          dist = x + y
          if ( dist - 3 <= 0) then
            p(i,j) = 1.0
          end if
        end do
      end do      
C     ====
C     evolution      
      do istep = 1, nsteps
        do 10 i = 1, Nx
          do 20 j = 1, Ny
C           boundary conditions
            jp = j + 1
            jm = j - 1
            ip = i + 1
            im = i - 1
            if (im == 0) then
              im = Nx
            end if      
            if (ip == (Nx+1)) then
              ip=1
            end if      
            if (jm == 0) then
              jm = Ny
            end if
            if (jp == (Ny+1)) then
              jp = 1
            end if
C           laplacian
            p_lap(i,j)=(p(ip,j)+p(im,j)+p(i,jm)+p(i,jp)-
     &      4.0*p(i,j))/(dx*dy)
            T_lap(i,j)=(T(ip,j)+T(im,j)+T(i,jm)+T(i,jp)-
     &      4.0*T(i,j))/(dx*dy)      
C           gradients and angle
            pdx(i,j) = (p(ip,j)-p(im,j))/(2.0*dx)
            pdy(i,j) = (p(i,jp)-p(i,jm))/(2.0*dy)
            theta =  atan2(pdy(i,j),pdx(i,j))      
C           anisotropy and its derivative
            a(i,j) = a0*(1.0+alpha*cos(aniso*(theta-theta0)))
            da(i,j)= -a0*aniso*alpha*sin(aniso*(theta-theta0))
            term1=(a(i,jp)*da(i,jp)*pdx(i,jp)-
     &      a(i,jm)*da(i,jm)*pdx(i,jm))/dy
            term2=-(a(ip,j)*da(ip,j)*pdy(ip,j)-
     &      a(im,j)*da(im,j)*pdy(im,j))/dx      
C           driving force
            f = -H*(T(i,j)-Tm)/Tm
C           phi and temperature
            pnew(i,j) = p(i,j)+(a(i,j)**2*p_lap(i,j)+term1+term2
     &      + 4.0*w*p(i,j)*(1.0-p(i,j))*(p(i,j)-0.5+15.0/(2.0*w)*f*
     &      p(i,j)*(1.0-p(i,j))))*dt*pmob
            Tnew(i,j) = T(i,j)+(Tk/Cp)*T_lap(i,j)*dt+30.0*p(i,j)*p(i,j)
     &      * (1.0-p(i,j))*(1.0-p(i,j))*(H/Cp)*(pnew(i,j)-p(i,j) )      
20        continue
10      continue      
C       update fields
        do i = 1, Ny
          do j = 1, Nx
            T(i,j) = Tnew(i,j)
            p(i,j) = pnew(i,j)
          end do
        end do
      end do      
C     get the end time
      call cpu_time (t2)
      print*, " time", (t2 - t1), " seconds"
      
      open ( 1, file = "phi.dat" )        
        do i = 1, Nx
            write( 1, 100 ) ( p(i,j),j = 1, Ny )
        end do
100   FORMAT(1000000F10.6)
        
      close( 1 )

      end