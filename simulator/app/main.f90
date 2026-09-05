program solidification_pure_Ni
  use precision_module
  use parameters_module
  use array_manage_module
  use error_status_module
  use initial_module
  use Anisotropy_Module
  use drivingForce_module
  use stencil_module
  use update_module
  use noise_module
  use output_module
  use timestamp_module
  implicit none
    
  call Allocate_arrays()
  
  ! phase field properties
  delta   = 4.0_dp*dx
  alpha   = 0.05_dp
  aniso   = 4.0_dp
  theta0  = 0.0_dp
  mu      = 2.0_dp
  lambda  = 0.1_dp
  b       = 2.0_dp*atanh(1.0_dp-2.0_dp*lambda)
  w       = 6.0_dp*sigma*b/delta
  a0      = sqrt(3.0_dp*delta*sigma/b)
  T0      = Tm-H/Cp*0.3_dp
  dt      = dx*dx/(5.0_dp*Tk/Cp)
  mobility= b*Tm*mu/(3.0_dp*delta*H)
  
  !===
  call Generate_initial_profiles (p,T)  
  print*," evolution ..."
  call system_clock (time_start,count_rate=time_rate)

  time_loop:do counter=1,iteration

    !--- gradient
#IFDEF DC
    dc_loop_gradient:do concurrent (j=1:Ny,i=1:Nx) default (none)    &
    shared (p,pdx,pdy,a,da)
#ELSE	
    !$omp parallel do collapse(2) default (none)                     &
    !$omp& private (i,j)                                             &
    !$omp& shared  (p,pdx,pdy,a,da)
    do j=1,Ny
      do i=1,Nx
#ENDIF	  
        call Compute_gradients (p,pdx,pdy,a,da,i,j)
#IFDEF DC
    end do dc_loop_gradient
#ELSE
      end do
    end do
    !$omp end parallel do
#ENDIF

    ! --- phase field & temperature	
#IFDEF DC
    dc_loop_fields:do concurrent (j=1:Ny,i=1:Nx) default (none)      &
    local (anisotropy)                                               &
	shared (p,pnew,T,Tnew,pdx,pdy,a,da,w,dt,mobility,Tk,Cp,H)
#ELSE	
    !$omp parallel do collapse(2) default (none)                     &
    !$omp& private (i,j,anisotropy)                                  &
    !$omp& shared  (p,pnew,T,Tnew,pdx,pdy,a,da,w,dt,mobility,Tk,Cp,H)
    do j=1,Ny
      do i=1,Nx
#ENDIF	  
        call Solve_anisotropy (pdx,pdy,a,da,i,j,anisotropy)

        ! phase field
        pnew(i,j) = p(i,j)+(a(i,j)**2*Laplace(p,i,j) + anisotropy    &
        + 4.0_dp*w*p(i,j)*(1.0_dp-p(i,j))*(p(i,j) - 0.5_dp           &
        + 15.0_dp/(2.0_dp*w) * DrivingForce(T,i,j)                   &
        * p(i,j)*(1.0_dp-p(i,j))))* dt*mobility

        ! temperature
        Tnew(i,j) = T(i,j) + (Tk/Cp)*Laplace(T,i,j)*dt               &
        + DH(p,i,j)*(H/Cp)*(pnew(i,j) - p(i,j) )

#IFDEF DC
    end do dc_loop_fields
#ELSE
      end do
    end do
    !$omp end parallel do
#ENDIF
    
    call Update_field (p,pnew,T,Tnew)
    call Side_branching (p,T)
    call Print_counter_frequency ()
  end do time_loop
  
  ! get the output
  call system_clock (time_end)
  call Output ()
  call Free_Memory ()
end program solidification_pure_Ni