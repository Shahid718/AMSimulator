!###########################################################
!  
!   Program :
!             Fortran code for plotting grids vs timing.
!             The program uses Dislin graphical library.
!
!   Programmer :
!                 Shahid Maqbool
!
!   Modified   :
!                    Jan 20, 2025
!
!   To compile and run with gfortran :
!                                     f90link -a grids
!
!###########################################################

program scaling_plot
  implicit none
  integer (4),parameter :: n=4
  real (4) , dimension (1:n) :: Nx,Time_dc,Time_omp
  
  ! data
  Nx          = [ 256 , 512  , 1024 , 2048   ]
  Time_dc     = [ 2.6 , 12.0 , 40.2 , 100.0  ]
  Time_omp    = [ 3.0 , 11.6 , 39.9 , 101.0  ]
  
  call plot ()
  
  contains
  
  subroutine plot( )
    use dislin
    implicit none
    
    character*80 :: cbuf
    integer (4)  :: IC
    
    ! initial setting of file
    call metafl ( 'cons' )
    call scrmod ( 'revers' )
    call disini
    
    ! for axis and graph setting 
    call hwfont
    call axspos ( 400, 1700 )
    call axslen ( 2400, 1400 )
    call labdig ( 1,'xy' )
    call ticks ( 0,'XY' )
    call height ( 60 )
    call hname ( 60 )
    
    ! for axis label
    call name ('Time (seconds)', 'Y' )
    call name ( 'No. of grids', 'X')
    call intax ()
    call graf ( 0.0, 2304.0, 0.0, 256.0, -40.0, 200.0,0.0 , 40.0 )
    call setrgb ( 0.8,0.8,0.8 )
    
    ! for curves
    call thkcrv ( 7 )
    call color ( 'fore' )
    call color ( 'red' )
    call curve (Nx,Time_dc,n)
    call color ( 'blue')
    call curve (Nx,Time_omp,n)
    call color ( 'fore' )
    
    ! for legend
    call legini (cbuf,2,20)
    call legtit ('')
    call leglin (cbuf, 'Do Concurrent', 1)    
    call leglin (cbuf, 'OpenMP', 2)
    call legend (cbuf,8)
    
    ! for data points on graph 
    call messag ("101.0", 2400, 770)
    call messag ("100.0", 2500, 930)
    call messag ("39.9", 1250, 1180)
    call messag ("40.2", 1400, 1280)
    call messag ("11.6", 770, 1300)
    call messag ("12.0", 900,1450)
    call messag ("3.0", 550, 1330)
    call messag ("2.6", 600, 1480)
    
    ! finalize dislin
    call disfin( )
    
  end subroutine plot
end program scaling_plot