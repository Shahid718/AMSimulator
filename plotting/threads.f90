!###########################################################
!  
!   Program :
!             Fortran code for plotting threads vs timing.
!             The program uses Dislin graphical library.
!
!   Programmer :
!                 Shahid Maqbool
!
!   Modified   :
!                    Jan 20, 2025
!
!   To compile and run with gfortran:
!                                    f90link -a grids
!
!###########################################################

program PlotThreads
  character*80 CBUF
  real,dimension(3) :: x1,y1,dc, x2,y2,omp
  
  !declare data
  x1 = [2.,4.,6.]
  y1 = [0.,0.,0.]
  dc = [25.77,15.0,4.02 ]  
  x2 = [1.,3.,5. ]
  y2 = [0.,0.,0. ]
  omp= [25.9,14.9,4.9 ] 
  
  ! dislin plotting routines
  call metafl('cons')
  call scrmod ( 'revers' )
  call page ( 2100,1500)
  call disini
  call hwfont
  call shdpat(16)
  call axslen(1500,1000)
  call axspos(400,1200)
  call barwth(0.4)
  call height (60)
  call ticks (0,'XY')
  call hname (60)
  call name ('Time (seconds)', 'Y')
  call labels ( "none" ,'X')
  call messag ("1", 700, 1250)
  call messag ("2", 1120, 1250)
  call messag ("4", 1550, 1250)
  call messag ("No. of cores", 800, 1350)
  call labels('second','bars')
  call labpos('outside','bars')
  call labclr(255,'bars')
  call graf(0.,7.,0.,1.,0.,45.,0.,5.)   ! x and y axis
  call color('red')
  call bars(x1,y1,dc,3)
  call shdpat(16)
  call color ('fore')
  call color ('blue')
  call bars(x2,y2,omp,3)
  call color ('fore')
  
  ! for legend
  call legini (cbuf,2,20)
  call legtit ('')
  call leglin (cbuf, 'Do Concurrent', 1)    
  call leglin (cbuf, 'OpenMP', 2)
  call legend (cbuf,7)
  
  ! finalize dislin
  call disfin
  
end program plotThreads