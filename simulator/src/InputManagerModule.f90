!    The module declare and the input parameters. 
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                October 5, 2024

module parameters_module
  use precision_module
  implicit none
  
  ! grid size in 2D    
  integer(sp),parameter :: Nx = 200      ! x direction
  integer(sp),parameter :: Ny = 200      ! y direction 
  real(dp) :: dx = 20e-9                 ! x direction
  real(dp) :: dy = 20e-9                 ! y direction
  
  ! time parameters    
  integer(sp) :: counter                  ! counter variable    
  integer(sp) :: iteration = 3000         ! no.of computing steps
  integer(sp) :: frequency = 1000         ! frequency of output steps
  integer(sp) :: time_start,time_end      ! computing time variables
  integer(sp) :: time_rate                ! computing time variables 
  real(dp)    :: time_computed            ! time taken for computation
  integer(sp) :: i,j
  
  ! error check
  integer(sp) :: istat
  
  ! metallurgical properties of Nickel
  real(dp) :: H     = 2.35e9       ! latent heat        
  real(dp) :: sigma = 0.37         ! interface energy
  real(dp) :: Tm    = 1728.0       ! Nickel melting temperature
  real(dp) :: Cp    = 5.42e6       ! Heat capactiy 
  real(dp) :: Tk    = 0.84e2       ! thermal conductivity
  
  ! phase field properties
  real(dp) :: delta        ! interface    
  real(dp) :: alpha        ! anisotropy of interface energy
  real(dp) :: aniso
  real(dp) :: theta0
  real(dp) :: mu
  real(dp) :: lambda     
  real(dp) :: b
  real(dp) :: w            ! height of double well potential           
  real(dp) :: a0           ! isotropic gradient energy coefficient                
  real(dp) :: T0           ! initial undercooling
  real(dp) :: dt           ! time interval
  real(dp) :: mobility     ! phase field mobility
  
  ! evolution variables
  real(dp) :: anisotropy
  
  ! 2D arrays of temperature and phi
  real(dp),dimension(:,:),allocatable :: T     ! temperature  
  real(dp),dimension(:,:),allocatable :: Tnew  ! temp. at new time step
  real(dp),dimension(:,:),allocatable :: p     ! phi
  real(dp),dimension(:,:),allocatable :: pnew  ! phi at new time step
  real(dp),dimension(:,:),allocatable :: pdx
  real(dp),dimension(:,:),allocatable :: pdy
  real(dp),dimension(:,:),allocatable :: a
  real(dp),dimension(:,:),allocatable :: da
  
end module parameters_module