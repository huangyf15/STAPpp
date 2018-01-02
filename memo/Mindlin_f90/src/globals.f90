module GLOBALS

   integer, parameter :: IELMNT=1	! Unit storing element data 
   integer, parameter :: ILOAD=2	! Unit storing load vectors
   integer, parameter :: IIN=5		! Unit used for input
   integer, parameter :: IOUT=6		! Unit used for output
    integer, parameter:: IPOST=7    ! Unit used for postprocess
   
   integer :: NUMNP		! Total number of nodal points
						! = 0 : Program stop
   integer :: NEQ		! Number of equations
   integer :: NWK		! Number of matrix elements
   integer :: MK		! Maximum half bandwidth

   integer :: IND		! Solution phase indicator
						!   1 - Read and generate element information
						!   2 - Assemble structure stiffness matrix
						!   3 - Stress calculations
   integer :: NPAR(10)	! Element group control data
						!   NPAR(1) - Element type
						!             1 : Truss element
						!   NPAR(2) - Number of elements
						!   NPAR(3) - Number of different sets of material and 
						!             cross-sectional  constants
   integer :: NUMEG		! Total number of element groups, > 0

   integer :: MODEX		! Solution mode: 0 - data check only;  1 -  execution                                   

   real :: TIM(5)		! Timing information
   character*80 :: HED	! Master heading information for use in labeling the output

   integer :: NFIRST
   integer :: NLAST
   integer :: MIDEST
   integer :: MAXEST

   integer :: NG
   
   real(8) :: GRIVATY(3)

   REAL(8)       :: GAUSS(4,4),GAUSSW(4,4)       ! Gaussian points and their weight values
   
   DATA GAUSS/0.D0,   0.D0,   0.D0,   0.D0, &
              -.5773502691896D0,     .5773502691896D0,   0.D0,   0.D0, &
              -.7745966692415D0,   0.D0,   .7745966692415D0,   0.D0, &
              -.8611363115941D0,   -.3399810435849D0,   .3399810435849D0,   .8611363115941D0 /
   
   DATA GAUSSW/2.D0,   0.D0,   0.D0,   0.D0, &
               1.D0,   1.D0,   0.D0,   0.D0, &
               .5555555555556D0,   .8888888888889D0,   .5555555555556D0,   0.D0, &
               .3478548451375D0,   .6521451548625D0,   .6521451548625D0,   .3478548451375D0 /
   
end module GLOBALS