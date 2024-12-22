MODULE mod_parameters
IMPLICIT NONE
SAVE

INTEGER*4, PARAMETER :: N_NEIGHBOURS = 3
INTEGER*4, PARAMETER :: N_NEXT_NEIGHBOURS = 6
INTEGER*4, PARAMETER :: N_ALL_NEIGHBOURS = N_NEIGHBOURS + N_NEXT_NEIGHBOURS
COMPLEX*16, PARAMETER :: imag = DCMPLX(0. , 1.)

REAL*8, PARAMETER :: meV2au = 1./27211.
REAL*8, PARAMETER :: nm2au = 1./0.05292
REAL*8, PARAMETER :: T2au = 4.254382E-6

REAL*8, PARAMETER :: k_B = 8.617333262 * 1e-5 * 1e3 * meV2au
REAL*8, PARAMETER :: A_TILDE =  SQRT(2./3.)*0.3905 * nm2au !length

REAL*8, PARAMETER :: PI = 4*ATAN(1.0d0)

! REAL*8, PARAMETER :: K1_MAX = (2. * PI * 2./3.)/A_TILDE !Full Brillouin zone to integrate over
! REAL*8, PARAMETER :: K2_MAX = (2. * PI * 2./3.)/A_TILDE
REAL*8, PARAMETER :: K1_MAX = 1. !Full Brillouin zone to integrate over
REAL*8, PARAMETER :: K2_MAX = 1.
REAL*8, PARAMETER :: KX_MAX = 4.0d0*PI/(3.0d0*SQRT(3.0d0)) !This defines maximum kx of hexagon corresponding to the first Brillouin Zone
REAL*8, PARAMETER :: KY_MAX = 2.0d0*PI/3.0d0 !This defines maximum ky of hexagon corresponding to the first Brillouin Zone
REAL*8, PARAMETER :: JACOBIAN = 8*PI**2/(3.*SQRT(3.0d0))



END MODULE mod_parameters