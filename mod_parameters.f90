MODULE mod_parameters
IMPLICIT NONE 
SAVE 

INTEGER*4, PARAMETER :: SUBLATTICES = 2
INTEGER*4, PARAMETER :: ORBITALS = 3
INTEGER*4, PARAMETER :: TBA_DIM = ORBITALS*SUBLATTICES
INTEGER*4, PARAMETER :: DIM_POSITIVE_K = TBA_DIM*2  !Hamiltonian for positive k i.e half of the Nambu space, *2 due to spin
INTEGER*4, PARAMETER :: DIM = DIM_POSITIVE_K*2    !*2 to transform to Nambu Space. 

INTEGER*4, PARAMETER :: N_NEIGHBOURS = 3
COMPLEX*16, PARAMETER :: imag = DCMPLX(0. , 1.)

REAL*8, PARAMETER :: meV2au = 1./27211.
REAL*8, PARAMETER :: nm2au = 1./0.05292

REAL*8, PARAMETER :: k_B = 8.617333262 * 1e-5 * 1e3 * meV2au
REAL*8, PARAMETER :: A_TILDE =  SQRT(2./3.)*0.3905 * nm2au !length

REAL*8, PARAMETER :: PI = 4*ATAN(1.0d0)

REAL*8, PARAMETER :: K1_MAX = (2. * PI * 2./3.)/A_TILDE !Full Brillouin zone to integrate over
REAL*8, PARAMETER :: K2_MAX = (2. * PI * 2./3.)/A_TILDE



END MODULE mod_parameters 