MODULE mod_reader
USE mod_parameters
IMPLICIT NONE
SAVE

!Default values, overwritten in get_input

!Physical parameters
REAL*8 :: T = 0.
REAL*8 :: t_D = 0.
REAL*8 :: t_I = 0.
REAL*8 :: lambda_SOC = 0.
REAL*8 :: DELTA_TRI = 0.
REAL*8 :: v = 0.
REAL*8 :: V_pdp = 0.
REAL*8 :: V_pds = 0.
REAL*8 :: J_SC = 0.
REAL*8 :: J_SC_PRIME = 0.
REAL*8 :: U_HUB = 0.
REAL*8 :: V_HUB = 0.
REAL*8 :: E_Fermi = 0.

!Discretization
INTEGER*4 :: k1_steps = 0
INTEGER*4 :: k2_steps = 0



!Self-consistency
INTEGER*4 :: max_sc_iter = 0
REAL*8 :: sc_alpha = 0.
REAL*8 :: eps_convergence = 0.


!Derived
REAL*8 eta_p
REAL*8 :: dk1, dk2, domega


NAMELIST /physical_params/  &
& T,                        &
& t_D,                      &
& t_I,                      &
& lambda_SOC,               &
& DELTA_TRI,                &
& v,                        &
& V_pdp,                    &
& V_pds,                    &
& J_SC,                     &
& J_SC_PRIME,               &
& U_HUB,                    &
& V_HUB,                    &
& E_Fermi

NAMELIST /discretization/ &
& k1_steps,               &
& k2_steps                

NAMELIST /self_consistency/ &
& max_sc_iter,              &
& sc_alpha,                 &
& eps_convergence


CONTAINS
SUBROUTINE GET_INPUT(nmlfile)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: nmlfile

    OPEN(unit = 9, FILE=nmlfile, FORM = "FORMATTED", ACTION = "READ", STATUS="OLD")

    !TODO: WRITE BETTER CHECKS!!!!!!!!!!!!!!!!!!
    READ(9,NML=physical_params)
    !Check input data
    IF (T < 0) STOP "Temperature in kelvins must be >= 0!"               
    
    !Change to atomic units
    t_D = t_D * meV2au           
    t_I = t_I * meV2au
    lambda_SOC = lambda_SOC * meV2au
    DELTA_TRI = DELTA_TRI * meV2au 
    v = v * meV2au 
    V_pdp = V_pdp * meV2au   
    V_pds = V_pds * meV2au
    J_SC = J_SC * meV2au
    J_SC_PRIME = J_SC_PRIME * meV2au
    U_HUB = U_HUB * meV2au
    V_HUB = V_HUB * meV2au
    E_Fermi = E_Fermi * meV2au

    READ(9,NML=discretization)
    IF ((k1_steps .LE. 0) .OR. (k2_steps .LE. 0)) STOP "k_steps must be > 0"

    READ(9,NML=self_consistency)
    eps_convergence = eps_convergence * meV2au
    !Calculating derived values
    dk1 = K1_MAX / k1_steps
    dk2 = K2_MAX / k2_steps
    domega = ABS(dk1*dk2*SIN(PI/3.))/(2*PI)**2
    eta_p = v * SQRT(3.) / 3.905 * nm2au

END SUBROUTINE GET_INPUT

END MODULE mod_reader
