PROGRAM Tests
    USE mod_hamiltonians
    USE mod_parameters
    USE mod_writers
    USE mod_reader

    IMPLICIT NONE

    COMPLEX*16, ALLOCATABLE :: Hamiltonian(:,:)
    REAL*8, ALLOCATABLE :: Charge_dens(:)

    ALLOCATE(Hamiltonian(DIM, DIM))
    ALLOCATE(Charge_dens(DIM_POSITIVE_K))

    Hamiltonian = DCMPLX(0. , 0.)
    Charge_dens = 1.

    CALL GET_INPUT("./test_input.nml")

    CALL COMPUTE_HUBBARD(Hamiltonian, Charge_dens)
    CALL PRINT_HAMILTONIAN(Hamiltonian)



    DEALLOCATE(Hamiltonian)
    DEALLOCATE(Charge_dens)

END PROGRAM Tests