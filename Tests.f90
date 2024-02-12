PROGRAM Tests
    USE mod_integrate

    IMPLICIT NONE

    COMPLEX*16, ALLOCATABLE :: X(:), Y(:)
    COMPLEX*16 :: dz, x_target, y_approx, dy
    COMPLEX*16 :: x_max, x_min, integral
    COMPLEX*16, EXTERNAL :: funky_square
    INTEGER*4 :: i,j, N, deg
    REAL*8 :: T_start, T_end

    N = 10
    deg = 8
    dz = DCMPLX(.5, .5)
    x_target = 10.2*dz
    ALLOCATE(X(N))
    ALLOCATE(Y(N))

    DO i = 1, N
        X(i) = dz*i
        Y(i) = funky_square(dz*i)
    END DO

    CALL POLINT(X, Y, N, x_target, y_approx, dy)
    PRINT*, "Approx value  = ", y_approx, "Error ~ ", dy
    PRINT*, "Exact Value = ", funky_square(x_target)
    PRINT*, "Relative error = ", ABS(funky_square(x_target) - y_approx)/ABS(funky_square(x_target))


    N = 8
    x_min = 0.
    x_max = 10.
    integral = 0.
    !Integrator testing
    !it generates 2**N + 1 points
    CALL CPU_TIME(T_start)
    DO i = 1, N + 1
        CALL TRAPZD(funky_square, x_min, x_max, integral, i)
    END DO
    CALL CPU_TIME(T_end)
    PRINT*, "Exec. time = ", T_end - T_start
    PRINT*, "Integration result = ", integral


    integral = 0.
    CALL CPU_TIME(T_start)
    CALL QROMB(funky_square, x_min, x_max, integral)
    CALL CPU_TIME(T_end)
    PRINT*, "Exec. time = ", T_end - T_start
    PRINT*, "QROMB result = ", integral


END PROGRAM Tests

COMPLEX*16 FUNCTION funky_square(z)
    COMPLEX*16 :: z
    funky_square = z**2 * EXP(DCMPLX(0., 1.)*ABS(z))
    RETURN
END FUNCTION funky_square