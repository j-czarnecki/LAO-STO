MODULE mod_logger
IMPLICIT NONE
SAVE
INTEGER*4, PARAMETER :: LOGGER_UNIT = 66
INTEGER*4, PARAMETER :: MAX_LOG_LEN = 1000
CONTAINS

SUBROUTINE INIT_LOGGER()
    
    OPEN(unit = LOGGER_UNIT, FILE = "./log.log", FORM = "FORMATTED", ACTION = "WRITE")
    WRITE(LOGGER_UNIT, *) "==== START ===="
    !Print some info about simulation

END SUBROUTINE

SUBROUTINE CLOSE_LOGGER()

    WRITE(LOGGER_UNIT, *) "==== END ===="
    CLOSE(LOGGER_UNIT)

END SUBROUTINE

SUBROUTINE LOG_STRING_DEBUG(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(I0, a, I0, a, I0, a, I0, a, I0, a, I0, a, I0, a)', ADVANCE = 'NO') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' '
    WRITE(LOGGER_UNIT, '(a)') logMsg
END SUBROUTINE


SUBROUTINE LOG_STRING_INFO(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(I0, a, I0, a, I0, a, I0, a, I0, a, I0, a, I0, a)', ADVANCE = 'NO') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' '
    WRITE(LOGGER_UNIT, '(a)') "INFO: " // logMsg
END SUBROUTINE

SUBROUTINE LOG_STRING_ABNORMAL(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(I0, a, I0, a, I0, a, I0, a, I0, a, I0, a, I0, a)', ADVANCE = 'NO') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' '
    WRITE(LOGGER_UNIT, '(a)') "ABNORMAL: " // logMsg
END SUBROUTINE

SUBROUTINE LOG_STRING_ERROR(logMsg)
    CHARACTER(LEN=*), INTENT(IN) :: logMsg
    INTEGER*4 :: Values(8)

    CALL DATE_AND_TIME(VALUES = Values)
    !Add time printing
    WRITE(LOGGER_UNIT, '(I0, a, I0, a, I0, a, I0, a, I0, a, I0, a, I0, a)', ADVANCE = 'NO') Values(1), '-', Values(2), '-', Values(3), ' ', Values(5), ':', Values(6), ':', Values(7), '.', Values(8), ' '
    WRITE(LOGGER_UNIT, '(a)') "ERROR: " // logMsg
END SUBROUTINE

END MODULE