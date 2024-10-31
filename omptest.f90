PROGRAM omptest
  USE omp_lib
  IMPLICIT NONE

  !OMP specific
  INTEGER*4 :: max_num_threads, used_threads


  max_num_threads = omp_get_max_threads()
  WRITE (*,*) "Max num threads", max_num_threads
  CALL omp_set_num_threads(max_num_threads)
  !$omp parallel
  PRINT*, "Hello from process", omp_get_thread_num()
  used_threads = omp_get_num_threads()
  CALL SLEEP(2)
  !$omp end parallel
  PRINT*, "Allocated threads", used_threads


END PROGRAM