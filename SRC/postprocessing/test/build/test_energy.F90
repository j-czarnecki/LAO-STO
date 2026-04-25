!! This file is part of LAO-STO.
!!
!! Copyright (C) 2025 Julian Czarnecki
!!
!! Unit tests for energy module

!@test
subroutine test_placeholder()
  use funit
#line 10 "/home/czarnecki/LAO-STO/SRC/postprocessing/test/test_energy.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_energy.pf', &
 & 10) )
  if (anyExceptions()) return
#line 11 "/home/czarnecki/LAO-STO/SRC/postprocessing/test/test_energy.pf"
end subroutine test_placeholder

module Wraptest_energy
   use FUnit
   implicit none
   private

contains


end module Wraptest_energy

function test_energy_suite() result(suite)
   use FUnit
   use Wraptest_energy
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   external test_placeholder


   suite = TestSuite('test_energy_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_placeholder', &
      test_placeholder))
   call suite%addTest(t)


end function test_energy_suite

