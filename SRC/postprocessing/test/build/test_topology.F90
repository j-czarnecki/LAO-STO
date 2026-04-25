!! This file is part of LAO-STO.
!!
!! Copyright (C) 2025 Julian Czarnecki
!!
!! Unit tests for topology module

!@test
subroutine test_placeholder()
  use funit
#line 10 "/home/czarnecki/LAO-STO/SRC/postprocessing/test/test_topology.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_topology.pf', &
 & 10) )
  if (anyExceptions()) return
#line 11 "/home/czarnecki/LAO-STO/SRC/postprocessing/test/test_topology.pf"
end subroutine test_placeholder

module Wraptest_topology
   use FUnit
   implicit none
   private

contains


end module Wraptest_topology

function test_topology_suite() result(suite)
   use FUnit
   use Wraptest_topology
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   external test_placeholder


   suite = TestSuite('test_topology_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_placeholder', &
      test_placeholder))
   call suite%addTest(t)


end function test_topology_suite

