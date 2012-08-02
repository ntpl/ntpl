  subroutine centre
!
!  Adjusts unit cell to handle centring operators
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2005
!
!  Julian Gale, NRI, Curtin University, July 2005
!
  use configurations
  use current
  implicit none
!
!  Local variables     
! 
  integer(i4)        :: i
  integer(i4)        :: j
  logical            :: ltransform
  real(dp)           :: cop(3,3)
!
!
  ltransform = .false.
!
!  Zero lattice vectors
!
  do i = 1,3
    do j = 1,3
      rv(j,i) = rvcfg(j,i,ncf)
      cop(j,i) = 0.0_dp
    enddo
  enddo
!
!  Create new lattice vectors based on reduced unit cell
!
  if (ncbl.eq.2) then
!
!  A-centring
!
    ltransform = .true.
    cop(1,1) = 1.0_dp
    cop(2,2) = 0.5_dp
    cop(2,3) = - 0.5_dp
    cop(3,2) = 0.5_dp
    cop(3,3) = 0.5_dp
  elseif (ncbl.eq.3) then
!
!  B-centring
!
    ltransform = .true.
    cop(1,1) = 0.5_dp
    cop(3,1) = - 0.5_dp
    cop(2,2) = 1.0_dp
    cop(1,3) = 0.5_dp
    cop(3,3) = 0.5_dp
  elseif (ncbl.eq.4) then
!
!  C-centring
!
    ltransform = .true.
    cop(1,1) = 0.5_dp
    cop(2,1) = 0.5_dp
    cop(1,2) = - 0.5_dp
    cop(2,2) = 0.5_dp
    cop(3,3) = 1.0_dp
  elseif (ncbl.eq.5) then
!
!  Face centring
!
    ltransform = .true.
    cop(2,1) = 0.5_dp
    cop(3,1) = 0.5_dp
    cop(1,2) = 0.5_dp
    cop(3,2) = 0.5_dp
    cop(1,3) = 0.5_dp
    cop(2,3) = 0.5_dp
  elseif (ncbl.eq.6) then
!
!  Body centring
!
    ltransform = .true.
    cop(1,1) = - 0.5_dp
    cop(2,1) = 0.5_dp
    cop(3,1) = 0.5_dp
    cop(1,2) = 0.5_dp
    cop(2,2) = - 0.5_dp
    cop(3,2) = 0.5_dp
    cop(1,3) = 0.5_dp
    cop(2,3) = 0.5_dp
    cop(3,3) = - 0.5_dp
!
!  Hexagonal/rhombohedral (in hexagonal setting)
!
  elseif (ncbl.eq.7) then
    ltransform = .true.
    cop(1,1) = 2.0_dp/3.0_dp
    cop(2,1) = 1.0_dp/3.0_dp
    cop(3,1) = 1.0_dp/3.0_dp
    cop(1,2) = - 1.0_dp/3.0_dp
    cop(2,2) = 1.0_dp/3.0_dp
    cop(3,2) = 1.0_dp/3.0_dp
    cop(1,3) = - 1.0_dp/3.0_dp
    cop(2,3) = - 2.0_dp/3.0_dp
    cop(3,3) = 1.0_dp/3.0_dp
  endif
!
!  Multiply matrices
!
  if (ltransform) then
    do i = 1,3
      do j = 1,3
        rvcfg(j,i,ncf) = cop(1,i)*rv(j,1) + cop(2,i)*rv(j,2) + cop(3,i)*rv(j,3)
      enddo
    enddo
  endif
!
  return
  end
