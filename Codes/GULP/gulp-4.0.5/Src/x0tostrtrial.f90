  subroutine x0tostrtrial(ntrialatom,nptrtrialatom)
!
!  Subroutine to convert linear structure array to main structure arrays for a subset of atoms
!
!   1/08 Created from x0tostr
!   6/09 Module name changed from three to m_three
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
!  Copyright Curtin University 2009
!
!  Julian Gale, NRI, Curtin University, June 2009
!
  use configurations
  use current
  use four
  use m_three
  use symmetry
  implicit none
!  
!  Passed variables
! 
  integer(i4), intent(in)                      :: ntrialatom
  integer(i4), intent(in)                      :: nptrtrialatom(ntrialatom)
!
!  Local variables
!
  integer(i4) :: i
  integer(i4) :: mvar
  integer(i4) :: nt
!
  if (ndim.eq.3) then
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    if (c.gt.1.0d-12) then
      recipc = 1.0_dp/c
    else
      recipc = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = rv(3,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = rv(3,2)
    r3x = rv(1,3)
    r3y = rv(2,3)
    r3z = rv(3,3)
    call rlist
  elseif (ndim.eq.2) then
    call uncell2D(rv,a,b,alpha)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    call rlist
  elseif (ndim.eq.1) then
    call uncell1D(rv,a)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    r1x = rv(1,1)
    call rlist
  endif
  mvar = 3*nasym + nstrains
!*************************************
!  Substitute parameters into place  *
!*************************************
!
!  Fractional coordinates
!
  do nt = 1,ntrialatom
    i = nptrtrialatom(nt)
    xafrac(i) = x0(3*i+nstrains-2)
    yafrac(i) = x0(3*i+nstrains-1)
    zafrac(i) = x0(3*i+nstrains)
  enddo
!
!  Radii
!
  if (nbsmat.gt.0) then
    do nt = 1,ntrialatom
      i = nptrtrialatom(nt)
      rada(i) = x0(mvar+i)
    enddo
  endif
!
!  Generate full coordinate set
!
  do nt = 1,ntrialatom
    i = nptrtrialatom(nt)
    xfrac(i) = xafrac(i)
    yfrac(i) = yafrac(i)
    zfrac(i) = zafrac(i)
  enddo
  if (nbsmat.gt.0) then
    do nt = 1,ntrialatom
      i = nptrtrialatom(nt)
      radf(i) = rada(i)
    enddo
  endif
!
!  Convert cell parameters and internal coordinates
!  into cartesian coordinates
!
  if (ndim.eq.3) then
    do nt = 1,ntrialatom
      i = nptrtrialatom(nt)
      xclat(i) = xfrac(i)*r1x + yfrac(i)*r2x + zfrac(i)*r3x
      yclat(i) = xfrac(i)*r1y + yfrac(i)*r2y + zfrac(i)*r3y
      zclat(i) = xfrac(i)*r1z + yfrac(i)*r2z + zfrac(i)*r3z
    enddo
  elseif (ndim.eq.2) then
    do nt = 1,ntrialatom
      i = nptrtrialatom(nt)
      xclat(i) = xfrac(i)*r1x + yfrac(i)*r2x
      yclat(i) = xfrac(i)*r1y + yfrac(i)*r2y
      zclat(i) = zfrac(i)
    enddo
  elseif (ndim.eq.1) then
    do nt = 1,ntrialatom
      i = nptrtrialatom(nt)
      xclat(i) = xfrac(i)*r1x
      yclat(i) = yfrac(i)
      zclat(i) = zfrac(i)
    enddo
  else
    do nt = 1,ntrialatom
      i = nptrtrialatom(nt)
      xclat(i) = xfrac(i)
      yclat(i) = yfrac(i)
      zclat(i) = zfrac(i)
    enddo
  endif
  do nt = 1,ntrialatom
    i = nptrtrialatom(nt)
    xalat(i) = xclat(i)
    yalat(i) = yclat(i)
    zalat(i) = zclat(i)
  enddo
!
  return
  end
