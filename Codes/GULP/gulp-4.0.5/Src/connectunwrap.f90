  subroutine connectunwrap
!
!  Uncorrects connect cell indices for atoms being wrapped
!  round by periodic boundary conditions. Aim is to try
!  to maintain continuity of images so that energy is
!  continuous. Strategy is to look for changes in fractionals
!  of more than half of a cell.
!
!   5/07 Created from connectwrap
!  11/11 Case when nconnectind = 0 trapped and excluded from unwrapping.
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
!  Copyright Curtin University 2011
!
!  Julian Gale, NRI, Curtin University, November 2011
!
  use current
  use molecule
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: ii
  integer(i4)                                  :: ixd
  integer(i4)                                  :: iyd
  integer(i4)                                  :: izd
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: kk
!
!  Loop over connections to find ones for this configuration
!
  do ic = 1,nconnect
    if (nconnectcfg(ic).eq.ncf.and.nconnectind(ic).gt.0) then
      i = n1connect(ic)
      j = n2connect(ic)
!
!  Find relative integer shift in vector
!
      ixd = ixshift(j) - ixshift(i)
      iyd = iyshift(j) - iyshift(i)
      izd = izshift(j) - izshift(i)
!
!  If any of the shifts are not zero then adjust index
!
      if ((abs(ixd)+abs(iyd)+abs(izd)).gt.0) then
        call mindtoijk(nconnectind(ic),ii,jj,kk)
        ii = ii - ixd
        jj = jj - iyd
        kk = kk - izd
        nconnectind(ic) = (ii+5) + 10*(jj+5) + 100*(kk+5)
      endif
    endif
  enddo
!
  return
  end
