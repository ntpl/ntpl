  subroutine changemaxnptstot
!
!  Alters the size of the arrays associated with maxnptstot
!
!   2/05 spxyzouter added
!  12/08 Migrated to version 3.5 and converted to f90 format
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
!  Copyright Curtin University 2008
!
!  Julian Gale, NRI, Curtin University, December 2008
!
  use cosmo
  use reallocate
  implicit none
!
  integer(i4) :: ierror
!
  call realloc(nar,maxnptstot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnptstot','nar')
  call realloc(nsetf,maxnptstot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnptstot','nsetf')
  call realloc(qonsas,maxnptstot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnptstot','qonsas')
  call realloc(spxyz,3_i4,maxnptstot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnptstot','spxyz')
  call realloc(spxyzouter,3_i4,maxnptstot,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnptstot','spxyzouter')
!
  return
  end
