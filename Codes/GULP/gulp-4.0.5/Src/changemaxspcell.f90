  subroutine changemaxspcell
!
!  Alters the size of the arrays associated with maxspcell 
!
!   6/07 Structure of arrays changed to make atom and cell
!        arrays 1-D and a pointer to where the cell data
!        starts added
!   9/10 Initialisations now performed in a subroutine
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
!  Copyright Curtin University 2010
!
!  Julian Gale, NRI, Curtin University, September 2010
!
  use spatial
  use reallocate
  implicit none
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxspcell = 0
!
!  Spatial decomposition data
!
  call realloc(nspcellat,maxspcell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxspcell','nspcellat')
  call realloc(nspcellat1ptr,maxspcell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxspcell','nspcellatptr')
!
!  Initialise defaults for new part of array
!
  if (maxspcell.gt.oldmaxspcell) then
    do i = oldmaxspcell+1,maxspcell
      call initmaxspcelldefaults(i)
    enddo
  endif
!
!  Save current value of maxspcell for next call
!
  oldmaxspcell = maxspcell
!
  return
  end
