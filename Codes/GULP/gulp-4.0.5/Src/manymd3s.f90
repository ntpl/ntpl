  subroutine manymd3s(emany,esregion12,esregion2,eattach,lgrad1)
!
!  Subroutine for calculating the many-body energy from the Sutton-Chen potential(s). 
!  Uses a spatial decomposition algorithm.
!  Requires real space routine to be called first to evaluate the repulsive contribution
!  and the rho parameters.
!
!  On entry the array scrho must contain the density at each atomic site.
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   7/03 Created from manymd3
!  10/03 Modifications for new spatial algorithm added
!  11/03 ndennat/ndentyp replaced
!  11/03 Alloy scaling added
!   7/05 Constant scaling added to sqrt and power functionals
!   7/05 Use of EAM species pointers introduced
!   9/05 rhoderv called to get density derivatives
!   9/05 Call to eamfnderv used to replace EAM function derivative code
!  10/05 Call to eamfnderv added for energy
!   4/06 Modified to handle species specific densities
!   7/06 Clean up for small amount of extra speed
!   3/07 Printing of EAM densities and energies added as an option
!   3/07 Calculation of emany parallelised
!   5/07 QM/MM schemes added
!   5/07 Call to rhoderv modified by adding rpot
!   6/07 Structure of arrays for storing spatial distribution changed to 1-D
!  12/07 Unused variables removed
!   4/08 Modified for variable domain size
!   4/08 xvec1cell replaced by xvec2cell etc for spatial algorithm
!  11/08 Call to rhoderv updated to include the density
!  11/08 Call to rhoderv modified to include x,y,z Cartesian components
!  11/08 rho arrays changed to 2-D to benefit MEAM
!  11/08 call to rhoderv replaced by calls to meamrho/eamrho according to MEAM vs EAM
!  12/08 rho switched back to 1-D array with condensed components
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!  12/08 MEAM calculation of density added as option
!   1/09 Total density now printed for MEAM case
!   1/09 Derivatives modified to accommodate MEAM : deriv -> deriv(3)
!   2/09 Derivatives for MEAM added
!   2/09 Expanding and contracting of maxdis arrays removed
!   3/09 Skipping when density is zero modified to ensure no terms are missed
!   3/09 Sign of return arguments from EAM/MEAM function routines corrected for
!   4/09 MEAM screening function derivatives added
!   5/09 MEAM third derivative arguments removed
!   6/09 Site energy and virial added
!  11/09 Region derivatives added
!   3/10 Extra virial contributions added
!   8/10 lgrad1 in calls to meamscreen set to true
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
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
!  Copyright Curtin University 2012
!
!  Julian Gale, NRI, Curtin University, May 2012
!
  use configurations, only : nregions, nregionno, nregiontype, QMMMmode
  use control
  use current
  use derivatives
  use eam
  use energies,       only : siteenergy
  use general
  use iochannels,     only : ioout
  use mdlogic
  use numbers,        only : third
  use optimisation
  use parallel
  use spatial
  use sutton
  use symmetry
  use times
  use two
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  real(dp),    intent(inout)                   :: eattach
  real(dp),    intent(out)                     :: emany
  real(dp),    intent(inout)                   :: esregion12
  real(dp),    intent(inout)                   :: esregion2
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ic
  integer(i4)                                  :: ii
  integer(i4)                                  :: imx
  integer(i4)                                  :: imy
  integer(i4)                                  :: imz
  integer(i4)                                  :: ind
  integer(i4)                                  :: ind2
  integer(i4)                                  :: indn
  integer(i4)                                  :: ixyz
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jc
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmx
  integer(i4)                                  :: jmy
  integer(i4)                                  :: jmz
  integer(i4)                                  :: jndn
  integer(i4)                                  :: k
  integer(i4)                                  :: kc
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: klp
  integer(i4)                                  :: m
  integer(i4)                                  :: maxx
  integer(i4)                                  :: maxxy
  integer(i4)                                  :: n
  integer(i4)                                  :: n1i
  integer(i4)                                  :: n1j
  integer(i4)                                  :: n1k
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: np
  integer(i4)                                  :: npartial
  integer(i4)                                  :: npot
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionk
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: nsplower(3)
  integer(i4)                                  :: nspupper(3)
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: status
  logical                                      :: lnonzeroSij
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopk
  logical                                      :: lpartial
  logical                                      :: lQMMMok
  logical                                      :: lreg2ok
  logical                                      :: lsg1
  logical                                      :: lvalidij
  real(dp)                                     :: cputime
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: deltarho(maxmeamcomponent)
  real(dp)                                     :: deriv(3)
  real(dp)                                     :: derivs(6)
  real(dp)                                     :: drhoij(3,maxmeamcomponent)
  real(dp)                                     :: drhoijs(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoij2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoji(3,maxmeamcomponent)
  real(dp)                                     :: drhojis(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoji2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhototij(3)
  real(dp)                                     :: drhototijs(6)
  real(dp)                                     :: drhototij2(6)
  real(dp)                                     :: drhototij2s(21)
  real(dp)                                     :: drhototij2m(6,3)
  real(dp)                                     :: drhototij3(10)
  real(dp)                                     :: drhototji(3)
  real(dp)                                     :: drhototjis(6)
  real(dp)                                     :: drhototji2(6)
  real(dp)                                     :: drhototji2s(21)
  real(dp)                                     :: drhototji2m(6,3)
  real(dp)                                     :: drhototji3(10)
  real(dp)                                     :: eeam
  real(dp)                                     :: emanytrm
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: r2ijmid
  real(dp)                                     :: r2ik
  real(dp)                                     :: r2jk
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: rhoi
  real(dp)                                     :: rhoj
  real(dp)                                     :: rhoij(maxmeamcomponent)
  real(dp)                                     :: rhoji(maxmeamcomponent)
  real(dp)                                     :: rhorr
  real(dp)                                     :: rp
  real(dp)                                     :: rscrhoi
  real(dp)                                     :: rscrhoi3
  real(dp)                                     :: rscrhoi5
  real(dp)                                     :: rscrhoj
  real(dp)                                     :: rscrhoj3
  real(dp)                                     :: rscrhoj5
  real(dp)                                     :: scmax
  real(dp)                                     :: Sij
  real(dp)                                     :: Sikj
  real(dp)                                     :: dSikjdr(3)
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xij0
  real(dp)                                     :: yij0
  real(dp)                                     :: zij0
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
  real(dp)                                     :: xki
  real(dp)                                     :: yki
  real(dp)                                     :: zki
  real(dp)                                     :: xkj
  real(dp)                                     :: ykj
  real(dp)                                     :: zkj
  type(screening_atoms)                        :: partial
!
  time1 = cputime()
!
!  For screened MEAM, set scale factor that determines searching cutoff based on which axis of the ellipse is largest.
!  Note: the cutoff is applied to the mid point of the i-j vector and so this is guaranteed to find all distances.
!
  if (lMEAMscreen) then
    rcutfactor = max(1.0_dp,0.25_dp*(1.0_dp + meam_Cmax))
  endif
!
!  Scale density
!
  call eamscalescrho(1_i4)
  if (lPrintEAM) then
    call mpbarrier
    if (ioproc) then
!
!  Openning banner for energy decomposition
!
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  EAM : Atom No.                Density                 Atom energy (eV) '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Energy calculation
!
  emany = 0.0_dp
  lreg2ok = (lseok.and.nregions(ncf).gt.1)
  do i = procid+1,numat,nprocs
    neamspeci = neamfnspecptr(i)
    rhoi = scrho(1,i)
    nregioni = nregionno(nsft+nrelat(i))
    nregiontypi = nregiontype(nregioni,ncf)
!
!  QM/MM handling : i is a QM atom => exclude
!
    lQMMMok = .true.
    if (QMMMmode(ncf).gt.0) then
      if (nregiontypi.eq.1) lQMMMok = .false.
    endif
    if (neamspeci.gt.0.and.rhoi.gt.1.0d-12.and.lQMMMok) then
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      else
        call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      endif
      emanytrm = occuf(i)*eeam
      if (lreg2ok.and.nregionno(nsft+nrelat(i)).gt.1) then
        esregion2 = esregion2 + emanytrm
      else
        emany = emany + emanytrm
      endif
      siteenergy(i) = siteenergy(i) + emanytrm
      if (lPrintEAM) then                      
        write(ioout,'(7x,i8,4x,f24.10,4x,f24.12)') i,rhoi,emanytrm
      endif
      eattach = eattach + emanytrm
      if (lMEAMfn) then
        deltarho(1:maxmeamcomponent) = scrho(1:maxmeamcomponent,i) - scrho12(1:maxmeamcomponent,i)
        call meamfnderv(neamfn,neamspeci,deltarho,rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      else
        rhorr = rhoi - scrho12(1,i)
        call eamfnderv(neamfn,neamspeci,rhorr,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      endif
      emanytrm = occuf(i)*eeam
      eattach = eattach - emanytrm
    endif
  enddo
  if (lPrintEAM) then
    call mpbarrier
    if (ioproc) then
!     
!  Closing banner for energy decomposition
!
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  If no forces are needed then we don't need to do loops over atoms, so just return
!
  if (.not.lgrad1) goto 1000
!
!  From here on we can assume that lgrad1 = .true.
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('manymd3s','npotl')
!
!  Initialise local variables
!
  lsg1 = (lstr.or.lmd)
!***************************
!  Set up local variables  *
!***************************
!
!  Find maximum cut-off radius
!
  scmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.19) then
      if (rpot(i).gt.scmax) scmax = rpot(i)
    endif
  enddo
!
!  Skip if real space terms are not to be done
!
  if (lnoreal) goto 999
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
!
!  Set up spatial decomposition terms
!
  maxxy = nspcell(1)*nspcell(2)
  maxx  = nspcell(1)
!        
!  Loop over all local spatial cells except buffer regions
!
  do ixyz = 1,ncellpernode
    ind = ncellnodeptr(ixyz)
    ind2 = ind - 1
    iz = ind2/maxxy
    ind2 = ind2 - maxxy*iz
    iy = ind2/maxx
    ix = ind2 - maxx*iy + 1
    iy = iy + 1
    iz = iz + 1
    if (.not.lbuffercell(ixyz)) then
!
!  Set cell search bounds
!  
      nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
      nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
      nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
      nsplower(1) = max(ix-ncellsearch(1),1)
      nsplower(2) = max(iy-ncellsearch(2),1)
      nsplower(3) = max(iz-ncellsearch(3),1)
!
!  Get number of atoms in this cell
!
      ni = nspcellat(ind)
      n1i = nspcellat1ptr(ind)
!
!  Outer loop over atoms within this cell
!
      iloop: do ii = 1,ni      
        i = nspcellatptr(n1i+ii)
        ic = nspcellatptrcell(n1i+ii)
!
!  Set coordinates of atom i
!        
        xi = xinbox(i) + xvec2cell(ic)
        yi = yinbox(i) + yvec2cell(ic)
        zi = zinbox(i) + zvec2cell(ic)
!
!  Set other properties of atom i
!
        nati = nat(i)       
        ntypi = nftype(i)
        nregioni = nregionno(nsft+nrelat(i))
        nregiontypi = nregiontype(nregioni,ncf)
        oci = occuf(i)
        lopi = (.not.lfreeze.or.lopf(nrelat(i)))
!     
!  Find EAM species for i
!  
        neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
        if (lMEAMfn) then
          call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.false.,.false.)
        else
          call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.false.,.false.)
          rhoi = scrho(1,i)
        endif
!
!  Loop over neighbouring cells
!
        do imz = nsplower(3),nspupper(3)
          do imy = nsplower(2),nspupper(2)
            do imx = nsplower(1),nspupper(1)
              indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!
!  Loop over atoms within neighbouring cells      
!
              nj = nspcellat(indn)
              n1j = nspcellat1ptr(indn)
              jloop: do jj = 1,nj
                j = nspcellatptr(n1j+jj)
!
!  Perform quick skip tests
!
                if (j.ge.i) cycle jloop
!
                lopj = (.not.lfreeze.or.lopf(nrelat(j)))
                if (.not.lopi.and..not.lopj) cycle jloop
!
                jc = nspcellatptrcell(n1j+jj)
!     
!  Find EAM species for j
!  
                neamspecj = neamfnspecptr(j)
!
!  Evaluate functional derivatives
!
                if (lMEAMfn) then
                  call meamfnderv(neamfn,neamspecj,scrho(1,j),rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.false.,.false.)
                else
                  call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.false.,.false.)
                  rhoj = scrho(1,j)
                endif
!
                if (rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp) cycle jloop
!
                natj = nat(j)
                ntypj = nftype(j)
                nregionj = nregionno(nsft+nrelat(j))
                nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are both QM atoms => no forces to compute
!
                if (QMMMmode(ncf).gt.0) then
                  if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
                endif
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
                if (nati.eq.natj) then
                  nat1 = nati
                  nat2 = natj
                  if (ntypi.lt.ntypj) then
                    ntyp1 = ntypi
                    ntyp2 = ntypj
                  else
                    ntyp1 = ntypj
                    ntyp2 = ntypi
                  endif
                elseif (nati.lt.natj) then
                  nat1 = nati
                  nat2 = nat(j)
                  ntyp1 = ntypi
                  ntyp2 = nftype(j)
                else
                  nat1 = nat(j)
                  nat2 = nati
                  ntyp1 = nftype(j)
                  ntyp2 = ntypi
                endif
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
                npots = 0
                rp = 0.0_dp
                do n = 1,npote
                  if (nptype(n).eq.19) then
                    if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
                      if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
                        npots = npots + 1
                        npotl(npots) = n
                        if (rpot(n).gt.rp) rp = rpot(n)
                      endif
                    endif
                  endif
                enddo
!
!  If no valid potentials then there is no need to continue with this pair, unless this is
!  a second derivative calculation, in which case there may be a contribution from triangles
!  of interactions.
!
                if (npots.eq.0) cycle jloop
!
                xji = xvec2cell(jc) + xinbox(j) - xi 
                yji = yvec2cell(jc) + yinbox(j) - yi
                zji = zvec2cell(jc) + zinbox(j) - zi
                ocj = occuf(j)
                ofct = oci*ocj
!
                cut2r = rp*rp
                if (cut2r.gt.cut2p) cut2r = cut2p
                cut2 = cut2r
                rp = sqrt(cut2)
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
                deriv(1:3) = 0.0_dp
                if (lsg1) derivs(1:6) = 0.0_dp
                r2 = xji*xji + yji*yji + zji*zji
                r = sqrt(r2)
!***************************************
!  Valid many-body potentials for i-j  *
!***************************************
                if (lMEAM) then
                  rhoij(1:maxmeamcomponent) = 0.0_dp
                  rhoji(1:maxmeamcomponent) = 0.0_dp
                  drhoij(1:3,1:maxmeamcomponent) = 0.0_dp
                  drhoji(1:3,1:maxmeamcomponent) = 0.0_dp
                  if (lsg1) then
                    drhoijs(1:6,1:maxmeamcomponent) = 0.0_dp
                    drhojis(1:6,1:maxmeamcomponent) = 0.0_dp
                  endif
                else
                  rhoij(1) = 0.0_dp
                  rhoji(1) = 0.0_dp
                endif
                drhototij(1:3) = 0.0_dp
                drhototji(1:3) = 0.0_dp
                if (lsg1) then
                  drhototijs(1:6) = 0.0_dp
                  drhototjis(1:6) = 0.0_dp
                endif
                lvalidij = .false.
                if (npots.gt.0) then
                  do m = 1,npots
                    npot = npotl(m)
                    if (r.gt.rpot2(npot).and.r.le.rpot(npot).and.r.le.rpmax) then
                      lvalidij = .true.
!**********************************
!  Calculate density derivatives  *
!**********************************
                      if (lMEAMden) then
                        call meamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xji,yji,zji,rhoij,rhoji,drhoij,drhoji, &
                                     drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                                     1.0_dp,1.0_dp,.true.,lsg1,.true.,.false.)
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
                        lnonzeroSij = .true.
                        Sij = 1.0_dp
!
                        if (lMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
                          rcut2 = rcutfactor*r2
!
!  Initialise counter to number of partial screening atoms
!
                          npartial = 0
!
!  Loop over neighbouring cells to search for atoms that may contribute to the screening
!
                          jmzloop: do jmz = nsplower(3),nspupper(3)
                            do jmy = nsplower(2),nspupper(2)
                              do jmx = nsplower(1),nspupper(1)
                                jndn = (jmz-1)*maxxy + (jmy-1)*maxx + jmx
!
!  Loop over atoms within neighbouring cells
!
                                nk = nspcellat(jndn)
                                n1k = nspcellat1ptr(jndn)
                                kloop: do kk = 1,nk
                                  k = nspcellatptr(n1k+kk)
                                  kc = nspcellatptrcell(n1k+kk)
!
!  Set vectors between atoms
!
                                  xki = xvec2cell(kc) + xinbox(k) - xi
                                  yki = yvec2cell(kc) + yinbox(k) - yi
                                  zki = zvec2cell(kc) + zinbox(k) - zi
                                  xkj = xki - xji
                                  ykj = yki - yji
                                  zkj = zki - zji
                                  xij0 = 0.5_dp*(xki + xkj)
                                  yij0 = 0.5_dp*(yki + ykj)
                                  zij0 = 0.5_dp*(zki + zkj)
!
!  Compute square of distance to i-j mid point
!
                                  r2ijmid = xij0*xij0 + yij0*yij0 + zij0*zij0
                                  if (r2ijmid.lt.rcut2) then
!
!  Complete distances
!
                                    r2ik = xki*xki + yki*yki + zki*zki
                                    r2jk = xkj*xkj + ykj*ykj + zkj*zkj
!
!  Compute screening function
!
                                    call meamscreen(r2,r2ik,r2jk,Sikj,dSikjdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                                    if (Sikj.eq.0.0_dp) then
                                      lnonzeroSij = .false.
                                      Sij = 0.0_dp
                                    else
!
!  Multiply total screening product
!
                                      Sij = Sij*Sikj
                                      if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                                        npartial = npartial + 1
                                        if (npartial.gt.partial%sa_maxdim) then
                                          call changemaxsa(partial,npartial)
                                        endif
                                        partial%sa_atom(npartial) = k
                                        partial%sa_rij(npartial) = sqrt(r2)
                                        partial%sa_rik(npartial) = sqrt(r2ik)
                                        partial%sa_rjk(npartial) = sqrt(r2jk)
                                        partial%sa_xik(npartial) = xki
                                        partial%sa_yik(npartial) = yki
                                        partial%sa_zik(npartial) = zki
                                        partial%sa_xjk(npartial) = xkj
                                        partial%sa_yjk(npartial) = ykj
                                        partial%sa_zjk(npartial) = zkj
                                        partial%sa_Sikj(npartial) = Sikj
                                        partial%sa_dSikjdr(1:3,npartial) = dSikjdr(1:3)
                                      endif
                                    endif
!
!  If screening factor has gone to zero then we can exit all loops over k
!
                                    if (.not.lnonzeroSij) exit jmzloop
                                  endif
                                enddo kloop
!
!  End of loops over cells for k
!
                              enddo
                            enddo
                          enddo jmzloop
                        endif
!
!  Only do remainder of work if the screening factor is non-zero
!
                        if (lnonzeroSij) then
                          if (lMEAMscreen) then
                            if (npartial.gt.0) then
!
!  Compute derivative contributions of the screening function
!
                              call meamtotalrhoscreenderv(neamspeci,npartial,partial,xji,yji,zji,scrho(1,i), &
                                                          rhoi,rscrhoi,rhoij,Sij,lsg1)
!
                              do np = 1,npartial
                                k = partial%sa_atom(np)
                                lopk = (.not.lfreeze.or.lopf(nrelat(k)))
                                nregionk = nregionno(nsft+nrelat(k))
!
!  i-j contribution
!
                                if (lopi) then
                                  xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np)*ofct
                                  ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np)*ofct
                                  zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np)*ofct
                                endif
                                if (lopj) then
                                  xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np)*ofct
                                  ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np)*ofct
                                  zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np)*ofct
                                endif
                                if (nregioni.ne.nregionj) then
                                  xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np)*ofct
                                  yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np)*ofct
                                  zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np)*ofct
                                  xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np)*ofct
                                  yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np)*ofct
                                  zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np)*ofct
                                endif
!
!  i-k contribution
!
                                if (lopi) then
                                  xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np)*ofct
                                  ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np)*ofct
                                  zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np)*ofct
                                endif
                                if (lopk) then
                                  xdrv(k) = xdrv(k) + partial%sa_drhototik(1,np)*ofct
                                  ydrv(k) = ydrv(k) + partial%sa_drhototik(2,np)*ofct
                                  zdrv(k) = zdrv(k) + partial%sa_drhototik(3,np)*ofct
                                endif
                                if (nregioni.ne.nregionk) then
                                  xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np)*ofct
                                  yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np)*ofct
                                  zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np)*ofct
                                  xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototik(1,np)*ofct
                                  yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototik(2,np)*ofct
                                  zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototik(3,np)*ofct
                                endif
!
!  j-k contribution
!
                                if (lopj) then
                                  xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np)*ofct
                                  ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np)*ofct
                                  zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np)*ofct
                                endif
                                if (lopk) then
                                  xdrv(k) = xdrv(k) + partial%sa_drhototjk(1,np)*ofct
                                  ydrv(k) = ydrv(k) + partial%sa_drhototjk(2,np)*ofct
                                  zdrv(k) = zdrv(k) + partial%sa_drhototjk(3,np)*ofct
                                endif
                                if (nregionj.ne.nregionk) then
                                  xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np)*ofct
                                  yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np)*ofct
                                  zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np)*ofct
                                  xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototjk(1,np)*ofct
                                  yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototjk(2,np)*ofct
                                  zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototjk(3,np)*ofct
                                endif
                                if (lstr) then
                                  do kl = 1,nstrains
                                    klp = nstrptr(kl)
                                    rstrd(kl) = rstrd(kl) + partial%sa_drhotots(klp,np)*ofct
                                  enddo
                                  if (latomicstress) then
                                    do kl = 1,nstrains
                                      klp = nstrptr(kl)
                                      atomicstress(kl,i) = atomicstress(kl,i) + third*partial%sa_drhotots(klp,np)*ofct
                                      atomicstress(kl,j) = atomicstress(kl,j) + third*partial%sa_drhotots(klp,np)*ofct
                                      atomicstress(kl,k) = atomicstress(kl,k) + third*partial%sa_drhotots(klp,np)*ofct
                                    enddo
                                  endif
                                endif
                              enddo

                              call meamtotalrhoscreenderv(neamspecj,npartial,partial,xji,yji,zji,scrho(1,j), &
                                                          rhoj,rscrhoj,rhoji,Sij,lsg1)
!
                              do np = 1,npartial
                                k = partial%sa_atom(np)
                                lopk = (.not.lfreeze.or.lopf(nrelat(k)))
                                nregionk = nregionno(nsft+nrelat(k))
!
!  i-j contribution
!  
                                if (lopi) then
                                  xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np)*ofct
                                  ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np)*ofct
                                  zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np)*ofct
                                endif
                                if (lopj) then
                                  xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np)*ofct
                                  ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np)*ofct
                                  zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np)*ofct
                                endif
                                if (nregioni.ne.nregionj) then
                                  xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np)*ofct
                                  yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np)*ofct
                                  zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np)*ofct
                                  xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np)*ofct
                                  yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np)*ofct
                                  zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np)*ofct
                                endif
!
!  i-k contribution
!
                                if (lopi) then
                                  xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np)*ofct
                                  ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np)*ofct
                                  zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np)*ofct
                                endif
                                if (lopk) then
                                  xdrv(k) = xdrv(k) + partial%sa_drhototik(1,np)*ofct
                                  ydrv(k) = ydrv(k) + partial%sa_drhototik(2,np)*ofct
                                  zdrv(k) = zdrv(k) + partial%sa_drhototik(3,np)*ofct
                                endif
                                if (nregioni.ne.nregionk) then
                                  xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np)*ofct
                                  yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np)*ofct
                                  zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np)*ofct
                                  xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototik(1,np)*ofct
                                  yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototik(2,np)*ofct
                                  zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototik(3,np)*ofct
                                endif
!
!  j-k contribution
!
                                if (lopj) then
                                  xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np)*ofct
                                  ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np)*ofct
                                  zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np)*ofct
                                endif
                                if (lopk) then
                                  xdrv(k) = xdrv(k) + partial%sa_drhototjk(1,np)*ofct
                                  ydrv(k) = ydrv(k) + partial%sa_drhototjk(2,np)*ofct
                                  zdrv(k) = zdrv(k) + partial%sa_drhototjk(3,np)*ofct
                                endif
                                if (nregionj.ne.nregionk) then
                                  xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np)*ofct
                                  yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np)*ofct
                                  zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np)*ofct
                                  xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototjk(1,np)*ofct
                                  yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototjk(2,np)*ofct
                                  zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototjk(3,np)*ofct
                                endif
                                if (lstr) then
                                  do kl = 1,nstrains
                                    klp = nstrptr(kl)
                                    rstrd(kl) = rstrd(kl) + partial%sa_drhotots(klp,np)*ofct
                                  enddo
                                  if (latomicstress) then
                                    do kl = 1,nstrains
                                      klp = nstrptr(kl)
                                      atomicstress(kl,i) = atomicstress(kl,i) + third*partial%sa_drhotots(klp,np)*ofct
                                      atomicstress(kl,j) = atomicstress(kl,j) + third*partial%sa_drhotots(klp,np)*ofct
                                      atomicstress(kl,k) = atomicstress(kl,k) + third*partial%sa_drhotots(klp,np)*ofct
                                    enddo
                                  endif
                                endif
                              enddo
                            endif
                          endif
!
!  Scale density and derivatives by screening factor
!
                          drhoij(1:3,1:maxmeamcomponent) = Sij*drhoij(1:3,1:maxmeamcomponent)
                          drhoji(1:3,1:maxmeamcomponent) = Sij*drhoji(1:3,1:maxmeamcomponent)
                          if (lsg1) then
                            drhoijs(1:6,1:maxmeamcomponent) = Sij*drhoijs(1:6,1:maxmeamcomponent)
                            drhojis(1:6,1:maxmeamcomponent) = Sij*drhojis(1:6,1:maxmeamcomponent)
                          endif
!
!  Compute total derivatives of MEAM density
!
                          call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                                drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                                lsg1,.false.)
                          call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                                drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                                lsg1,.false.)
                        endif
                      else
                        call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xji,yji,zji,rhoij,rhoji,drhototij,drhototji, &
                                    drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                                    drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,lsg1,.true.,.false.,.false.)
                      endif
                    endif
                  enddo
!
!  Combine derivative terms
!
                  if (QMMMmode(ncf).gt.0) then
                    if (nregiontypi.ne.1) then
                      deriv(1:3) = deriv(1:3) + rscrhoi*drhototij(1:3)*ofct
                      if (lsg1) then
                        derivs(1:6) = derivs(1:6) + rscrhoi*drhototijs(1:6)*ofct
                      endif
                    endif
                    if (nregiontypj.ne.1) then
                      deriv(1:3) = deriv(1:3) + rscrhoj*drhototji(1:3)*ofct
                      if (lsg1) then
                        derivs(1:6) = derivs(1:6) + rscrhoj*drhototjis(1:6)*ofct
                      endif
                    endif
                  else
                    deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
                    if (lsg1) then
                      derivs(1:6) = derivs(1:6) + (rscrhoi*drhototijs(1:6) + rscrhoj*drhototjis(1:6))*ofct
                    endif
                  endif
                endif
                if (lvalidij) then
!******************************
!  Internal first derivatives *
!******************************
                  if (lMEAM) then
                    if (lopi) then
                      xdrvnr(i) = xdrvnr(i) - deriv(1)
                      ydrvnr(i) = ydrvnr(i) - deriv(2)
                      zdrvnr(i) = zdrvnr(i) - deriv(3)
                    endif
                    if (lopj) then
                      xdrvnr(j) = xdrvnr(j) + deriv(1)
                      ydrvnr(j) = ydrvnr(j) + deriv(2)
                      zdrvnr(j) = zdrvnr(j) + deriv(3)
                    endif
                  else
                    if (lopi) then
                      xdrv(i) = xdrv(i) - deriv(1)
                      ydrv(i) = ydrv(i) - deriv(2)
                      zdrv(i) = zdrv(i) - deriv(3)
                    endif
                    if (lopj) then
                      xdrv(j) = xdrv(j) + deriv(1)
                      ydrv(j) = ydrv(j) + deriv(2)
                      zdrv(j) = zdrv(j) + deriv(3)
                    endif
                  endif
                  if (nregioni.ne.nregionj) then
                    xregdrv(nregioni) = xregdrv(nregioni) - deriv(1)
                    yregdrv(nregioni) = yregdrv(nregioni) - deriv(2)
                    zregdrv(nregioni) = zregdrv(nregioni) - deriv(3)
                    xregdrv(nregionj) = xregdrv(nregionj) + deriv(1)
                    yregdrv(nregionj) = yregdrv(nregionj) + deriv(2)
                    zregdrv(nregionj) = zregdrv(nregionj) + deriv(3)
                  endif
!*****************************
!  Strain first derivatives  *
!*****************************
                  if (lstr) then
                    if (lMEAM) then
                      do kl = 1,nstrains
                        klp = nstrptr(kl)
                        rstrdnr(kl) = rstrdnr(kl) + derivs(klp)
                      enddo
                    else
                      do kl = 1,nstrains
                        klp = nstrptr(kl)
                        rstrd(kl) = rstrd(kl) + derivs(klp)
                      enddo
                    endif
                    if (latomicstress) then
                      do kl = 1,nstrains
                        klp = nstrptr(kl)
                        atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*derivs(klp)
                        atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*derivs(klp)
                      enddo
                    endif
                  endif
                endif
!
!  End of loop over atoms within neighbouring cells
!
              enddo jloop
!
!  End of loops over neighbouring cells
!
            enddo
          enddo
        enddo
!
!  End loop over atoms within cell
!
      enddo iloop
!  
!  End if for non-buffer cell
!
    endif
!
!  End loop over spatial decomposition cells
!
  enddo
!
!  End of real space part - perform general tasks
!
999 continue
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('manymd3s','npotl')
!
!  Exit point
!
1000 continue
!
!  Unscale density
!
  call eamscalescrho(-1_i4)
!
!  Timing
!
  time2 = cputime()
  tmany = tmany + time2 - time1
!
  return
  end
