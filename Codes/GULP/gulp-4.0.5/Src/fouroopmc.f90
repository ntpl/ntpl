  subroutine fouroopmc(eoop,ntrialatom,nptrtrialatom,ltrialatom)
!
!  Subroutine for four-body energy from out of plane potentials - subset version
!
!   1/08 Created from fouroopmd
!   4/08 Minimum cutoff added for out of plane potentials
!   5/08 UFFoop potential added
!   5/08 only3 check added as an option
!   7/08 New type checking algorithm introduced
!   7/08 Handling of xangle potential added w.r.t. to type swapping
!  11/08 loutofplane array introduced to indicate whether a potential is out-of-plane 
!        type or not to replace potential number checking
!  11/08 Option to output energy terms added
!   4/10 Code modified to increase speed for bonded potentials by only searching over
!        bonded atoms
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
!  Julian Gale, NRI, Curtin University, April 2010
!
  use current
  use four
  use iochannels,     only : ioout
  use molecule
  use parallel
  use symmetry
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: ntrialatom
  integer(i4), intent(in)                      :: nptrtrialatom(ntrialatom)
  logical,     intent(in)                      :: ltrialatom(numat)
  real(dp),   intent(inout)                    :: eoop
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: isgn
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmk
  integer(i4)                                  :: indml
  integer(i4)                                  :: ixi
  integer(i4)                                  :: iyi
  integer(i4)                                  :: izi
  integer(i4)                                  :: ixj
  integer(i4)                                  :: iyj
  integer(i4)                                  :: izj
  integer(i4)                                  :: ixk
  integer(i4)                                  :: iyk
  integer(i4)                                  :: izk
  integer(i4)                                  :: ixl
  integer(i4)                                  :: iyl
  integer(i4)                                  :: izl
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjmax
  integer(i4)                                  :: jloop
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: k
  integer(i4)                                  :: kmax
  integer(i4)                                  :: kxx
  integer(i4)                                  :: kyy
  integer(i4)                                  :: kzz
  integer(i4)                                  :: l
  integer(i4)                                  :: l4
  integer(i4)                                  :: lj
  integer(i4)                                  :: lk
  integer(i4)                                  :: ll
  integer(i4)                                  :: lu
  integer(i4)                                  :: llmax
  integer(i4)                                  :: lmax
  integer(i4)                                  :: n
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeik
  integer(i4)                                  :: nbtypeil
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nbtypeik2
  integer(i4)                                  :: nbtypeil2
  integer(i4)                                  :: nfortype
  integer(i4)                                  :: ni
  integer(i4)                                  :: nj
  integer(i4)                                  :: nk
  integer(i4)                                  :: nl 
  integer(i4)                                  :: nmi 
  integer(i4)                                  :: nmj 
  integer(i4)                                  :: nmk 
  integer(i4)                                  :: nml 
  integer(i4)                                  :: nout
  integer(i4)                                  :: nouterloop
  integer(i4)                                  :: nt1
  integer(i4)                                  :: nt2
  integer(i4)                                  :: nt3
  integer(i4)                                  :: nt4
  integer(i4)                                  :: ntp2
  integer(i4)                                  :: ntp3
  integer(i4)                                  :: ntp4
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntyp3
  integer(i4)                                  :: ntyp4
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: ntypl
  integer(i4)                                  :: nunique
  integer(i4), dimension(:), allocatable       :: nuniqueptr
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: lbonded
  logical                                      :: lbtyp
  logical                                      :: linter_only
  logical                                      :: lintra_only
  logical                                      :: ljbond
  logical                                      :: lmatch
  logical                                      :: lmatch2
  logical                                      :: lmatch3
  logical                                      :: lmatchanyof2
  logical                                      :: lmatchanyof3
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lsamemol
  logical                                      :: lunique
  real(dp)                                     :: e1d(6)
  real(dp)                                     :: e2d(21)
  real(dp)                                     :: e3d(1)
  real(dp)                                     :: eterm
  real(dp)                                     :: fpoly(5)
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ocl 
  real(dp)                                     :: ofct 
  real(dp)                                     :: phi0 
  real(dp)                                     :: r21
  real(dp)                                     :: r212
  real(dp)                                     :: r31
  real(dp)                                     :: r312
  real(dp)                                     :: r32
  real(dp)                                     :: r322
  real(dp)                                     :: r41
  real(dp)                                     :: r412
  real(dp)                                     :: r42
  real(dp)                                     :: r422
  real(dp)                                     :: r43
  real(dp)                                     :: r432
  real(dp)                                     :: rkfor
  real(dp)                                     :: rkfor4
  real(dp)                                     :: rko
  real(dp)                                     :: rn
  real(dp)                                     :: tr1
  real(dp)                                     :: tr2
  real(dp)                                     :: tr3
  real(dp)                                     :: tr1min
  real(dp)                                     :: tr2min
  real(dp)                                     :: tr3min
  real(dp)                                     :: x21
  real(dp)                                     :: y21
  real(dp)                                     :: z21
  real(dp)                                     :: x21t
  real(dp)                                     :: y21t
  real(dp)                                     :: z21t
  real(dp)                                     :: x31
  real(dp)                                     :: y31
  real(dp)                                     :: z31
  real(dp)                                     :: x31t
  real(dp)                                     :: y31t
  real(dp)                                     :: z31t
  real(dp)                                     :: x41
  real(dp)                                     :: y41
  real(dp)                                     :: z41
  real(dp)                                     :: x41t
  real(dp)                                     :: y41t
  real(dp)                                     :: z41t
  real(dp)                                     :: x32
  real(dp)                                     :: y32
  real(dp)                                     :: z32
  real(dp)                                     :: x42
  real(dp)                                     :: y42
  real(dp)                                     :: z42
  real(dp)                                     :: x43
  real(dp)                                     :: y43
  real(dp)                                     :: z43
  real(dp)                                     :: xc1
  real(dp)                                     :: yc1
  real(dp)                                     :: zc1
  real(dp)                                     :: xc2t
  real(dp)                                     :: yc2t
  real(dp)                                     :: zc2t
  real(dp)                                     :: xc3t
  real(dp)                                     :: yc3t
  real(dp)                                     :: zc3t
  real(dp)                                     :: xc4t
  real(dp)                                     :: yc4t
  real(dp)                                     :: zc4t
!
!  Allocate local memory
!
  allocate(nuniqueptr(numat),stat=status)
  if (status/=0) call outofmemory('fouroopmc','nuniqueptr')
!
!  Openning banner for energy decomposition
!
  if (lPrintFour) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  OOP  : Atom No. 1  Atom No. 2  Atom No. 3  Atom No. 4  OutOfPlane energy (eV) '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!*************************
!  Loop over potentials  *
!*************************
!
!  Convolute outer two loops for parallelisation
!
  nouterloop = nfor*numat
  outerloop: do nout = procid+1,nouterloop,nprocs
    n = ((nout-1)/numat) + 1
!orig do 5 n=1,nfor
    nfortype = nforty(n)
    if (.not.loutofplane(n)) cycle outerloop
    ntp2 = 2
    ntp3 = 3
    ntp4 = 4
    nt1 = nfspec1(n)
    nt2 = nfspec2(n)
    nt3 = nfspec3(n)
    nt4 = nfspec4(n)
    ntyp1 = nfptyp1(n)
    ntyp2 = nfptyp2(n)
    ntyp3 = nfptyp3(n)
    ntyp4 = nfptyp4(n)
    tr1 = for1(n)**2
    tr2 = for2(n)**2
    tr3 = for3(n)**2
    tr1min = for1min(n)**2
    tr2min = for2min(n)**2
    tr3min = for3min(n)**2
    lbtyp = (mmfexc(n).eq.1)
    rkfor = fork(n)
    rkfor4 = forpoly(1,n)
    lintra_only = (lfintra(n).and..not.lfinter(n))
    linter_only = (lfinter(n).and..not.lfintra(n))
    lneedmol = (lintra_only.or.linter_only.or.lbtyp)
!********************************
!  Loop over middle site 1 / i  *
!********************************
    i = nout - (n-1)*numat
    ni = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
!
!  Check i is allowed for n
!
    if (.not.lmatch(ni,ntypi,nt1,ntyp1,.true.)) cycle outerloop
!
!  Only 3 bonds check
!
    if (lbtyp.and.lonly3oop(n).and.nbonds(i).ne.3) cycle outerloop
!
!  Set loop range for other atoms
!
    if (lbtyp) then
      ljbond = .true.
      if (nbonds(i).gt.0) then
        nunique = 1
        nuniqueptr(1) = nbonded(1,i)
        do jloop = 2,nbonds(i)
          lunique = .true.
          do lu = 1,nunique
            if (nbonded(jloop,i).eq.nuniqueptr(lu)) lunique = .false.
          enddo
          if (lunique) then
            nunique = nunique + 1
            nuniqueptr(nunique) = nbonded(jloop,i)
          endif
        enddo
        jloop = nunique
      else
        jloop = 0
      endif
    else
      ljbond = .false.
      jloop = numat
    endif
!
!  Skip if jloop is zero
!
    if (jloop.eq.0) cycle outerloop
!
!  i has been accepted
!
    xc1 = xclat(i)
    yc1 = yclat(i)
    zc1 = zclat(i)
!
!  Molecule handling
!
    if (lmol.and.lneedmol) then
      nmi = natmol(i)
      if (ndim.gt.0) then
        indm = nmolind(i)
        call mindtoijk(indm,ixi,iyi,izi)
      endif
    endif
!***********************************
!  Loop over first end site 2 / j  *
!***********************************
    ljloop: do lj = 1,jloop
      if (ljbond) then
        j = nuniqueptr(lj)
      else
        j = lj
      endif
      nj = nat(j)
      ntypj = nftype(j)
      ocj = occuf(j)
!
!  Check j is allowed for n
!
      lmatch3 = lmatchanyof3(nj,ntypj,ntp2,nt2,ntyp2,tr1,tr1min,ntp3,nt3,ntyp3,tr2,tr2min,ntp4,nt4,ntyp4,tr3,tr3min)
      if (.not.lmatch3) cycle ljloop
!
      if (lmol.and.lneedmol) then
!
!  Molecule handling
!
        nmj = natmol(j)
        if (ndim.gt.0) then
          indmj = nmolind(j)
          call mindtoijk(indmj,ixj,iyj,izj)
          ixj = ixj - ixi
          iyj = iyj - iyi
          izj = izj - izi
        endif
        lmolok = (nmi.eq.nmj.and.nmi.gt.0)
      else
        lmolok = .false.
      endif
!
!  Check for intra and but not in same molecule
!
      if (lintra_only.and..not.lmolok) cycle ljloop
      if (lbtyp.and..not.lmolok) cycle ljloop
!
      xc2t = xclat(j)
      yc2t = yclat(j)
      zc2t = zclat(j)
      x21t = xc2t - xc1
      y21t = yc2t - yc1
      z21t = zc2t - zc1
!
!  Check r21 is OK
!  Loop over cell vectors
!
      iiloop: do ii = 1,iimax
        r212 = (xvec1cell(ii)+x21t)**2 + (yvec1cell(ii)+y21t)**2 + (zvec1cell(ii)+z21t)**2
        if (r212.lt.1d-12) cycle iiloop
!
!  Molecule checking
!
        lbonded = .false.
        if (lmolok) then
          if (ndim.eq.0) then
            if (linter_only) cycle iiloop
            if (lbtyp) then
              call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,0_i4,0_i4,0_i4)
              if (.not.lbonded) cycle iiloop
            endif
          else
            call lintoijk(ixx,iyy,izz,ii,imaxl,jmaxl,kmaxl)
            if (lbtyp) then
              call bonded(lbonded,l2bonds,nbtypeij,nbtypeij2,i,j,ixx,iyy,izz)
              if (.not.lbonded) cycle iiloop
              lsamemol = (lbonded.or.l2bonds)
            else
              lsamemol = .false.
            endif
            if (.not.lsamemol) then
              call samemol(lsamemol,nmi,ixx,iyy,izz,ixj,iyj,izj)
            endif
            if (lintra_only.and..not.lsamemol) cycle iiloop
            if (linter_only.and.lsamemol) cycle iiloop
          endif
        endif
!
!  Distance checking
!
        if ((r212.gt.tr1.or.r212.lt.tr1min).and.(.not.lbtyp.or..not.lbonded)) cycle iiloop
        r21 = sqrt(r212)
        x21 = x21t + xvec1cell(ii)
        y21 = y21t + yvec1cell(ii)
        z21 = z21t + zvec1cell(ii)
        if (ndim.eq.0) then
          kmax = lj - 1
        else
          kmax = lj
        endif
!************************************
!  Loop over second end site 3 / k  *
!************************************
        lkloop: do lk = 1,kmax
          if (ljbond) then
            k = nuniqueptr(lk)
          else
            k = lk
          endif
          nk = nat(k)
          ntypk = nftype(k)
          ock = occuf(k)
!
!  Check k is allowed for n
!
          lmatch2 = lmatchanyof2(nk,ntypk,ntp3,nt3,ntyp3,tr2,tr2min,ntp4,nt4,ntyp4,tr3,tr3min)
          if (.not.lmatch2) cycle lkloop
!
          if (lmol.and.lneedmol) then
!
!  Molecule handling
!
            nmk = natmol(k)
            if (ndim.gt.0) then
              indmk = nmolind(k)
              call mindtoijk(indmk,ixk,iyk,izk)
              ixk = ixk - ixi
              iyk = iyk - iyi
              izk = izk - izi
            endif
            lmolok = (nmi.eq.nmk.and.nmi.gt.0)
          else
            lmolok = .false.
          endif
!
!  Check for intra and but not in same molecule
!
          if (lintra_only.and..not.lmolok) cycle lkloop
          if (lbtyp.and..not.lmolok) cycle lkloop
          xc3t = xclat(k)
          yc3t = yclat(k)
          zc3t = zclat(k)
          x31t = xc3t - xc1
          y31t = yc3t - yc1
          z31t = zc3t - zc1
          if (j.eq.k) then
            jjmax = ii - 1
          else
            jjmax = iimax
          endif
!
!  Check r31 is OK
!  Loop over cell vectors
!
          jjloop: do jj = 1,jjmax
            r312 = (xvec1cell(jj)+x31t)**2 + (yvec1cell(jj)+y31t)**2 + (zvec1cell(jj)+z31t)**2
            if (r312.lt.1d-12) cycle jjloop
!
!  Prevent atoms i and k being the same atom
!           
            if (k.eq.i.and.(jj.eq.iimid)) cycle jjloop
!  
!  Prevent atoms j and k being the same atom
!             
            if (k.eq.j.and.jj.eq.ii) cycle jjloop
!
!  Molecule checking
!
            lbonded = .false.
            if (lmolok) then
              if (ndim.eq.0) then
                if (linter_only) cycle jjloop
                if (lbtyp) then
                  call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,0_i4,0_i4,0_i4)
                  if (.not.lbonded) cycle jjloop
                endif
              else
                call lintoijk(jxx,jyy,jzz,jj,imaxl,jmaxl,kmaxl)
                if (lbtyp) then
                  call bonded(lbonded,l2bonds,nbtypeik,nbtypeik2,i,k,jxx,jyy,jzz)
                  if (.not.lbonded) cycle jjloop
                  lsamemol = (lbonded.or.l2bonds)
                else
                  lsamemol = .false.
                endif
                if (.not.lsamemol) then
                  call samemol(lsamemol,nmi,jxx,jyy,jzz,ixk,iyk,izk)
                endif
                if (lintra_only.and..not.lsamemol) cycle jjloop
                if (linter_only.and.lsamemol) cycle jjloop
              endif
            endif
!
!  Distance checking
!
            if ((r312.gt.tr2.or.r312.lt.tr2min).and.(.not.lbtyp.or..not.lbonded)) cycle jjloop
!
            r31 = sqrt(r312)
            x31 = x31t + xvec1cell(jj)
            y31 = y31t + yvec1cell(jj)
            z31 = z31t + zvec1cell(jj)
!
            if (ndim.eq.0) then
              lmax = lk - 1
            else
              lmax = lk
            endif
!**********************************
!  Loop over last end site 4 / l  *
!**********************************
            l4loop: do l4 = 1,lmax
              if (ljbond) then
                l = nuniqueptr(l4)
              else
                l = l4
              endif
              nl = nat(l)
              ntypl = nftype(l)
              ocl = occuf(l)
!
!  Check l is allowed for n
!
              if (.not.lmatch(nl,ntypl,nt4,ntyp4,.true.)) cycle l4loop
!
!  If no atoms are from trial set then skip
!
              if (.not.ltrialatom(i).and..not.ltrialatom(j).and..not.ltrialatom(k).and..not.ltrialatom(l)) cycle l4loop
!
              if (lmol.and.lneedmol) then
!
!  Molecule handling
!
                nml = natmol(l)
                if (ndim.gt.0) then
                  indml = nmolind(l)
                  call mindtoijk(indml,ixl,iyl,izl)
                  ixl = ixl - ixi
                  iyl = iyl - iyi
                  izl = izl - izi
                endif
                lmolok = (nmi.eq.nml.and.nmi.gt.0)
              else
                lmolok = .false.
              endif
!
!  Check for intra and but not in same molecule
!
              if (lintra_only.and..not.lmolok) cycle l4loop
              if (lbtyp.and..not.lmolok) cycle l4loop
!
              xc4t = xclat(l)
              yc4t = yclat(l)
              zc4t = zclat(l)
              x41t = xc4t - xc1
              y41t = yc4t - yc1
              z41t = zc4t - zc1
              if (k.eq.l) then
                llmax = jj - 1
              else
                llmax = iimax
              endif
!
!  Check r41 is OK
!  Loop over cell vectors
!
              llloop: do ll = 1,llmax
                r412 = (xvec1cell(ll)+x41t)**2 + (yvec1cell(ll)+y41t)**2 + (zvec1cell(ll)+z41t)**2
                if (r412.lt.1d-12) cycle llloop
!
!  Prevent atoms i and l being the same atom
!
                if (l.eq.i.and.ll.eq.iimid) cycle llloop
!
!  Prevent atoms j and l being the same atom
!
                if (l.eq.j.and.ll.eq.ii) cycle llloop
!
!  Prevent atoms k and l being the same atom
!
                if (l.eq.k.and.ll.eq.jj) cycle llloop
!
!  Molecule checking
!
                lbonded = .false.
                if (lmolok) then
                  if (ndim.eq.0) then
                    if (linter_only) cycle llloop
                    if (lbtyp) then
                      call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,0_i4,0_i4,0_i4)
                      if (.not.lbonded) cycle llloop
                    endif
                  else
                    call lintoijk(kxx,kyy,kzz,ll,imaxl,jmaxl,kmaxl)
                    if (lbtyp) then
                      call bonded(lbonded,l2bonds,nbtypeil,nbtypeil2,i,l,kxx,kyy,kzz)
                      if (.not.lbonded) cycle llloop
                      lsamemol = (lbonded.or.l2bonds)
                    else
                      lsamemol = .false.
                    endif
                    if (.not.lsamemol) then
                      call samemol(lsamemol,nmi,kxx,kyy,kzz,ixl,iyl,izl)
                    endif
                    if (lintra_only.and..not.lsamemol) cycle llloop
                    if (linter_only.and.lsamemol) cycle llloop
                  endif
                endif
!
!  Distance checking
!
                if ((r412.gt.tr3.or.r412.lt.tr3min).and.(.not.lbtyp.or..not.lbonded)) cycle llloop
!
                x41 = x41t + xvec1cell(ll)
                y41 = y41t + yvec1cell(ll)
                z41 = z41t + zvec1cell(ll)
!
!  Calculate other vectors needed
!
                x32 = x31 - x21
                y32 = y31 - y21
                z32 = z31 - z21
                r322 = x32*x32 + y32*y32 + z32*z32
                r32 = sqrt(r322)
                x43 = x41 - x31
                y43 = y41 - y31
                z43 = z41 - z31
!*********************************
!  Valid four-body term located  *
!*********************************
!
!  Calculate remaining distances
!
                r41 = sqrt(r412)
                x42 = x43 + x32
                y42 = y43 + y32
                z42 = z43 + z32
                r422 = x42*x42 + y42*y42 + z42*z42
                r432 = x43*x43 + y43*y43 + z43*z43
                r42 = sqrt(r422)
                r43 = sqrt(r432)
!
                ofct = oci*ocj*ock*ocl
                rko = rkfor*ofct
                if (nfortype.eq.14.or.nfortype.eq.16) then
                  if (ntp2.eq.2.and.ntp3.eq.3) then
                    fpoly(1) = forpoly(1,n)*ofct
                    fpoly(2) = forpoly(2,n)*ofct
                    fpoly(3) = forpoly(3,n)
                    fpoly(4) = forpoly(4,n)
                    fpoly(5) = forpoly(5,n)
                  elseif (ntp2.eq.2.and.ntp3.eq.4) then
                    fpoly(1) = rko
                    rko = forpoly(1,n)*ofct
                    fpoly(2) = forpoly(2,n)*ofct
                    fpoly(3) = forpoly(4,n)
                    fpoly(4) = forpoly(3,n)
                    fpoly(5) = forpoly(5,n)
                  elseif (ntp2.eq.3.and.ntp3.eq.2) then
                    fpoly(1) = forpoly(2,n)*ofct
                    fpoly(2) = forpoly(1,n)*ofct
                    fpoly(3) = forpoly(3,n)
                    fpoly(4) = forpoly(5,n)
                    fpoly(5) = forpoly(4,n)
                  elseif (ntp2.eq.3.and.ntp3.eq.4) then
                    fpoly(1) = rko
                    fpoly(2) = forpoly(1,n)*ofct
                    rko = forpoly(2,n)*ofct
                    fpoly(3) = forpoly(5,n)
                    fpoly(4) = forpoly(3,n)
                    fpoly(5) = forpoly(4,n)
                  elseif (ntp2.eq.4.and.ntp3.eq.2) then
                    fpoly(2) = rko
                    rko = forpoly(1,n)*ofct
                    fpoly(1) = forpoly(2,n)*ofct
                    fpoly(3) = forpoly(4,n)
                    fpoly(4) = forpoly(5,n)
                    fpoly(5) = forpoly(3,n)
                  elseif (ntp2.eq.4.and.ntp3.eq.3) then
                    fpoly(2) = rko
                    rko = forpoly(2,n)*ofct
                    fpoly(1) = forpoly(1,n)*ofct
                    fpoly(3) = forpoly(5,n)
                    fpoly(4) = forpoly(4,n)
                    fpoly(5) = forpoly(3,n)
                  endif
                elseif (nfortype.eq.15) then
                  fpoly(1) = forpoly(1,n)
                  fpoly(2) = forpoly(2,n)
                  fpoly(3) = forpoly(3,n)
                else
                  fpoly(1) = rkfor4*ofct
                  fpoly(2) = forpoly(2,n)
                endif
!
                call fourbody(n,nfortype,r21,r31,r41,r32,r42,r43,eterm,e1d,e2d,e3d,rko,rn,phi0,isgn,fpoly, &
                  .false.,.false.,.false.)
                eoop = eoop + eterm
!
!  Output energy contribution
!
                if (lPrintFour) then
                  write(ioout,'(4x,4i12,1x,f22.10)') i,j,k,l,eterm
                endif
!
!  End of inner loops over atoms and cell vectors
!
              enddo llloop
            enddo l4loop
          enddo jjloop
        enddo lkloop
      enddo iiloop
    enddo ljloop
!
!  End of outer loops
!
  enddo outerloop
!
!  Closing banner for energy decomposition
!
  if (lPrintFour) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Free local memory
!
  deallocate(nuniqueptr,stat=status)
  if (status/=0) call deallocate_error('fouroopmc','nuniqueptr')
!
!  All tidying up of derivatives is handled by four so we can just return here
!
  return
  end
