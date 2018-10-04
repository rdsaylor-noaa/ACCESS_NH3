!======================================================================================================================!
!                                                                                                                      !
!     Program:      ACCESS-NH3                                                                                         !
!                   Atmospheric Chemistry and Canopy Exchange Simulation System                                        !
!                   Ammonia version                                                                                    !
!                                                                                                                      !
!     Version:      3.0.0                                                                                              !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Initiated:    May 2014                                                                                           !
!     Last Update:  July 2018                                                                                           !
!                                                                                                                      !
!     Contact:      Rick D. Saylor, PhD                                                                                !
!                   Physical Scientist                                                                                 !
!                   U. S. Department of Commerce                                                                       !
!                   National Oceanic and Atmospheric Administration                                                    !
!                   Air Resources Laboratory                                                                           !
!                   Atmospheric Turbulence and Diffusion Division                                                      !
!                   456 S. Illinois Ave                                                                                !
!                   Oak Ridge, TN 37830                                                                                !
!                   email: Rick.Saylor@noaa.gov                                                                        !
!                                                                                                                      !
!**********************************************************************************************************************!
!                   NOTE: See Legal Notice in Main.f90                                                                 !
!**********************************************************************************************************************!
!                                                                                                                      ! 
!     Module:       Output                                                                                             !
!                                                                                                                      ! 
!     Description:  contains routines related to output of results                                                     !
!                                                                                                                      !
!======================================================================================================================!
module Output
  use GlobalData
  use Utils
  implicit none

  private PrinttoSTDOUT, SaveResults, CalculateVertFluxes, CalculateBudgetRates, &
          IntegrateCanopyFluxes
  public OutputResult, PrintFinaltoFile, PrintCPUtime
contains

!**********************************************************************************************************************!
! OutputResult - prints results to STDOUT and stores for later
!**********************************************************************************************************************!
subroutine OutputResult()

  ! save results for printing at end of simulation
  call SaveResults()

  ! print defined results to STDOUT and simulation runtime file
  call PrinttoSTDOUT()

  nt = nt + 1

  return
end subroutine OutputResult

!**********************************************************************************************************************!
! PrinttoSTDOUT - print selected results to STDOUT
!                 and simultaneously to simulation runtime file
!**********************************************************************************************************************!
subroutine PrinttoSTDOUT()
  integer(kind=i4)  :: i
  integer(kind=i4)  :: l
  character(len=2)  :: ncols
  character(len=7)  :: f0
  character(len=30) :: f1
  character(len=30) :: f2
  character(len=30) :: f3
  character(len=4)  :: col1
  real(kind=dp)     :: xfnet, xflxg, xflxh, xflxcs, xflxcw

  ! dynamically create format string
  write(ncols,'(i2)') nstdsp
  f0='(1x, a)'
  f1='(4x,a4,1x,' // ncols // '(5x,a4,5x))'
  f2='(f8.2,' // ncols // 'e14.6)'
  f3='(a8,' // ncols // 'e14.6)'
  col1='z(m)'

  ! array stdsp maps species indices for STDOUT output
  write(*,1001) (t-tstart)/3600.0_dp
  write(*,fmt=f0) sdtout(nt)

  write(*,1000) zadeg, dcos(zarad)
  write(*,1022) sradzref, ppfdzref
  write(*,1023) ppfd_direct, ppfd_diffus, nir_direct, nir_diffus
! write(*,fmt=f1) col1, (sspc(stdsp(l)),l=1,nstdsp)

  write(URUN,1001) (t-tstart)/3600.0_dp
  write(URUN,fmt=f0) sdtout(nt)
  write(URUN,fmt=f1) col1, (sspc(stdsp(l)),l=1,nstdsp)

  ! -ConvertOneToOutputFormat determines the appropriate output units
  !  from data in CTRL file
  do i=npts,1,-1
    write(*,f2) 0.01*z(i), (ConvertOneToOutputFormat(cint(i,stdsp(l)),stdsp(l),i),l=1,nstdsp) 

    write(URUN,f2) 0.01*z(i), (ConvertOneToOutputFormat(cint(i,stdsp(l)),stdsp(l),i),l=1,nstdsp) 
  end do

  write(*,*)
  write(URUN,*)
  write(*,1005) gammaso, gammast(1,1)          ! Only print out the value at the lowest level, for now
  write(URUN,1005) gammaso, gammast(1,1)
  write(*,1002) rakv(nt)
  write(URUN,1002) rakv(nt)

  if (ANY(outppb == iNH3)) then
    write(*,*) 'Fluxes (ppbv-cm/s)'
    write(URUN,*) 'Fluxes (ppbv-cm/s)'
  else
    write(*,*) 'Fluxes (ng/m2-s)'
    write(URUN,*) 'Fluxes (ng/m2-s)'
  end if

  xfnet=ConvertFluxToOutputFormat(fnet(1,nt),1,ncnpy)
  write(*,1010) xfnet
  write(URUN,1010) xfnet
  xflxg=ConvertFluxToOutputFormat(flxg(1,nt),1,1)
  write(*,1011) xflxg
  write(URUN,1011) xflxg
  xflxh=ConvertFluxToOutputFormat(flxh(1,nt),1,int(ncnpy/2))
  write(*,1012) xflxh
  write(URUN,1012) xflxh
  xflxcs=ConvertFluxToOutputFormat(flxcs(1,nt),1,int(ncnpy/2))
  write(*,1013) xflxcs
  write(URUN,1013) xflxcs
  xflxcw=ConvertFluxToOutputFormat(flxcw(1,nt),1,int(ncnpy/2))
  write(*,1014) xflxcw
  write(URUN,1014) xflxcw

  write(*,1015) xflxh*100.0/xflxg, xfnet*100.0/xflxg
  write(URUN,1015) xflxh*100.0/xflxg, xfnet*100.0/xflxg
  write(*,1016) xflxcs*100.0/xflxh, xflxcw*100.0/xflxh
  

1000 format(5x,'Zenith angle (deg)       = ', f7.1/ &
            5x,'Cosine(zenith angle)     = ', f8.2)
1022 format(5x,'Meas. Rg (W/m2)          = ', f7.1/ &
            5x,'Meas. PPFD = (umol/m2-s) = ', f7.1)
1023 format(5x,'PPFD direct (umol/m2-s)  = ', f7.1/ &
            5x,'PPFD diffuse (umol/m2-s) = ', f7.1/ &
            5x,'NIR direct (W/m2)        = ', f7.1/ &
            5x,'NIR diffuse (W/m2)       = ', f7.1)
1001 format(/' t-t0 = ', f6.2, ' hrs')
1002 format(' Rakv (s/cm) = ', e14.6)
1005 format(' Soil emission potential     = ', e14.6/' Stomatal emission potential = ', e14.6/)
1010 format(' Net flux             = ', e14.6)
1011 format(' Surface flux         = ', e14.6)
1012 format(' Total canopy flux    = ', e14.6)
1013 format(' Total stomatal flux  = ', e14.6)
1014 format(' Total cuticular flux = ', e14.6)
1015 format(/' Canopy % of surface flux = ', f6.1/' Net % of surface flux =    ', f6.1/)
1016 format(' Stomatal % of canopy flux =  ', f6.1/' Cuticular % of canopy flux = ', f6.1/)

  return
end subroutine PrinttoSTDOUT

!**********************************************************************************************************************!
! SaveResults - save all results to storage array
!**********************************************************************************************************************!
subroutine SaveResults()
  integer(kind=i4) :: l, i

  do l=1,ninteg
    do i=1,npts
      cout(i,l,nt) = cint(i,l)
      vdout(i,l,nt) = vd(i,l)
    end do
  end do

  do i=1,npts
    ! met data
    tkout(i,nt)    = tk(i)
    pmbout(i,nt)   = pmb(i)
    qhout(i,nt)    = qh(i)
    ubarout(i,nt)  = ubar(i)
    kvout(i,nt)    = kv(i)
    ppfdout(i,nt)  = ppfd(i)
    cairout(i,nt)  = cair(i)
    h2oout(i,nt)   = h2o(i) 
    rhout(i,nt)    = RelativeHumidity(tk(i), pmb(i), qh(i))

    ! canopy physics data
    ppfdsunout(i,nt) = ppfd_sun(i)
    ppfdshdout(i,nt) = ppfd_shd(i)
    ppfdwgtout(i,nt) = ppfd_wgt(i)
    nirsunout(i,nt)  = nir_sun(i)
    nirshdout(i,nt)  = nir_shd(i)
    nirwgtout(i,nt)  = nir_wgt(i)
    rabssunout(i,nt) = rabs_sun(i)
    rabsshdout(i,nt) = rabs_shd(i)
    rabswgtout(i,nt) = rabs_wgt(i)
    rtsunout(i,nt)   = rt_sun(i)
    rtshdout(i,nt)   = rt_shd(i)
    rtwgtout(i,nt)   = rt_wgt(i)
    lwupout(i,nt)    = lw_up(i)
    lwdnout(i,nt)    = lw_dn(i)
    tlsunout(i,nt)   = tl_sun(i)
    tlshdout(i,nt)   = tl_shd(i)
    tlwgtout(i,nt)   = tl_wgt(i)
    gssunout(i,nt)   = gs_sun(i)
    gsshdout(i,nt)   = gs_shd(i)
    gswgtout(i,nt)   = gs_wgt(i)
    rssunout(i,nt)   = rs_sun(i)
    rsshdout(i,nt)   = rs_shd(i)
    rswgtout(i,nt)   = rs_wgt(i)
    anetsunout(i,nt) = anet_sun(i)
    anetshdout(i,nt) = anet_shd(i)
    anetwgtout(i,nt) = anet_wgt(i)
    fsunout(i,nt)    = fsun(i)
    fshdout(i,nt)    = fshd(i)
  end do
  ppfddirout(nt) = ppfd_direct
  ppfddifout(nt) = ppfd_diffus
  nirdirout(nt)  = nir_direct
  nirdifout(nt)  = nir_diffus

  ! soil physics data
  vsh2oout(nt) = vsh2o
  qsoilout(nt) = qsoil
  effrhsoilout(nt) = effrhsoil
  rbgout(nt) = rbg
  gbgout(nt) = gbg
  rsoilout(nt) = rsoil

  call CalculateVertFluxes()
  call IntegrateCanopyFluxes()

  timeout(nt)=(t-tstart)/3600.0_dp     ! convert from seconds to hrs

  return
end subroutine SaveResults

!**********************************************************************************************************************!
! PrintFinaltoFile - prints selected results to individual species files
!                    output species are defined in CTRL file
!**********************************************************************************************************************!
subroutine PrintFinaltoFile()
  integer(kind=i4) :: i, l, m, j
  character(len=4)  :: ncols
  character(len=30) :: f1, f2, f3, f4, f5, f6
  character(len=4)  :: col1
  character(len=8)  :: srxn
  character(len=45) :: ofname
  character(len=40) :: hdr

  ! dynamically create format strings for output
  write(ncols,'(i4)') ntout+1
  f1='(6x,a4,3x,' // ncols // '(f8.1,5x))'
  f2='(f10.2,' // ncols // '(e13.5))'
  f3='(a)'
  f4='(a,i5.5)'
  f5='(6a10)'
  f6='(2i10,4f10.4)'
  col1='z(m)'

  ! output grid info
  ofname='./out/' // trim(simname) // '/grid/grid.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f5) 'npts','ncnpy','dz(m)','z0(m)','H(m)','hc(m)'
  write(UOUT,fmt=f6) npts, ncnpy, dzhc*0.01, z0*0.01, zi*0.01, hc   ! all heights in meters
  close(UOUT)

  ! output all species profiles
  ofname='./out/' // trim(simname) // '/sp/sp.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertOneToOutputFormat(cout(i,l,m),l,i),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output deposition velocities
  ofname='./out/' // trim(simname) // '/vd/vd.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(vdout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT) 

  ! output source/sink
  ofname='./out/' // trim(simname) // '/s/s.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertSourceToOutputFormat(s(i,l,m),l,i),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output vertical fluxes
  ofname='./out/' // trim(simname) // '/vflux/vflux.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertFluxToOutputFormat(vfout(i,l,m),l,i),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output stomatal fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxs.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertFluxToOutputFormat(flxs(i,l,m),l,i),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output cuticular fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxw.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertFluxToOutputFormat(flxw(i,l,m),l,i),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output canopy fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxc.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertFluxToOutputFormat(flxc(i,l,m),l,i),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output soil fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr        flxg(ng/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ConvertFluxToOutputFormat(flxg(1,m),1,1)
  end do
  close(UOUT) 

  ! output integrated canopy fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxh.dat'
  open(UOUT,file=ofname)
  hdr = ' hr        flxh(ng/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ConvertFluxToOutputFormat(flxh(1,m),1,1)
  end do
  close(UOUT) 

  ! output integrated stomatal fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxcs.dat'
  open(UOUT,file=ofname)
  hdr = ' hr        flxcs(ng/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ConvertFluxToOutputFormat(flxcs(1,m),1,1)
  end do
  close(UOUT) 

  ! output integrated cuticular fluxes
  ofname='./out/' // trim(simname) // '/cflux/flxcw.dat'
  open(UOUT,file=ofname)
  hdr = ' hr        flxcw(ng/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ConvertFluxToOutputFormat(flxcw(1,m),1,1)
  end do
  close(UOUT) 

  ! output net fluxes
  ofname='./out/' // trim(simname) // '/cflux/fnet.dat'
  open(UOUT,file=ofname)
  hdr = ' hr        fnet(ng/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ConvertFluxToOutputFormat(fnet(1,m),1,1)
  end do
  close(UOUT) 

  ! output canopy compensation points over simulation
  ofname='./out/' // trim(simname) // '/gp/gp.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertOneToOutputFormat(gpout(i,l,m),l,i),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output stomatal compensation points over simulation
  ofname='./out/' // trim(simname) // '/gp/gpst.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertOneToOutputFormat(gpstout(i,l,m),l,i),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output soil compensation points over simulation
  ofname='./out/' // trim(simname) // '/gp/gpg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         gpg(ug/m3)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ConvertOneToOutputFormat(gpgout(1,m),1,1)
  end do
  close(UOUT) 

  ! output soil exchange velocities over simulation
  ofname='./out/' // trim(simname) // '/vs/vs.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          vs(cm/s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), vsout(1,m)
  end do
  close(UOUT) 

  ! output leaf boundary layer resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rb.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rbout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output stomatal resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rs.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rsout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output cuticular resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rw.dat'
  open(UOUT,file=ofname)
  do l=1,ninteg
  ! write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rwout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output ground boundary layer resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rbg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         rbg(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), rbgout(1,m)
  end do
  close(UOUT) 

  ! output aerodynanic resistances over simulation
  ofname='./out/' // trim(simname) // '/r/ra.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          ra(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ra(m)
  end do
  close(UOUT) 

  ! output soil resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         rsoil(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), rsoilout(1,m)
  end do
  close(UOUT) 

  ! output temperature profiles over simulation
  ofname='./out/' // trim(simname) // '/met/tk.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tkout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output pressure profiles over simulation
  ofname='./out/' // trim(simname) // '/met/pmb.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(pmbout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output specific humidity profiles over simulation
  ofname='./out/' // trim(simname) // '/met/qh.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(qhout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output mean wind speed profiles over simulation
  ofname='./out/' // trim(simname) // '/met/ubar.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ubarout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output eddy diffusivity profiles over simulation
  ofname='./out/' // trim(simname) // '/met/kv.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(kvout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output PPFD profiles over simulation
  ofname='./out/' // trim(simname) // '/met/ppfd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output radiation attenuation profiles over simulation
  ofname='./out/' // trim(simname) // '/met/fj.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(fjout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output cair profiles over simulation
  ofname='./out/' // trim(simname) // '/met/cair.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(cairout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output h2o profiles over simulation
  ofname='./out/' // trim(simname) // '/met/h2o.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(h2oout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output relative humidity profiles over simulation
  ofname='./out/' // trim(simname) // '/met/rh.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rhout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output (z-d)/L over simulation
  ofname='./out/' // trim(simname) // '/met/zol.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          zol()'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), zolout(m)
  end do
  close(UOUT)

  ! output aerodynamic conductances over simulation
  ofname='./out/' // trim(simname) // '/met/gaero.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          gaero(mol/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), gaero(m)
  end do
  close(UOUT)

  ! output Ra's over simulation
  ofname='./out/' // trim(simname) // '/met/ra.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         ra(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ra(m)
  end do
  close(UOUT)


  ! output ppfd_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/ppfdsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output ppfd_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/ppfdshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output ppfd_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/ppfdwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output nir_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/nirsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(nirsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output nir_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/nirshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(nirshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output nir_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/nirwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(nirwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output lw_up profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/lwup.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(lwupout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output lw_dn profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/lwdn.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(lwdnout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rt_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rtsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rtsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rt_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rtshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rtshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rt_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rtwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rtwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rabs_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rabssun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rabssunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rabs_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rabsshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rabsshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rabs_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rabswgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rabswgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output tl_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/tlsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tlsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output tl_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/tlshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tlshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output tl_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/tlwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tlwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output gs_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/gssun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(gssunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output gs_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/gsshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(gsshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output gs_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/gswgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(gswgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rs_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rssun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rssunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rs_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rsshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rsshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rs_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rswgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rswgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output anet_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/anetsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(anetsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output anet_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/anetshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(anetshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output anet_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/anetwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(anetwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output fsun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/fsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(fsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output fshd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/fshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(fshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output clai and lai profiles
  ofname='./out/' // trim(simname) // '/canopy/laiprof.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt='(a10,2a13)') col1, 'lai', 'clai'
  do i=1,npts
    write(UOUT,fmt='(f10.1, 2e13.5)') 0.01*z(i), lai(i), clai(i)
  end do
  close(UOUT)

  ! output vsh2o exchange coeffients over simulation
  ofname='./out/' // trim(simname) // '/soil/vsh2o.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         vsh2o(cm/s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), vsh2oout(m)
  end do
  close(UOUT)

  ! output qsoil values over simulation
  ofname='./out/' // trim(simname) // '/soil/qsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         qsoil(mol/cm3)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), qsoilout(m)
  end do
  close(UOUT)

  ! output effrhsoil over simulation
  ofname='./out/' // trim(simname) // '/soil/effrhsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         effrhsoil'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), effrhsoilout(m)
  end do
  close(UOUT)

  ! output rbg over simulation
  ofname='./out/' // trim(simname) // '/soil/rbg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         rbg(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), rbgout(m)
  end do
  close(UOUT)

  ! output gbg over simulation
  ofname='./out/' // trim(simname) // '/soil/gbg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         gbg(mol/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), gbgout(m)
  end do
  close(UOUT)

  ! output rsoil over simulation
  ofname='./out/' // trim(simname) // '/soil/rsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         rsoil(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), rsoilout(m)
  end do
  close(UOUT)

  ! output tsoilk over simulation
  ofname='./out/' // trim(simname) // '/soil/tsoilk.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         tsoilk(K)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), tsoilkout(m)
  end do
  close(UOUT)

  ! output tk0 over simulation
  ofname='./out/' // trim(simname) // '/soil/tk0.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         tk0(K)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), tk0out(m)
  end do
  close(UOUT)

  !==========================================================
  ! output elapsed hour/datetime key file
  ofname='./out/' // trim(simname) // '/ACCESS_timekey.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          datetime'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1002) timeout(m), sdtout(m)
  end do
  close(UOUT) 
  

1000 format(4x, a)
1001 format(f8.1, 5x, e13.5)
1002 format(f8.1, 5x, a)

  return
end subroutine PrintFinaltoFile

!**********************************************************************************************************************!
! subroutine PrintCPUtime - prints simulation CPU time to STDOUT and summary
!**********************************************************************************************************************!
subroutine PrintCPUtime(tcpu)
   real(kind=4) :: tcpu

   ! to STDOUT
   write(*,1001) tcpu/3600.,tcpu

   ! to USUMM for archive
   open(USUMM,file=simsummary,status='old',action='write',position='append')
   write(USUMM,1001) tcpu/3600.,tcpu
   close(USUMM)

1001 format('Total CPU time: ', f12.3,' hrs'/'                ',f12.3,' sec')

end subroutine PrintCPUtime

!**********************************************************************************************************************!
! subroutine CalculateVertFluxes - calculate vertical fluxes of all species 
!                                  at all grid points
!**********************************************************************************************************************!
subroutine CalculateVertFluxes()
  integer(kind=i4) :: i, l
  real(kind=dp)    :: kvm, grad

  do l=1,ninteg
    do i=1,npts-1
      kvm = 0.5*(kv(i+1)+kv(i))
      grad  = (cint(i+1,l)-cint(i,l))/(z(i+1)-z(i))
      vfout(i,l,nt) = - kvm*grad
    end do
    kvm = kv(npts)
    grad = (caloft(l)-cint(npts,l))/(z(npts)-z(npts-1))
    vfout(npts,l,nt) = - kvm*grad
  end do
    
  return
end subroutine CalculateVertFluxes

!**********************************************************************************************************************!
! subroutine IntegrateCanopyFluxes - integrate component fluxes over the
!                                  the entire canopy at time nt 
!**********************************************************************************************************************!
subroutine IntegrateCanopyFluxes()
  integer(kind=i4) :: i, l
  real(kind=dp)    :: canopy_flxs, canopy_flxw, canopy_flxh

  do l=1,ninteg
    canopy_flxs = 0.0_dp
    canopy_flxw = 0.0_dp
    canopy_flxh = 0.0_dp
    do i=1,npts-1
      canopy_flxs = canopy_flxs + flxs(i,l,nt)*lad(i)*dzhc 
      canopy_flxw = canopy_flxw + flxw(i,l,nt)*lad(i)*dzhc 
      canopy_flxh = canopy_flxh + flxc(i,l,nt)*lad(i)*dzhc 
    end do
    fnet(l,nt) = flxg(l,nt) + canopy_flxs + canopy_flxw
    flxh(l,nt) = canopy_flxh 
    flxcs(l,nt) = canopy_flxs
    flxcw(l,nt) = canopy_flxw
  end do

  return
end subroutine IntegrateCanopyFluxes

!**********************************************************************************************************************!
! subroutine CalculateBudgetRates - calculate budget rates for all processes 
!                                   (ppbv/hr)
!**********************************************************************************************************************!
subroutine CalculateBudgetRates()
  integer(kind=i4) :: i, l

  ! xxtout accumulations are in ppbv
  ! dividing by dtout gives the rate in ppbv/s
  ! multiplying by 3600.0 converts to ppbv/hr
  do l=1,ninteg
    do i=1,npts
      emtout(i,l,nt)=emtout(i,l,nt)*3600.0/dtout
      adtout(i,l,nt)=adtout(i,l,nt)*3600.0/dtout
      dptout(i,l,nt)=dptout(i,l,nt)*3600.0/dtout
      vttout(i,l,nt)=vttout(i,l,nt)*3600.0/dtout
      chtout(i,l,nt)=chtout(i,l,nt)*3600.0/dtout
    end do
  end do 

  return
end subroutine CalculateBudgetRates

end module Output
