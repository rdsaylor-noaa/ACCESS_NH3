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
!     Last Update:  July 2018                                                                                          !
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
!     Module:       Utils                                                                                              !
!                                                                                                                      ! 
!     Description:  contains general utility routines                                                                  !
!                                                                                                                      !
!======================================================================================================================!
module Utils
  use GlobalData
  use PhysChemData
  implicit none

  private 
  public SolveTridiag, CleanUp, ConvertPPBToMolecCC, ConvertMolecCCToPPB, Convertugm3ToMolecCC, &
         ConvertMolecCCTougm3, ConvertOneToOutputFormat, ConvertSourceToOutputFormat, &
         ConvertFluxToOutputFormat, ConvertEmissionRate, SolarZenithAngle, &
         TimeString, RelativeHumidity, SpecificHumidity, rtog, gtor, DateTime2SimTime, &
         Convert_h2o_to_qh, Convert_qh_to_h2o

contains

!**********************************************************************************************************************!
!  SolveTridiag - solves tridiagonal system
!**********************************************************************************************************************!
subroutine SolveTridiag(a, b, c, r, x, n)
  ! a = sub-diagonal (i.e., diagonal below the main diagonal)
  ! b = main diagonal
  ! c = super-diagonal (i.e., diagonal above the main diagonal)
  ! r = right hand side vector
  ! x = solution vector
  ! n = number of equations
  integer(kind=i4), intent(in) :: n
  real(kind=dp), dimension(n), intent(in) :: a,b,c,r
  real(kind=dp), dimension(n), intent(out) :: x
  real(kind=dp), dimension(n) :: beta, gam
  integer(kind=i4) :: i

  if (n .eq. 1) then
    x(1) = r(1)/b(1)
    return
  end if

  beta(1) = b(1)
  gam(1) = r(1)/beta(1)
  do i = 2,n
    beta(i) = b(i) - a(i)*c(i-1)/beta(i-1)
    gam(i) = (r(i) - a(i)*gam(i-1))/beta(i)
  end do 
   
  x(n) = gam(n)
  do i = n-1, 1, -1
    x(i) = gam(i) - c(i)*x(i+1)/beta(i) 
  end do

  return

end subroutine SolveTridiag

!**********************************************************************************************************************!
!  CleanUp - cleans up any loose ends before simulation is normally halted
!**********************************************************************************************************************!
subroutine CleanUp()

  ! deallocate arrays
  deallocate(cout)
  deallocate(vdout)
  deallocate(vsout)
  deallocate(gpout)
  deallocate(gpstout)
  deallocate(gpgout)
  deallocate(rbout)
  deallocate(rsout)
  deallocate(rwout)
  deallocate(rbgout)
  deallocate(rsoilout)
  deallocate(flxs)
  deallocate(flxw)
  deallocate(flxc)
  deallocate(fnet)
  deallocate(flxg)
  deallocate(flxh)
  deallocate(flxcs)
  deallocate(flxcw)
  deallocate(s)

  deallocate(ra)
  deallocate(rakv)
  deallocate(zolout)
  deallocate(gaero)

  deallocate(tkout)
  deallocate(pmbout)
  deallocate(qhout)
  deallocate(ubarout)
  deallocate(kvout)
  deallocate(ppfdout)
  deallocate(fjout)
  deallocate(cairout)
  deallocate(h2oout)
  deallocate(rhout)

  deallocate(qout)
  deallocate(vfout)
  deallocate(emtout)
  deallocate(adtout)
  deallocate(dptout)
  deallocate(vttout)
  deallocate(chtout)
  deallocate(timeout)
  deallocate(sdtout)
  deallocate(stdsp)
  if (ADVECTION .eqv. .TRUE.) then
    deallocate(advsp)
    deallocate(cadv)
  end if

  deallocate(ppfddirout)
  deallocate(ppfddifout)
  deallocate(nirdirout)
  deallocate(nirdifout)
  deallocate(ppfdsunout)
  deallocate(ppfdshdout)
  deallocate(ppfdwgtout)
  deallocate(nirsunout)
  deallocate(nirshdout)
  deallocate(nirwgtout)
  deallocate(lwupout)
  deallocate(lwdnout)
  deallocate(rtsunout)
  deallocate(rtshdout)
  deallocate(rtwgtout)
  deallocate(rabssunout)
  deallocate(rabsshdout)
  deallocate(rabswgtout)
  deallocate(fsunout)
  deallocate(fshdout)
  deallocate(tlsunout)
  deallocate(tlshdout)
  deallocate(tlwgtout)
  deallocate(gssunout)
  deallocate(gsshdout)
  deallocate(gswgtout)
  deallocate(rssunout)
  deallocate(rsshdout)
  deallocate(rswgtout)
  deallocate(anetsunout)
  deallocate(anetshdout)
  deallocate(anetwgtout)

  ! close metfile
  close(UENV)

  ! close simulation runtime file
  close(URUN)

  ! All done
  write(6,200) trim(simname)
200 format(/'ACCESS simulation ',a,' finished!')

  return
end subroutine CleanUp

!**********************************************************************************************************************!
! function rtog - converts resistance (s/m) to conductance (mol/m2-s)
!                 CAUTION: Assumes non-zero value for rz!
!
!**********************************************************************************************************************!
function rtog(rz, pmbi, tki)
  real(kind=dp), intent(in)  :: rz                ! resistance, s/m
  real(kind=dp), intent(in)  :: pmbi              ! air pressure, mb
  real(kind=dp), intent(in)  :: tki               ! temperature, K
  real(kind=dp)              :: rm2smol           ! resistance, m2-s/mol
  real(kind=dp)              :: rtog              ! conductance, mol/m2-s
  real(kind=dp), parameter   :: rgas=8.205D-05    ! ideal gas constant, m3-atm/K-mol

  rm2smol = rz*(pmbi/1013.)/(rgas*tki)
  rtog = 1.0/rm2smol

  return

end function rtog

!**********************************************************************************************************************!
! function gtor - converts conductance (mol/m2-s) to resistance (s/m)
!                 CAUTION: Assumes non-zero value for gz!
!
!**********************************************************************************************************************!
function gtor(gz, pmbi, tki)
  real(kind=dp), intent(in)  :: gz                ! conductance, mol/m2-s
  real(kind=dp), intent(in)  :: pmbi              ! air pressure, mb
  real(kind=dp), intent(in)  :: tki               ! temperature, K
  real(kind=dp)              :: gtor              ! resistance, s/m
  real(kind=dp)              :: gms               ! conductance, m/s
  real(kind=dp), parameter   :: rgas=8.205D-05    ! ideal gas constant, m3-atm/K-mol

  gms = gz*rgas*tki/(pmbi/1013.)
  gtor = 1.0/gms

  return

end function gtor

!**********************************************************************************************************************!
! DateTime2SimTime ... convert date string and time string to simulation time
!                      for calculation of solar zenith angle
!
!**********************************************************************************************************************!
function DateTime2SimTime(sdate, stime)
  character(len=10), intent(in)  :: sdate              ! date string "YYYY-MM-DD"
  character(len=8), intent(in)   :: stime              ! time string "hh:mm:ss"
  real(kind=dp)                  :: DateTime2SimTime
  real(kind=dp)    :: delta, leap, mjd, hrlocal, hrutc, time21
  real(kind=dp)    :: doy
  integer(kind=i4) :: m, year, month, daymonth, hz, mz, sz
  integer(kind=i4), dimension(12) :: months

  data months /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

  read(sdate(1:4), *) year
  read(sdate(6:7), *) month
  read(sdate(9:10), *) daymonth
  read(stime(1:2), *) hz
  read(stime(4:5), *) mz
  read(stime(7:8), *) sz
  
  ! calculate doy (day-of-year)
  doy = 0.0
  do m=1,month-1
    doy = doy + dble(months(m))
  end do
  doy = doy + dble(daymonth)

  ! calculate decimal UTC hour
  hrlocal = dble(hz) + (dble(mz)/60.) + (dble(sz)/3600.)
  hrutc = hrlocal - tzone

  ! years from 1949
  delta = dble(year)-1949.

  ! leap days over this period (note correction below)
  leap = aint(delta/4.)

  ! modified Julian Day
  mjd=32916.5+delta*365.+leap+dble(doy)+hrutc/24.

  ! the last year of century is not a leap year unless divisible by 400
  if (dmod(dble(year),100.0D+00) .eq. 0.0 .and. dmod(dble(year),400.0D+00) .ne. 0.0) then
    mjd=mjd-1.
  end if

  ! time in days since Jan 1, 2000 Greenwich Noon
  time21 = mjd-51545.0 
 
  ! time in seconds since Jan 1, 2000 Greenwich Noon
  ! simulation time is kept in seconds
! tstart=time21*24.*3600.

  DateTime2SimTime = time21

  return
end function DateTime2SimTime

!**********************************************************************************************************************!
! function ConvertPPBToMolecCC - converts ppbv to molecules/cm3
!**********************************************************************************************************************!
function ConvertPPBToMolecCC(cppb,izpt)
  integer(kind=i4) :: izpt
  real(kind=dp)    :: cppb
  real(kind=dp)    :: ConvertPPBToMolecCC

  ConvertPPBToMolecCC = cppb*cair(izpt)/1.0D+09

  return
end function ConvertPPBToMolecCC

!**********************************************************************************************************************!
! function ConvertMolecCCToPPB - converts molecules/cm3 to ppbv
!**********************************************************************************************************************!
function ConvertMolecCCToPPB(cmcc,izpt)
  integer(kind=i4) :: izpt
  real(kind=dp)    :: cmcc
  real(kind=dp)    :: ConvertMolecCCToPPB

  ConvertMolecCCToPPB = cmcc*1.0D+09/cair(izpt)

  return
end function ConvertMolecCCToPPB


!**********************************************************************************************************************!
! function ConvertMolecCCTougm3 - converts molecules/cm3 to micrograms/m3
!**********************************************************************************************************************!
function ConvertMolecCCTougm3(cmcc, isp)
  integer(kind=i4) :: isp
  real(kind=dp)    :: cmcc
  real(kind=dp)    :: ConvertMolecCCTougm3

  ConvertMolecCCTougm3 = cmcc*molecmass(isp)*1.0D+12/navo

  return 
end function ConvertMolecCCTougm3

!**********************************************************************************************************************!
! function Convertugm3ToMolecCC - converts micrograms/m3 to molecules/cm3
!**********************************************************************************************************************!
function Convertugm3ToMolecCC(cugm3, isp)
  integer(kind=i4) :: isp
  real(kind=dp)    :: cugm3
  real(kind=dp)    :: Convertugm3ToMolecCC

  Convertugm3ToMolecCC = cugm3*navo*1.0D-12/molecmass(isp)

end function Convertugm3ToMolecCC

!**********************************************************************************************************************!
! function ConvertEmissionRate - converts emission rate from ppb/hr 
!                                to molecules/s 
!**********************************************************************************************************************!
function ConvertEmissionRate(qppb, dzi, izpt)
  integer(kind=i4) :: izpt
  real(kind=dp)    :: qppb, dzi
  real(kind=dp)    :: ConvertEmissionRate

  ConvertEmissionRate = qppb*cair(izpt)*dzi/(3600.0D+09)

  return
end function ConvertEmissionRate

!**********************************************************************************************************************!
! function ConvertSourceToOutputFormat - based on data provided in the
!                                     simulation control file (via outppb)
!                                     either convert to ppbv/s or ng/m3-s
!**********************************************************************************************************************!
function ConvertSourceToOutputFormat(smcc, isp, izpt)
  integer(kind=i4) :: isp
  integer(kind=i4) :: izpt
  real(kind=dp)    :: smcc
  real(kind=dp)    :: ConvertSourceToOutputFormat

  ! if species isp is contained in the array outppb, 
  ! then convert source to ppb/s, otherwise convert to nanograms/m3-s
  if (ANY(outppb == isp)) then
    ConvertSourceToOutputFormat = smcc*1.0D+09/cair(izpt)
  else
    ConvertSourceToOutputFormat = smcc*molecmass(isp)*1.0D+15/navo
  end if
  
  return
end function ConvertSourceToOutputFormat

!**********************************************************************************************************************!
! function ConvertFluxToOutputFormat - based on data provided in the
!                                     simulation control file (via outppb)
!                                     either convert to ppbv-cm/s or ng/m2-s
!**********************************************************************************************************************!
function ConvertFluxToOutputFormat(fmcc, isp, izpt)
  integer(kind=i4) :: isp
  integer(kind=i4) :: izpt
  real(kind=dp)    :: fmcc
  real(kind=dp)    :: ConvertFluxToOutputFormat

  ! if species isp is contained in the array outppb, 
  ! then convert source to ppbv-cm/s, otherwise convert to nanograms/m2-s
  if (ANY(outppb == isp)) then
    ConvertFluxToOutputFormat = fmcc*1.0D+09/cair(izpt)
  else
    ConvertFluxToOutputFormat = fmcc*molecmass(isp)*1.0D+13/navo
  end if

  return
end function ConvertFluxToOutputFormat

!**********************************************************************************************************************!
! function ConvertOneToOutputFormat - based on data provided in the
!                                     simulation control file (via outppb)
!                                     either convert to ppbv or microg/m3
!**********************************************************************************************************************!
function ConvertOneToOutputFormat(cmcc, isp, izpt)
  integer(kind=i4) :: isp
  integer(kind=i4) :: izpt
  real(kind=dp)    :: cmcc
  real(kind=dp)    :: ConvertOneToOutputFormat

  ! if species isp is contained in the array outppb, 
  ! then convert it to PPB, otherwise convert to micrograms/m3
  if (ANY(outppb == isp)) then
    ConvertOneToOutputFormat = ConvertMolecCCToPPB(cmcc, izpt)
  else
    ConvertOneToOutputFormat = ConvertMolecCCTougm3(cmcc, isp) 
  end if

  return
end function ConvertOneToOutputFormat

!**********************************************************************************************************************!
! function RelativeHumidity - calculates the relative humidity from the
!                             supplied specific humidity, temperature and
!                             pressure
!**********************************************************************************************************************!
function RelativeHumidity(tki, pmbi, qhi)
  real(kind=dp), intent(in) :: tki
  real(kind=dp), intent(in) :: pmbi
  real(kind=dp), intent(in) :: qhi
  real(kind=dp)             :: RelativeHumidity
  real(kind=dp)             :: e, es, qhd, pkpa, rhi
  real(kind=dp), parameter  :: rhmin=0.1
  real(kind=dp), parameter  :: rhmax=99.0
  qhd=0.001*qhi
  pkpa=0.1*pmbi
  e = pkpa*qhd/(0.622+qhd)
  es = esat(tki)
  rhi = max(rhmin, min(rhmax, 100.0*e/es))   ! bound RH to (rhmin, rhmax)
  RelativeHumidity = rhi
  return 
end function RelativeHumidity

!**********************************************************************************************************************!
! function SpecificHumidity - calculates specific humidity (g/kg) from the supplied relative humidity, 
!                             temperature, and pressure
!**********************************************************************************************************************!
function SpecificHumidity(rhi, tki, pmbi)
  real(kind=dp), intent(in) :: rhi                 ! relative humidity (%)
  real(kind=dp), intent(in) :: tki                 ! air temperature (K)
  real(kind=dp), intent(in) :: pmbi                ! air pressure (mb)
  real(kind=dp)             :: SpecificHumidity    ! specific humidity (g/kg)
  real(kind=dp)             :: es                  ! saturation vapor pressure at tki (mb)
  real(kind=dp)             :: e                   ! ambient vapor pressure (mb)

  es = esat(tki)*10.0            ! kPa -> mb
  e  = es*rhi*0.01               ! mb

  SpecificHumidity = 622.0*e/(pmbi-0.378*e)

  return

end function SpecificHumidity

!**********************************************************************************************************************!
! function Convert_qh_to_h2o - convert qh (g/kg) to h2o (molecs/cm3)
!       
!**********************************************************************************************************************!
function Convert_qh_to_h2o(qhi, cairi)
  real(kind=dp), intent(in) :: qhi                   ! specific humidity at i (g/kg)
  real(kind=dp), intent(in) :: cairi                 ! air concentration at i (molecs/cm3)
  real(kind=dp)             :: Convert_qh_to_h2o     ! h2o concentration (molecs/cm3)

  Convert_qh_to_h2o = 0.001611*qhi*cairi
  
  return
end function Convert_qh_to_h2o

!**********************************************************************************************************************!
! function Convert_h2o_to_qh - convert h2o (molecs/cm3) to qh (g/kg)
!       
!**********************************************************************************************************************!
function Convert_h2o_to_qh(h2oi, cairi)
  real(kind=dp), intent(in) :: h2oi                 ! water vapor concentration at i (molecs/cm3)
  real(kind=dp), intent(in) :: cairi                ! air concentration at i (molecs/cm3)
  real(kind=dp)             :: Convert_h2o_to_qh    ! specific humidity at i (g/kg)

  Convert_h2o_to_qh = h2oi/(0.001611*cairi)

  return
end function Convert_h2o_to_qh

!**********************************************************************************************************************!
! function TimeString - construct and return a datetime string from the
!                       current simulation time
!       
!**********************************************************************************************************************!
function TimeString(tc)
  real(kind=dp), intent(in) :: tc    ! current simulation time (secs from 
                                     ! Jan 1, 2000 at Greenwich Noon
  character(len=19)         :: TimeString
  character(len=19)         :: tstr
  real(kind=dp)             :: dec_elapsedhours, dec_elapsedmins, dec_elapsedsecs
  real(kind=dp)             :: dec_elapseddays
  integer(kind=i4)          :: int_elapsedhours, int_elapsedmins, int_elapsedsecs
  integer(kind=i4)          :: int_elapseddays, int_elapsedmonths, int_elpasedyears
  integer(kind=i4)          :: yr, mn, dy, hr, mi, ss
  integer(kind=i4), dimension(12) :: ndays, ndaysleap
  integer(kind=i4)          :: ndaysmn

  data ndays/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data ndaysleap/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

  ! elapsed hours since simulation start
  dec_elapsedhours = (tc-tstart)/3600.
  int_elapsedhours = aint(dec_elapsedhours) 

  ! remainder minutes and seconds
  dec_elapsedmins  = (dec_elapsedhours - float(int_elapsedhours))*60.0
  int_elapsedmins  = aint(dec_elapsedmins)
  dec_elapsedsecs  = (dec_elapsedmins - float(int_elapsedmins))*60.0
  int_elapsedsecs  = aint(dec_elapsedsecs)
  mi = int_elapsedmins
  ss = int_elapsedsecs

  ! elapsed days since simulation start
  dec_elapseddays  = float(int_elapsedhours)/24.0
  int_elapseddays  = aint(dec_elapseddays)

  ! current hour (24hr clock) is 
  !  starting hour + total elapsed hours - number of hours in total elapsed days
  hr = hz + int_elapsedhours - int_elapseddays*24
  
  ! initial value for ndaysmn (number of days in current month)
  !  taking leap years into account
  mn = month
  yr = year
  if ( ((mod(yr,4) .eq. 0) .and. (mod(yr,100) .ne. 0)) .or. (mod(yr,400) .eq. 0) ) then 
    ! leap year
    ndaysmn = ndaysleap(mn)
  else
    ! non leap year
    ndaysmn = ndays(mn)
  end if

  ! total number of elapsed days + initial day of the month
  dy = daymonth + int_elapseddays

  ! keep going until dy is less than number of days in the current month
  do while (dy .gt. ndaysmn) 
    dy = dy - ndaysmn
    mn = mn+1
    if (mn .gt. 12) then
      ! if December has passed, reset mn and increment yr
      mn=1
      yr=yr+1
    end if
    ! determine number of days in new current month
    if ( ((mod(yr,4) .eq. 0) .and. (mod(yr,100) .ne. 0)) .or. (mod(yr,400) .eq. 0) ) then 
      ndaysmn = ndaysleap(mn)
    else
      ndaysmn = ndays(mn)
    end if
  end do

  !write(*,*) yr, mn, dy, hr, mi, ss

  ! construct datetime string
  write(tstr, 999) yr, mn, dy, hr, mi, ss

  TimeString = tstr

999 format(i4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)

  return
end function TimeString

!**********************************************************************************************************************!
! SolarZenithAngle ... calculates the topocentric (i.e., local) zenith angle                                           !
!                          (in degrees) based on the algorithm of ...                                                  !
!                                                                                                                      ! 
!      Michalsky, J.J., 1988. The Astronomical Almanac’s algorithm for                                                 !
!      approximate solar position (1950–2050). Solar Energy 40 (3), 227– 235.                                          !
!                                                                                                                      ! 
!      For elevation refraction correction, algorithm assumes                                                          !
!      surface T=288K and P=1013.25mb                                                                                  !
!                                                                                                                      !
!    input parameters                                                                                                  ! 
!      time21=modified Julian Day from Noon on Jan 1 2000                                                              !
!      slat=latitude in degrees (north is positive)                                                                    !
!      slon=longitude in degrees (east is positive)                    	                                               !
!                                                                                                                      !
!    internal parameters                                                                                               !
!	el= sun elevation angle (degs)                                 	                                               !
!	ha= solar hour angle						                                               !
!	dec = declination						                                               !
!                                                                                                                      !
!     returns                                                                                                          !
!	za = zenith angle (degs)					                                               !
!**********************************************************************************************************************!
function SolarZenithAngle(time21)
  real(kind=dp), intent(in) :: time21

  real(kind=dp) :: SolarZenithAngle

  real(kind=dp) :: twopi, pi, deg2rad
  real(kind=dp) :: mnlong, mnanom, eclong, oblqec, num, den, ra
  real(kind=dp) :: gmst, lmst, latrad, elc, refrac, soldia
  real(kind=dp) :: el, dec, ha

  pi = 2.0*acos(0.0)
  twopi = 2.0*pi
  deg2rad = pi/180.0
  
  ! force mean longitude between 0 and 360 degs
  mnlong=280.460+.9856474*time21
  mnlong=mod(mnlong,360.)
  if (mnlong.lt.0.) then
    mnlong=mnlong+360.
  end if

  ! mean anomaly in radians between 0 and 2*pi
  mnanom=357.528+.9856003*time21
  mnanom=mod(mnanom,360.)
  if(mnanom.lt.0.) then
    mnanom=mnanom+360.
  end if
  mnanom=mnanom*deg2rad

  ! compute the ecliptic longitude and obliquity of ecliptic in radians
  eclong=mnlong+1.915*sin(mnanom)+.020*sin(2.*mnanom)
  eclong=mod(eclong,360.)
  if (eclong.lt.0.) then
    eclong=eclong+360.
  end if
  oblqec=23.439-.0000004*time21
  eclong=eclong*deg2rad
  oblqec=oblqec*deg2rad

  ! calculate right ascension and declination
  num=cos(oblqec)*sin(eclong)
  den=cos(eclong)
  ra=atan(num/den)

  ! force ra between 0 and 2*pi
  if (den.lt.0) then
    ra=ra+pi
  elseif (num.lt.0) then
    ra=ra+twopi
  endif

  ! dec in radians
  dec=asin(sin(oblqec)*sin(eclong))

  ! calculate Greenwich mean sidereal time in hours
  ! Substitute the approximate gmst formula from the U. S. Naval Observatory web site
  ! http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time
  gmst = 18.697374558 + 24.06570982441908 * time21

  gmst = mod(gmst,24.)
  if (gmst.lt.0.) then
    gmst=gmst+24.
  end if

  ! calculate local mean sidereal time in radians 
  lmst=gmst+slon/15.
  lmst=mod(lmst,24.)
  if (lmst.lt.0.) then
    lmst=lmst+24.
  end if
  lmst=lmst*15.*deg2rad

  ! calculate hour angle in radians between -pi and pi
  ha=lmst-ra
  if (ha.lt.-pi) then
    ha=ha+twopi
  end if
  if (ha.gt.pi) then
    ha=ha-twopi
  end if

  ! change latitude to radians
  latrad=slat*deg2rad

  ! calculate elevation
  el=asin(sin(dec)*sin(latrad)+cos(dec)*cos(latrad)*cos(ha))

  ! calculate refraction correction for US stan. atmosphere
  ! need to have el in degs before calculating correction
  el=el/deg2rad

  ! note that 3.51823=1013.25 mb/288 K
  if(el.ge.19.225) then 
    refrac=.00452*3.51823/tan(el*deg2rad)
  else if (el.gt.-.766.and.el.lt.19.225) then
    refrac=3.51823*(.1594+.0196*el+.00002*el**2)/(1.+.505*el+.0845*el**2)
  else if (el.le.-.766) then
    refrac=0.0
  end if

  ! elevation in degs
  el=el+refrac

  ! zenith angle in degs
  SolarZenithAngle = 90.0 - el

end function SolarZenithAngle

end module Utils
