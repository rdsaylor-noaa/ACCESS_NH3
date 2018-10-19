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
!     Module:       GasChem                                                                                            !
!                   Not currently used in ACCESS-NH3                                                                   !
!                                                                                                                      ! 
!     Description:  Contains all routines related to gas-phase chemistry                                               !
!                                                                                                                      !
!======================================================================================================================!
module GasChem
  use GlobalData
  use Utils
  implicit none

  private GasPhotoRateCoeffs, GasThermalRateCoeffs,   &
          CalcJPhoto, fjo3a, fjo3b, fjh2o2, fjno2, fjno3a, fjno3b, fjhono,     &
          fjhno3, fjho2no2, fjhchoa, fjhchob, fjch3cho, fjc2h5cho, fjc3h7choa, &
          fjc3h7chob, fjic3h7cho, fjmacra, fjmacrb, fjacet, fjmek, fjmvka,     &
          fjmvkb, fjglyoxa, fjglyoxb, fjglyoxc, fjmglyox, fjbiacet, fjch3ooh,  &
          fjch3no3, fjc2h5no3, fjnc3h7no3, fjic3h7no3, fjtc4h9no3, fjnoaa,     &
          fjnoab, fjhoch2cho, fjacetol, fjpaa, fjpyra, fjhpald,                &
          ARR, ARR2, TERM, TJPL, TERCO, KTYP2, KTYP3, KTYPHO2
  public fgaschem, RxnRates, GasRateCoeffs
contains

!**********************************************************************************************************************!
! fgaschem - used by ODE solver to define the chemical system
!**********************************************************************************************************************!
subroutine fgaschem(neqn, ts, y, yp )
  integer(kind=i4) :: neqn
  real(kind=dp) ts
  real(kind=dp), dimension(neqn) :: y
  real(kind=dp), dimension(neqn) :: yp
  real(kind=dp), dimension(nrxn) :: k       ! reaction rate coefficients
  real(kind=dp), dimension(nrxn) :: r       ! reaction rates

  ! get reaction rate coefficients
  call GasRateCoeffs(ts, k)

  ! calculate reaction rates
  call RxnRates(neqn, k, y, r)

   ! Dummy routine for ACCESS_NH3

  return
end subroutine fgaschem

!**********************************************************************************************************************!
! RxnRates - calculates instantaneous values of the reaction rates
!            (molecules/cm3-s)
!**********************************************************************************************************************!
subroutine RxnRates(neqn,k,y,r)
  integer(kind=i4) :: neqn
  real(kind=dp), dimension(nrxn) :: k
  real(kind=dp), dimension(nrxn) :: r
  real(kind=dp), dimension(neqn) :: y
  real(kind=dp) :: H2, O2, H2Oz

   ! non-integrated species concentrations
   H2 = 500.0e-9*cair(zpt)        ! molecular hydrogen
   O2 = 0.21*cair(zpt)            ! diatomic oxygen
   H2Oz = h2o(zpt)                ! water vapor

   ! Dummy routine for ACCESS_NH3

  return
end subroutine rxnrates

!**********************************************************************************************************************!
! GasRateCoeffs - calculates gas-phase reaction rate coefficients
!                 at ts
!**********************************************************************************************************************!
subroutine GasRateCoeffs(ts, k)
  real(kind=dp) :: ts
  real(kind=dp), dimension(nrxn) :: k

  ! get photolysis rate coefficients
  call GasPhotoRateCoeffs(ts, k)

  ! get thermal rate coefficients
  call GasThermalRateCoeffs(ts, k)

  return
end subroutine GasRateCoeffs

!**********************************************************************************************************************!
! GasPhotoRateCoeffs - handles calculation of photolysis rate
!                      coefficients at ts
!**********************************************************************************************************************!
subroutine GasPhotoRateCoeffs(ts, k)
  real(kind=dp), dimension(nrxn) :: k
  real(kind=dp) :: ts, cz, za, pi, time21, zalt_m
  real(kind=dp), parameter :: mincz=0.017452406   ! cos(89deg)
  real(kind=dp), dimension(65) :: J    ! MCM J-values

  ! set current z value
  zalt_m=zcurrent

  ! calculate current zenith angle
  if (OPT_SIMTYPE .eq. SSTATE) then
    ! for steady-state, solar time doesn't change
    ! convert tstart to days
    time21 = tstart/(24.*3600.)
  else
    ! for dynamic, use current simulation time
    ! convert ts to days
    time21 = ts/(24.*3600.)
  end if
  za = SolarZenithAngle(time21)
  pi=2.0*acos(0.0)
  cz = dcos(za*pi/180.) 
  if (cz < mincz) then
    ! nighttime, so don't bother calculating anything
  else
    ! daytime, calculate the Js
    call CalcJPhoto(cz, zalt_m)

  end if

  return
end subroutine GasPhotoRateCoeffs

!**********************************************************************************************************************!
! CalcJPhoto - calculates photolysis frequencies (1/s) for a given
!              altitude and cosine(zenith angle).
!              Photolysis frequencies were generated from TUV 5.0
!              (Madronich, S. and S. Flocke, Theoretical estimation of 
!              biologically effective UV radiation at the Earth's surface, 
!              in Solar Ultraviolet Radiation - Modeling, Measurements and 
!              Effects (Zerefos, C., ed.). NATO ASI Series Vol. I52, 
!              Springer-Verlag, Berlin, 1997.) as a function of altitude 
!              above sea level and zenith angle.  For each photolysis reaction 
!              polynomial fits were created as a function of cosine(zenith 
!              angle) at seven altitudes.  At a given altitude, the 
!              photolysis frequency for each reaction is determined 
!              by interpolation between bounding altitudes. Frequencies are
!              valid from 0 to 12 km above mean sea level.
!**********************************************************************************************************************!
subroutine CalcJPhoto(cz, zalt_m)
  integer(kind=i4) :: i, j, l, lm1, jp
  real(kind=dp)    :: cz, zalt_m
  real(kind=dp)    :: alpha, zf, zf1, cz2, cz3, cz4, cz5, cz6

  ! calculate powers of cz
  cz2=cz*cz
  cz3=cz2*cz
  cz4=cz3*cz
  cz5=cz4*cz
  cz6=cz5*cz

  ! determine levels for interpolation
  if (zalt_m .ge. 12000.) then    ! just use values at 12km
    alpha = 1.0
    l = 7
    lm1 = 7
  else if (zalt_m .lt. 12000. .and. zalt_m .ge. 10000.) then
    alpha = (12000. - zalt_m)/2000.
    l = 7
    lm1 = 6
  else if (zalt_m .lt. 10000. .and. zalt_m .ge. 8000.) then
    alpha = (10000. - zalt_m)/2000.
    l = 6
    lm1 = 5
  else if (zalt_m .lt. 8000. .and. zalt_m .ge. 6000.) then
    alpha = (8000. - zalt_m)/2000.
    l = 5
    lm1 = 4
  else if (zalt_m .lt. 6000. .and. zalt_m .ge. 4000.) then
    alpha = (6000. - zalt_m)/2000.
    l = 4
    lm1 = 3
  else if (zalt_m .lt. 4000. .and. zalt_m .ge. 2000.) then
    alpha = (4000. - zalt_m)/2000.
    l = 3
    lm1 = 2
  else if (zalt_m .lt. 2000.) then
    alpha = (2000. - zalt_m)/2000.
    l = 2
    lm1 = 1
  end if

  ! calculate interpolation parameters
  zf  = alpha
  zf1 = 1.0 - alpha

  ! calculate photolysis frequencies
  kphoto(jo3a) = fjo3a(cz, cz2, cz3, cz4, l, lm1, zf, zf1)
  kphoto(jo3b) = fjo3b(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  kphoto(jh2o2) = fjh2o2(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jno2) = fjno2(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  kphoto(jno3a) = fjno3a(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  kphoto(jno3b) = fjno3b(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  kphoto(jhono) = fjhono(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jhno3) = fjhno3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jho2no2) = fjho2no2(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jhchoa) = fjhchoa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jhchob) = fjhchob(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jch3cho) = fjch3cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jc2h5cho) = fjc2h5cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jc3h7choa) = fjc3h7choa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jc3h7chob) = fjc3h7chob(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jic3h7cho) = fjic3h7cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jmacra) = fjmacra(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jmacrb) = fjmacrb(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jacet) = fjacet(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jmek) = fjmek(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jmvka) = fjmvka(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jmvkb) = fjmvkb(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jglyoxa) = fjglyoxa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jglyoxb) = fjglyoxb(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jglyoxc) = fjglyoxc(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jmglyox) = fjmglyox(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jbiacet) = fjbiacet(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  kphoto(jch3ooh) = fjch3ooh(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jch3no3) = fjch3no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jc2h5no3) = fjc2h5no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jnc3h7no3) = fjnc3h7no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jic3h7no3) = fjic3h7no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jtc4h9no3) = fjtc4h9no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jnoaa) = fjnoaa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jnoab) = fjnoab(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jhoch2cho) = fjhoch2cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jacetol) = fjacetol(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jpaa) = fjpaa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jpyra) = fjpyra(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  kphoto(jhpald) = fjhpald(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)

  ! apply fj factor to account for reduction of actinic flux
  ! within the canopy
  do jp=1,nphoto
    kphoto(jp)=kphoto(jp)*fj(zpt)
  end do

  return
end subroutine CalcJPhoto

!**********************************************************************************************************************!
! fjo3a - function to calculate photolysis frequency (1/s) for the reaction
!
!         O3 + hv --> O(1D) + O2
!         ozone
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjo3a(cz, cz2, cz3, cz4, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjo3a
  
  data a0/+4.651D-08,+8.749D-08,+1.055D-07,+1.133D-07,+1.151D-07,+1.183D-07,+1.272D-07/
  data a1/+1.471D-06,-5.579D-07,-1.934D-06,-2.890D-06,-3.460D-06,-3.884D-06,-4.293D-06/
  data a2/-1.785D-06,+1.137D-05,+2.178D-05,+3.069D-05,+3.851D-05,+4.690D-05,+5.839D-05/
  data a3/+7.576D-05,+6.905D-05,+5.626D-05,+4.228D-05,+2.885D-05,+1.512D-05,-2.091D-07/
  data a4/-3.260D-05,-3.490D-05,-3.152D-05,-2.661D-05,-2.164D-05,-1.637D-05,-1.054D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1) + a1(lm1)*cz + a2(lm1)*cz2 + a3(lm1)*cz3 + a4(lm1)*cz4
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjo3a = zf1*jphi + zf*jplo 

end function fjo3a

!**********************************************************************************************************************!
! fjo3b - function to calculate photolysis frequency (1/s) for the reaction
!
!         O3 + hv --> O(3P) + O2
!         ozone
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjo3b(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5, cz6
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5, a6
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjo3b
  
  data a0/+4.516D-06,+7.547D-06,+1.859D-05,+2.901D-05,+3.846D-05,+4.970D-05,+6.275D-05/
  data a1/+0.0002296,+0.0011360,+0.0015160,+0.0017150,+0.0018930,+0.0020350,+0.0021490/
  data a2/+0.0048210,-0.0004056,-0.0031760,-0.0047490,-0.0062100,-0.0074570,-0.0085580/
  data a3/-0.0181000,-0.0048620,+0.0032900,+0.0080470,+0.0125900,+0.0164900,+0.0200600/
  data a4/+0.0290600,+0.0114500,-0.0006084,-0.0077740,-0.0147600,-0.0207100,-0.0262900/
  data a5/-0.0220600,-0.0101800,-0.0013790,+0.0039210,+0.0091710,+0.0136000,+0.0178300/
  data a6/+0.0064610,+0.0032550,+0.0007384,-0.0007947,-0.0023330,-0.0036150,-0.0048620/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5 + a6(l)*cz6
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5+a6(lm1)*cz6
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjo3b = zf1*jphi + zf*jplo 

end function fjo3b

!**********************************************************************************************************************!
! fjh2o2 -function to calculate photolysis frequency (1/s) for the reaction
!
!         H2O2 + hv --> 2OH
!         hydrogen peroxide
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjh2o2(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjh2o2
  
  data a0/+5.089D-08,+6.023D-08,+5.978D-08,+5.443D-08,+4.627D-08,+3.871D-08,+3.673D-08/
  data a1/+1.201D-06,+5.192D-07,+2.428D-07,+4.996D-07,+1.497D-06,+3.430D-06,+6.449D-06/
  data a2/+1.137D-05,+2.353D-05,+3.276D-05,+3.898D-05,+4.103D-05,+3.737D-05,+2.728D-05/
  data a3/-8.381D-07,-2.751D-05,-5.306D-05,-7.435D-05,-8.808D-05,-8.974D-05,-7.671D-05/
  data a4/-8.771D-06,+1.309D-05,+3.825D-05,+6.138D-05,+7.863D-05,+8.497D-05,+7.717D-05/
  data a5/+4.551D-06,-1.956D-06,-1.065D-05,-1.915D-05,-2.595D-05,-2.914D-05,-2.739D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjh2o2 = zf1*jphi + zf*jplo 

end function fjh2o2

!**********************************************************************************************************************!
! fjno2 - function to calculate photolysis frequency (1/s) for the reaction
!
!         NO2 + hv --> NO + O(3P)
!         nitrogen dioxide
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjno2(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5, cz6
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5, a6
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjno2
  
  data a0/+0.0001472,+0.0001613,+0.0001628,+0.0001622,+0.0001649,+0.0001825,+0.0002437/
  data a1/+0.0028890,+0.0013940,+0.0023140,+0.0054820,+0.0113900,+0.0204400,+0.0321600/
  data a2/+0.0341300,+0.0841500,+0.1090000,+0.1127000,+0.0920100,+0.0416500,-0.0359800/
  data a3/-0.0451700,-0.2353000,-0.3583000,-0.4164000,-0.3967000,-0.2779000,-0.0626700/
  data a4/+0.0050470,+0.3056000,+0.5269000,+0.6547000,+0.6644000,+0.5174000,+0.2125000/
  data a5/+0.0254200,-0.1974000,-0.3759000,-0.4889000,-0.5149000,-0.4221000,-0.2063000/
  data a6/-0.0127800,+0.0509600,+0.1052000,+0.1414000,+0.1527000,+0.1291000,+0.0686500/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5 + a6(l)*cz6
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5+a6(lm1)*cz6
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjno2 = zf1*jphi + zf*jplo 

end function fjno2

!**********************************************************************************************************************!
! fjno3a - function to calculate photolysis frequency (1/s) for the reaction
!
!         NO3 + hv --> NO + O2
!         nitrate radical
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjno3a(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5, cz6
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5, a6
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjno3a
  
  data a0/+0.0002113,+0.0003791,+0.0010420,+0.0016620,+0.0021950,+0.0027890,+0.0034410/
  data a1/+0.0127900,+0.0685600,+0.0898500,+0.0999900,+0.1084000,+0.1154000,+0.1216000/
  data a2/+0.2924000,-0.0422400,-0.2031000,-0.2895000,-0.3624000,-0.4298000,-0.4928000/
  data a3/-1.1290000,-0.2637000,+0.2125000,+0.4797000,+0.7074000,+0.9254000,+1.1330000/
  data a4/+1.8310000,+0.6624000,-0.0407600,-0.4453000,-0.7923000,-1.1310000,-1.4550000/
  data a5/-1.3960000,-0.5982000,-0.0862600,+0.2130000,+0.4710000,+0.7260000,+0.9697000/
  data a6/+0.4097000,+0.1925000,+0.0462600,-0.0401200,-0.1150000,-0.1896000,-0.2607000/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5 + a6(l)*cz6
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5+a6(lm1)*cz6
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjno3a = zf1*jphi + zf*jplo 

end function fjno3a

!**********************************************************************************************************************!
! fjno3b - function to calculate photolysis frequency (1/s) for the reaction
!
!         NO3 + hv --> NO2 + O(3P)
!         nitrate radical
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjno3b(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5, cz6
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5, a6
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjno3b
  
  data a0/+0.0020070,+0.0026260,+0.0057260,+0.0092890,+0.0129700,+0.0176800,+0.0235300/
  data a1/+0.0571400,+0.4115000,+0.6146000,+0.7388000,+0.8471000,+0.9384000,+1.0090000/
  data a2/+2.3310000,+0.5458000,-0.8554000,-1.7840000,-2.6550000,-3.4460000,-4.1340000/
  data a3/-8.4190000,-4.5620000,-0.5793000,+2.1500000,+4.8380000,+7.3330000,+9.5830000/
  data a4/+13.180000,+8.7900000,+3.0430000,-0.9543000,-5.0410000,-8.8670000,-12.390000/
  data a5/-9.8160000,-7.2660000,-3.1420000,-0.2579000,+2.7830000,+5.6390000,+8.3050000/
  data a6/+2.8340000,+2.2370000,+1.0700000,+0.2534000,-0.6305000,-1.4610000,-2.2440000/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5 + a6(l)*cz6
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5+a6(lm1)*cz6
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjno3b = zf1*jphi + zf*jplo 

end function fjno3b

!**********************************************************************************************************************!
! fjhono -function to calculate photolysis frequency (1/s) for the reaction
!
!         HONO + hv --> NO + OH
!         nitrous acid
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjhono(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjhono
  
  data a0/+2.977D-05,+3.129D-05,+2.922D-05,+2.659D-05,+2.507D-05,+2.717D-05,+3.790D-05/
  data a1/+0.0005187,+0.0005138,+0.0009115,+0.0016810,+0.0028910,+0.0046170,+0.0067360/
  data a2/+0.0068190,+0.0120700,+0.0138000,+0.0128200,+0.0088960,+0.0012700,-0.0095280/
  data a3/-0.0103500,-0.0258100,-0.0337600,-0.0357100,-0.0311000,-0.0177400,+0.0034860/
  data a4/+0.0064020,+0.0222400,+0.0317300,+0.0357200,+0.0335400,+0.0227000,+0.0038820/
  data a5/-0.0014490,-0.0070690,-0.0107400,-0.0125800,-0.0123200,-0.0089600,-0.0027230/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjhono = zf1*jphi + zf*jplo 

end function fjhono

!**********************************************************************************************************************!
! fjhno3 -function to calculate photolysis frequency (1/s) for the reaction
!
!         HNO3 + hv --> NO2 + OH
!         nitric acid
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjhno3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjhno3
  
  data a0/+1.754D-09,+2.216D-09,+2.304D-09,+2.186D-09,+1.905D-09,+1.555D-09,+1.284D-09/
  data a1/+6.207D-08,+2.199D-08,-1.322D-08,-3.606D-08,-3.789D-08,-1.576D-08,+3.566D-08/
  data a2/+3.135D-07,+7.865D-07,+1.251D-06,+1.669D-06,+1.975D-06,+2.141D-06,+2.219D-06/
  data a3/+1.480D-06,+8.807D-07,-1.622D-07,-1.315D-06,-2.336D-06,-3.109D-06,-3.621D-06/
  data a4/-1.713D-06,-1.578D-06,-7.441D-07,+3.377D-07,+1.377D-06,+2.243D-06,+2.843D-06/
  data a5/+5.908D-07,+6.380D-07,+3.952D-07,+3.581D-08,-3.276D-07,-6.481D-07,-8.732D-07/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjhno3 = zf1*jphi + zf*jplo 

end function fjhno3

!**********************************************************************************************************************!
! fjho2no2 -function to calculate photolysis frequency (1/s) for the reaction
!
!         HO2NO2 + hv --> NO2 + HO2
!         peroxynitric acid
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjho2no2(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjho2no2
  
  data a0/+1.077D-08,+1.596D-08,+1.883D-08,+1.953D-08,+1.828D-08,+1.518D-08,+1.130D-08/
  data a1/+3.754D-07,+5.299D-08,-3.052D-07,-5.758D-07,-6.780D-07,-4.707D-07,+1.485D-07/
  data a2/+2.651D-06,+7.015D-06,+1.231D-05,+1.782D-05,+2.316D-05,+2.730D-05,+2.950D-05/
  data a3/+1.070D-05,+5.948D-06,-4.574D-06,-1.772D-05,-3.227D-05,-4.550D-05,-5.497D-05/
  data a4/-1.318D-05,-1.270D-05,-4.641D-06,+7.175D-06,+2.146D-05,+3.534D-05,+4.611D-05/
  data a5/+4.683D-06,+5.393D-06,+3.123D-06,-7.115D-07,-5.656D-06,-1.066D-05,-1.470D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjho2no2 = zf1*jphi + zf*jplo 

end function fjho2no2

!**********************************************************************************************************************!
! fjhchoa -function to calculate photolysis frequency (1/s) for the reaction
!
!         HCHO + hv --> HCO + H
!         formaldehyde
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjhchoa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjhchoa
  
  data a0/+1.412D-07,+1.931D-07,+2.137D-07,+2.093D-07,+1.841D-07,+1.449D-07,+1.084D-07/
  data a1/+3.562D-06,+4.664D-07,-2.063D-06,-2.984D-06,-1.134D-06,+4.994D-06,+1.618D-05/
  data a2/+4.027D-05,+9.071D-05,+0.0001419,+0.0001895,+0.0002266,+0.0002427,+0.0002268/
  data a3/+2.702D-05,-6.359D-05,-0.0001846,-0.0003149,-0.0004364,-0.0005222,-0.0005401/
  data a4/-6.610D-05,-6.181D-06,+0.0001025,+0.0002315,+0.0003619,+0.0004661,+0.0005084/
  data a5/+2.836D-05,+1.440D-05,-2.070D-05,-6.559D-05,-0.0001132,-0.0001537,-0.0001736/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjhchoa = zf1*jphi + zf*jplo 

end function fjhchoa

!**********************************************************************************************************************!
! fjhchob -function to calculate photolysis frequency (1/s) for the reaction
!
!         HCHO + hv --> CO + H2
!         formaldehyde
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjhchob(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjhchob
  
  data a0/+4.517D-07,+5.656D-07,+5.925D-07,+5.692D-07,+5.128D-07,+4.635D-07,+5.013D-07/
  data a1/+1.030D-05,+5.606D-06,+5.209D-06,+1.105D-05,+2.716D-05,+5.819D-05,+0.0001064/
  data a2/+0.0001000,+0.0002175,+0.0003129,+0.0003823,+0.0004014,+0.0003365,+0.0001600/
  data a3/-6.849D-05,-0.0003420,-0.0006091,-0.0008481,-0.0009996,-0.0009793,-0.0007096/
  data a4/-1.357D-05,+0.0002268,+0.0004952,+0.0007589,+0.0009558,+0.0009989,+0.0008038/
  data a5/+1.939D-05,-5.669D-05,-0.0001509,-0.0002491,-0.0003282,-0.0003565,-0.0003019/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjhchob = zf1*jphi + zf*jplo 

end function fjhchob

!**********************************************************************************************************************!
! fjch3cho -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3CHO + hv --> CH3 + HCO
!         acetaldehyde
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjch3cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjch3cho
  
  data a0/+7.448D-09,+1.429D-08,+2.167D-08,+2.877D-08,+3.462D-08,+3.746D-08,+3.692D-08/
  data a1/+3.715D-07,+1.288D-07,-3.725D-07,-1.043D-06,-1.802D-06,-2.335D-06,-2.195D-06/
  data a2/+6.601D-07,+4.398D-06,+1.177D-05,+2.305D-05,+3.922D-05,+5.962D-05,+8.454D-05/
  data a3/+1.673D-05,+1.966D-05,+1.250D-05,-4.520D-06,-3.488D-05,-7.804D-05,-0.0001334/
  data a4/-1.762D-05,-2.649D-05,-2.524D-05,-1.456D-05,+9.714D-06,+4.745D-05,+9.704D-05/
  data a5/+5.770D-06,+9.924D-06,+1.063D-05,+8.171D-06,+8.189D-07,-1.141D-05,-2.766D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjch3cho = zf1*jphi + zf*jplo 

end function fjch3cho

!**********************************************************************************************************************!
! fjc2h5cho -function to calculate photolysis frequency (1/s) for the reaction
!
!         C2H5CHO + hv --> C2H5 + HCO
!         propionaldehyde
!         (propanal)
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjc2h5cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjc2h5cho
  
  data a0/+4.819D-08,+7.950D-08,+1.040D-07,+1.200D-07,+1.232D-07,+1.125D-07,+9.759D-08/
  data a1/+1.598D-06,+3.400D-07,-1.303D-06,-2.712D-06,-2.840D-06,-1.236D-07,+7.272D-06/
  data a2/+1.079D-05,+3.318D-05,+6.501D-05,+0.0001047,+0.0001476,+0.0001847,+0.0002079/
  data a3/+4.057D-05,+1.875D-05,-3.860D-05,-0.0001255,-0.0002329,-0.0003418,-0.0004333/
  data a4/-5.098D-05,-5.029D-05,-8.666D-06,+6.710D-05,+0.0001690,+0.0002802,+0.0003823/
  data a5/+1.808D-05,+2.210D-05,+1.083D-05,-1.342D-05,-4.804D-05,-8.758D-05,-0.0001257/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjc2h5cho = zf1*jphi + zf*jplo 

end function fjc2h5cho

!**********************************************************************************************************************!
! fjc3h7choa -function to calculate photolysis frequency (1/s) for the reaction
!
!         n-C3H7CHO + hv --> n-C3H7 + HCO
!         n-butyraldehyde
!         (n-butanal)
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjc3h7choa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjc3h7choa
  
  data a0/+6.362D-08,+8.398D-08,+9.039D-08,+8.724D-08,+7.663D-08,+6.238D-08,+5.170D-08/
  data a1/+1.620D-06,+3.794D-07,-4.969D-07,-6.757D-07,+2.742D-07,+2.872D-06,+7.450D-06/
  data a2/+1.624D-05,+3.674D-05,+5.602D-05,+7.306D-05,+8.532D-05,+8.897D-05,+8.094D-05/
  data a3/+1.112D-05,-2.710D-05,-7.339D-05,-0.0001213,-0.0001639,-0.0001910,-0.0001932/
  data a4/-2.675D-05,-3.826D-08,+4.166D-05,+8.947D-05,+0.0001359,+0.0001704,+0.0001819/
  data a5/+1.146D-05,+4.686D-06,-8.779D-06,-2.545D-05,-4.258D-05,-5.624D-05,-6.216D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjc3h7choa = zf1*jphi + zf*jplo 

end function fjc3h7choa

!**********************************************************************************************************************!
! fjc3h7chob -function to calculate photolysis frequency (1/s) for the reaction
!
!         n-C3H7CHO + hv --> CH3CHO + CH2=CH2
!         n-butyraldehyde
!         (n-butanal)
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjc3h7chob(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjc3h7chob
  
  data a0/+3.028D-08,+3.983D-08,+4.299D-08,+4.166D-08,+3.664D-08,+2.985D-08,+2.463D-08/
  data a1/+7.698D-07,+1.974D-07,-2.312D-07,-3.335D-07,+1.157D-07,+1.353D-06,+3.546D-06/
  data a2/+7.752D-06,+1.733D-05,+2.663D-05,+3.491D-05,+4.077D-05,+4.250D-05,+3.856D-05/
  data a3/+5.230D-06,-1.240D-05,-3.480D-05,-5.812D-05,-7.847D-05,-9.136D-05,-9.206D-05/
  data a4/-1.265D-05,-6.451D-07,+1.968D-05,+4.309D-05,+6.529D-05,+8.160D-05,+8.667D-05/
  data a5/+5.415D-06,+2.501D-06,-4.119D-06,-1.233D-05,-2.052D-05,-2.696D-05,-2.960D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjc3h7chob = zf1*jphi + zf*jplo 

end function fjc3h7chob

!**********************************************************************************************************************!
! fjic3h7cho -function to calculate photolysis frequency (1/s) for the reaction
!
!         i-C3H7CHO + hv --> i-C3H7 + HCO
!         isobutyraldehyde
!         (2-methyl propanal)
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjic3h7cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjic3h7cho
  
  data a0/+2.541D-07,+3.452D-07,+3.809D-07,+3.751D-07,+3.298D-07,+2.602D-07,+1.958D-07/
  data a1/+6.469D-06,+7.109D-07,-4.251D-06,-6.580D-06,-4.178D-06,+5.271D-06,+2.331D-05/
  data a2/+7.232D-05,+0.0001620,+0.0002543,+0.0003406,+0.0004091,+0.0004424,+0.0004273/
  data a3/+6.852D-05,-8.565D-05,-0.0003006,-0.0005341,-0.0007541,-0.0009143,-0.0009741/
  data a4/-0.0001404,-4.573D-05,+0.0001441,+0.0003732,+0.0006073,+0.0007970,+0.0008971/
  data a5/+5.814D-05,+3.855D-05,-2.188D-05,-0.0001009,-0.0001860,-0.0002589,-0.0003022/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjic3h7cho = zf1*jphi + zf*jplo 

end function fjic3h7cho

!**********************************************************************************************************************!
! fjmacra -function to calculate photolysis frequency (1/s) for the reaction
!
!         MACR + hv --> CH2=C(CH3) + HCO
!         methacrolein
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjmacra(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjmacra
  
  data a0/+1.313D-08,+1.477D-08,+1.434D-08,+1.313D-08,+1.176D-08,+1.137D-08,+1.379D-08/
  data a1/+2.596D-07,+1.895D-07,+2.677D-07,+5.081D-07,+9.673D-07,+1.684D-06,+2.648D-06/
  data a2/+2.874D-06,+5.569D-06,+7.099D-06,+7.589D-06,+6.704D-06,+4.131D-06,-2.571D-07/
  data a3/-3.249D-06,-1.030D-05,-1.541D-05,-1.848D-05,-1.864D-05,-1.506D-05,-7.153D-06/
  data a4/+1.250D-06,+7.973D-06,+1.352D-05,+1.741D-05,+1.863D-05,+1.635D-05,+9.760D-06/
  data a5/-6.028D-08,-2.325D-06,-4.369D-06,-5.916D-06,-6.563D-06,-6.016D-06,-3.932D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjmacra = zf1*jphi + zf*jplo 

end function fjmacra

!**********************************************************************************************************************!
! fjmacrb -function to calculate photolysis frequency (1/s) for the reaction
!
!         MACR + hv --> CH2=C(CH3)CO + H
!         methacrolein
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjmacrb(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjmacrb
  
  data a0/+1.313D-08,+1.477D-08,+1.434D-08,+1.313D-08,+1.176D-08,+1.137D-08,+1.379D-08/
  data a1/+2.596D-07,+1.895D-07,+2.677D-07,+5.081D-07,+9.673D-07,+1.684D-06,+2.648D-06/
  data a2/+2.874D-06,+5.569D-06,+7.099D-06,+7.589D-06,+6.704D-06,+4.131D-06,-2.571D-07/
  data a3/-3.249D-06,-1.030D-05,-1.541D-05,-1.848D-05,-1.864D-05,-1.506D-05,-7.153D-06/
  data a4/+1.250D-06,+7.973D-06,+1.352D-05,+1.741D-05,+1.863D-05,+1.635D-05,+9.760D-06/
  data a5/-6.028D-08,-2.325D-06,-4.369D-06,-5.916D-06,-6.563D-06,-6.016D-06,-3.932D-06/
  data a0/+1.313D-08,+1.477D-08,+1.434D-08,+1.313D-08,+1.176D-08,+1.137D-08,+1.379D-08/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjmacrb = zf1*jphi + zf*jplo 

end function fjmacrb
  
!**********************************************************************************************************************!
! fjacet -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3COCH3 + hv --> CH3CO + CH3
!         acetone
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjacet(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjacet
  
  data a0/+8.831D-11,+1.412D-10,+1.994D-10,+2.483D-10,+2.882D-10,+3.339D-10,+4.428D-10/
  data a1/+3.321D-08,+3.899D-08,+3.642D-08,+3.192D-08,+2.809D-08,+2.553D-08,+2.534D-08/
  data a2/-2.003D-07,-2.946D-07,-3.008D-07,-2.779D-07,-2.549D-07,-2.425D-07,-2.346D-07/
  data a3/+1.275D-06,+1.904D-06,+2.172D-06,+2.329D-06,+2.509D-06,+2.788D-06,+3.598D-06/
  data a4/-6.661D-07,-1.383D-06,-1.728D-06,-1.943D-06,-2.180D-06,-2.501D-06,-3.479D-06/
  data a5/+9.562D-08,+3.419D-07,+4.707D-07,+5.532D-07,+6.420D-07,+7.532D-07,+1.117D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjacet = zf1*jphi + zf*jplo 

end function fjacet

!**********************************************************************************************************************!
! fjmek  -function to calculate photolysis frequency (1/s) for the reaction
!
!         MEK + hv --> CH3CO + CH2CH3
!         methyl ethyl ketone
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjmek(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjmek
  
  data a0/+1.789D-08,+2.900D-08,+3.756D-08,+4.337D-08,+4.513D-08,+4.159D-08,+3.669D-08/
  data a1/+5.801D-07,+1.167D-07,-4.939D-07,-1.083D-06,-1.343D-06,-7.477D-07,+1.058D-06/
  data a2/+4.065D-06,+1.205D-05,+2.324D-05,+3.739D-05,+5.324D-05,+6.748D-05,+7.910D-05/
  data a3/+1.533D-05,+8.566D-06,-9.999D-06,-3.910D-05,-7.630D-05,-0.0001141,-0.0001497/
  data a4/-1.955D-05,-2.049D-05,-7.915D-06,+1.674D-05,+5.123D-05,+8.825D-05,+0.0001255/
  data a5/+7.013D-06,+8.869D-06,+5.722D-06,-2.037D-06,-1.363D-05,-2.649D-05,-3.997D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjmek = zf1*jphi + zf*jplo 

end function fjmek

!**********************************************************************************************************************!
! fjmvka  -function to calculate photolysis frequency (1/s) for the reaction
!
!         MVK + hv --> CH3CH=CH2 + CO
!         methyl vinyl ketone
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjmvka(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjmvka
  
  data a0/+1.826D-08,+2.205D-08,+2.255D-08,+2.108D-08,+1.850D-08,+1.610D-08,+1.628D-08/
  data a1/+4.065D-07,+2.043D-07,+1.444D-07,+3.082D-07,+7.905D-07,+1.693D-06,+3.043D-06/
  data a2/+4.191D-06,+8.727D-06,+1.237D-05,+1.483D-05,+1.556D-05,+1.376D-05,+9.037D-06/
  data a3/-1.659D-06,-1.181D-05,-2.184D-05,-3.021D-05,-3.545D-05,-3.548D-05,-2.890D-05/
  data a4/-1.789D-06,+6.771D-06,+1.673D-05,+2.582D-05,+3.249D-05,+3.449D-05,+3.019D-05/
  data a5/+1.192D-06,-1.427D-06,-4.894D-06,-8.238D-06,-1.089D-05,-1.199D-05,-1.091D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjmvka = zf1*jphi + zf*jplo 

end function fjmvka

!**********************************************************************************************************************!
! fjmvkb  -function to calculate photolysis frequency (1/s) for the reaction
!
!         MVK + hv --> CH3CO + CH2=CH2
!         methyl vinyl ketone
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjmvkb(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjmvkb
  
  data a0/+1.826D-08,+2.205D-08,+2.255D-08,+2.108D-08,+1.850D-08,+1.610D-08,+1.628D-08/
  data a1/+4.065D-07,+2.043D-07,+1.444D-07,+3.082D-07,+7.905D-07,+1.693D-06,+3.043D-06/
  data a2/+4.191D-06,+8.727D-06,+1.237D-05,+1.483D-05,+1.556D-05,+1.376D-05,+9.037D-06/
  data a3/-1.659D-06,-1.181D-05,-2.184D-05,-3.021D-05,-3.545D-05,-3.548D-05,-2.890D-05/
  data a4/-1.789D-06,+6.771D-06,+1.673D-05,+2.582D-05,+3.249D-05,+3.449D-05,+3.019D-05/
  data a5/+1.192D-06,-1.427D-06,-4.894D-06,-8.238D-06,-1.089D-05,-1.199D-05,-1.091D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjmvkb = zf1*jphi + zf*jplo 

end function fjmvkb

!**********************************************************************************************************************!
! fjglyoxa  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CHOCHO + hv --> 2CO + H2
!         glyoxal
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjglyoxa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjglyoxa
  
  data a0/+3.772D-08,+4.850D-08,+5.171D-08,+4.976D-08,+4.393D-08,+3.658D-08,+3.226D-08/
  data a1/+9.435D-07,+3.161D-07,-1.007D-07,-9.776D-08,+5.753D-07,+2.203D-06,+4.926D-06/
  data a2/+9.010D-06,+2.011D-05,+3.042D-05,+3.916D-05,+4.481D-05,+4.524D-05,+3.904D-05/
  data a3/+4.647D-06,-1.660D-05,-4.202D-05,-6.731D-05,-8.835D-05,-9.947D-05,-9.616D-05/
  data a4/-1.296D-05,+2.298D-06,+2.575D-05,+5.142D-05,+7.485D-05,+9.002D-05,+9.177D-05/
  data a5/+5.671D-06,+1.699D-06,-6.046D-06,-1.512D-05,-2.385D-05,-3.001D-05,-3.160D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjglyoxa = zf1*jphi + zf*jplo 

end function fjglyoxa

!**********************************************************************************************************************!
! fjglyoxb  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CHOCHO + hv --> HCHO + CO
!         glyoxal
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjglyoxb(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjglyoxb
  
  data a0/+2.399D-07,+2.851D-07,+2.899D-07,+2.758D-07,+2.501D-07,+2.329D-07,+2.597D-07/
  data a1/+5.193D-06,+3.060D-06,+3.133D-06,+5.781D-06,+1.245D-05,+2.434D-05,+4.159D-05/
  data a2/+5.992D-05,+0.0001186,+0.0001614,+0.0001910,+0.0001985,+0.0001748,+0.0001159/
  data a3/-2.728D-05,-0.0001594,-0.0002785,-0.0003812,-0.0004442,-0.0004429,-0.0003623/
  data a4/-1.874D-05,+9.276D-05,+0.0002101,+0.0003220,+0.0004023,+0.0004256,+0.0003733/
  data a5/+1.386D-05,-2.032D-05,-6.089D-05,-0.0001021,-0.0001338,-0.0001469,-0.0001338/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjglyoxb = zf1*jphi + zf*jplo 

end function fjglyoxb

!**********************************************************************************************************************!
! fjglyoxc  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CHOCHO + hv --> HCO + HCO
!         glyoxal
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjglyoxc(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjglyoxc
  
  data a0/+9.520D-07,+9.648D-07,+9.245D-07,+9.023D-07,+9.475D-07,+1.175D-06,+1.795D-06/
  data a1/+1.699D-05,+2.543D-05,+4.539D-05,+7.443D-05,+0.0001155,+0.0001696,+0.0002328/
  data a2/+0.0002703,+0.0004060,+0.0004215,+0.0003647,+0.0002258,-6.151D-06,-0.0003170/
  data a3/-0.0004279,-0.0008613,-0.0010210,-0.0010180,-0.0008326,-0.0004325,+0.0001659/
  data a4/+0.0002965,+0.0007498,+0.0009568,+0.0010120,+0.0008964,+0.0005728,+4.740D-05/
  data a5/-7.923D-05,-0.0002421,-0.0003243,-0.0003555,-0.0003279,-0.0002276,-5.423D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjglyoxc = zf1*jphi + zf*jplo 

end function fjglyoxc

!**********************************************************************************************************************!
! fjmglyox  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3COCHO + hv --> CH3CO + HCO
!         methyl glyoxal
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjmglyox(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjmglyox
  
  data a0/+1.692D-06,+1.784D-06,+1.710D-06,+1.676D-06,+1.790D-06,+2.378D-06,+4.203D-06/
  data a1/+2.801D-05,+3.444D-05,+7.028D-05,+0.0001337,+0.0002380,+0.0003912,+0.0006030/
  data a2/+0.0004289,+0.0007404,+0.0008221,+0.0007389,+0.0004187,-0.0001883,-0.0011310/
  data a3/-0.0006435,-0.0015680,-0.0020040,-0.0021000,-0.0016820,-0.0006362,+0.0011600/
  data a4/+0.0004140,+0.0013660,+0.0018940,+0.0021230,+0.0018670,+0.0010150,-0.0005633/
  data a5/-0.0001004,-0.0004404,-0.0006456,-0.0007539,-0.0006952,-0.0004292,+9.303D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjmglyox = zf1*jphi + zf*jplo 

end function fjmglyox

!**********************************************************************************************************************!
! fjbiacet  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3COCOCH3 + hv --> 2CH3CO
!         biacetyl
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjbiacet(cz, cz2, cz3, cz4, cz5, cz6, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5, cz6
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5, a6
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjbiacet
  
  data a0/+5.075D-06,+5.195D-06,+5.306D-06,+5.783D-06,+7.026D-06,+1.005D-05,+1.643D-05/
  data a1/+4.451D-05,+0.0001043,+0.0002929,+0.0005467,+0.0008658,+0.0012380,+0.0016060/
  data a2/+0.0022410,+0.0035970,+0.0032100,+0.0019680,+1.457D-05,-0.0025860,-0.0054320/
  data a3/-0.0058280,-0.0126300,-0.0132300,-0.0106300,-0.0055040,+0.0020070,+0.0107100/
  data a4/+0.0070570,+0.0193300,+0.0218300,+0.0189800,+0.0120300,+0.0010860,-0.0120700/
  data a5/-0.0042580,-0.0141700,-0.0167800,-0.0152000,-0.0104400,-0.0025430,+0.0072110/
  data a6/+0.0010290,+0.0040440,+0.0049480,+0.0045960,+0.0033010,+0.0010590,-0.0017670/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5 + a6(l)*cz6
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5+a6(lm1)*cz6
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjbiacet = zf1*jphi + zf*jplo 

end function fjbiacet

!**********************************************************************************************************************!
! fjch3ooh  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3OOH + hv --> CH3O + OH
!         methyl hydroperoxide
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjch3ooh(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjch3ooh
  
  data a0/+4.346D-08,+5.243D-08,+5.337D-08,+4.990D-08,+4.396D-08,+3.882D-08,+3.985D-08/
  data a1/+9.772D-07,+4.985D-07,+3.891D-07,+8.013D-07,+1.964D-06,+4.102D-06,+7.289D-06/
  data a2/+9.647D-06,+2.035D-05,+2.866D-05,+3.424D-05,+3.572D-05,+3.134D-05,+2.002D-05/
  data a3/-3.267D-06,-2.719D-05,-5.005D-05,-6.913D-05,-8.081D-05,-8.036D-05,-6.406D-05/
  data a4/-4.562D-06,+1.567D-05,+3.828D-05,+5.905D-05,+7.402D-05,+7.812D-05,+6.695D-05/
  data a5/+2.866D-06,-3.356D-06,-1.121D-05,-1.886D-05,-2.481D-05,-2.718D-05,-2.421D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjch3ooh = zf1*jphi + zf*jplo 

end function fjch3ooh

!**********************************************************************************************************************!
! fjch3no3  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3NO3 + hv --> CH3O + NO2
!         methyl nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjch3no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjch3no3
  
  data a0/+2.286D-09,+2.839D-09,+2.907D-09,+2.693D-09,+2.302D-09,+1.821D-09,+1.463D-09/
  data a1/+8.414D-08,+2.690D-08,-2.458D-08,-5.720D-08,-6.521D-08,-4.499D-08,+4.077D-09/
  data a2/+4.113D-07,+1.014D-06,+1.607D-06,+2.105D-06,+2.461D-06,+2.639D-06,+2.779D-06/
  data a3/+2.225D-06,+1.468D-06,+5.243D-08,-1.416D-06,-2.707D-06,-3.651D-06,-4.396D-06/
  data a4/-2.547D-06,-2.446D-06,-1.292D-06,+9.875D-08,+1.433D-06,+2.508D-06,+3.359D-06/
  data a5/+8.724D-07,+9.756D-07,+6.356D-07,+1.739D-07,-2.947D-07,-6.944D-07,-1.009D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjch3no3 = zf1*jphi + zf*jplo 

end function fjch3no3

!**********************************************************************************************************************!
! fjc2h5no3  -function to calculate photolysis frequency (1/s) for the reaction
!
!         C2H5NO3 + hv --> C2H5O + NO2
!         ethyl nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjc2h5no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjc2h5no3
  
  data a0/+3.996D-09,+4.867D-09,+4.987D-09,+4.517D-09,+3.808D-09,+2.949D-09,+2.388D-09/
  data a1/+1.405D-07,+4.557D-08,-4.460D-08,-8.914D-08,-9.693D-08,-5.359D-08,+3.146D-08/
  data a2/+7.566D-07,+1.772D-06,+2.809D-06,+3.550D-06,+4.093D-06,+4.282D-06,+4.478D-06/
  data a3/+3.476D-06,+2.108D-06,-4.739D-07,-2.722D-06,-4.799D-06,-6.098D-06,-7.300D-06/
  data a4/-4.095D-06,-3.789D-06,-1.565D-06,+5.563D-07,+2.755D-06,+4.256D-06,+5.689D-06/
  data a5/+1.423D-06,+1.542D-06,+8.471D-07,+1.534D-07,-6.332D-07,-1.185D-06,-1.735D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjc2h5no3 = zf1*jphi + zf*jplo 

end function fjc2h5no3

!**********************************************************************************************************************!
! fjnc3h7no3  -function to calculate photolysis frequency (1/s) for the reaction
!
!         n-C3H7NO3 + hv --> n-C3H7O + NO2
!         n-propyl nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjnc3h7no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjnc3h7no3
  
  data a0/+7.852D-09,+1.041D-08,+1.137D-08,+1.125D-08,+1.017D-08,+8.660D-09,+7.643D-09/
  data a1/+2.370D-07,+9.642D-08,-1.791D-08,-5.960D-08,+3.815D-08,+3.325D-07,+8.602D-07/
  data a2/+1.591D-06,+3.923D-06,+6.297D-06,+8.527D-06,+1.018D-05,+1.085D-05,+1.020D-05/
  data a3/+3.145D-06,-5.498D-07,-5.861D-06,-1.176D-05,-1.701D-05,-2.056D-05,-2.131D-05/
  data a4/-4.461D-06,-2.378D-06,+2.177D-06,+7.866D-06,+1.334D-05,+1.752D-05,+1.916D-05/
  data a5/+1.672D-06,+1.275D-06,-1.534D-07,-2.103D-06,-4.063D-06,-5.658D-06,-6.410D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjnc3h7no3 = zf1*jphi + zf*jplo 

end function fjnc3h7no3

!**********************************************************************************************************************!
! fjic3h7no3  -function to calculate photolysis frequency (1/s) for the reaction
!
!         i-C3H7NO3 + hv --> i-C3H7O + NO2
!         isopropyl nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjic3h7no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjic3h7no3
  
  data a0/+7.592D-09,+9.173D-09,+9.096D-09,+8.186D-09,+6.764D-09,+5.255D-09,+4.135D-09/
  data a1/+2.505D-07,+6.836D-08,-7.116D-08,-1.479D-07,-1.423D-07,-5.650D-08,+1.300D-07/
  data a2/+1.540D-06,+3.516D-06,+5.188D-06,+6.504D-06,+7.290D-06,+7.532D-06,+7.535D-06/
  data a3/+5.474D-06,+2.333D-06,-1.928D-06,-6.144D-06,-9.477D-06,-1.170D-05,-1.297D-05/
  data a4/-6.834D-06,-5.438D-06,-1.801D-06,+2.394D-06,+6.015D-06,+8.760D-06,+1.041D-05/
  data a5/+2.435D-06,+2.321D-06,+1.214D-06,-2.279D-07,-1.536D-06,-2.596D-06,-3.231D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjic3h7no3 = zf1*jphi + zf*jplo 

end function fjic3h7no3

!**********************************************************************************************************************!
! fjtc4h9no3  -function to calculate photolysis frequency (1/s) for the reaction
!
!         t-C4H9NO3 + hv --> t-C4H9O + NO2
!         t-butyl nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjtc4h9no3(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjtc4h9no3
  
  data a0/+1.894D-08,+2.730D-08,+3.145D-08,+3.208D-08,+2.942D-08,+2.407D-08,+1.805D-08/
  data a1/+5.860D-07,+6.299D-08,-4.678D-07,-8.393D-07,-8.849D-07,-3.961D-07,+7.858D-07/
  data a2/+4.954D-06,+1.231D-05,+2.067D-05,+2.924D-05,+3.700D-05,+4.251D-05,+4.449D-05/
  data a3/+1.365D-05,+3.946D-06,-1.358D-05,-3.486D-05,-5.697D-05,-7.608D-05,-8.822D-05/
  data a4/-1.864D-05,-1.535D-05,-1.258D-06,+1.850D-05,+4.066D-05,+6.134D-05,+7.604D-05/
  data a5/+6.886D-06,+7.104D-06,+2.952D-06,-3.614D-06,-1.138D-05,-1.896D-05,-2.463D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjtc4h9no3 = zf1*jphi + zf*jplo 

end function fjtc4h9no3

!**********************************************************************************************************************!
! fjnoaa  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3C(O)CH2NO3 + hv --> CH3C(O)CH2O + NO2
!         1-hydroxy-2-propanone nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjnoaa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjnoaa
  
  data a0/+1.837D-08,+2.530D-08,+2.833D-08,+2.844D-08,+2.590D-08,+2.158D-08,+1.753D-08/
  data a1/+5.392D-07,+1.260D-07,-2.556D-07,-4.712D-07,-3.687D-07,+2.203D-07,+1.408D-06/
  data a2/+4.343D-06,+1.058D-05,+1.731D-05,+2.390D-05,+2.942D-05,+3.269D-05,+3.276D-05/
  data a3/+9.562D-06,+4.236D-07,-1.420D-05,-3.107D-05,-4.758D-05,-6.061D-05,-6.709D-05/
  data a4/-1.371D-05,-9.472D-06,+2.703D-06,+1.864D-05,+3.551D-05,+5.009D-05,+5.891D-05/
  data a5/+5.136D-06,+4.700D-06,+9.973D-07,-4.365D-06,-1.034D-05,-1.578D-05,-1.934D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjnoaa = zf1*jphi + zf*jplo 

end function fjnoaa

!**********************************************************************************************************************!
! fjnoab  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3C(O)CH2NO3 + hv --> CH3CO + HCHO + NO2
!         1-hydroxy-2-propanone nitrate
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjnoab(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjnoab
  
  data a0/+5.205D-09,+7.487D-09,+8.887D-09,+9.317D-09,+9.013D-09,+8.114D-09,+6.975D-09/
  data a1/+2.126D-07,+9.487D-08,-5.808D-08,-1.696D-07,-2.109D-07,-1.226D-07,+1.404D-07/
  data a2/+7.017D-07,+2.343D-06,+4.590D-06,+6.925D-06,+9.215D-06,+1.104D-05,+1.208D-05/
  data a3/+6.389D-06,+5.452D-06,+1.336D-06,-3.902D-06,-9.835D-06,-1.531D-05,-1.926D-05/
  data a4/-6.803D-06,-7.907D-06,-5.020D-06,-5.569D-07,+4.999D-06,+1.048D-05,+1.473D-05/
  data a5/+2.239D-06,+3.022D-06,+2.271D-06,+8.760D-07,-9.802D-07,-2.896D-06,-4.440D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjnoab = zf1*jphi + zf*jplo 

end function fjnoab

!**********************************************************************************************************************!
! fjhoch2cho  -function to calculate photolysis frequency (1/s) for the reaction
!
!         HOCH2CHO + hv --> HCO + HCHO + HO2
!         glycolaldehyde
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjhoch2cho(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjhoch2cho
  
  data a0/+1.258D-08,+1.810D-08,+2.114D-08,+2.204D-08,+2.080D-08,+1.787D-08,+1.448D-08/
  data a1/+4.319D-07,+1.089D-07,-2.510D-07,-5.272D-07,-5.993D-07,-3.359D-07,+3.536D-07/
  data a2/+2.591D-06,+7.080D-06,+1.255D-05,+1.836D-05,+2.375D-05,+2.781D-05,+2.990D-05/
  data a3/+1.207D-05,+7.491D-06,-3.255D-06,-1.704D-05,-3.174D-05,-4.484D-05,-5.418D-05/
  data a4/-1.456D-05,-1.459D-05,-6.506D-06,+5.829D-06,+2.012D-05,+3.380D-05,+4.440D-05/
  data a5/+5.097D-06,+6.055D-06,+3.815D-06,-1.723D-07,-5.075D-06,-9.985D-06,-1.395D-05/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjhoch2cho = zf1*jphi + zf*jplo 

end function fjhoch2cho

!**********************************************************************************************************************!
! fjacetol  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH2(OH)COCH3 + hv --> CH3CO + HCHO + HO2
!         hydroxyacetone
!         (acetol)
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjacetol(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjacetol
  
  data a0/+3.901D-09,+5.227D-09,+5.750D-09,+5.636D-09,+5.070D-09,+4.203D-09,+3.485D-09/
  data a1/+1.178D-07,+4.340D-08,-2.088D-08,-4.239D-08,-1.423D-09,+1.449D-07,+4.182D-07/
  data a2/+8.724D-07,+2.088D-06,+3.355D-06,+4.500D-06,+5.440D-06,+5.845D-06,+5.518D-06/
  data a3/+1.490D-06,-4.314D-07,-3.284D-06,-6.276D-06,-9.241D-06,-1.125D-05,-1.166D-05/
  data a4/-2.151D-06,-1.067D-06,+1.411D-06,+4.261D-06,+7.375D-06,+9.735D-06,+1.061D-05/
  data a5/+8.111D-07,+6.041D-07,-1.853D-07,-1.148D-06,-2.276D-06,-3.181D-06,-3.575D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjacetol = zf1*jphi + zf*jplo 

end function fjacetol

!**********************************************************************************************************************!
! fjpaa  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3C(O)OOH + hv --> CH3O2 + OH
!         peroxyacetic acid
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjpaa(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjpaa
  
  data a0/+4.011D-09,+5.221D-09,+5.620D-09,+5.436D-09,+4.805D-09,+3.977D-09,+3.399D-09/
  data a1/+1.047D-07,+3.535D-08,-1.631D-08,-2.385D-08,+3.932D-08,+2.053D-07,+4.931D-07/
  data a2/+9.502D-07,+2.146D-06,+3.309D-06,+4.320D-06,+5.022D-06,+5.173D-06,+4.600D-06/
  data a3/+7.238D-07,-1.462D-06,-4.262D-06,-7.109D-06,-9.580D-06,-1.103D-05,-1.091D-05/
  data a4/-1.580D-06,-9.187D-08,+2.452D-06,+5.304D-06,+8.006D-06,+9.870D-06,+1.029D-05/
  data a5/+6.626D-07,+2.974D-07,-5.339D-07,-1.533D-06,-2.531D-06,-3.273D-06,-3.526D-06/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjpaa = zf1*jphi + zf*jplo 

end function fjpaa

!**********************************************************************************************************************!
! fjpyra  -function to calculate photolysis frequency (1/s) for the reaction
!
!         CH3COCOOH + hv --> CH3O3 + HO2
!         pyruvic acid
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjpyra(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjpyra
  
  data a0/+4.807D-06,+5.286D-06,+5.046D-06,+4.582D-06,+4.163D-06,+4.175D-06,+5.319D-06/
  data a1/+9.071D-05,+7.254D-05,+0.0001133,+0.0002153,+0.0003931,+0.0006618,+0.0010120/
  data a2/+0.0010350,+0.0019660,+0.0024170,+0.0024500,+0.0019970,+0.0009242,-0.0007661/
  data a3/-0.0013520,-0.0039060,-0.0055580,-0.0063240,-0.0060470,-0.0043530,-0.0011570/
  data a4/+0.0006754,+0.0031900,+0.0050480,+0.0061250,+0.0062440,+0.0049960,+0.0022380/
  data a5/-0.0001042,-0.0009717,-0.0016710,-0.0021160,-0.0022410,-0.0018880,-0.0009910/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjpyra = zf1*jphi + zf*jplo 

end function fjpyra

!**********************************************************************************************************************!
! fjhpald  -function to calculate photolysis frequency (1/s) for the reaction
!
!         HPALD + hv --> Products
!         hydroperoxyaldehydes
!
!         for a given cosine(zenith angle) and altitude above mean sea level.
!         fourth-order polynomials of cosine(zenith angle) at each of 7 levels
!         1=0km, 2=2km, 3=4km, 4=6km, 5=8km, 6=10km, 7=12km
!         Derived from TUV5.0, http://http://cprm.acd.ucar.edu/Models/TUV/
!         (Madronich and Flocke, 1997)
!**********************************************************************************************************************!
function fjhpald(cz, cz2, cz3, cz4, cz5, l, lm1, zf, zf1)
  integer(kind=i4) :: l, lm1
  real(kind=dp)    :: cz, cz2, cz3, cz4, cz5
  real(kind=dp)    :: zf, zf1
  integer(kind=i4), parameter :: nalts=7
  real(kind=dp), dimension(nalts) :: a0, a1, a2, a3, a4, a5
  real(kind=dp)    :: jphi, jplo
  real(kind=dp)    :: fjhpald
  
  data a0/+6.571D-06,+7.390D-06,+7.174D-06,+6.558D-06,+5.893D-06,+5.685D-06,+6.906D-06/
  data a1/+0.0001295,+9.449D-05,+0.0001333,+0.0002545,+0.0004821,+0.0008420,+0.0013230/
  data a2/+0.0014390,+0.0027870,+0.0035560,+0.0037920,+0.0033670,+0.0020660,-0.0001210/
  data a3/-0.0016270,-0.0051620,-0.0077230,-0.0092330,-0.0093680,-0.0075380,-0.0035960/
  data a4/+0.0006247,+0.0040030,+0.0067830,+0.0086940,+0.0093740,+0.0081840,+0.0049020/
  data a5/-2.929D-05,-0.0011700,-0.0021940,-0.0029530,-0.0033070,-0.0030130,-0.0019750/

  ! polynomial value at l (highest altitude)
  jphi = a0(l) + a1(l)*cz + a2(l)*cz2 + a3(l)*cz3 + a4(l)*cz4 + a5(l)*cz5
  jphi = max(0.0D0, jphi)   ! no negative values

  ! polynomial value at lm1 (lowest altitude)
  jplo = a0(lm1)+a1(lm1)*cz+a2(lm1)*cz2+a3(lm1)*cz3+a4(lm1)*cz4+a5(lm1)*cz5
  jplo = max(0.0D0, jplo)   ! no negative values

  ! interpolate to required altitude
  fjhpald = zf1*jphi + zf*jplo 

end function fjhpald

!**********************************************************************************************************************!
! GasThermalRateCoeffs - calculates thermal rate coefficients
!**********************************************************************************************************************!
subroutine GasThermalRateCoeffs(ts, k)
  real(kind=dp), dimension(nrxn) :: k
  real(kind=dp) :: ts  ! not used yet
  real(kind=dp) :: TEMP, AIR, N2, O2, H2Oz, RO2

  ! thermal reaction rate coefficients (units: molecules-cm3-s)
  ! Generated by chemgen.py from mechanism anthro_isop

  TEMP = tk(zpt)
  AIR  = cair(zpt)
  N2 = 0.79*cair(zpt)
  O2 = 0.21*cair(zpt)
  H2Oz = h2o(zpt) 

   ! Dummy routine for ACCESS_NH3

  return
end subroutine GasThermalRateCoeffs

!**********************************************************************************************************************!
!  Reaction rate coefficient functions
!**********************************************************************************************************************!
!  Arrhenius equation with temperature dependence of A0
function ARR(A0,E0,B0,TEMP)
   real(kind=dp), intent(in) :: A0, E0, B0, TEMP  ! temp (K)
   real(kind=dp) :: ARR
   ARR = A0*exp(E0/TEMP)*(TEMP/300.0_dp)**B0
   return
end function ARR

!  Arrhenius equation
function ARR2(A0,E0,TEMP)
   real(kind=dp), intent(in) :: A0, E0, TEMP      ! temp (K)
   real(kind=dp) :: ARR2
   ARR2 = A0*exp(E0/TEMP)
   return
end function ARR2

!  full termolecular rate expression
function TERM(A0,B0,E0,Ai,Bi,Ei,F,TEMP,AIR)
   real(kind=dp), intent(in) :: A0, B0, E0        ! Arrhenius parameters for the low
                                                  ! pressure limit
   real(kind=dp), intent(in) :: Ai, Bi, Ei        ! Arrhenius parameters for the high
                                                  ! pressure limit
   real(kind=dp), intent(in) :: F                 ! base and exponent parameters
   real(kind=dp), intent(in) :: TEMP, AIR         ! temp (K) and cair (molec/cc)
   real(kind=dp) :: TERM
   real(kind=dp)    :: k0, ki, rk, logrk, g

   k0 = ARR(A0,E0,B0,TEMP)
   ki = ARR(Ai,Ei,Bi,TEMP)
   rk = k0*AIR/ki
   logrk = DLOG10(rk)
   g = 1.0_dp/(1.0_dp + (logrk*logrk))
   TERM = (k0*AIR/(1.0_dp + rk))*F**g
   return
end function TERM

!  JPL-style termolecular rate expression
function TJPL(A0,B0,Ai,Bi,TEMP,AIR)
   real(kind=dp), intent(in) :: A0, B0        ! Arrhenius parameters for the low
                                              ! pressure limit
   real(kind=dp), intent(in) :: Ai, Bi        ! Arrhenius parameters for the high
                                              ! pressure limit
   real(kind=dp), intent(in) :: TEMP, AIR     ! temp (K) and cair (molec/cc)
   real(kind=dp)             :: F             ! base and exponent parameters
   real(kind=dp)             :: TJPL
   real(kind=dp)             :: k0, ki, rk, logrk, g
   F = 0.6_dp
   k0 = ARR(A0,0.0_dp,B0,TEMP)
   ki = ARR(Ai,0.0_dp,Bi,TEMP)
   rk = k0*AIR/ki
   logrk = DLOG10(rk)
   g = 1.0_dp/(1.0_dp + (logrk*logrk))
   TJPL = (k0*AIR/(1.0_dp + rk))*F**g
   return
end function TJPL

!  JPL-style termolecular rate expression for CO + OH = H + CO2
function TERCO(A0,B0,Ai,Bi,TEMP,AIR)
   real(kind=dp), intent(in) :: A0, B0       ! Arrhenius parameters for the low
                                             ! pressure limit
   real(kind=dp), intent(in) :: Ai, Bi       ! Arrhenius parameters for the high
                                             ! pressure limit
   real(kind=dp), intent(in) :: TEMP, AIR    ! temp (K) and cair (molec/cc)
   real(kind=dp)             :: F            ! base and exponent parameters
   real(kind=dp)             :: TERCO
   real(kind=dp)             :: k0, ki, rk, logrk, g
   F = 0.6_dp
   k0 = ARR(A0,0.0_dp,B0,TEMP)
   ki = ARR(Ai,0.0_dp,Bi,TEMP)
   rk = k0/(ki/AIR)
   logrk = DLOG10(rk)
   g = 1.0_dp/(1.0_dp + (logrk*logrk))
   TERCO = (k0/(1.0_dp + rk))*F**g
   return
end function TERCO

!  special rate expression function type 2
function KTYP2(A0,E0,A2,E2,A3,E3,TEMP,AIR)
   real(kind=dp), intent(in) :: A0, E0
   real(kind=dp), intent(in) :: A2, E2
   real(kind=dp), intent(in) :: A3, E3
   real(kind=dp), intent(in) :: TEMP, AIR          ! temp (K) and cair (molec/cc)
   real(kind=dp) :: k0, k2, k3, k1
   real(kind=dp) :: KTYP2
   k0 = ARR(A0,E0,0.0_dp,TEMP)
   k2 = ARR(A2,E2,0.0_dp,TEMP)
   k3 = ARR(A3,E3,0.0_dp,TEMP)
   k1 = k3*AIR/(1.0_dp+(k3*AIR/k2)) 
   KTYP2 = k0 + k1
   return
end function KTYP2

!  special rate expression function type 3
function KTYP3(A1,E1,A2,E2,TEMP,AIR)
   real(kind=dp), intent(in) :: A1, E1
   real(kind=dp), intent(in) :: A2, E2
   real(kind=dp), intent(in) :: TEMP, AIR         ! temp (K) and cair (molec/cc)
   real(kind=dp) :: k1, k2
   real(kind=dp) :: KTYP3
   k1 = ARR(A1,E1,0.0_dp,TEMP)
   k2 = ARR(A2,E2,0.0_dp,TEMP)
   KTYP3 = k1 + k2*AIR
   return
end function KTYP3

!  special rate expression function for HO2 + HO2 = H2O2
function KTYPHO2(A0,E0,A1,E1,A2,E2,TEMP,AIR,H2O)
   real(kind=dp), intent(in) :: A0, E0
   real(kind=dp), intent(in) :: A1, E1
   real(kind=dp), intent(in) :: A2, E2
   real(kind=dp), intent(in) :: TEMP, AIR, H2O       ! temp (K), cair (molec/cc) & h2o (molec/cc)
   real(kind=dp) :: ka, kb, kh2o
   real(kind=dp) :: KTYPHO2
   ka = A0*EXP(E0/TEMP)
   kb = A1*AIR*EXP(E1/TEMP)
   kh2o = A2*H2O*EXP(E2/TEMP)
   KTYPHO2 = (ka + kb)*(1.0_dp+kh2o)
   return
end function KTYPHO2

!  end of reaction rate coefficient functions
!**********************************************************************************************************************!

end module GasChem
