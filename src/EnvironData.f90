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
!     Module:       EnvironData                                                                                        !
!                                                                                                                      !
!     Description:  contains algorithms related to reading and analyzing environmental data                            !
!                                                                                                                      !
!======================================================================================================================!
module EnvironData
  use GlobalData
  use PhysChemData
  use CanopyPhysics
  use Utils
  implicit none

  public ReadEnvironData, AeroConductanceHeat, AeroConductanceMomen, AeroConductanceWater

  private CalcTemp, CalcPressure, CalcMeanWindSpeed, CalcCair, CalcEddyDiffRa, CalcEddyDiffStull, &
          Psi_m, Psi_h, Psi_w, CalcRiB, rav, MapRiBtoZoL, Phi

contains

!**********************************************************************************************************************!
! ReadEnvironData ... read one time step of environmental data
!
!**********************************************************************************************************************!
  subroutine ReadEnvironData()
    integer(kind=i4)          :: j, n
    integer(kind=i4), parameter :: nmhlines=17
    character(len=19)         :: sdt      ! time slice datetime as a string
    character(len=10)         :: sdate    ! time slice date as a string
    character(len=8)          :: stime    ! time slice time as a string
    real(kind=dp)             :: xaugm3   ! NH3 measured conc at zref (ug/m3)
    real(kind=dp)             :: tac      ! air temperature at zref (C)
    real(kind=dp)             :: tsc      ! surface temperature (C)
    real(kind=dp)             :: tak      ! air temperature at zref (K)
    real(kind=dp)             :: ts05k    ! soil temperature at 5 cm depth (K)
    real(kind=dp)             :: ts55k    ! soil temperature at 55 cm depth (K)
    real(kind=dp)             :: takm     ! air temperature at measurement height (K)
    real(kind=dp)             :: pmbm     ! air pressure at measurement height (mb)
    real(kind=dp)             :: umsm     ! wind speed at measurement height (m/s)
    real(kind=dp)             :: rhm      ! relative humidity at measurement height (%)
    real(kind=dp)             :: ppfdm    ! PPFD at measurement height (umol/m2-s)
    real(kind=dp)             :: sradm    ! shortwave solar radiation at measurement height (W/m2)
    real(kind=dp)             :: ustrms   ! friction velocity (m/s)
    real(kind=dp)             :: lewm2    ! latent heat flux (W/m2)
    real(kind=dp)             :: hwm2     ! sensible heat flux (W/m2)
    real(kind=dp)             :: gwm2     ! soil heat flux (W/m2)
    real(kind=dp)             :: lmo      ! Monin-Obukhov length scale, L (m)
    real(kind=dp)             :: zol      ! Monin-Obukhov stability parameter
    real(kind=dp)             :: qhkgkg   ! specific humidity in units of kg H2O/kg air
    real(kind=dp)             :: rhomolm3 ! air density (mol air/m3)
    real(kind=dp), parameter  :: umsmin=0.01   ! minimum wind speed (m/s)
    real(kind=dp)             :: time21        ! simulation time (seconds since midnight Jan 1, 2000)
    logical                   :: initcall
    data initcall /.TRUE./
    save initcall

    KVSELECT=KVRA

    ! initial call
    if (initcall) then
      open(unit=UENV, file=('./data/' // envfile))
      do j=1,nmhlines
        read(UENV,*)
      end do
      initcall = .FALSE.
      zref = zi          ! For Coweeta, domain top is the height of the topmost measurements
    end if

    ! read environmental data for the next time slice
    read(UENV,*) sdate, stime, xaugm3, takm, ts05k, ts55k, pmbm, umsm, rhm, ppfdm, sradm, &
                  ustrms, lewm2, hwm2, gwm2, stheta

    sdt = sdate // ' ' // stime
    time21 = DateTime2SimTime(sdate, stime)
    zadeg = SolarZenithAngle(time21)
    zarad = zadeg*pi/180.

    sdtout(nt) = sdt

    ! adjust measurement heights to zref
    ! For Coweeta, zref == domain top, so no adjustment necessary
    tak      = takm
    pmbzref  = pmbm
    umszref  = max(umsm, umsmin)  ! assume no height adjustment necessary
    rhzref   = rhm                ! assume no height adjustment necessary
    ppfdzref = ppfdm              ! assume no height adjustment necessary
    sradzref = sradm              ! assume no height adjustment necessary

    ! calculate water vapor pressure above the canopy, kPa
    eatm = esat(tak)*rhm/100.

    ! unit conversions
    tac    = tak - 273.15
    tsc    = ts05k - 273.15
    ubzref = umszref*100.0

    ! calculate Monin-Obukhov length scale
    rhomolm3 = pmbm/(0.08314*tak)
    qhkgkg = SpecificHumidity(rhm, takm, pmbm)*0.001
    lmo = MOLength(ustrms, takm, pmbm, rhomolm3, qhkgkg, hwm2, lewm2)

    ! calculate MO stability parameter
    zol = (zref-d)/(lmo*100.)     ! convert lmo from m to cm
    zolout(nt) = zol

    ! calculate Ra from MO stability parameter
    ra(nt) = ramol(zref*0.01, zol, umszref)*0.01     ! convert Ra from s/m to s/cm

    ! calculate aerodynamic conductance (mol/m2-s)
    gaero(nt) = AeroConductanceHeat(zref*0.01, zol, rhomolm3, umszref)

    do n=npts,1,-1
       ! Initial guesses for tk and tl (linear interpolation)
       tk(n)        = CalcTemp(z(n), tac, tsc)
       tl_sun(n)    = tk(n) 
       tl_shd(n)    = tk(n)

       ! Calculate pmb
       pmb(n)       = CalcPressure(pmbzref, z(n), tak, ts05k)

       ! Initial guess for qh based on rhzref and tk
       qh(n)        = SpecificHumidity(rhzref, tk(npts), pmb(npts))

       ! Calculate mean wind speed (Meyers et al., 1998)
       ubar(n)      = CalcMeanWindSpeed(z(n), ubzref)

       select case(KVSELECT)
         case(KVRA)
           ! Kv based on Ra
           kv(n)        = CalcEddyDiffRa(z(n), ubzref, ubar(n), ra(nt))
         case(KVSTULL)
           ! Kv from Stull (1998)
           kv(n)        = CalcEddyDiffStull(z(n), ubar(n), zol)
         case default
           kv(n)        = CalcEddyDiffRa(z(n), ubzref, ubar(n), ra(nt))
       end select

       ! Calculate air density based on pmb and tk
       ! (updated as tk is updated during integration)
       cair(n)      = CalcCair(pmb(n), tk(n))

       ! Convert specific humidity to molecules/cm3
       h2o(n)       = Convert_qh_to_h2o(qh(n), cair(n))

       caloft(iNH3) = Convertugm3ToMolecCC(xaugm3, iNH3)     ! ug/m3 --> molec/cm3
    end do 

    qzref  = h2o(npts)/navo          ! moles/cm3
    tkzref = tak
    tsoilk = ts05k                   ! 5 cm temperature
    dtdzsoil = (tsoilk-ts55k)/0.5    ! soil temperature gradient (K/m)

    ! Canopy Physics
    !  Partition measured solar radiation and PPFD into direct and diffuse components
    call PartitionRAD()

    !  Calculate radiation profiles throughout entire domain
    call CalcRadProfiles()

    !  Calculate canopy physics, including temperature profiles, photsynthetic assimilation rates,
    !  stomatal conductances, and integrate water vapor and air temperature profiles, if desired
    call CalcCanopyPhysics()

    !  Calculate sun/shade weighted canopy profiles
    call CalcWeightedProfiles()

    return
  end subroutine ReadEnvironData

!**********************************************************************************************************************!
! MOLength ... compute Monin-Obukhov length scale (m)
!
!**********************************************************************************************************************!
  function MOLength(ustr, ta, p, rho, qhd, h, le)
    real(kind=dp), intent(in)  :: ustr        ! friction velocity (m/s)
    real(kind=dp), intent(in)  :: ta          ! air temperature (K)
    real(kind=dp), intent(in)  :: p           ! air pressure (mb)
    real(kind=dp), intent(in)  :: rho         ! air density (mol/m3)
    real(kind=dp), intent(in)  :: qhd         ! specific humidity (kg/kg)
    real(kind=dp), intent(in)  :: h           ! sensible heat flux (W/m2)
    real(kind=dp), intent(in)  :: le          ! latent heat flux (W/m2)
    real(kind=dp)              :: pta         ! potential temperature (K)
    real(kind=dp)              :: vta         ! virtual temperature (K)
    real(kind=dp)              :: eh2o        ! water vapor flux (mol/m2-s)
    real(kind=dp)              :: hv          ! virtual sensible heat flux (W/m2)
    real(kind=dp), parameter   :: g=9.822     ! gravitational acceleration (m/s2)
    real(kind=dp)              :: MOLength    ! Monin-Obukhov length scale (m)

    pta = ta*(1000./p)**0.286
    vta = pta*(1.0 + 0.61*qhd)
    hv  = h + 0.61*cpair*pta*(le/lambda(ta))
 
    MOLength = - (ustr*ustr*ustr*cpair*rho*vta)/(0.4*g*hv)

    return

  end function MOLength

!**********************************************************************************************************************!
! ramol ... compute the aerodynamic resistance based on the Monin-Obukhov length
!            scale, L
!
!**********************************************************************************************************************!
  function ramol(zm, zol, ums)
    real(kind=dp), intent(in)  :: zm           ! height (m)
    real(kind=dp), intent(in)  :: zol          ! (z-d)/L, MO-stability parameter
    real(kind=dp), intent(in)  :: ums          ! mean wind velocity at height zm (m/s)
    real(kind=dp), parameter   :: vk=0.4       ! von Karman constant
    real(kind=dp)              :: zcm          ! height in cm
    real(kind=dp)              :: lnz0m, lnz0h, fm, fh
    real(kind=dp)              :: ramol        ! aerodynamic resistance (s/m)
    real(kind=dp), parameter   :: ramin=0.01   ! minimum Ra (s/m)
    real(kind=dp), parameter   :: ramax=1000.  ! maximum Ra (s/m)

    zcm = 100.0*zm   ! convert to cm since d and z0m, z0h are in cm
    lnz0m = log((zcm-d)/z0m)
    lnz0h = log((zcm-d)/z0h)
    fm = lnz0m - Psi_m(zol)
    fh = lnz0h - Psi_h(zol)
    ramol = fm*fh/(vk*vk*ums)
    ramol = min(max(ramol, ramin), ramax)

    return

  end function ramol
!**********************************************************************************************************************!
! AeroConductanceHeat ... compute the aerodynamic conductance for heat (mol/m2-s)
!                         Bonan (2016) p. 215
!
!**********************************************************************************************************************!
  function AeroConductanceHeat(zm, zol, rhoair, ums)
    real(kind=dp), intent(in)  :: zm           ! height (m)
    real(kind=dp), intent(in)  :: zol          ! (z-d)/L, MO-stability parameter
    real(kind=dp), intent(in)  :: rhoair       ! air density at height zm (mol/m3)
    real(kind=dp), intent(in)  :: ums          ! mean wind velocity at height zm (m/s)
    real(kind=dp), parameter   :: vk=0.4       ! von Karman constant
    real(kind=dp), parameter   :: gahmin=0.5   ! minimum value of AeroConductanceHeat
    real(kind=dp), parameter   :: gahmax=10.0  ! maximum value of AeroDonductanceHeat
    real(kind=dp)              :: zcm          ! height in cm
    real(kind=dp)              :: gah          ! AeroConductanceHeat (mol/m2-s)
    real(kind=dp)              :: lnz0m, lnz0h, fm, fh
    real(kind=dp)              :: AeroConductanceHeat  ! (mol/m2-s)

    zcm = 100.0*zm  ! convert to cm since d and z0m, z0h are in cm
    lnz0m = log((zcm-d)/z0m)
    lnz0h = log((zcm-d)/z0h)
    fm = lnz0m - Psi_m(zol)
    fh = lnz0h - Psi_h(zol)
    gah = rhoair*vk*vk*ums/(fm*fh)
    AeroConductanceHeat = min(max(gah, gahmin), gahmax)

    return

  end function AeroConductanceHeat

!**********************************************************************************************************************!
! AeroConductanceMomen ... compute the aerodynamic conductance for momentum (mol/m2-s)
!                          Bonan (2016) p. 215
!
!**********************************************************************************************************************!
  function AeroConductanceMomen(zm, zol, rhoair, ums)
    real(kind=dp), intent(in)  :: zm           ! height (m)
    real(kind=dp), intent(in)  :: zol          ! (z-d)/L, MO-stability parameter
    real(kind=dp), intent(in)  :: rhoair       ! air density at height zm (mol/m3)
    real(kind=dp), intent(in)  :: ums          ! mean wind velocity at height zm (m/s)
    real(kind=dp), parameter   :: vk=0.4       ! von Karman constant
    real(kind=dp), parameter   :: gammin=0.5   ! minimum value of AeroConductanceMomen
    real(kind=dp), parameter   :: gammax=10.0  ! maximum value of AeroDonductanceMomen
    real(kind=dp)              :: zcm          ! height in cm
    real(kind=dp)              :: gam          ! AeroConductanceMomen (mol/m2-s)
    real(kind=dp)              :: lnz0m, fm
    real(kind=dp)              :: AeroConductanceMomen  ! (mol/m2-s)

    zcm = 100.0*zm  ! convert to cm since d and z0m, z0h are in cm
    lnz0m = log((zcm-d)/z0m)
    fm = lnz0m - Psi_m(zol)
    gam = rhoair*vk*vk*ums/(fm*fm)
    AeroConductanceMomen = min(max(gam, gammin), gammax)

    return

  end function AeroConductanceMomen

!**********************************************************************************************************************!
! AeroConductanceWater ... compute the aerodynamic conductance for water vapor 
!                          (mol/m2-s)
!                          Bonan (2016) p. 215
!
!**********************************************************************************************************************!
  function AeroConductanceWater(zm, zol, rhoair, ums)
    real(kind=dp), intent(in)  :: zm           ! height (m)
    real(kind=dp), intent(in)  :: zol          ! (z-d)/L, MO-stability parameter
    real(kind=dp), intent(in)  :: rhoair       ! air density at height zm (mol/m3)
    real(kind=dp), intent(in)  :: ums          ! mean wind velocity at height zm (m/s)
    real(kind=dp), parameter   :: vk=0.4       ! von Karman constant
    real(kind=dp), parameter   :: gawmin=0.5   ! minimum value of AeroConductanceWater
    real(kind=dp), parameter   :: gawmax=10.0  ! maximum value of AeroDonductanceWater
    real(kind=dp)              :: zcm          ! height in cm
    real(kind=dp)              :: gaw          ! AeroConductanceWater (mol/m2-s)
    real(kind=dp)              :: lnz0m, lnz0w, fm, fw
    real(kind=dp)              :: AeroConductanceWater  ! (mol/m2-s)

    zcm = 100.0*zm  ! convert to cm since d and z0m, z0h are in cm
    lnz0m = log((zcm-d)/z0m)
    lnz0w = log((zcm-d)/z0h)
    fm = lnz0m - Psi_m(zol)
    fw = lnz0w - Psi_w(zol)
    gaw = rhoair*vk*vk*ums/(fm*fw)
    AeroConductanceWater = min(max(gaw, gawmin), gawmax)
    
    return 

  end function AeroConductanceWater

!**********************************************************************************************************************!
! Psi_m ... stability function for momentum
!            Bonan (2016) p. 215
!
!**********************************************************************************************************************!
  function Psi_m(zol)
    real(kind=dp), intent(in)  :: zol        ! (z-d)/L, MO-stability parameter
    real(kind=dp)              :: x, p1, p2, p3
    real(kind=dp)              :: Psi_m      ! momentum stability function
    
    if (zol < 0.0) then
      x = (1.0 - 16.0*zol)**(0.25)
      p1 = 2.0*log(0.5*(1.0+x))
      p2 = log(0.5*(1.0+x*x))
      p3 = 2.0*atan(x)
      Psi_m = p1 + p2 - p3 + 0.5*pi
    else
      Psi_m = -5.0*zol
    end if

    return

  end function Psi_m

!**********************************************************************************************************************!
! Psi_h ... stability function for heat
!            Bonan (2016) p. 215
!
!**********************************************************************************************************************!
  function Psi_h(zol)
    real(kind=dp), intent(in)  :: zol        ! (z-d)/L, MO-stability parameter
    real(kind=dp)              :: x
    real(kind=dp)              :: Psi_h      ! heat stability function

    if (zol < 0.0) then
      x = (1.0 - 16.0*zol)**(0.25)
      Psi_h = 2.0*log(0.5*(1.0+x*x))
    else
      Psi_h = -5.0*zol
    end if

    return

  end function Psi_h

!**********************************************************************************************************************!
! Psi_w ... stability function for water vapor
!            Bonan (2016) p. 215
!
!**********************************************************************************************************************!
  function Psi_w(zol)
    real(kind=dp), intent(in)  :: zol        ! (z-d)/L, MO-stability parameter
    real(kind=dp)              :: x
    real(kind=dp)              :: Psi_w      ! water vapor stability function

    if (zol < 0.0) then
      x = (1.0 - 16.0*zol)**(0.25)
      Psi_w = 2.0*log(0.5*(1.0+x*x))
    else
      Psi_w = -5.0*zol
    end if
  
    return

  end function Psi_w

!**********************************************************************************************************************!
! CalcTemp ... compute air temperature (K) at z
!
!**********************************************************************************************************************!
  function CalcTemp(zn, tac, tsc)
    real(kind=dp), intent(in)  :: zn        ! current z (cm)
    real(kind=dp), intent(in)  :: tac       ! air temperature at zref (C)
    real(kind=dp), intent(in)  :: tsc       ! soil temperature (C)
    real(kind=dp)              :: dtmp      ! temperature gradient (C/cm)
    real(kind=dp)              :: CalcTemp  ! air temp (K) at z 

    ! zref and zn must be same units!!!

    dtmp = (tac-tsc)/zref
    CalcTemp = (tsc + dtmp*zn) + 273.15  

    return
  end function CalcTemp

!**********************************************************************************************************************!
! CalcPressure ... compute air pressure (mb) at z
!
!**********************************************************************************************************************!
  function CalcPressure(pmbzref, zk, ti, t0)
    real(kind=dp), intent(in)  :: pmbzref   ! air pressure at zref (mb)
    real(kind=dp), intent(in)  :: zk        ! current z (cm)
    real(kind=dp), intent(in)  :: ti        ! temperature at zref (K)
    real(kind=dp), intent(in)  :: t0        ! temperature at z0 (K)
    real(kind=dp), parameter   :: a0=3.42D-04   ! Mg/R (K/cm)
    real(kind=dp)              :: CalcPressure  ! air pressure at current z (mb)

    CalcPressure = pmbzref*exp(-a0*(zk-zi)/(0.5*(ti+t0)))

    return
  end function CalcPressure

!**********************************************************************************************************************!
! CalcMeanWindSpeed ... compute mean wind speed (cm/s) at z
!
!**********************************************************************************************************************!
  function CalcMeanWindSpeed(zk, ubzref)
    real(kind=dp), intent(in)  :: zk         ! current z (cm)
    real(kind=dp), intent(in)  :: ubzref     ! mean wind speed at zref (cm/s)
    real(kind=dp)              :: gamlai     ! coefficient of wind speed attentuation
                                             !  function of LAI
    real(kind=dp)              :: blad       ! exponent of wind speed function
                                             !  function of LAD profile type
    real(kind=dp)              :: pladtype   ! LAD profile type = 1, 2 or 3
    real(kind=dp)   :: CalcMeanWindSpeed     ! mean wind speed at current z (cm/s)
    ! Meyers et al. (1998) A multilayer model for inferring dyr deposition using standard meteorological
    !  measurements, JGR, 103, D17, 22645-22661.

    pladtype=3  ! appropriate for forest canopy
    blad = 1.0 - (pladtype-1.0)*0.25 
    gamlai = 0.65*laitot
    if (gamlai > 4.0) then
      gamlai = 4.0
    end if

    if (zk <= hccm) then
      CalcMeanWindSpeed=ubzref*exp(-gamlai*((1.0 - (zk/hccm))**blad))
    else
      CalcMeanWindSpeed=ubzref
    end if

  end function CalcMeanWindSpeed

!**********************************************************************************************************************!
! CalcCair ... compute cair at z
!
!**********************************************************************************************************************!
  function CalcCair(pmbi, tki)
    real(kind=dp), intent(in)  :: pmbi       ! air pressure at current z (mb)
    real(kind=dp), intent(in)  :: tki        ! air temperature at current z (K)
    real(kind=dp)              :: CalcCair   ! air concentration (molec/cm3)

    CalcCair  = pmbi*7.2428D+18/tki
    return
  end function CalcCair

!**********************************************************************************************************************!
! CalcEddyDiffRa ... compute a Kv value at current z using calculated Ra value
!
!   Formulation based on ...
!     Surface layer Kv:  uses Kvsl = (zref-hc)/Ra
!     Canopy Kv: Raupach, M. R. (1989) A practical Lagrangian method for relating
!                 scalar concentrations to source distributions in vegetation
!                 canopies, Q. J. R. Meteorol. Soc., 115, 609-632.
!     Continuity of Kv enforced at z=hc via "fc"
!
!**********************************************************************************************************************!
  function CalcEddyDiffRa(zk, ubarzref, ubari, razref)
    real(kind=dp), intent(in)  :: zk        ! current z (cm)
    real(kind=dp), intent(in)  :: ubarzref  ! mean wind speed at zref (cm/s)
    real(kind=dp), intent(in)  :: ubari     ! mean wind speed at current z (cm/s)
    real(kind=dp), intent(in)  :: razref    ! aerodynamic resistance at zref (s/cm)
    real(kind=dp)              :: fc        ! Kv continuity factor
    real(kind=dp)              :: CalcEddyDiffRa  ! Kv at current z (cm2/s) 
    
    fc  = 15.1*(zref - hccm)/(ubarzref*hccm*razref)

    if (zk <= hccm) then
      CalcEddyDiffRa = fc*0.06643*ubari*hccm
    else
      CalcEddyDiffRa = (zref - hccm)/razref
    end if

    return
  end function CalcEddyDiffRa

!**********************************************************************************************************************!
! CalcEddyDiffStull ... compute a Kv value at current z using Stull's surface layer
!                          Kv approximation
!
!   Formulation based on ...
!     Surface layer Kv: Stull (1988) An Introduction to Boundary Layer Meteorology,
!                        Kluwer Academic Publishers, Dordrecht, The Netherlands.
!     Canopy Kv: Raupach, M. R. (1989) A practical Lagrangian method for relating
!                 scalar concentrations to source distributions in vegetation
!                 canopies, Q. J. R. Meteorol. Soc., 115, 609-632.
!     Continuity of Kv enforced at z=hc via "fc"
!
!**********************************************************************************************************************!
  function CalcEddyDiffStull(zk, ubari, zol)
    real(kind=dp), intent(in)  :: zk        ! current z (cm)
    real(kind=dp), intent(in)  :: ubari     ! mean wind speed at current z (cm/s)
    real(kind=dp), intent(in)  :: zol       ! (zk-d)/L  MO stability parameter
    real(kind=dp), parameter   :: vkc=0.4   ! von Karman constant
    real(kind=dp)              :: fc        ! Kv continuity factor
    real(kind=dp)              :: CalcEddyDiffStull  ! Kv at current z (cm2/s) 
    
    fc  = 0.843/Phi(zol)

    if (zk <= hccm) then
      CalcEddyDiffStull = fc*0.06643*ubari*hccm
    else
      CalcEddyDiffStull = vkc*ubari*0.14*zk/Phi(zol)
    end if

    return
  end function CalcEddyDiffStull

!**********************************************************************************************************************!
! CalcRiB ... compute the bulk Richardson number
!
!**********************************************************************************************************************!
  function CalcRiB(tak, tsk, ubari)
    real(kind=dp), intent(in) :: tak            ! air temperature at zref (K)
    real(kind=dp), intent(in) :: tsk            ! surface temperature (K)
    real(kind=dp), intent(in) :: ubari          ! mean wind speed at zref (cm/s)
    real(kind=dp), parameter  :: g=982.2        ! gravitational acceleration (cm/s2)
    real(kind=dp)             :: CalcRiB        ! bulk Richardson number ()
     
    CalcRiB = ((zref-d)*g*(tak-tsk))/(ubari*ubari*tak)
    return
  end function CalcRiB

!**********************************************************************************************************************!
! subroutine rav - calculate aerodynamic resistance (ra)
!
! Uses formulation of ...
!  Viney (1991) BLM, 56, 381-393.
!
!**********************************************************************************************************************!
function rav(ubar, zref, d, hc, rib)
  real(kind=dp)              :: rav     ! aerodynamic resistance, ra (s/cm)
  real(kind=dp), intent(in)  :: ubar    ! wind speed at reference height (cm/s)
  real(kind=dp), intent(in)  :: zref    ! reference height (m or cm)
  real(kind=dp), intent(in)  :: d       ! displacement height (m or cm)
  real(kind=dp), intent(in)  :: hc      ! height of canopy (m or cm)
                                        ! zref, d, and hc must all be same units!
  real(kind=dp), intent(out) :: rib     ! bulk Richardson number
  real(kind=dp), parameter   :: kv=0.4  ! von Karman constant
  real(kind=dp)              :: gmm     ! momentum gamma
  real(kind=dp)              :: gmh     ! heat gamma
  real(kind=dp)              :: z0m     ! momentum roughness height (m or cm)
  real(kind=dp)              :: z0h     ! heat roughness height (m or cm)
  real(kind=dp)              :: rapr    ! ra neutral (s/cm) 
  real(kind=dp)              :: phim    ! momentum stability function
  real(kind=dp)              :: phih    ! heat stability function
  real(kind=dp)              :: a, b, c ! Viney empirical constants for unstable 
  real(kind=dp)              :: a1, a2, a3

  z0m = 0.13*hc
  z0h = z0m/7.0

  gmm = log((zref-d)/z0m)
  gmh = log((zref-d)/z0h)

  rapr = (gmm*gmh)/(kv*kv*ubar)

  if (rib > 0.0) then
    ! stable
    phim = (gmh+10.0*rib*gmm-((gmh*gmh+20.0*rib*gmm*(gmh-gmm))**0.5))/(2.0-10.0*rib)
    phim = max(-5.0, phim)
    phih = phim
    rav = ((gmh-phih)*(gmm-phim))/(kv*kv*ubar)
  else if (rib < 0.0) then
    ! unstable
    a = 1.0591 - 0.0552*log(1.72+(4.03-gmm)**2)
    b = 1.9117 - 0.2237*log(1.86+(2.12-gmm)**2)
    c = 0.8437 - 0.1243*log(3.49+(2.79-gmm)**2)
    rav = rapr/(a + b*(-rib)**c)
  else
    ! neutral
    rav = rapr
  end if

  return
end function rav

!**********************************************************************************************************************!
! MapRiBtoZoL ... maps bulk Richardson number to Monin-Obukhov stability, z/L
!
!         Uses ...
!           Li, Y., Gao, Z., Lenschow, D. H., Chen, F. (2010) An improved approach
!            for parameterizing surface-layer turbulent transfer coefficients in
!            numerical models, BLM, 137, 153-165.
!
!**********************************************************************************************************************!
  function MapRiBtoZoL(rib)
    real(kind=dp), intent(in)  :: rib           ! bulk Richardson number ()
    real(kind=dp)              :: alph          ! ln(ratio of zref to z0m)
    real(kind=dp)              :: beta          ! ln(ratio of z0m to z0h)
    real(kind=dp)              :: a, b          ! fitting parameters
    real(kind=dp)              :: MapRiBtoZoL   ! Monin-Obukhov stability  parameter ()

    alph = log(zref/z0m)
    beta = log(z0m/z0h) 

    if (rib > 0.2) then
      ! strongly stable
      a = 0.7529*alph+14.94
      b = 0.1569*alph-0.3091*beta-1.303
      MapRibtoZoL = a*rib + b

    else if (rib > 0.0 .and. rib <= 0.2) then
      ! weakly stable
      a = (0.5738*beta-0.4399)*alph + (-4.901*beta+52.50)
      b = (-0.0539*beta+1.540)*alph + (-0.6690*beta-3.282)
      MapRibtoZoL = a*rib*rib + b*rib

    else
      ! unstable
      a = 0.0450*alph
      b  = (0.003*beta+0.0059)*alph*alph + (-0.0828*beta+0.8845)*alph   &
           + (0.1739*beta*beta-0.9213*beta-0.1057)
      MapRiBtoZoL = a*rib*rib + b*rib
    end if

    return
  end function MapRiBtoZoL

!**********************************************************************************************************************!
! Phi ... calculates Businger phi function
!
!         Businger et al. (1971) Flux profile relationships in the
!          atmospheric surface layer, J. Atmos. Sci., 28, 181-189.
!
!**********************************************************************************************************************!
  function Phi(zol)
    real(kind=dp), intent(in) :: zol  ! Monin-Obukhov stability parameter, z/L
    real(kind=dp)             :: Phi  ! phi parameter

    if (zol > 0.0) then
      ! stable
      Phi = 1.0 + 4.7*zol
    else if (zol < 0.0) then
      ! unstable
      Phi = (1.0 - 15*zol)**(-0.5)
    else
      ! neutral
      Phi = 1.0
    end if

    return
  end function Phi

end module EnvironData
