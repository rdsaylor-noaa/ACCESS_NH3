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
!     Program:      CANACC                                                                                             !
!                   Canopy physics model                                                                               !
!                                                                                                                      !
!     Version:      2.0.0                                                                                              !
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
!     Module:       CanopyPhysics                                                                                      !
!                                                                                                                      !
!     Description:  contains algorithms related to canopy physics                                                      !
!                                                                                                                      !
!======================================================================================================================!
module CanopyPhysics
  use GlobalData
  use PhysChemData
  use VertTransport
  use Utils
  implicit none

  private CalcLeafEnergyBal, IntegrateWaterVapor, IntegrateTempAir, CanopyPhysConverged, CalcSoilExchangeParams, &
          SurfaceAirTemp, CalcOneLevelLEB, LeafTemp, LeafTemp1, LeafTemp2, LeafBalEq, gh_cn, gv_cn, gb_cn, &
          gs_medlyn, NetPhoto, kbe
  public CalcCanopyPhysics, PartitionRAD, CalcRadProfiles, CalcWeightedProfiles, &
         rbgl, rsoill, rbl, rw_massad, rw_flechard, rs_zhang_nh3, rs_zhang, rs_rds, &
         rs_jarvis1, rs_jarvis2, rs_jim, rs_nmj, rincl, gpl, gpstl, gpgl

contains

!**********************************************************************************************************************!
! subroutine CalcCanopyPhysics - calculates leaf temperatures, net photosynthesis assimilation rates and
!                                stomatal conductances for all levels in the sunlit and shaded canopy
!
!**********************************************************************************************************************!
subroutine CalcCanopyPhysics()
  integer(kind=i4)            :: i
  integer(kind=i4)            :: iter
  integer(kind=i4), parameter :: maxiter=20

  ! Calculate ksoil, qsoil, vsh2o
  call CalcSoilExchangeParams()

  do iter=1,maxiter

    ! save previous iteration results
    do i=1,npts
      tkp(i) = tk(i)
      tlp_sun(i) = tl_sun(i)
      tlp_shd(i) = tl_shd(i)
      qhp(i) = qh(i)
    end do

    ! Calculate Tleaf(z)
    call CalcLeafEnergyBal()

    if (INTGWVAP .eqv. .true.) then
       ! Calculate qh(z)
       call IntegrateWaterVapor(t, t+dtout)
    end if

    if (INTGTAIR .eqv. .true.) then
       ! Calculate Ta(z)
       call IntegrateTempAir(t, t+dtout)
    end if

    if (CanopyPhysConverged()) then
!     print *, 'CalcCanopyPhysics: iter = ', iter
      return
    end if

  end do
  
  write(*,*) 'Canopy physics did not converge!'

  return

end subroutine CalcCanopyPhysics

!**********************************************************************************************************************!
! subroutine CalcLeafEnergyBal - calculates leaf energy balance over the entire canopy
!                                for sunlit and shaded fractions
!
!**********************************************************************************************************************!
subroutine CalcLeafEnergyBal()
  integer(kind=i4)            :: i
  real(kind=dp)               :: relhumi            ! relative humidity, %
  real(kind=dp)               :: ubarms             ! mean wind speed, m/s

  do i=ncnpy,1,-1

    relhumi = RelativeHumidity(tk(i), pmb(i), qh(i))
    ubarms = ubar(i)*0.01

    ! sunlit canopy fraction
    call CalcOneLevelLEB(rabs_sun(i), ppfd_sun(i), tk(i), pmb(i), ubarms, relhumi, gs_sun(i), gl_sun(i), &
                         gb(i), anet_sun(i), tl_sun(i))

    ! shaded canopy fraction
    call CalcOneLevelLEB(rabs_shd(i), ppfd_shd(i), tk(i), pmb(i), ubarms, relhumi, gs_shd(i), gl_shd(i), &
                         gb(i), anet_shd(i), tl_shd(i))

  end do

  return

end subroutine CalcLeafEnergyBal

!**********************************************************************************************************************!
! subroutine IntegrateWaterVapor - integrates water vapor continuity equation throughout domain
!
!**********************************************************************************************************************!
subroutine IntegrateWaterVapor(tin, tout)
  real(kind=dp), intent(in) :: tin, tout
  real(kind=dp)             :: ts         ! current integration time [tin, tout]                  
  real(kind=dp)             :: dtphys     ! integration time step (secs)
  integer(kind=i4)     :: i
  real(kind=dp), dimension(npts) :: q       ! local water vapor concentration, moles/cm3
  real(kind=dp), dimension(npts) :: kvm     ! Kv (cm2/s)
  real(kind=dp), dimension(npts) :: rho     ! air density (moles/cm3)
  real(kind=dp), dimension(npts) :: phiq    ! new solution array
  real(kind=dp), dimension(npts) :: phiqp   ! previous solution array
  real(kind=dp), dimension(npts) :: vl_sun  ! gl_sun in units of cm/s
  real(kind=dp), dimension(npts) :: vl_shd  ! gl_shd in units of cm/s

  ts = tin
  dtphys = 5.0_dp

  ! initialize q with water vapor concentration from last full time step
  do i=1,npts
    q(i) = h2o(i)/navo     ! convert from molec/cm3 to moles/cm3
  end do

  ! integration over tin -> tout
  do
    ! Integrate source term (sun/shade)
    ! gl_sun, gl_shd, fsun calculated in CalcLeafEnergyBal
    do i=1,npts
      vl_sun(i) = 0.08314*gl_sun(i)*tl_sun(i)*100./pmb(i)    ! convert gl from mol/m2-s to cm/s
      vl_shd(i) = 0.08314*gl_shd(i)*tl_shd(i)*100./pmb(i)
      esun(i) = vl_sun(i)*(qsat(tl_sun(i)) - q(i))
      eshd(i) = vl_shd(i)*(qsat(tl_shd(i)) - q(i))
      q(i) = q(i) + dtphys*(esun(i)*fsun(i) + eshd(i)*(1.0-fsun(i)))*lad(i)
    end do  

    ! Integrate vert transport
    do i=1,npts
      phiq(i) = 0.0_dp
      phiqp(i) = q(i)
      kvm(i) = kv(i)
      rho(i) = cair(i)/navo     ! convert from molec/cm3 to moles/cm3
    end do

    call SubIntegrateVertTransportSplitBC(phiq, phiqp, kvm, qzref, qsoil, vsh2o, rho, dtphys)

    do i=1,npts
      q(i) = phiq(i)
    end do

    ts = ts + dtphys

    if (ts >= tout) exit

  end do

  ! fill qh and h2o arrays with final result
  do i=1,npts
    h2o(i) = q(i)*navo                     ! convert from moles/cm3 to molec/cm3
    qh(i) = Convert_h2o_to_qh(h2o(i), cair(i))    ! convert molec/cm3 to g/kg
  end do 

  return 

end subroutine IntegrateWaterVapor

!**********************************************************************************************************************!
! subroutine IntegrateTempAir - integrates sensible heat equation throughout domain, producing
!                               air temperature profile 
!
!**********************************************************************************************************************!
subroutine IntegrateTempAir(tin, tout)
  real(kind=dp), intent(in) :: tin, tout
  real(kind=dp)             :: ts         ! current integration time [tin, tout]                  
  real(kind=dp)             :: dtphys     ! integration time step (sec)
  real(kind=dp)             :: tk0        ! lower boundary condition for Tair derived from surface energy balance (K)
  real(kind=dp)             :: rcp        ! rho*cpair (moles/m3)*(J/mole-K)
  real(kind=dp), dimension(npts) :: phitk    ! new solution array
  real(kind=dp), dimension(npts) :: phitkp   ! previous solution array
  real(kind=dp), dimension(npts) :: kvm      ! Kv (m2/s)
  real(kind=dp), dimension(npts) :: rho      ! molar air density (moles/m3)
  integer(kind=i4)     :: i

  ts = tin
  dtphys = 1.0_dp

  tk0 = SurfaceAirTemp()

  ! integration over tin -> tout
  do
    ! Integrate source term (sun/shade)
    ! gb, fsun calculated in CalcLeafEnergybal
    do i=1,npts
      rho(i) = cair(i)/navo          ! convert from molec/cm3 to moles/cm3
      rcp = rho(i)*cpair
      hsun(i) = 2.0*cpair*gb(i)*1.0E-04*(tl_sun(i) - tk(i))  ! convert gb from mol/m2-s to mol/cm2-s
      hshd(i) = 2.0*cpair*gb(i)*1.0E-04*(tl_shd(i) - tk(i))  ! convert gb from mol/m2-s to mol/cm2-s
      tk(i) = tk(i) + dtphys*(hsun(i)*fsun(i) + hshd(i)*(1.0-fsun(i)))*lad(i)/rcp
    end do  

    ! Integrate vert transport
    do i=1,npts
      phitk(i) = 0.0_dp
      phitkp(i) = tk(i)
      kvm(i) = kv(i)
    end do

    call SubIntegrateVertTransportConstBC(phitk, phitkp, kvm, tkzref, tk0, rho, dtphys)

    do i=1,npts
      tk(i) = phitk(i)
      cair(i) = pmb(i)*7.2428D+18/tk(i)    ! update cair with new Tair profile (molec/cm3)
    end do

    ts = ts + dtphys

    if (ts >= tout) exit

  end do

  return 

end subroutine IntegrateTempAir

!**********************************************************************************************************************!
! function CanopyPhysConverged - checks for convergence among air temperature, water vapor and leaf temperature 
!                                variables (the infinity norm of the |relative error| of each variable must be
!                                less than a prescribed tolerance (ctol) 
!
!**********************************************************************************************************************!
function CanopyPhysConverged()
  integer(kind=i4)       :: i
  logical                :: CanopyPhysConverged
  real(kind=dp)          :: mtk, mtlsun, mtlshd, mqh
  real(kind=dp)          :: tkrd, tlsunrd, tlshdrd, qhrd

  ! have tk, tl and q all converged?
  mtk=0.0
  mtlsun=0.0
  mtlshd=0.0
  mqh=0.0
  do i=1,npts 
    tkrd = abs(tk(i)-tkp(i))/tkp(i)
    tlsunrd = abs(tl_sun(i)-tlp_sun(i))/tlp_sun(i)
    tlshdrd = abs(tl_shd(i)-tlp_shd(i))/tlp_shd(i) 
    qhrd = abs(qh(i)-qhp(i))/qhp(i)
    if (tkrd > mtk) mtk=tkrd
    if (tlsunrd > mtlsun) mtlsun=tlsunrd
    if (tlshdrd > mtlshd) mtlshd=tlshdrd
    if (qhrd > mqh) mqh=qhrd 
  end do

  CanopyPhysConverged = .false.

  if ( (mtk <= ctol) .and. (mtlsun <= ctol) .and. (mtlshd <= ctol) .and. (mqh <= ctol) ) then
    CanopyPhysConverged = .true.
  end if

  return

end function CanopyPhysConverged

!**********************************************************************************************************************!
! subroutine CalcSoilExchangeParams - calculations parameters for exchange of water vapor with soil surface
!
! Uses formulation of ...
!    Pleim et al. (2013) JGR, 118, 3794-3806.
!    Schuepp (1977) BLM, 12, 171-186.
!
! With data from ...
!    Rawls et al. (1982) Trans. ASAE, 25, 1316-1320.
!    Clapp and Hornberger (1978) Water Resources Res., 14, 601-604.
!
!**********************************************************************************************************************!
subroutine CalcSoilExchangeParams()
  real(kind=dp)              :: fsat           ! fractional soil saturation
  real(kind=dp)              :: mdiff          ! molecular diffusivity of water vapor (cm2/s)
  real(kind=dp), parameter   :: viscair=0.155  ! dynamic viscosity of air (cm2/s)
  real(kind=dp)              :: ugstar, del0 
  real(kind=dp)              :: ldry, mdiffp, xe, rsh2o

  ! soil thermal conductivity
  fsat = (stheta-rtheta)/(sattheta-rtheta)
  ksoil = ksoildry(isoiltype) + (ksoilwet(isoiltype)-ksoildry(isoiltype))*fsat

  ! surface molar water vapor concentration (moles/cm3)
  qsoil = hsoil*qsat(tsoilk)

  ! ground boundary layer resistance (s/cm)
  mdiff=mdiffh2o(tsoilk, pmb(1))
  ugstar = 0.05*ubar(ncnpy+1)+0.00005*ubar(ncnpy+1)*ubar(ncnpy+1) 
  del0 = viscair/(0.4*ugstar)
  rbg = ((viscair/mdiff)-dlog(del0/10.0))/(0.4*ugstar)

  ! Note:  These are the correct units needed in SurfaceAirTemp!
  ! ground aerodynamic conductance (mol/m2-s) as in Bonan (2016)
  ! gaero is calculated in EnvironData
  gbg = gaero(nt)

  ! soil diffusion resistance (s/cm)
  xe =(1.0-(stheta/sattheta))**5.0
  ldry = dsoil*(dexp(xe)-1.0)/1.7183
  mdiffp = mdiff*sattheta*sattheta*(1.0-(rtheta/sattheta))**(2.0+3.0/sbcoef)
  rsh2o = ldry/mdiffp

  ! deposition velocity to ground surface (cm/s)
  vsh2o = 1.0/(ra(nt)+rsh2o)       

  return

end subroutine CalcSoilExchangeParams

!**********************************************************************************************************************!
! function SurfaceAirTemp - provides the lower (i.e., ground surface) boundary condition for integration of 
!                           air temperature through domain.  Uses an energy balance at the surface to calculate
!                           the effective air temperature right at the ground surface.
!
!**********************************************************************************************************************!
function SurfaceAirTemp()
  real(kind=dp)            :: SurfaceAirTemp    ! air temperature at the ground surface (K)
  real(kind=dp)            :: q0                ! water vapor concentration at the ground surface (mol/cm3)
  real(kind=dp)            :: rn                ! net radiative flux at the surface (W/m2)
  real(kind=dp)            :: le                ! latent heat flux at the surface (W/m2)
  real(kind=dp)            :: g                 ! ground heat flux (W/m2)
  real(kind=dp)            :: delt              ! temperature difference from measured soil temperature (K)
  real(kind=dp), parameter :: rhovis=0.15
  real(kind=dp), parameter :: rhonir=0.25

  ! water vapor density at the surface
  q0 = h2o(1)/navo            ! convert from molecs/cm3 to mol/cm3

  ! net radiaive flux
  rn = (ppfd_wgt(1)*(1.0-rhovis)/4.6) + nir_wgt(1)*(1.0-rhonir) + (lw_dn(1) - lw_up(1))

  ! latent heat flux
  le = lambda(tsoilk)*(qsoil-q0)*vsh2o*1.0E+04   ! convert from W/cm2 to W/m2

  ! ground heat flux 
  g  = ksoil*dtdzsoil 

  ! temperature difference between soil and air based on surface energy balance, Rn = H + LE + G
  ! where H is the sensible heat flux, H = cpair*gbg*(Tsoil-Tair)
  delt = (rn - le - g)/(cpair*gbg)

  SurfaceAirTemp = tsoilk - delt

end function SurfaceAirTemp

!**********************************************************************************************************************!
! subroutine CalcOneLevelLEB - calculates leaf energy balance and leaf temperature, net photosynthesis assimilation 
!                               rate and stomatal conductance for one level in the canopy
!
! iterates through equations until stomatal conductance converges to specified tolerance 
!
!**********************************************************************************************************************!
subroutine CalcOneLevelLEB(rabs, ppfdi, tki, pmbi, ubarms, relhumi, gs, gv, gb, anet, tleaf)
  real(kind=dp), intent(in)   :: rabs                 !
  real(kind=dp), intent(in)   :: ppfdi                !
  real(kind=dp), intent(in)   :: tki                  !
  real(kind=dp), intent(in)   :: pmbi                 !
  real(kind=dp), intent(in)   :: ubarms               !
  real(kind=dp), intent(in)   :: relhumi              !
  real(kind=dp), intent(out)  :: gs                   ! leaf stomatal conductance, mol/m2-s
  real(kind=dp), intent(out)  :: gb                   ! boundary layer conductance, mol/m2-s
  real(kind=dp), intent(out)  :: gv                   ! overall leaf conductance, mol/m2-s
  real(kind=dp), intent(out)  :: anet                 ! net photosynthesis assimilation rate, umol/m2-s
  real(kind=dp), intent(out)  :: tleaf                ! leaf temperature, K
  real(kind=dp)               :: cci                  ! intercellular CO2 mixing ratio, umol/mol
  real(kind=dp)               :: delgs                ! absolute difference in gs between iterations
  real(kind=dp)               :: gh                   ! heat conductance, mol/m2-s
  real(kind=dp),parameter     :: gr=0.15              ! radiative conductance, mol/m2-s (Rr=300 s/m)
  real(kind=dp)               :: gsnew                ! updated value of gs, mol/m2-s
  integer(kind=i4), parameter :: maxiter=10
  integer(kind=i4)            :: iter

  ! initial guesses
  cci = 200.
  gs = 0.2

  ! heat conductance w/ radiative conductance
  gh = gh_cn(ubarms, dleaf, gr)

  ! boundary layer conductance
  gb = gb_cn(ubarms, dleaf)

  delgs=gs
  iter=0
  do while (delgs/gs > ltol .and. iter <= maxiter)

    ! overall leaf conductance
    gv = gv_cn(gb, gs)
  
    ! calculate leaf temperature
    tleaf = LeafTemp(tki, pmbi, relhumi, rabs, gh, gv)

    ! calculate photosynthesis assimilation
    anet = NetPhoto(tleaf, cci, ppfdi)

    ! calculate new stomatal conductance
    gsnew = gs_medlyn(tleaf, relhumi, anet, cci)

    delgs = abs(gsnew-gs)

    gs = gsnew

    iter=iter+1

  end do

  return

end subroutine CalcOneLevelLEB

!**********************************************************************************************************************!
! function LeafTemp - calculates the leaf temperature
!
! Uses simple Newton-Raphson method to compute tleaf based on code from
!  Press et al. (1992) Numerical Recipes in Fortran 77, Chap. 9. Root Finding and Nonlinear Sets of Equations
! No approximations are made in this formulation, which should produce the most accurate estimate of leaf temperature.
!
!**********************************************************************************************************************!
function LeafTemp(tki, pmbi, relhumi, rabs, gh, gv)
  real(kind=dp), intent(in)  :: tki               ! air temperature, K
  real(kind=dp), intent(in)  :: pmbi              ! air pressure, mb
  real(kind=dp), intent(in)  :: relhumi           ! relative humidity, %
  real(kind=dp), intent(in)  :: rabs              ! absorbed radiation, W/m2
  real(kind=dp), intent(in)  :: gh                ! heat conductance for leaf-atmosphere
  real(kind=dp), intent(in)  :: gv                ! leaf water vapor conductance 
  real(kind=dp)              :: LeafTemp          ! leaf temperature, K
  integer(kind=i4)           :: i
  integer(kind=i4), parameter :: maxiter=10       ! maximum number of iterations
  real(kind=dp), parameter    :: rltol=0.1        ! convergence tolerance
  real(kind=dp)               :: tair, tleaf
  real(kind=dp)               :: dt, f, df

  tleaf=tki
  tair=tki

  do i=1,maxiter
    call LeafBalEq(tleaf, tair, pmbi, relhumi, rabs, gh, gv, f, df)
    dt = f/df
    tleaf = tleaf - dt
    if (abs(dt) < rltol) then
      LeafTemp = tleaf
      return 
    end if
  end do
  write(*,*) 'LeafTemp non-convergence!'
  LeafTemp=tki
  return

end function LeafTemp

!**********************************************************************************************************************!
! subroutine LeafBalEq - computes value for leaf energy balance equation and its first derivative
!
! Used for full non-linear version of the solution employing the Newton_Raphson method.
!
!**********************************************************************************************************************!
subroutine LeafBalEq(tleaf, tair, pmbi, relhumi, rabs, gh, gv, f, df)
  real(kind=dp), intent(in)  :: tleaf             ! leaf temperature, K
  real(kind=dp), intent(in)  :: tair              ! air temperature, K
  real(kind=dp), intent(in)  :: pmbi              ! air pressure, mb
  real(kind=dp), intent(in)  :: relhumi           ! relative humidity, %
  real(kind=dp), intent(in)  :: rabs              ! absorbed radiation, W/m2
  real(kind=dp), intent(in)  :: gh                ! heat conductance for leaf-atmosphere
  real(kind=dp), intent(in)  :: gv                ! leaf water vapor conductance 
  real(kind=dp)              :: vpd               ! vapor pressure deficit, kPa
  real(kind=dp)              :: patm              ! air pressure, kPa
  real(kind=dp)              :: s                 ! first derivative of esat parameter
  real(kind=dp)              :: f
  real(kind=dp)              :: df 

  vpd = esat(tleaf) - esat(tair)*relhumi/100.
  patm = pmbi*0.1   ! convert mb to kPa
  s = desdt(tleaf)/patm

  f = -rabs + 2.0*epsleaf*sbsig*tleaf**4.0 + cpair*gh*(tleaf-tair) + lambda(tleaf)*gv*vpd/patm
  df = 8.0*epsleaf*sbsig*tleaf**3.0 + cpair*gh + lambda(tleaf)*gv*s

end subroutine LeafBalEq

!**********************************************************************************************************************!
! function LeafTemp2 - calculates the leaf temperature
!
! Uses a formulation of the leaf energy balance equation where 2nd-order Taylor Series expansions are used to
! approximate the leaf's longwave emission and surface water vapor saturation pressure, producing a quadratic
! equation for leaf temprature. 2nd most accurate estimate of the leaf temperature is produced.
!
!**********************************************************************************************************************!
function LeafTemp2(tki, pmbi, relhumi, rabs, gh, gv)
  real(kind=dp), intent(in)  :: tki               ! air temperature, K
  real(kind=dp), intent(in)  :: pmbi              ! air pressure, mb
  real(kind=dp), intent(in)  :: relhumi           ! relative humidity, %
  real(kind=dp), intent(in)  :: rabs              ! absorbed radiation, W/m2
  real(kind=dp), intent(in)  :: gh                ! heat conductance for leaf-atmosphere
  real(kind=dp), intent(in)  :: gv                ! leaf water vapor conductance 
  real(kind=dp)              :: LeafTemp2         ! leaf temperature, K
  real(kind=dp)              :: vpd               ! vapor pressure deficit, kPa
  real(kind=dp)              :: patm              ! air pressure, kPa
  real(kind=dp)              :: s                 ! first derivative of esat parameter
  real(kind=dp)              :: sp                ! second derivative of esat parameter
  real(kind=dp)              :: a, b, c

  vpd = esat(tki)*(1.0-relhumi/100.)
  patm = pmbi*0.1   ! convert mb to kPa
  s = desdt(tki)/patm
  sp = des2dt2(tki)/patm

  a = 12.0*epsleaf*sbsig*tki*tki + lambda(tki)*gv*sp*0.5
  b = 8.0*epsleaf*sbsig*tki**3.0 + 2.0*cpair*gh + lambda(tki)*gv*s
  c = -rabs + 2.0*epsleaf*sbsig*tki**4.0 + lambda(tki)*gv*vpd/patm

  LeafTemp2 = tki + (- b + ((b*b-4.0*a*c)**0.5))/(2.0*a)

  return

end function LeafTemp2

!**********************************************************************************************************************!
! function LeafTemp1 - calculates the leaf temperature
!
! Uses a formulation of the leaf energy balance equation where 1st-order Taylor Series expansions are used to
! approximate the leaf's longwave emission and surface water vapor saturation pressure, prooducing a linear
! equation for leaf temperature. Least accurate estimate of the leaf temperature is produced.
!
!**********************************************************************************************************************!
function LeafTemp1(tki, pmbi, relhumi, rabs, gh, gv)
  real(kind=dp), intent(in)  :: tki               ! air temperature, K
  real(kind=dp), intent(in)  :: pmbi              ! air pressure, mb
  real(kind=dp), intent(in)  :: relhumi           ! relative humidity, %
  real(kind=dp), intent(in)  :: rabs              ! absorbed radiation, W/m2
  real(kind=dp), intent(in)  :: gh                ! heat conductance for leaf-atmosphere
  real(kind=dp), intent(in)  :: gv                ! leaf water vapor conductance 
  real(kind=dp)              :: LeafTemp1         ! leaf temperature, K
  real(kind=dp)              :: vpd               ! vapor pressure deficit, kPa
  real(kind=dp)              :: patm              ! air pressure, kPa
  real(kind=dp)              :: s                 ! first derivative of esat parameter

  vpd = esat(tki)*(1.0-relhumi/100.)
  patm = pmbi*0.1   ! convert mb to kPa
  s = desdt(tki)/patm
  LeafTemp1 = tki + ((rabs - 2.0*epsleaf*sbsig*tki*tki*tki*tki - lambda(tki)*gv*vpd/patm)/ &
                    (8.0*epsleaf*sbsig*tki*tki*tki + 2.0*cpair*gh + lambda(tki)*s*gv))

  return

end function LeafTemp1

!**********************************************************************************************************************!
! function gh_cn - calculates the heat conductance using the forced convection formulation
!                  of Campbell & Norman, Table 7.6, p. 109, plus an assumed value for the radiative conductance.
!
!  Campbell, G. S., and Norman, J. M. (1998) An Introduction to Environmental Biophysics, 2nd Ed., Springer, New York.
!
!**********************************************************************************************************************!
function gh_cn(ubar, dleaf, gr)
  real(kind=dp), intent(in) :: ubar              ! mean wind speed, m/s
  real(kind=dp), intent(in) :: dleaf             ! leaf dimension, m
  real(kind=dp), intent(in) :: gr                ! radiative conductance, mol/m2-s
  real(kind=dp)             :: gh_cn             ! total heat conductance, mol/m2-s

  gh_cn = 1.4*0.135*sqrt(ubar/dleaf)
  gh_cn = gh_cn*gr/(gh_cn+gr)

  return

end function gh_cn

!**********************************************************************************************************************!
! function gb_cn - calculates the boundary layer conductance using the forced convection formulation
!                  of Campbell & Norman, Table 7.6, p. 109

!  Campbell, G. S., and Norman, J. M. (1998) An Introduction to Environmental Biophysics, 2nd Ed., Springer, New York.
!
!**********************************************************************************************************************!
function gb_cn(ubar, dleaf)
  real(kind=dp), intent(in) :: ubar              ! mean wind speed, m/s
  real(kind=dp), intent(in) :: dleaf             ! leaf dimension, m
  real(kind=dp)             :: gb_cn             ! boundary layer conductance, mol/m2-s

  gb_cn = 1.4*0.147*sqrt(ubar/dleaf)

  return

end function gb_cn

!**********************************************************************************************************************!
! function gv_cn - calculates the overall leaf conductance using approximate formulation
!                  of Campbell & Norman, Table 7.6, p. 109

!  Campbell, G. S., and Norman, J. M. (1998) An Introduction to Environmental Biophysics, 2nd Ed., Springer, New York.
!
!**********************************************************************************************************************!
function gv_cn(gb, gs)
  real(kind=dp), intent(in) :: gb                ! boundary layer conductance, mol/m2-s
  real(kind=dp), intent(in) :: gs                ! stomatal conductance, mol/m2-s
  real(kind=dp)             :: gv_cn             ! overall leaf conductance, mol/m2-s

  gv_cn = gb*gs/(gb+gs)

  return

end function gv_cn

!**********************************************************************************************************************!
! function NetPhoto - calculates net photosynthesis assimilation rate
!
! Uses the formulation of ...
!  Collatz et al. (1991) Physiological and environmental regulation of stomatal conductance, photosynthesis, and
!  transpiration: a model that includes a laminar boundary layer, Agric. For. Meteorol., 54, 107-136.
!
! as expressed in ...
!  Campbell, G. S., and Norman, J. M. (1998) An Introduction to Environmental Biophysics, 2nd Ed., Springer, New York.
!
!**********************************************************************************************************************!
function NetPhoto(tleaf, cci, ppfd)
  implicit none
  real(kind=dp), intent(in)  :: tleaf            ! leaf temperature, K
  real(kind=dp), intent(in)  :: cci              ! intercellular CO2 concentration, umol/mol
  real(kind=dp), intent(in)  :: ppfd             ! photosynthetic photon flux density, umol/m2-s
  real(kind=dp)              :: NetPhoto         ! net photosynthesis assimilation rate, umol/m2-s
  real(kind=dp)              :: tlc              ! leaf temperature, C
  real(kind=dp)              :: agross           ! gross photosynthesis assimilation rate, umol/m2-s
  real(kind=dp)              :: rd               ! leaf respiration rate, umol/m2-s
  real(kind=dp)              :: jp               ! tmp variable, colimitation assimilation, umol/m2-s
  real(kind=dp)              :: js               ! photosynthesis assimilation rate imposed by sucrose 
                                                 ! synthesis, umol/m2-s
  real(kind=dp)              :: jc               ! Rubisco-limited photosynthesis assimilation rate, umol/m2-s
  real(kind=dp)              :: je               ! light-limited photosynthesis assimilation rate, umol/m2-s
  real(kind=dp)              :: gamstr           ! light compensation point - CO2 conc at which assimilation is 
                                                 ! zero, umol/mol
  real(kind=dp)              :: tau              ! CO2/O2 specificity ratio, mmol/umol
  real(kind=dp)              :: vm               ! max Rubisco capacity per unit leaf area, umol/m2-s
  real(kind=dp)              :: kc               ! Michaelis constant for CO2, umol/mol
  real(kind=dp)              :: k0               ! inhibition constant for O2, umol/mol
  real(kind=dp), parameter   :: tau25=2600.      ! CO2/O2 specificity ratio @ 25degC, mmol/mmol
  real(kind=dp), parameter   :: vm25=100.0       ! max Rubisco capacity per unit leaf area @25degC, umol/m2-s
  real(kind=dp), parameter   :: kc25=300.        ! Michaelis constant for CO2 @ 25degC, umol/mol
  real(kind=dp), parameter   :: k025=300000.     ! inhibition constant for O2 @ 25degC, umol/mol
  real(kind=dp), parameter   :: coa=210000.      ! oxygen mole fraction, umol/mol          
  real(kind=dp), parameter   :: alfp=0.8         ! leaf absorptivity for PPFD
  real(kind=dp), parameter   :: em=0.08          ! maximum quantum efficiency, mol/mol
  real(kind=dp), parameter   :: beta=0.98        ! transition parameter between Jp and Js
  real(kind=dp), parameter   :: theta=0.95       ! transition parameter between Jp and Js
  real(kind=dp), parameter   :: rd25=1.5         ! leaf respiration rate @ 25degC, umol/m2-s

  ! equations require leaf temperature in degC
  tlc = tleaf - 273.15

  ! ratio describing partitioning of RuBP to the carboxylase and oxygenase 
  ! reactions of Rubisco
  tau = tau25*exp(-0.056*(tlc-25.0))

  ! CO2 conc where assimilation is zero
  gamstr = coa/(2.0*tau)

  ! light-limited photosynthesis assimilation rate
  je = alfp*em*ppfd*(cci - gamstr)/(cci+2.0*gamstr)

  ! max Rubisco capacity per unit leaf area
  vm = vm25*exp(0.088*(tlc-25.0))/(1.0+exp(0.29*(tlc-41.0)))

  ! Michaelis constant for CO2 
  kc = kc25*exp(0.074*(tlc-25.0))

  ! inhibition constant for O2
  k0 = k025*exp(0.018*(tlc-25.0))

  ! Rubisco-limited photosynthesis assimilation rate
  jc = vm*(cci - gamstr)/(cci + kc*(1.0+coa/k0))

  ! colimitation between Je and Jc
  jp = (je + jc - sqrt((je+jc)*(je+jc) - 4.0*theta*je*jc))/(2.0*theta)

  ! sucrose-limited photosynthesis assimilation rate
  js = vm*0.5

  ! gross photosynthesis assimilation rate
  agross = (jp + js - sqrt((jp+js)*(jp+js) - 4.0*beta*jp*js))/(2.0*beta)

  ! leaf respiration rate
  rd = rd25*exp(0.069*(tlc-25.0))/(1.0+exp(1.3*(tlc-55.0)))
! rd = 0.015*vm   ! used by Collatz et al.

  ! net photosynthesis assimilation rate
  NetPhoto = max(0.0, agross - rd)

end function NetPhoto

!**********************************************************************************************************************!
! function gs_medlyn - calculates the stomatal conductance for water vapor, mol/m2-s
!
! Uses the formulation of ...
!  Medlyn, B. E., et al. (2011) Reconciling the optimal and empirical approaches to modelling stomatal conductance,
!  Global Change Biology, 17, 2134-2144, doi: 10.1111/j.1365-2486.2010.02375.x.
!
!**********************************************************************************************************************!
function gs_medlyn(tki, relhumi, anet, cci)
  real(kind=dp), intent(in)  :: tki              ! air temperature, K
  real(kind=dp), intent(in)  :: relhumi          ! relative humidity
  real(kind=dp), intent(in)  :: anet             ! net photosynthesis assimilation rate, umol/m2-s
  real(kind=dp), intent(out) :: cci              ! intercellular CO2 concentration, umol/mol 
  real(kind=dp)              :: gs_medlyn        ! stomatal conductance for water vapor, mol/m2-s
  real(kind=dp)              :: vpd              ! vapor pressure deficit
  real(kind=dp)              :: gs

  vpd = esat(tki)*(1.0-relhumi/100.)
  if (vpd > 0.0) then
     cci = cca*(g1/(g1 + sqrt(vpd)))
     gs = g0 + 1.6*anet*(1.0 + g1/sqrt(vpd))/cca
  else
     cci = cca
     gs = gsmax
  end if

  gs_medlyn = gs
  return

end function gs_medlyn

!**********************************************************************************************************************!
! subroutine PartitionRAD - partitions measured global and incoming PAR into diffuse and direct beams
!
! Uses the algorithm of ...
!   Weiss & Norman (1985) Ag. & Forest Meteor., 34, 205-213.
!
! with corrections included in CANOAK by Baldocchi
!
!   Eq (3) of Weiss & Norman corrected to
!      Rdv = 0.4*(600*cos(theta) - RDV)
!
!   Eq (5) of Weiss & Norman corrected to
!      RdN = 0.6*(720 - RDN/cos(theta) - w)*cos(theta)
!
!   and adjusting all of the equations as necessary for a solar constant value of
!   1373 W/m2 instead of the 1320 W/m2 used by Weiss & Norman
!
!**********************************************************************************************************************!
subroutine PartitionRAD()
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer, parameter :: i4=kind(1)
  real(kind=dp)              :: cosza           ! cosine of zenith angle
  real(kind=dp)              :: rdpar           ! potential direct beam PAR (W/m2)
  real(kind=dp)              :: rspar           ! potential diffuse PAR (W/m2)
  real(kind=dp)              :: rdnir           ! potential direct beam NIR (W/m2)
  real(kind=dp)              :: rsnir           ! potential diffuse NIR (W/m2)
  real(kind=dp)              :: rtpar           ! potential total PAR (W/m2)
  real(kind=dp)              :: rtnir           ! potential total NIR (W/m2)
  real(kind=dp)              :: wa              ! water absorption in NIR
  real(kind=dp)              :: radratio        ! ratio of actual to potential total radiation
  real(kind=dp)              :: rr              ! temp variable for radratio
  real(kind=dp)              :: nirx            ! IR inferred from sradzref and ppfdzref (W/m2)
  real(kind=dp)              :: rfac            ! computation variable
  real(kind=dp)              :: fdpar           ! calculated direct beam fraction of PAR
  real(kind=dp)              :: fspar           ! calculated diffuse fraction of PAR
  real(kind=dp)              :: fdnir           ! calculated direct beam fraction of NIR 
  real(kind=dp)              :: fsnir           ! calculated diffuse fraction of NIR 
  real(kind=dp)              :: pc              ! pressure and optical factor
  real(kind=dp), parameter   :: coszamin=0.017  ! cos(89deg)

  cosza=cos(zarad)
  if (cosza < coszamin) then
    ppfd_direct=0.0
    ppfd_diffus=0.0
    nir_direct=0.0
    nir_diffus=0.0
    return
  end if
 
  pc = (pmbzref/1013.25)/cosza

  ! potential direct beam PAR 
  rdpar = 624.0*exp(-0.185*pc)*cosza

  ! potential diffuse PAR
  rspar = 0.4*(624.0*cosza - rdpar)

  ! water absorption in NIR for 10 mm precipitable water
  ! (expression used here is that used in CANOAK, but is different than that given
  ! in Weiss & Norman)
  wa = 1373.0*0.077*(2.0*pc)**0.30

  ! potential direct beam NIR
  rdnir = (748.0*exp(-0.06*pc) - wa)*cosza
  
  ! potential diffuse NIR
  rsnir = 0.6*(748.0 - (rdnir/cosza) - wa)*cosza

  ! total potential PAR & NIR
  rtpar = rdpar + rspar
  rtnir = rdnir + rsnir

  ! ratio of observed and potential total radiation
  radratio = sradzref/(rtpar+rtnir)

  ! measured IR
  nirx = sradzref - (ppfdzref/4.6)

  ! PAR direct and diffuse (umol/m2-s)
  if (radratio > 0.9) then
    rr = 0.9
  else
    rr = radratio
  endif 
  rfac = 1.0 - ((0.9 - rr)/0.7)**0.667
  fdpar = (rdpar/rtpar)*rfac

  if (fdpar < 0.0) fdpar=0.0
  if (fdpar > 1.0) fdpar=1.0

  fspar = 1.0 - fdpar

  ppfd_direct = fdpar*ppfdzref
  ppfd_diffus = fspar*ppfdzref

  if (ppfd_direct < 0.0) then
    ppfd_direct = 0.0
    ppfd_diffus = ppfdzref
  end if

  if (ppfdzref == 0.0) then
    ppfd_direct = 0.001
    ppfd_diffus = 0.001
  end if

  ! NIR direct and diffuse (W/m2)
  if (radratio > 0.88) then
    rr = 0.88
  else
    rr = radratio
  end if
  rfac = 1.0 - ((0.88 - rr)/0.68)**0.667
  fdnir = (rdnir/rtnir)*rfac

  if (fdnir < 0.0) fdnir=0.0
  if (fdnir > 1.0) fdnir=1.0

  fsnir = 1.0 - fdnir

  nir_direct = fdnir*nirx
  nir_diffus = fsnir*nirx

  if (nir_direct < 0.0) then
    nir_direct = 0.0
    nir_diffus = nirx
  endif

  if (nirx <= 0.0) then
    nir_direct = 0.1
    nir_diffus = 0.1
  endif

  return 

end subroutine PartitionRAD

!**********************************************************************************************************************!
! function kbe - canopy extinction coefficient
!
! Uses the ellipsoidal leaf distribution function of ...
!   Campbell and Norman (1998) An Introduction to Environmental Biophysics, 
!    2nd Ed., Springer, New York. (p. 251)
!
!**********************************************************************************************************************!
function kbe(x, za)
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer, parameter :: i4=kind(1)
  real(kind=dp), intent(in)  :: x             ! leaf angle distribution parameter (canopy specific)
  real(kind=dp), intent(in)  :: za            ! zenith angle (radians)
  real(kind=dp)              :: kbe           ! canopy extinction coefficient
  real(kind=dp)              :: numer, denom  ! tmps

  numer = (x*x + dtan(za)*dtan(za))**0.5
  denom = x + 1.774*((x + 1.182)**(-0.733))

  kbe = numer/denom

  return
end function kbe

!**********************************************************************************************************************!
! subroutine CalcRadProfiles - compute radiation profiles within the canopy
!
! Uses algorithms of Bodin & Franklin (2012) Efficient modeling of sun/shade canopy radiation dynamics
!                    explicitly accounting for scattering, Geoscientific Model Development, 5, 535-541. 
!   ... with elliptical kb from Campbell & Norman (1998) pp. 247-259
!   ... clear-sky effective emissivity of Prata (1996) Q.J.R. Meteorol. Soc., 122, 1127-1151.
!   ... cloudy-sky correction algorithm of Crawford and Duchon (1999) J. Appl. Meteorol., 48, 474-480.
!   ... based on evaluations of Flerchinger et al. (2009) Water Resour. Res., 45, W03423, doi:10.1029/2008WR007394.
!
!**********************************************************************************************************************!
subroutine CalcRadProfiles()
  integer(kind=i4)           :: i
  real(kind=dp)              :: kbexza       ! canopy extinction coefficient for given x and za values
  real(kind=dp)              :: sky_ir       ! effective sky thermal radiation (W/m2)
  real(kind=dp)              :: grnd_ir      ! ground surface thermal radiation (W/m2)
  real(kind=dp)              :: epssky       ! effective emissivity for sky thermal radiation
  real(kind=dp)              :: w                           ! precipitable water (cm)
  real(kind=dp), parameter   :: rfl_par=0.11                ! leaf reflectance coefficient for PAR
  real(kind=dp), parameter   :: trns_par=0.16               ! leaf transmission coefficient for PAR
  real(kind=dp), parameter   :: sigma_par=rfl_par+trns_par  ! leaf scattering coefficient for PAR
  real(kind=dp), parameter   :: rfl_nir=0.43                ! leaf reflectance coefficient for NIR
  real(kind=dp), parameter   :: trns_nir=0.26               ! leaf transmission coefficient for NIR
  real(kind=dp), parameter   :: sigma_nir=rfl_nir+trns_nir  ! leaf scattering coefficient for NIR
  real(kind=dp), parameter   :: alfpar=0.85                 ! leaf absorptance for PAR
  real(kind=dp), parameter   :: alfnir=0.25                 ! leaf absorptance for NIR
  real(kind=dp), parameter   :: alfir=1.0                   ! leaf absortpance for IR
  real(kind=dp)              :: taudt_ir, s    
  real(kind=dp)              :: kdsig_par, kdrfl_par, kdtrs_par, rho_par
  real(kind=dp)              :: kdsig_nir, kdrfl_nir, kdtrs_nir, rho_nir
  real(kind=dp)              :: rdpar, rdnir, rsuppar, rsupnir, rsdnpar, rsdnnir

  ! canopy extinction coefficient for ellipitical leaf angle distribution
  kbexza = kbe(x, zarad)

  ! effective downwelling sky thermal radiation  
  w = 4.65*eatm/tk(npts)
  epssky = 1.0 - (1 + w)*dexp(-(1.2 + 3.*w)**0.5)
  s = (((ppfd_direct+ppfd_diffus)/4.6) + nir_direct + nir_diffus)/1000.0
  epssky = (1.0-s) + s*epssky
  sky_ir = epssky*sbsig*tk(npts)**4.0

  ! longwave upwelling radiation
  grnd_ir = epsgrnd*sbsig*tsoilk**4.0
! grnd_ir = epsgrnd*sbsig*tk(1)**4.0
  lw_up(1) = grnd_ir
  do i=2,npts
    if (lai(i) > 0.0) then
      taudt_ir = dexp(-kd*lai(i))
      lw_up(i) = lw_up(i-1)*taudt_ir + (1.0-taudt_ir)*epsleaf*sbsig*tk(i-1)**4.0
    else
      lw_up(i) = lw_up(i-1)
    end if
  end do

  ! diffuse extinction coefficients and canopy reflectances
  kdsig_par = kd/dsqrt(1.0-sigma_par)
  kdrfl_par = kd/dsqrt(1.0-rfl_par)
  kdtrs_par = kd/dsqrt(1.0-trns_par)
  rho_par = ((1.0-dsqrt(1.0-sigma_par))/(1.0+dsqrt(1.0-sigma_par)))*(2.0/(1.0+1.6*dcos(zarad)))

  kdsig_nir = kd/dsqrt(1.0-sigma_nir)
  kdrfl_nir = kd/dsqrt(1.0-rfl_nir)
  kdtrs_nir = kd/dsqrt(1.0-trns_nir)
  rho_nir = ((1.0-dsqrt(1.0-sigma_nir))/(1.0+dsqrt(1.0-sigma_nir)))*(2.0/(1.0+1.6*dcos(zarad)))

  ! loop over the domain from the top
  do i=npts,1,-1

    ! within or above canopy?
    if (i > ncnpy) then                             ! above canopy is all sunny
      ppfd_sun(i) = ppfd_direct + ppfd_diffus
      ppfd_shd(i) = 0.0
      nir_sun(i) = nir_direct + nir_diffus
      nir_shd(i) = 0.0
      rt_sun(i)  = ppfd_sun(i)/4.6 + nir_sun(i)
      rt_shd(i)  = 0.0

      ! day or night?
      if (rt_sun(i) > 0.0) then
        fsun(i) = 1.0
        fshd(i) = 0.0
      else
        fsun(i) = 0.0
        fshd(i) = 1.0
      end if

      ! no canopy, so no absorption
      rabs_sun(i)= 0.0
      rabs_shd(i)= 0.0

      lw_dn(i) = sky_ir
    else                                             ! within canopy
      
      ! day or night?
      if (rt_sun(ncnpy+1) > 0.0) then
        fsun(i) = dexp(-kbexza*clai(i))
        fshd(i) = 1.0 - fsun(i)
      else
        fsun(i) = 0.0
        fshd(i) = 1.0
      end if

      ! downwelling radiation
      if (lai(i) > 0.0) then
        taudt_ir = dexp(-kd*lai(i))
        lw_dn(i) = lw_dn(i+1)*taudt_ir + (1.0-taudt_ir)*epsleaf*sbsig*tk(i+1)**4.0
      else
        lw_dn(i) = lw_dn(i+1)
      end if

      ! diffuse radiation
      rdpar = ppfd_diffus*(1.0 - rho_par)*dexp(-kd*clai(i))
      rdnir = nir_diffus*(1.0 - rho_nir)*dexp(-kd*clai(i))

      ! upwelling scattered radiation
      rsuppar = ppfd_direct*rfl_par*(dexp(-kbexza*clai(i)) - dexp(kd*clai(i)-(kbexza+kd)*clai(1)))/(kd+kbexza)
      rsupnir = nir_direct*rfl_nir*(dexp(-kbexza*clai(i)) - dexp(kd*clai(i)-(kbexza+kd)*clai(1)))/(kd+kbexza)

      ! downwelling scattered radiation
      rsdnpar = ppfd_direct*trns_par*(dexp(-kbexza*clai(i)) - dexp(-kd*clai(i)))/(kd-kbexza)
      rsdnnir = nir_direct*trns_nir*(dexp(-kbexza*clai(i)) - dexp(-kd*clai(i)))/(kd-kbexza)

      ! sunlit leaves
      ppfd_sun(i) = kdsig_par*rdpar + kdrfl_par*rsuppar + kdtrs_par*rsdnpar + kbexza*ppfd_direct
      nir_sun(i)  = kdsig_nir*rdnir + kdrfl_nir*rsupnir + kdtrs_nir*rsdnnir + kbexza*nir_direct
      rt_sun(i)   = alfpar*ppfd_sun(i)/4.6 + alfnir*nir_sun(i)
      rabs_sun(i) = rt_sun(i) + alfir*(lw_up(i) + lw_dn(i))
      
      ! shaded leaves
      ppfd_shd(i) = kdsig_par*rdpar + kdrfl_par*rsuppar + kdtrs_par*rsdnpar
      nir_shd(i)  = kdsig_nir*rdnir + kdrfl_nir*rsupnir + kdtrs_nir*rsdnnir
      rt_shd(i)   = alfpar*ppfd_shd(i)/4.6 + alfnir*nir_shd(i)
      rabs_shd(i) = rt_shd(i) + alfir*(lw_up(i) + lw_dn(i))
    end if
  end do

  return 
end subroutine CalcRadProfiles

!**********************************************************************************************************************!
! CalcWeightedProfiles - calculate sun/shade weighted canopy profiles
!
!**********************************************************************************************************************!
subroutine CalcWeightedProfiles()
  integer(kind=i4) :: i

  do i=npts,1,-1
    ppfd_wgt(i) = fsun(i)*ppfd_sun(i) + fshd(i)*ppfd_shd(i)
    nir_wgt(i)  = fsun(i)*nir_sun(i)  + fshd(i)*nir_shd(i)
    rt_wgt(i)   = fsun(i)*rt_sun(i)   + fshd(i)*rt_shd(i)
    rabs_wgt(i) = fsun(i)*rabs_sun(i) + fshd(i)*rabs_shd(i)
    tl_wgt(i)   = fsun(i)*tl_sun(i)   + fshd(i)*tl_shd(i)
    gs_wgt(i)   = fsun(i)*gs_sun(i)   + fshd(i)*gs_shd(i)
    if (i <= ncnpy) then
      rs_sun(i)   = gtor(gs_sun(i), pmb(i), tk(i))
      rs_shd(i)   = gtor(gs_shd(i), pmb(i), tk(i))
    end if
    rs_wgt(i)   = fsun(i)*rs_sun(i)   + fshd(i)*rs_shd(i)
    anet_wgt(i) = fsun(i)*anet_sun(i) + fshd(i)*anet_shd(i)
  end do

  return
end subroutine CalcWeightedProfiles

!**********************************************************************************************************************!
! function rbgl - calculate surface boundary layer resistance
!
! Uses formulation of ...
!    Schuepp (1977) BLM, 12, 171-186.
!
!**********************************************************************************************************************!
function rbgl(mdiffl, ubarh)
  real(kind=dp)             :: rbgl           ! surface boundary layer resistance
                                              ! (s/cm)
  real(kind=dp), intent(in) :: mdiffl         ! molecular diffusivity of species 
                                              ! in air (cm2/s)
  real(kind=dp), intent(in) :: ubarh          ! mean wind speed at canopy top
                                              ! (cm/s)
  real(kind=dp), parameter  :: viscair=0.155  ! dynamic viscosity of air (cm2/s)
  real(kind=dp)             :: ugstar, del0 

  ugstar = 0.05*ubarh+0.00005*ubarh*ubarh 
  del0 = viscair/(0.4*ugstar)
  rbgl = ((viscair/mdiffl)-dlog(del0/10.0))/(0.4*ugstar)
  !print *, 'ubarh, mdiffl, rbgl=', ubarh, mdiffl, rbgl
  return
end function rbgl

!**********************************************************************************************************************!
! function rsoill - calculate resistance to diffusion through the soil
!
! Uses formulation of ...
!    Pleim et al. (2013) JGR, 118, 3794-3806.
!
! With data from ...
!    Rawls et al. (1982) Trans. ASAE, 25, 1316-1320.
!    Clapp and Hornberger (1978) Water Resources Res., 14, 601-604.
!
!**********************************************************************************************************************!
function rsoill(mdiffl, stheta, sattheta, rtheta, sbcoef, dsoil)
  real(kind=dp)             :: rsoill     ! resistance to diffusion 
                                          ! through the soil (s/cm)
  real(kind=dp), intent(in) :: mdiffl     ! molecular diffusivity of species in air (cm2/s)
  real(kind=dp), intent(in) :: stheta     ! volumetric soil water content (m3/m3) 
  real(kind=dp), intent(in) :: sattheta   ! saturation volumetric soil water content (m3/m3)
  real(kind=dp), intent(in) :: rtheta     ! residual volumetric soil water content (m3/m3)
  real(kind=dp), intent(in) :: sbcoef     ! Clapp and Hornberger exponent
  real(kind=dp), intent(in) :: dsoil      ! depth of topsoil (cm)
  real(kind=dp)             :: ldry
  real(kind=dp)             :: mdiffp
  real(kind=dp)             :: xe

  xe =(1.0-(stheta/sattheta))**5.0
  ldry = dsoil*(dexp(xe)-1.0)/1.7183
  mdiffp = mdiffl*sattheta*sattheta*(1.0-(rtheta/sattheta))**(2.0+3.0/sbcoef)
  rsoill = ldry/mdiffp
  return
end function rsoill

!**********************************************************************************************************************!
! function rbl - calculate leaf boundary resistance for trace species
!
! Uses the rb formulation of ...
!   Wu, Y, et al. (2003) JGR, 108, D1, 4013, doi:10.1029/2002JD002293.
! and approximates the friction velocity, ustar as 0.14*ubar as suggested in
!   Weber, R. (1999) Boundary-Layer Met., 93, 197-209.
!
!**********************************************************************************************************************!
function rbl(mdiffl, ubari)
  real(kind=dp)             :: rbl      ! leaf boundary layer resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl   ! molecular diffusivity of species in air (cm2/s)
  real(kind=dp), intent(in) :: ubari    ! mean wind speed (cm/s)
  rbl = 10.53_dp/((mdiffl**0.666667)*ubari)
  return
end function rbl

!**********************************************************************************************************************!
! function rw_massad - calculate cuticular resistance for ammonia
!
! Uses formulation of ...
!   Massad et al. (2010) Atmos. Chem. Phys., 10, 10359-10386.
!
!**********************************************************************************************************************!
function rw_massad(relhumi, tki, ar, laii)
  real(kind=dp)             :: rw_massad  ! ammonia cuticular resistance (s/cm)
  real(kind=dp), intent(in) :: relhumi    ! relative humidity (%)
  real(kind=dp), intent(in) :: tki        ! temperature
  real(kind=dp), intent(in) :: ar         ! acid ratio = [2SO2+HNO3+HCl]/[NH3]
                                          ! fixed at neutral for now
  real(kind=dp), intent(in) :: laii       ! leaf area index (either total or layer)
  real(kind=dp), parameter  :: a = 0.148  ! exp coeff for RH dependence for 
                                          ! arable land
  real(kind=dp)             :: rwmin
  real(kind=dp)             :: rwuncorr

  rwmin = 0.315/ar
  rwuncorr = rwmin*dexp(a*(100.0_dp-relhumi))
  rw_massad = rwuncorr*dexp(0.15*(tki-273.15))/(laii**0.5)
  return
end function rw_massad

!**********************************************************************************************************************!
! function rw_flechard - calculate cuticular resistance for ammonia
!
! Uses formulation of ...
!   Flechard et al. (2010) Biogeosciences, 7, 537-556.
!
!**********************************************************************************************************************!
function rw_flechard(relhumi, tki)
  real(kind=dp)             :: rw_flechard ! ammonia cuticular resistance (s/cm)
  real(kind=dp), intent(in) :: relhumi     ! relative humidity (%)
  real(kind=dp), intent(in) :: tki         ! temperature
  real(kind=dp), parameter  :: a = 0.11    ! exp coeff for RH dependence
  real(kind=dp), parameter  :: b = 0.15    ! exp coeff for T dependence (degC-1)
  real(kind=dp), parameter  :: rwmin=0.10  ! minimum rw (s/cm)
  real(kind=dp), parameter  :: rwmax=12.0  ! maximum rw (s/cm)
  real(kind=dp)             :: rwf

  rwf = rwmin*dexp(a*(100.0_dp-relhumi))*dexp(b*abs(tki-273.15))
  rw_flechard = min(rwmax, rwf)

  return
end function rw_flechard

!**********************************************************************************************************************!
! function rs_zhang_nh3 - calculate stomatal resistance for ammonia
!
! Uses formulation of ...
!   Zhang, L. et al. (2003) Atmos. Chem. Phys., 3, 2067-2082.
!   Zhang, L. et al. (2002) Atmos. Environ., 36, 537-560.
!
!**********************************************************************************************************************!
function rs_zhang_nh3(mdiffl, tki, pmbi, ppfdi, srad, relhumi, rsmin, rsnight)
  real(kind=dp)             :: rs_zhang_nh3 ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl       ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki          ! temperature
  real(kind=dp), intent(in) :: pmbi         ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi        ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: srad         ! solar irradiation (W/m2)
  real(kind=dp), intent(in) :: relhumi      ! relative humidity (%)
  real(kind=dp), intent(in) :: rsmin        ! minimum leaf stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: rsnight      ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), parameter  :: brsp=297.0   ! empirical constant (umol/m2-s) for maize canopy
  real(kind=dp), parameter  :: tmin=5.0     ! temperature correction parameter-maize canopy
  real(kind=dp), parameter  :: tmax=45.0    ! temperature correction parameter-maize canopy
  real(kind=dp), parameter  :: topt=25.0    ! temperature correction parameter-maize canopy
  real(kind=dp), parameter  :: bvpd=0.0     ! empirical constant for VPD correction-maize canopy
  real(kind=dp), parameter  :: phic1=-1.5   ! empirical constant for water stress correction-maize canopy
  real(kind=dp), parameter  :: phic2=-2.5   ! empirical constant for water stress correction-maize canopy
  real(kind=dp), parameter  :: sradmin=0.0  ! nighttime threshold for srad (W/m2)
  real(kind=dp)             :: cfppfd
  real(kind=dp)             :: cft
  real(kind=dp)             :: cfvpd
  real(kind=dp)             :: cfphi
  real(kind=dp)             :: cfmdiff
  real(kind=dp)             :: tcel, ft1, ft2, et
  real(kind=dp)             :: vpd
  real(kind=dp)             :: phi
  real(kind=dp)             :: rsx

  if (srad .le. sradmin) then

     rsx = rsnight

  else

     ! PPFD correction
     cfppfd = 1.0_dp/(1.0_dp + brsp/ppfdi)

     ! temperature correction
     tcel=tki-273.15
     et=(tmax-topt)/(topt-tmin)
     ft1=(tcel-tmin)/(topt-tmin)
     ft2=(tmax-tcel)/(tmax-topt)
     cft=ft1*(ft2**et)

     ! water vapor pressure deficit correction
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     cfvpd = 1.0_dp - bvpd*vpd

     ! water stress correction
     phi = -0.72_dp - 0.0013_dp*srad

     if (phi > phic1) then
        cfphi=1.0_dp
     else
        cfphi = (phi-phic2)/(phic1-phic2) 
     end if

     ! species diffusivity correction
     cfmdiff=mdiffl/mdiffh2o(tki,pmbi)

     rsx = rsmin/(cfppfd*cft*cfvpd*cfphi*cfmdiff)
  end if

  rs_zhang_nh3 = rsx

  return
end function rs_zhang_nh3

!**********************************************************************************************************************!
! function rs_zhang - calculate stomatal resistance - general version
!
! Uses formulation of ...
!   Zhang, L. et al. (2003) Atmos. Chem. Phys., 3, 2067-2082.
!   Zhang, L. et al. (2002) Atmos. Environ., 36, 537-560.
!
!**********************************************************************************************************************!
function rs_zhang(mdiffl, tki, pmbi, ppfdi, srad, relhumi, rsmin, rsnight, &
  brsp, tmin, tmax, topt, bvpd, phic1, phic2)
  real(kind=dp)             :: rs_zhang     ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl       ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki          ! temperature
  real(kind=dp), intent(in) :: pmbi         ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi        ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: srad         ! solar irradiation (W/m2)
  real(kind=dp), intent(in) :: relhumi      ! relative humidity (%)
  real(kind=dp), intent(in) :: rsmin        ! minimum leaf stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: rsnight      ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: brsp         ! empirical constant (umol/m2-s) for PPFD correction
  real(kind=dp), intent(in) :: tmin         ! temperature correction parameter
  real(kind=dp), intent(in) :: tmax         ! temperature correction parameter
  real(kind=dp), intent(in) :: topt         ! temperature correction parameter
  real(kind=dp), intent(in) :: bvpd         ! empirical constant for VPD correction
  real(kind=dp), intent(in) :: phic1        ! empirical constant for water stress correction
  real(kind=dp), intent(in) :: phic2        ! empirical constant for water stress correction
  real(kind=dp), parameter  :: sradmin=0.0  ! nighttime threshold for srad (W/m2)
  real(kind=dp)             :: cfppfd
  real(kind=dp)             :: cft
  real(kind=dp)             :: cfvpd
  real(kind=dp)             :: cfphi
  real(kind=dp)             :: cfmdiff
  real(kind=dp)             :: tcel, ft1, ft2, et
  real(kind=dp)             :: vpd
  real(kind=dp)             :: phi
  real(kind=dp)             :: rsx

  if (srad .le. sradmin) then

     rsx = rsnight

  else

     ! PPFD correction
     cfppfd = 1.0_dp/(1.0_dp + brsp/ppfdi)

     ! temperature correction
     tcel=tki-273.15
     et=(tmax-topt)/(topt-tmin)
     ft1=(tcel-tmin)/(topt-tmin)
     ft2=(tmax-tcel)/(tmax-topt)
     cft=ft1*(ft2**et)

     ! water vapor pressure deficit correction
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     cfvpd = 1.0_dp - bvpd*vpd

     ! water stress correction
     phi = -0.72_dp - 0.0013_dp*srad

     if (phi > phic1) then
        cfphi=1.0_dp
     else
        cfphi = (phi-phic2)/(phic1-phic2) 
     end if

     ! species diffusivity correction
     cfmdiff=mdiffl/mdiffh2o(tki,pmbi)

     rsx = rsmin/(cfppfd*cft*cfvpd*cfphi*cfmdiff)

  endif

  rs_zhang = rsx

  return
end function rs_zhang

!**********************************************************************************************************************!
! function rs_rds - calculate stomatal resistance for ammonia
!
! Uses formulation of ...
!   Zhang, L. et al. (2003) Atmos. Chem. Phys., 3, 2067-2082.
!   Zhang, L. et al. (2002) Atmos. Environ., 36, 537-560.
!
!   but, with cfppfd based on the data of ...
!   Irmak and Mutiibwa (2009) Transactions of the ASABE, 52, 1923-1939.
!   Fig. 2.
!
!**********************************************************************************************************************!
function rs_rds(mdiffl, tki, pmbi, ppfdi, srad, relhumi, rsnight)
  real(kind=dp)             :: rs_rds      ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl      ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki         ! temperature
  real(kind=dp), intent(in) :: pmbi        ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi       ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: srad        ! solar irradiation (W/m2)
  real(kind=dp), intent(in) :: relhumi     ! relative humidity (%)
  real(kind=dp), intent(in) :: rsnight     ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), parameter  :: a0=33.00    ! 
  real(kind=dp), parameter  :: a1=-0.5     ! 
  real(kind=dp), parameter  :: tmin=5.0    ! temperature correction parameter-maize canopy
  real(kind=dp), parameter  :: tmax=45.0   ! temperature correction parameter-maize canopy
  real(kind=dp), parameter  :: topt=25.0   ! temperature correction parameter-maize canopy
  real(kind=dp), parameter  :: bvpd=0.0    ! empirical constant for VPD correction-maize canopy
  real(kind=dp), parameter  :: phic1=-1.5  ! empirical constant for water stress correction-maize canopy
  real(kind=dp), parameter  :: phic2=-2.5  ! empirical constant for water stress correction-maize canopy
  real(kind=dp), parameter  :: sradmin=0.0 ! nighttime threshold for srad (W/m2)
  real(kind=dp)             :: cfppfd
  real(kind=dp)             :: cfmdiff
  real(kind=dp)             :: cft
  real(kind=dp)             :: cfvpd
  real(kind=dp)             :: cfphi
  real(kind=dp)             :: tcel, ft1, ft2, et
  real(kind=dp)             :: vpd
  real(kind=dp)             :: phi
  real(kind=dp)             :: rsx

  if (srad .le. sradmin) then

     rsx = rsnight

  else

     ! PPFD correction
     cfppfd = 1.0_dp/(a0*ppfdi**a1)

     ! temperature correction
     tcel=tki-273.15
     et=(tmax-topt)/(topt-tmin)
     ft1=(tcel-tmin)/(topt-tmin)
     ft2=(tmax-tcel)/(tmax-topt)
     cft=ft1*(ft2**et)

     ! water vapor pressure deficit correction
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     cfvpd = 1.0_dp - bvpd*vpd

     ! water stress correction
     phi = -0.72_dp - 0.0013_dp*srad
     if (phi > phic1) then
        cfphi=1.0_dp
     else
        cfphi = (phi-phic2)/(phic1-phic2) 
     end if

     ! species diffusivity correction
     cfmdiff=mdiffl/mdiffh2o(tki,pmbi)
  
     rsx = 1.0_dp/(cfmdiff*cfppfd*cft*cfvpd*cfphi)
  end if

  rs_rds = rsx

  return
end function rs_rds

!**********************************************************************************************************************!
! function rs_jarvis1 - calculate stomatal resistance
!
! Uses formulation of ...
!  Jarvis (1976) with data set I parameters
!
!**********************************************************************************************************************!
function rs_jarvis1(mdiffl, tki, pmbi, ppfdi, relhumi, phi1, rsnight)
  real(kind=dp)             :: rs_jarvis1   ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl       ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki          ! temperature
  real(kind=dp), intent(in) :: pmbi         ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi        ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: relhumi      ! relative humidity (%)
  real(kind=dp), intent(in) :: phi1         ! water leaf potential (MPa)
  real(kind=dp), intent(in) :: rsnight      ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), parameter  :: b1=0.43      ! cm/s
  real(kind=dp), parameter  :: b2=0.0296    ! cm s-1/umol m-2 s-1
  real(kind=dp), parameter  :: b10=0.0      ! cm/s
  real(kind=dp), parameter  :: t1=-5.0      ! degC
  real(kind=dp), parameter  :: t0=9.0       ! degC
  real(kind=dp), parameter  :: th=35.0      ! degC
  real(kind=dp), parameter  :: b5=0.26      ! kPa-1
  real(kind=dp), parameter  :: b6=40.0      ! MPa-1
  real(kind=dp), parameter  :: phim=-2.4    ! MPa-1
  real(kind=dp), parameter  :: ppfdmin=0.0  ! nighttime threshold for ppfd (umol/m2-s)
  real(kind=dp)             :: gs1, gs2, gs3, gs4, gs5
  real(kind=dp)             :: q
  real(kind=dp)             :: tc
  real(kind=dp)             :: et, ft1, ft2
  real(kind=dp)             :: vpd
  real(kind=dp)             :: delphi
  real(kind=dp)             :: gs
  real(kind=dp)             :: rsx

  if (ppfdi .le. ppfdmin) then

     rsx = rsnight

  else

     q = b10/b1

     ! PPFD
     gs1 = (b1*b2*(ppfdi-q))/(b1 + b2*(ppfdi-q))

     ! temperature correction
     tc=min(th,tki-273.15)
     et=(th-t0)/(th-t1)
     ft1=(tc-t1)/(t0-t1)
     ft2=(th-tc)/(th-t0)
     gs2=ft1*(ft2**et)

     ! water vapor pressure deficit correction
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     gs3 = max(0.0_dp, 1.0_dp - b5*vpd)

     ! leaf potential
     delphi = phi1 - phim
     gs4 = 1.0_dp - dexp(-b6*delphi)

     ! diffusivity correction
     gs5 = mdiffl/mdiffh2o(tki, pmbi)

     gs = gs1*gs2*gs3*gs4*gs5

     rsx = 1.0_dp/gs

  end if

  rs_jarvis1 = rsx

  return
end function rs_jarvis1

!**********************************************************************************************************************!
! function rs_jarvis2 - calculate stomatal resistance
!
! Uses formulation of ...
!  Jarvis (1976) with data set II parameters
!
!**********************************************************************************************************************!
function rs_jarvis2(mdiffl, tki, pmbi, ppfdi, relhumi, phi1, rsnight)
  real(kind=dp)             :: rs_jarvis2   ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl       ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki          ! temperature
  real(kind=dp), intent(in) :: pmbi         ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi        ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: relhumi      ! relative humidity (%)
  real(kind=dp), intent(in) :: phi1         ! water leaf potential (MPa)
  real(kind=dp), intent(in) :: rsnight      ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), parameter  :: b1=0.30      ! cm/s
  real(kind=dp), parameter  :: b2=0.0135    ! cm s-1/umol m-2 s-1
  real(kind=dp), parameter  :: b10=0.01     ! cm/s
  real(kind=dp), parameter  :: t1=-5.0      ! degC
  real(kind=dp), parameter  :: t0=20.0      ! degC
  real(kind=dp), parameter  :: th=45.0      ! degC
  real(kind=dp), parameter  :: b5=0.12      ! kPa-1
  real(kind=dp), parameter  :: b6=0.93      ! MPa-1
  real(kind=dp), parameter  :: phim=-2.4    ! MPa-1
  real(kind=dp), parameter  :: ppfdmin=0.0  ! nighttime threshold for ppfd (umol/m2-s)
  real(kind=dp)             :: gs1, gs2, gs3, gs4, gs5
  real(kind=dp)             :: q
  real(kind=dp)             :: tc
  real(kind=dp)             :: et, ft1, ft2
  real(kind=dp)             :: vpd
  real(kind=dp)             :: delphi
  real(kind=dp)             :: gs
  real(kind=dp)             :: rsx

  if (ppfdi .le. ppfdmin) then

     rsx = rsnight

  else

     q = b10/b1

     ! PPFD
     gs1 = (b1*b2*(ppfdi-q))/(b1 + b2*(ppfdi-q))

     ! temperature correction
     tc=tki-273.15
     et=(th-t0)/(th-t1)
     ft1=(tc-t1)/(t0-t1)
     ft2=(th-tc)/(th-t0)
     gs2=ft1*(ft2**et)

     ! water vapor pressure deficit correction
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     gs3 = 1.0_dp - b5*vpd

     ! leaf potential
     delphi = phi1 - phim
     gs4 = 1.0_dp - dexp(-b6*delphi)

     ! diffusivity correction
     gs5 = mdiffl/mdiffh2o(tki,pmbi)

     gs = gs1*gs2*gs3*gs4*gs5
  
     rsx = 1.0_dp/gs

  end if

  rs_jarvis2 = rsx

  return
end function rs_jarvis2

!**********************************************************************************************************************!
! function rs_jim - calculate stomatal resistance for ammonia
!
! Uses formulation of ...
!   Jarvis (1976) [as presented by Irmak & Mutiibwa (2009)]
!    with parameters optimized for maize by
!   Irmak and Mutiibwa (2009) Transactions of the ASABE, 52, 1923-1939.
!
!**********************************************************************************************************************!
function rs_jim(mdiffl, tki, pmbi, ppfdi, relhumi, rsnight)
  real(kind=dp)             :: rs_jim       ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl       ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki          ! temperature
  real(kind=dp), intent(in) :: pmbi         ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi        ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: relhumi      ! relative humidity (%)
  real(kind=dp), intent(in) :: rsnight      ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), parameter  :: b1=0.492     ! cm/s
  real(kind=dp), parameter  :: b2=0.0024    ! cm s-1/umol m-2 s-1
  real(kind=dp), parameter  :: b3=0.036     ! kPa-1
  real(kind=dp), parameter  :: b4=1.321     !
  real(kind=dp), parameter  :: b5= -0.0352  ! K
  real(kind=dp), parameter  :: tcoef= 2.0   !
  real(kind=dp), parameter  :: q=0.0        ! =b10/b1
  real(kind=dp), parameter  :: ppfdmin=0.0  ! nighttime threshold for ppfd (umol/m2-s)
  real(kind=dp)             :: f1           ! PPFD factor
  real(kind=dp)             :: f2           ! VPD factor
  real(kind=dp)             :: f3           ! T factor
  real(kind=dp)             :: f4           ! soil moisture factor
  real(kind=dp)             :: f5           ! diffusivity factor
  real(kind=dp)             :: vpd          ! vapor pressure deficit, kPa
  real(kind=dp)             :: rsx

  if (ppfdi .le. ppfdmin) then

     rsx = rsnight

  else

     ! PPFD
     f1 = (b1*b2*(ppfdi-q))/(b1 + b2*(ppfdi-q))

     ! VPD
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     f2 = (1.0_dp - b3*vpd)**b4

     ! T
     f3 = 1.0_dp - b5*((298.0_dp-tki)**tcoef)

     ! soil moisture
     f4 = 1.0_dp

     ! diffusivity
     f5 = mdiffl/mdiffh2o(tki,pmbi)

     rsx = 1.0_dp/(f1*f2*f3*f4*f5)

  end if

  rs_jim = rsx

  return
end function rs_jim

!**********************************************************************************************************************!
! function rs_nmj - calculate stomatal resistance for ammonia
!
! Uses Jarvis-type new modified formulation of ...
!   Irmak and Mutiibwa (2009) Transactions of the ASABE, 52, 1923-1939.
!    with parameters optimized for maize
!
!**********************************************************************************************************************!
function rs_nmj(mdiffl, tki, pmbi, ppfdi, relhumi, laitot, rsnight)
  real(kind=dp)             :: rs_nmj       ! stomatal resistance (s/cm)
  real(kind=dp), intent(in) :: mdiffl       ! molecular diffusivity of trace species in air (cm2/s)
  real(kind=dp), intent(in) :: tki          ! temperature
  real(kind=dp), intent(in) :: pmbi         ! air pressure (mb)
  real(kind=dp), intent(in) :: ppfdi        ! photosynthetic photon flux (umol/m2-s)
  real(kind=dp), intent(in) :: relhumi      ! relative humidity (%)
  real(kind=dp), intent(in) :: laitot       ! total canopy LAI (m2/m2)
  real(kind=dp), intent(in) :: rsnight      ! nighttime leaf stomatal resistance (s/cm)
  real(kind=dp), parameter  :: a0=32.56  
  real(kind=dp), parameter  :: a1=-0.343
  real(kind=dp), parameter  :: a2=0.035
  real(kind=dp), parameter  :: a3=2.8
  real(kind=dp), parameter  :: a4=0.025
  real(kind=dp), parameter  :: a5=-6.4
  real(kind=dp), parameter  :: rsmin=0.33   ! s/cm 
  real(kind=dp), parameter  :: ppfdmin=0.0  ! nighttime threshold for ppfd (umol/m2-s)
  real(kind=dp)             :: h1           ! PPFD factor
  real(kind=dp)             :: h2           ! VPD factor
  real(kind=dp)             :: h3           ! T factor
  real(kind=dp)             :: h4           ! diffusivity factor
  real(kind=dp)             :: h5           ! LAI factor
  real(kind=dp)             :: vpd          ! vapor pressure deficit, kPa
  real(kind=dp)             :: rsx

  if (ppfdi .le. ppfdmin) then

     rsx = rsnight

  else

     ! PPFD
     h1 = a0*(ppfdi**a1)

     ! VPD
     vpd = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
     h2 = (1.0_dp - a2*vpd)**a3

     ! T
     h3 = (1.0_dp - a4*(298.0_dp-tki))**a5
 
     ! diffusivity
     h4 = mdiffh2o(tki,pmbi)/mdiffl

     ! LAI correction
     h5 = rsmin*dexp(-laitot)

     rsx = (h1*h2*h3*h4)+h5 

  end if

  rs_nmj = rsx

end function rs_nmj

!**********************************************************************************************************************!
! function rincl - in-canopy resistance (s/cm)
!
!   Based on the parameterization of ...
!     Erisman et al. (1994) Atmos. Environ., 28, 2595-2607.
!
!**********************************************************************************************************************!
function rincl(ustar, hc, laitot)
  real(kind=dp)              :: rincl   ! in-canopy resistance (Pleim et al., 2013;
                                        !  Erisman et al., 1994)  (s/m)
  real(kind=dp), intent(in)  :: ustar   ! friction velocity (cm/s)
  real(kind=dp), intent(in)  :: hc      ! height of canopy (cm)
  real(kind=dp), intent(in)  :: laitot  ! total leaf area index (m3/m3)
  real(kind=dp), parameter   :: b=0.14  ! empirical constant (cm-1)

  rincl = b*laitot*hc/ustar

  return
end function rincl

!**********************************************************************************************************************!
! function gpstl - calculate stomatal compensation point for ammonia - ACCESS-NH3
!
!**********************************************************************************************************************!
!function gpstl(tkl,l)
!  real(kind=dp)             :: gpstl         ! stomatal comp point (molec/cm3)
!  real(kind=dp), intent(in) :: tkl           ! leaf temperature (K)
!  integer(kind=i4)          :: l
!
!  gpstl = 161500.0_dp*dexp(-10380.0_dp/tkl)*gammast(l)*navo/(tkl*1000.0_dp) 
!  return
!end function gpstl

!**********************************************************************************************************************!
! function gpl - calculate canopy compensation point for ammonia - ACCESS-NH3
!
!**********************************************************************************************************************!
function gpl(xi,xs,rs,rb,rw)
  real(kind=dp), intent(in) :: xi      ! ambient NH3 conc (molec/cm3)
  real(kind=dp), intent(in) :: xs      ! stomatal compensation point (molec/cm3)
  real(kind=dp), intent(in) :: rs      ! stomatal resistance (cm/s) 
  real(kind=dp), intent(in) :: rb      ! quasi-laminar boundary layer resistance (cm/s)
  real(kind=dp), intent(in) :: rw      ! cuticular resistance (cm/s)
  real(kind=dp)             :: gpl     ! canopy compensation point for ammonia (molec/cm3)
  real(kind=dp)             :: denom

  denom = rs*rw+rb*rw+rb*rs
  gpl = (rs*rw*xi + rb*rw*xs)/denom
  return
end function gpl

!**********************************************************************************************************************!
! function gpla - calculate stomatal compensation point for ACCESS-ORG trace species
!
!**********************************************************************************************************************!
function gpla(ilev, lsp)
  integer(kind=i4), intent(in) :: ilev
  integer(kind=i4), intent(in) :: lsp
  real(kind=dp)                :: gpla
  ! dummy stub for now
  gpla=0.0_dp
  return
end function gpla

!**********************************************************************************************************************!
! function gpstl - calculate stomatal compensation point for ammonia - AMM2L
!
!**********************************************************************************************************************!
function gpstl(tkl, gammastom)
  real(kind=dp)             :: gpstl           ! stomatal comp point (molec/cm3)
  real(kind=dp), intent(in) :: tkl             ! leaf temperature (K)
  real(kind=dp), intent(in) :: gammastom       ! stomatal emission potential

  gpstl = 161500.0_dp*dexp(-10380.0_dp/tkl)*gammastom*navo/(tkl*1000.0_dp) 
  return
end function gpstl

!**********************************************************************************************************************!
! function gpgl - calculate soil compensation point for ammonia - AMM2L
!
!**********************************************************************************************************************!
function gpgl(tks, gammasoil)
  real(kind=dp)             :: gpgl            ! soil comp point (molec/cm3)
  real(kind=dp), intent(in) :: tks             ! soil temperature (K)
  real(kind=dp), intent(in) :: gammasoil       ! soil emission potential

  gpgl = 161500.0_dp*dexp(-10380.0_dp/tks)*gammasoil*navo/(tks*1000.0_dp) 
  return
end function gpgl

end module CanopyPhysics
