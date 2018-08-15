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
!     Module:       DryDep                                                                                             !
!                                                                                                                      !
!     Description:  contains dry deposition velocity algorithms                                                        !
!                                                                                                                      !
!======================================================================================================================!
module DryDep
  use GlobalData
  use CanopyPhysics
  use PhysChemData
  use Utils
  implicit none

  private
  public GetDryDepExCoeffs, GetSoilDryDepExCoeffs

contains

!**********************************************************************************************************************!
! subroutine GetDryDepExCoeffs - calculate leaf-scale dry deposition velocities
!                                and compensation points
!**********************************************************************************************************************!
subroutine GetDryDepExCoeffs()
  integer(kind=i4)          :: i, l
  real(kind=dp)             :: mdiffl       ! molecular diffusivity of species l in air (cm2/s)
  real(kind=dp)             :: hstarl       ! effective Henry's Law coefficient (M/atm)
  real(kind=dp)             :: f0l          ! Wesley's reactivity parameter (dimensionless)
  real(kind=dp)             :: srad         ! solar radiation at canopy top - W/m2
  real(kind=dp)             :: relhumi      ! relative humidity (%)
  real(kind=dp), parameter  :: rsmin=0.4    ! minimum stomatal resistance (s/cm)
  real(kind=dp), parameter  :: rsnight=25.  ! nighttime stomatal resistance (s/cm)

  do l=1,ninteg
    do i=1,npts
      if(lad(i) .gt. 0.0) then               
        ! calculate vd and gp in canopy

        ! molecular diffusivity (cm2/s)
        mdiffl=MolecDiff(l, tk(i), pmb(i))

        ! relative humidity at level i
        relhumi=RelativeHumidity(tk(i),pmb(i),qh(i))

        ! leaf boundary layer resistance (s/cm)
        rb(i,l)=rbl(mdiffl, ubar(i)) 

        ! leaf cuticular resistance (s/cm)
        rw(i,l)=rw_flechard(relhumi, tk(i))

        ! leaf stomatal resistance (s/cm)
!       rs(i,l)=rs_zhang_nh3(mdiffl,tk(i),pmb(i),ppfd(i),sradzref,relhumi,rsmin,rsnight)
        rs(i,l)=(mdiffh2o(tk(i),pmb(i))/mdiffl)*rs_wgt(i)*0.01     ! convert rs_wgt from s/m to s/cm

        ! canopy exchange velocity
        vd(i,l) = 1.0_dp/rb(i,l)

        ! stomatal compensation point (molecules/cm3)
        gpst(i,l) = gpstl(tk(i),gammast(i,l))

        ! canopy compensation point (molecules/cm3)
        gp(i,l) = gpl(cint(i,l),gpst(i,l),rs(i,l),rb(i,l),rw(i,l))

      else                                   
        ! out of the canopy
        rb(i,l) = 0.0_dp
        rw(i,l) = 0.0_dp
        rs(i,l) = 0.0_dp
        vd(i,l) = 0.0_dp
        gpst(i,l) = 0.0_dp
        gp(i,l) = 0.0_dp
      end if

    end do
  end do

  return
end subroutine GetDryDepExCoeffs

!**********************************************************************************************************************!
! subroutine GetSoilDryDepExCoeffs - calculate dry deposition velocities for
!                                    deposition to the ground surface
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
subroutine GetSoilDryDepExCoeffs()
  integer(kind=i4)          :: l
  real(kind=dp)             :: tks            ! soil temperature (K)
  real(kind=dp)             :: mdiffl         ! molecular diffusivity (cm2/s)
  real(kind=dp), parameter  :: viscair=0.155  ! dynamic viscosity of air (cm2/s)
  real(kind=dp)             :: ugstar, del0 
  real(kind=dp)             :: ldry, mdiffp, xe

  tks = tk(1)  ! for now, maybe better later?

  ! soil exchange coeffs and soil compensation points
  do l=1,ninteg

    ! ground boundary layer resistance (s/cm)
    mdiffl=MolecDiff(l, tks, pmb(1))
    ugstar = 0.05*ubar(ncnpy+1)+0.00005*ubar(ncnpy+1)*ubar(ncnpy+1) 
    del0 = viscair/(0.4*ugstar)
    rbg(l) = ((viscair/mdiffl)-dlog(del0/10.0))/(0.4*ugstar)

    ! soil diffusion resistance (s/cm)
    xe =(1.0-(stheta/sattheta))**5.0
    ldry = dsoil*(dexp(xe)-1.0)/1.7183
    mdiffp = mdiffl*sattheta*sattheta*(1.0-(rtheta/sattheta))**(2.0+3.0/sbcoef)
    rsoil(l) = ldry/mdiffp

    ! deposition velocity to ground surface (cm/s)
    vs(l)=1.0/(rbg(l)+rsoil(l))       

    ! soil compensation point (molecules/cm3)
    csoil(l) = 161500.0_dp*dexp(-10380.0_dp/tks)*gammaso(l)*navo/(tks*1000.0_dp) 
  end do

  return
end subroutine GetSoilDryDepExCoeffs

end module DryDep
