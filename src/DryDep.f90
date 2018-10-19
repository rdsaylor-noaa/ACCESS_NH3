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
! Uses formulation as derived from ...
!    Sakaguchi & Zeng (2009) JGR, D01107, doi: 10.1029/2008JD010834.
!    Schuepp (1977) BLM, 12, 171-186.
!    Philip (1957) J. of Met., 14, 354-366.
!
! With data from ...
!    Rawls et al. (1982) Trans. ASAE, 25, 1316-1320.
!    Clapp and Hornberger (1978) Water Resources Res., 14, 601-604.
!
!**********************************************************************************************************************!
subroutine GetSoilDryDepExCoeffs()
  integer(kind=i4)          :: l
  real(kind=dp)             :: mdiffl         ! molecular diffusivity (cm2/s)

  ! soil exchange coeffs and soil compensation points
  do l=1,ninteg

    ! ground boundary layer resistance (Rbg, s/cm) calculated in CanopyPhysics
    ! and assumed to be invariant over species

    ! soil diffusion resistance (s/cm)
    mdiffl=MolecDiff(l, tsoilk, pmb(1))
    rsoill(l) = SoilResist(mdiffl)

    ! deposition velocity to ground surface (cm/s)
    vs(l)=1.0/(rbg+rsoill(l))       

    ! soil compensation point (molecules/cm3)
    csoil(l) = 161500.0_dp*dexp(-10380.0_dp/tsoilk)*gammaso(l)*navo/(tsoilk*1000.0_dp) 
  end do

  return
end subroutine GetSoilDryDepExCoeffs

end module DryDep
