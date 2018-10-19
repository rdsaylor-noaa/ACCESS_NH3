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
!     Module:       Emissions                                                                                          !
!                   Not used in ACCESS-NH3                                                                             !
!                                                                                                                      !
!     Description:  Contains emissions algorithms                                                                      !
!                                                                                                                      !
!======================================================================================================================!
module Emissions
  use GlobalData
  use Utils
  implicit none

  public GetEmissions

contains

!**********************************************************************************************************************!
! subroutine GetEmissions - calculate emission rates for all emitted species
!                           q - molecules species/cm3-s
!**********************************************************************************************************************!
subroutine GetEmissions()
  integer(kind=i4) :: i, l

  ! NH3 emissions are handled via a soil compensation point

  do l=1,ninteg
    do i=1,npts
      q(i,l)=0.0_dp
    end do
  end do

  return
end subroutine GetEmissions

end module Emissions
