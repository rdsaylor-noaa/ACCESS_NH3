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
!     Module:       Chem                                                                                               !
!                   Not currently used in ACCESS-NH3                                                                   !
!                                                                                                                      !
!     Description:  top-level module to handle chemistry                                                               !
!                                                                                                                      !
!======================================================================================================================!
module Chem
  use GlobalData
  use GasChem
  use DVODE_F90_M
  implicit none

  private fchemodes, jvode
  public IntegChemOneLevel

contains

!**********************************************************************************************************************!
! IntegChemOneLevel - integrates chemistry at one level from tin to tout
!                     using the DVODE software package (see documentation
!                     in MODULE dvode_f90_m)
!**********************************************************************************************************************!
subroutine IntegChemOneLevel(neqn, y, tin, tout, istate)
  integer(kind=i4) :: neqn
  real(kind=dp), dimension(neqn) :: y, w
  integer(kind=i4) :: itask, istate
  real(kind=dp)    :: tin, tout
  real(kind=dp)    :: tlin, tlout
  integer(kind=i4) :: i, n

  type (VODE_OPTS) :: options

  ! use local values of tin and tout, because these get changed
  ! within dvode_f90
  tlin=tin
  tlout=tout

  itask  = 1          ! normal computation of y at t=tout
  istate = 1          ! first call

  ! dense Jacobian = true
  ! absolute error tolerance = atol (defined in MODULE GlobalData)
  ! relative error tolerance = rtol (defined in MODULE GlobalData)
  ! MXSTEP and NZSWAG set to high values for stiff chemical system
  ! no user-supplied Jacobian is provided (tests showed that the
  ! internally generated Jacobian was much faster)
  options = SET_OPTS(DENSE_J=.TRUE.,ABSERR=atol, &
            RELERR=rtol,MXSTEP=500000,NZSWAG=50000, &
            USER_SUPPLIED_JACOBIAN=.FALSE.)

  call dvode_f90(fchemodes, neqn, y, tlin, tlout, itask, istate, options, j_fcn=jvode)

  if(istate /= 2) then
    write(*,*)
    write(*,*) 'DVODE - Fatal error!'
    write(*,'(a,i8)') '  DVODE returned istate = ', istate
    write(*,*)
    return
  end if

  return
end subroutine IntegChemOneLevel

!**********************************************************************************************************************!
! fchemodes - sets up the chemistry ODE system for VODE
!**********************************************************************************************************************!
subroutine fchemodes(neqn, ts, y, yp)
  integer(kind=i4) :: neqn
  real(kind=dp)    :: ts
  real(kind=dp), dimension(neqn) :: y, yp
  intent(in)  :: neqn, ts, y
  intent(out) :: yp

  ! TODO: When aqueous- and/or aerosol-phases are implemented, the coupled 
  !       chemistry will be handled here (somehow)
  !       For now, only gas-phase chemistry ODEs are set up

  call fgaschem(neqn, ts, y, yp)

  return
end subroutine fchemodes

!**********************************************************************************************************************!
! jvode - dummy routine for Jacobian calculation of ODE system
!          NOT USED, but has to be here for VODE
!**********************************************************************************************************************!
subroutine jvode (neqn, ts, y, ml, mu, pd, nrowpd)
  integer(kind=i4) :: neqn, ml, mu, nrowpd
  real(kind=dp)    :: ts
  real(kind=dp), dimension(neqn) :: y
  real(kind=dp), dimension(nrxn) :: k
  real(kind=dp), dimension(nrowpd, neqn) :: pd
  intent(in)    :: neqn, ts, y, ml, mu, nrowpd
  intent(inout) :: pd

  return
end subroutine jvode

end module Chem
