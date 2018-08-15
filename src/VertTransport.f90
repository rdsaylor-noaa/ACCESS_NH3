!======================================================================================================================!
!                                                                                                                      !
!     Program:      ACCESS                                                                                             !
!                   Atmospheric Chemistry and Canopy Exchange Simulation System                                        !
!                   Full BVOC chemistry version                                                                        !
!                                                                                                                      !
!     Version:      3.0.0                                                                                              !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Last Update:  Dec 2017                                                                                           !
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
!     Module:       VertTransport                                                                                      !
!                                                                                                                      !
!     Description:  Contains routines for integrating the diffusion equation in the vertical domain                    !
!                                                                                                                      !
!======================================================================================================================!
module VertTransport
  use GlobalData
  use Utils
  implicit none

  public SubIntegrateVertTransportConstBC, SubIntegrateVertTransportSplitBC

contains

!**********************************************************************************************************************!
! subroutine  SubIntegrateVertTransportConstBC - integrates the diffusion equation for dttot
!                                                using nvxtot timesteps of dtstep timestep
!                                                Constant BCs at surface and domain top
!
!     Follows the recommendation of Venkatram (1993) to use the mass conservative
!           form of the vertical transport equation.
!
!     Venkatram, A. (1993) The parameterization of the vertical dispersion of a 
!      scalar in the atmospheric boundary layer, Atmos. Environ., 27A, 1963-1966. 
!**********************************************************************************************************************!
subroutine SubIntegrateVertTransportConstBC(phi, phip, kv, ca, cs, rho, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip, rho
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: ca, cs
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: rho2m1, rho2p1, rhoim1, rhoip1, rhohm1, rhohp1
  real(kind=dp)                    :: k2m1, k2p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dz2m1, dz2p1, dzim1, dzip1, dzhm1, dzhp1
  real(kind=dp)                    :: dz2, dzi, dznm1
  real(kind=dp), dimension(npts-2) :: a, b, c, r
  real(kind=dp), dimension(npts-2) :: x
  integer(kind=i4) :: nvx, nvxtot, i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !  Surface
    rho2m1 = 0.5*(rho(2)+rho(1))
    rho2p1 = 0.5*(rho(3)+rho(2))
    k2m1   = 0.5*(kv(2)+kv(1))
    k2p1   = 0.5*(kv(3)+kv(2))
    dz2m1  = z(2)-z(1)   
    dz2p1  = z(3)-z(2)
    dz2    = 0.5*(z(3)-z(1))
    a(1) = 0.0
    b(1) = (0.5*rho2m1*k2m1/(rho(2)*dz2m1)) + (dz2/dtstep) + ( 0.5*rho2p1*k2p1/(rho(2)*dz2p1))
    c(1) = -0.5*rho2p1*k2p1/(rho(3)*dz2p1)
    r(1) = rho2m1*k2m1*cs/(rho(1)*dz2m1) +  &
             (-0.5*rho2m1*k2m1/(rho(2)*dz2m1) + (dz2/dtstep) -0.5*rho2p1*k2p1/(rho(2)*dz2p1))*phip(2) &
           + 0.5*rho2p1*k2p1*phip(3)/(rho(3)*dz2p1)

    !  Interior points
    do k=2,npts-3
      i = k+1
      rhoim1 = 0.5*(rho(i)+rho(i-1))
      rhoip1 = 0.5*(rho(i+1)+rho(i))
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*rhoim1*kim1/(rho(i-1)*dzim1)
      b(k) = 0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) + 0.5*rhoip1*kip1/(rho(i)*dzip1)
      c(k) = -0.5*rhoip1*kip1/(rho(i+1)*dzip1)
      r(k) = 0.5*rhoim1*kim1*phip(i-1)/(rho(i-1)*dzim1) +  &
             (-0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) - 0.5*rhoip1*kip1/(rho(i)*dzip1))*phip(i)  &
             + 0.5*rhoip1*kip1*phip(i+1)/(rho(i+1)*dzip1)
    end do 

    !  Domain top
    rhohm1 = 0.5*(rho(npts-1)+rho(npts-2))
    rhohp1 = 0.5*(rho(npts)+rho(npts-1))
    khm1   = 0.5*(kv(npts-1)+kv(npts-2))
    khp1   = 0.5*(kv(npts)+kv(npts-1))
    dzhm1  = z(npts-1)-z(npts-2)
    dzhp1  = z(npts)-z(npts-1)
    dznm1  = 0.5*(z(npts)-z(npts-2))
    a(npts-2) = -0.5*rhohm1*khm1/(rho(npts-2)*dzhm1)
    b(npts-2) = 0.5*rhohm1*khm1/(rho(npts-1)*dzhm1) + (dznm1/dtstep) + 0.5*rhohp1*khp1/(rho(npts-1)*dzhp1)
    c(npts-2) = 0.0
    r(npts-2) = 0.5*rhohm1*khm1*phip(npts-2)/(rho(npts-2)*dzhm1) +  &
              (-0.5*rhohm1*khm1/(rho(npts-1)*dzhm1)+(dznm1/dtstep)-0.5*rhohp1*khp1/(rho(npts-1)*dzhp1))*phip(npts-1) &
              + rhohp1*khp1*ca/(rho(npts)*dzhp1)

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-2)

    do k=1,npts-2
      phi(k+1) = x(k)
    end do

    phi(1)    = cs
    phi(npts) = ca 

    phip = phi

  end do

  return
end subroutine SubIntegrateVertTransportConstBC


!**********************************************************************************************************************!
! subroutine  SubIntegrateVertTransportSplitBC - integrates the diffusion equation for dttot
!                                                using nvxtot timesteps of dtstep timestep
!                                                Constant BC at domain top, Flux BC at surface
!
!     Follows the recommendation of Venkatram (1993) to use the mass conservative
!           form of the vertical transport equation.
!
!     Venkatram, A. (1993) The parameterization of the vertical dispersion of a 
!      scalar in the atmospheric boundary layer, Atmos. Environ., 27A, 1963-1966. 
!**********************************************************************************************************************!
subroutine SubIntegrateVertTransportSplitBC(phi, phip, kv, ca, cs, vs, rho, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip, rho
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: ca, cs, vs, alf
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: rho1p1, rhoim1, rhoip1, rhohm1, rhohp1
  real(kind=dp)                    :: k1p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dzim1, dzip1, dzhm1, dzhp1, dznm1
  real(kind=dp)                    :: dz0, dz02, dz2, dzi
  real(kind=dp), dimension(npts-1) :: a, b, c, r
  real(kind=dp), dimension(npts-1) :: x
  integer(kind=i4)                 :: nvx, nvxtot
  integer(kind=i4)                 :: i, k

  dtstep = 1.0      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !  Surface
    rho1p1 = 0.5*(rho(2)+rho(1))
    k1p1   = 0.5*(kv(2)+kv(1))
    dz0    = z(2)-z(1)
    dz02   = 2*dz0
    alf    = vs*dz02/(rho(1)*kv(1))
    a(1) = 0.0
    b(1) = (0.5*kv(1)/dz0)*(1.0+alf*rho(1)) + (2.0*dz0/dtstep) + (0.5*rho1p1*k1p1/(rho(1)*dz0))
    c(1) = -0.5*rho(1)*kv(1)/(rho(2)*dz0) - 0.5*rho1p1*k1p1/(rho(2)*dz0)
    r(1) = rho(1)*kv(1)*alf*cs/dz0 +  &
             ((-0.5*kv(1)/dz0)*(1.0+alf*rho(1)) + (2.0*dz0/dtstep) - 0.5*rho1p1*k1p1/(rho(1)*dz0))*phip(1) &
           + (0.5*rho(1)*kv(1)/(rho(2)*dz0) + 0.5*rho1p1*k1p1/(rho(2)*dz0))*phip(2)

    !  Interior points
    do k=2,npts-2
      i = k
      rhoim1 = 0.5*(rho(i)+rho(i-1))
      rhoip1 = 0.5*(rho(i+1)+rho(i))
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*rhoim1*kim1/(rho(i-1)*dzim1)
      b(k) = 0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) + 0.5*rhoip1*kip1/(rho(i)*dzip1)
      c(k) = -0.5*rhoip1*kip1/(rho(i+1)*dzip1)
      r(k) = 0.5*rhoim1*kim1*phip(i-1)/(rho(i-1)*dzim1) +  &
             (-0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) - 0.5*rhoip1*kip1/(rho(i)*dzip1))*phip(i)  &
             + 0.5*rhoip1*kip1*phip(i+1)/(rho(i+1)*dzip1)
    end do 

    !  Domain top
    rhohm1 = 0.5*(rho(npts-1)+rho(npts-2))
    rhohp1 = 0.5*(rho(npts)+rho(npts-1))
    khm1   = 0.5*(kv(npts-1)+kv(npts-2))
    khp1   = 0.5*(kv(npts)+kv(npts-1))
    dzhm1  = z(npts-1)-z(npts-2)
    dzhp1  = z(npts)-z(npts-1)
    dznm1  = 0.5*(z(npts)-z(npts-2))
    a(npts-1) = -0.5*rhohm1*khm1/(rho(npts-2)*dzhm1)
    b(npts-1) = 0.5*rhohm1*khm1/(rho(npts-1)*dzhm1) + (dznm1/dtstep) + 0.5*rhohp1*khp1/(rho(npts-1)*dzhp1)
    c(npts-1) = 0.0
    r(npts-1) = 0.5*rhohm1*khm1*phip(npts-2)/(rho(npts-2)*dzhm1) +  &
              (-0.5*rhohm1*khm1/(rho(npts-1)*dzhm1)+(dznm1/dtstep)-0.5*rhohp1*khp1/(rho(npts-1)*dzhp1))*phip(npts-1) &
              + rhohp1*khp1*ca/(rho(npts)*dzhp1)

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-1)

    do k=1,npts-1
      phi(k) = x(k)
    end do

    phi(npts) = ca 

    phip = phi

  end do

  return
end subroutine SubIntegrateVertTransportSplitBC

end module VertTransport
!======================================================================================================================!
