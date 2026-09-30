!>\file module_mp_aero_emissions.F90
!! This module provides aerosol emission for aerosol-aware microphysics.
!! Haiqin.Li@noaa.gov & jordan.schnell@noaa.gov 03/2026

module module_mp_aero_emissions

   use mpas_kind_types
   use mpas_constants,        only: cp
   use mpas_timer, only : mpas_timer_start, mpas_timer_stop
   use mpas_smoke_config

   implicit none

   private

   public :: mp_aero_emission

contains

  subroutine mp_aero_emission(em_dust,em_fire_oc,em_antho_oc,em_seas,       &
        dt, land, nwfa2d, nifa2d, rri, dz8w,                                &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte                                           )

    integer,intent(in):: ids,ide,jds,jde,kds,kde,                           &
                         ims,ime,jms,jme,kms,kme,                           &
                         its,ite,jts,jte,kts,kte

    real(RKIND),intent(in) :: dt
    real(RKIND),intent(in), dimension(ims:ime, jms:jme):: em_dust,em_fire_oc,em_antho_oc,em_seas,land
    real(RKIND),intent(in), dimension(ims:ime, kms:kme, jms:jme):: rri,dz8w
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme):: nwfa2d,nifa2d

    real(RKIND),dimension(ims:ime, jms:jme):: daero_emis_wfa,daero_emis_ifa

    real(RKIND), parameter :: pi  = 3.1415926
    ! -- aerosol diameter to caluclate wfa & ifa (m)
    real(RKIND), parameter :: mean_diameter1= 4.E-8, sigma1=1.8
    real(RKIND), parameter :: mean_diameter2= 1.E-6, sigma2=1.8
    !-- aerosol density to caluclate wfa & ifa (kg/m3)
    real(RKIND), parameter :: density_dust= 2.6e+3, density_sulfate=1.8e+3
    real(RKIND), parameter :: density_oc  = 1.4e+3, density_seasalt=2.2e+3

    real(RKIND) :: fact_wfa, fact_ifa

    integer i, j

    nwfa2d(:,:) = 0.
    nifa2d(:,:) = 0.
    fact_wfa = 1.e-9*6.0/pi*exp(4.5*log(sigma1)**2)/mean_diameter1**3
    fact_ifa = 1.e-9*6.0/pi*exp(4.5*log(sigma2)**2)/mean_diameter2**3

    do i=its, ite
    do j=jts, jte
     daero_emis_wfa(i,j) = (em_fire_oc(i,j)+em_antho_oc(i,j))/density_oc + em_seas(i,j)/density_seasalt
     daero_emis_ifa(i,j) = em_dust(i,j)/density_dust

     daero_emis_wfa(i,j) = daero_emis_wfa(i,j)*fact_wfa*rri(i,1,j)/dz8w(i,1,j)
     daero_emis_ifa(i,j) = daero_emis_ifa(i,j)*fact_ifa*rri(i,1,j)/dz8w(i,1,j)

     nwfa2d(i,j)=daero_emis_wfa(i,j)
     nifa2d(i,j)=daero_emis_ifa(i,j)

     !-- mimicking dry deposition
     if(land(i,j).eq.1)then
      nwfa2d(i,j)    = nwfa2d(i,j)*(1. - 0.10*dt/86400.)
      nifa2d(i,j)    = nifa2d(i,j)*(1. - 0.10*dt/86400.)
     else
      nwfa2d(i,j)    = nwfa2d(i,j)*(1. - 0.05*dt/86400.)
      nifa2d(i,j)    = nifa2d(i,j)*(1. - 0.05*dt/86400.)
     endif
     nwfa2d(i,j)     = MIN(9999.E6,nwfa2d(i,j))
     nifa2d(i,j)     = MIN(9999.E6,nifa2d(i,j))
    enddo
    enddo

  end subroutine mp_aero_emission

!> @}
  end module module_mp_aero_emissions
