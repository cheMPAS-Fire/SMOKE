!! Simple CO-proxy SOA production for cheMPAS-Fire
!! This code has been written by Minsu Choi, CIRES/NOAA GSL
!! For questions, please contact smoke.gsl@noaa.gov, minsu.choi@noaa.gov
!! This module forms SOA from smoke CO using
!! a prescribed OH surrogate + precursors and the Hodzic and Jimenez CO method.
!! Minsu Choi, CIRES/NOAA GSL, May 19, 2026
module module_simple_soa

  use mpas_kind_types
  use mpas_smoke_init
  use mpas_smoke_config

  implicit none

contains

  subroutine simple_soa(dtstep, chem, num_chem, swdown, coszen, &
                               p_co, p_soa, &
                               ids, ide, jds, jde, kds, kde,           &
                               ims, ime, jms, jme, kms, kme,           &
                               its, ite, jts, jte, kts, kte)

    implicit none


    integer, intent(in) :: ids, ide, jds, jde, kds, kde
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    integer, intent(in) :: its, ite, jts, jte, kts, kte
    integer, intent(in) :: p_co
    integer, intent(in) :: p_soa
    integer, intent(in) :: num_chem
    real(RKIND), intent(in) :: dtstep
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem), intent(inout) :: chem
    real(RKIND), dimension(ims:ime, jms:jme), intent(in) :: swdown, coszen

    integer :: i, j, k

    real(RKIND), parameter :: oh_ref = 1.5e6_RKIND  ! molecules/cm3
    real(RKIND), parameter :: oh_night = 1.3e4_RKIND  ! molecules/cm3, Y Lu et al. 1992. Chemical and Physical Meteorology
    real(RKIND), parameter :: koh = 1.25e-11_RKIND

    ! Simple diagnostic SOA yield using excess CO as proxy
    real(RKIND), parameter :: soa_co_yld = 0.08_RKIND

    real(RKIND), parameter :: sw_ref = 800.0_RKIND
    real(RKIND), parameter :: sw_min = 10.0_RKIND

    ! Background CO threshold for now. 50 ppbv = 50e-9 mol/mol
    ! Zero for now, b/c of anthro CO
    real(RKIND), parameter :: co_bg = 0.0
    !real(RKIND), parameter :: co_bg = 50.0e-9_RKIND

    ! Safety cap for one timestep SOA production, unit is ug/kg
    real(RKIND), parameter :: soa_cap = 20.0_RKIND

    ! Maximum bbSOA mass mixing ratio, unit is ug/kg
    real(RKIND), parameter :: soa_max = 5.0e3_RKIND

    real(RKIND), parameter :: mw_co = 28.0_RKIND
    real(RKIND), parameter :: mw_air = 28.97_RKIND

    real(RKIND) :: light_frac
    real(RKIND) :: oh_eff
    real(RKIND) :: aging_frac
    real(RKIND) :: co_excess
    real(RKIND) :: co_aged
    real(RKIND) :: soa_add

    if (p_co < 1 .or. p_co > num_chem) return
    if (p_soa < 1 .or. p_soa > num_chem) return

    do j = jts, jte
      do i = its, ite
        ! using swdown, suppressing OH radical at night
        ! simple method, Minsu Choi CIRES/NOAA GSL
        ! during daytime, oh_eff ranges from ~1.9e5 to 1.5e6 molecules cm-3
        ! Applying nighttime OH concentration
        light_frac = max(0.0_RKIND, &
                min(1.0_RKIND, swdown(i,j) / sw_ref))

        if (coszen(i,j) > 0.0_RKIND) then
          oh_eff = max(oh_night, oh_ref * light_frac)
        else
          oh_eff = oh_night
        endif
! Old version, keep it as reference        
!        if (swdown(i,j) > sw_min .and. coszen(i,j) > 0.0_RKIND) then
!          oh_eff = oh_ref * min(1.0_RKIND, swdown(i,j) / sw_ref)
!        else
!          oh_eff = 0.0_RKIND
!        endif

        ! From Hodzic and Jimenez, 2011
        aging_frac = 1.0_RKIND - exp(-koh * oh_eff * dtstep)

        if (aging_frac <= 0.0_RKIND) cycle

        do k = kts, kte

          ! Use only CO above background. Do not modify real CO.
          co_excess = max(chem(i,k,j,p_co) - co_bg, 0.0_RKIND)

          if (co_excess <= 0.0_RKIND) cycle

          co_aged = aging_frac * co_excess
          ! mol/mol CO -> kg/kg CO -> ug/kg SOA proxy
          soa_add = soa_co_yld * co_aged * (mw_co / mw_air) * 1.0e9_RKIND
          soa_add = min(soa_add, soa_cap)
          chem(i,k,j,p_soa) = chem(i,k,j,p_soa) + soa_add
          chem(i,k,j,p_soa) = min(max(chem(i,k,j,p_soa), epsilc), soa_max)

        enddo
      enddo
    enddo

  end subroutine simple_soa

  ! This version of SOA subroutine is directly reading the ebu, ebb carbon monoxide
  ! and calculate SOA formation. 
  ! Mainly designed for separately tracking Anthropogenic SOA, and Biomass Burning SOA
  ! This subroutine needs more version update as AC branch at NOAA GSL adding explicit VOCs from emission
  subroutine simple_soa_voc(dtstep, chem, num_chem, swdown, coszen, &
                          rho_phy, dz8w,                         &
                          e_bb_co, e_ant_co,                     &
                          p_bbvoc, p_antvoc,                     &
                          p_bbsoa, p_antsoa,                     &
                          ids, ide, jds, jde, kds, kde,          &
                          ims, ime, jms, jme, kms, kme,          &
                          its, ite, jts, jte, kts, kte)

  implicit none

  integer, intent(in) :: ids, ide, jds, jde, kds, kde
  integer, intent(in) :: ims, ime, jms, jme, kms, kme
  integer, intent(in) :: its, ite, jts, jte, kts, kte
  integer, intent(in) :: num_chem
  integer, intent(in) :: p_bbvoc, p_antvoc
  integer, intent(in) :: p_bbsoa, p_antsoa

  real(RKIND), intent(in) :: dtstep
  real(RKIND), dimension(ims:ime,kms:kme,jms:jme,1:num_chem), intent(inout) :: chem
  real(RKIND), dimension(ims:ime,jms:jme), intent(in) :: swdown, coszen
  real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: rho_phy, dz8w

  ! BB CO emission flux: mol m-2 s-1
  real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_bb_co
  ! Anthropogenic CO emission increment from e_ant_out: ppmv
  real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_ant_co

  integer :: i, j, k

  real(RKIND), parameter :: oh_ref      = 1.5e6_RKIND
  real(RKIND), parameter :: oh_night    = 1.3e4_RKIND  ! molecules/cm3, Y Lu et al. 1992. Chemical and Physical Meteorology
  real(RKIND), parameter :: koh         = 1.25e-11_RKIND
  real(RKIND), parameter :: soa_co_yld  = 0.08_RKIND
  real(RKIND), parameter :: sw_ref      = 800.0_RKIND
  real(RKIND), parameter :: sw_min      = 10.0_RKIND
  real(RKIND), parameter :: voc_max     = 5.0e3_RKIND
  real(RKIND), parameter :: soa_max     = 5.0e3_RKIND
  real(RKIND), parameter :: mw_co       = 28.0_RKIND      ! g mol-1
  real(RKIND), parameter :: mw_air      = 0.02897_RKIND   ! kg mol-1

  real(RKIND) :: light_frac
  real(RKIND) :: oh_eff
  real(RKIND) :: aging_frac
  real(RKIND) :: co_add
  real(RKIND) :: voc_add
  real(RKIND) :: soa_add

  if (p_bbvoc  < 1 .or. p_bbvoc  > num_chem) return
  if (p_bbsoa  < 1 .or. p_bbsoa  > num_chem) return

  do j = jts, jte
  do i = its, ite

     light_frac = max(0.0_RKIND, &
             min(1.0_RKIND, swdown(i,j) / sw_ref))

     if (coszen(i,j) > 0.0_RKIND) then
       oh_eff = max(oh_night, oh_ref * light_frac)
     else
       oh_eff = oh_night
     endif
     aging_frac = 1.0_RKIND - exp(-koh * oh_eff * dtstep)
!     if (swdown(i,j) > sw_min .and. coszen(i,j) > 0.0_RKIND) then
!        oh_eff = oh_ref * min(1.0_RKIND, swdown(i,j) / sw_ref)
!     else
!        oh_eff = 0.0_RKIND
!     endif
     do k = kts, kte
        ! Wildfire CO emission -> transported bbVOC precursor
        if (e_bb_co(i,k,j) > 0.0_RKIND) then
           co_add = dtstep * mw_air * e_bb_co(i,k,j) / &
                    (rho_phy(i,k,j) * dz8w(i,k,j))

           voc_add = soa_co_yld * co_add * &
                     ((mw_co * 1.0e-3_RKIND) / mw_air) * 1.0e9_RKIND

           chem(i,k,j,p_bbvoc) = chem(i,k,j,p_bbvoc) + voc_add
           chem(i,k,j,p_bbvoc) = min(max(chem(i,k,j,p_bbvoc), epsilc), voc_max)
        endif

        ! Transported bbVOC precursor -> bbSOA
        if (aging_frac > 0.0_RKIND .and. chem(i,k,j,p_bbvoc) > epsilc) then
           soa_add = aging_frac * chem(i,k,j,p_bbvoc)

           chem(i,k,j,p_bbvoc) = chem(i,k,j,p_bbvoc) - soa_add
           chem(i,k,j,p_bbvoc) = max(chem(i,k,j,p_bbvoc), epsilc)

           chem(i,k,j,p_bbsoa) = chem(i,k,j,p_bbsoa) + soa_add
           chem(i,k,j,p_bbsoa) = min(max(chem(i,k,j,p_bbsoa), epsilc), soa_max)
        endif

        ! Anthropogenic CO increment is already in ppmv
        if (p_antvoc > 0 .and. p_antvoc <= num_chem .and. &
            e_ant_co(i,k,j) > 0.0_RKIND) then
        
           ! ppmv -> mol/mol
           co_add = e_ant_co(i,k,j) * 1.0e-6_RKIND
        
           ! mol/mol CO -> ug/kg CO-equivalent SOA precursor
           voc_add = soa_co_yld * co_add * &
                     ((mw_co * 1.0e-3_RKIND) / mw_air) * 1.0e9_RKIND
        
           chem(i,k,j,p_antvoc) = chem(i,k,j,p_antvoc) + voc_add
           chem(i,k,j,p_antvoc) = min(max(chem(i,k,j,p_antvoc), epsilc), voc_max)
        
        endif

        ! Transported antVOC precursor -> antSOA
        if (p_antvoc > 0 .and. p_antvoc <= num_chem .and. &
            p_antsoa > 0 .and. p_antsoa <= num_chem) then

           if (aging_frac > 0.0_RKIND .and. chem(i,k,j,p_antvoc) > epsilc) then
              soa_add = aging_frac * chem(i,k,j,p_antvoc)

              chem(i,k,j,p_antvoc) = chem(i,k,j,p_antvoc) - soa_add
              chem(i,k,j,p_antvoc) = max(chem(i,k,j,p_antvoc), epsilc)

              chem(i,k,j,p_antsoa) = chem(i,k,j,p_antsoa) + soa_add
              chem(i,k,j,p_antsoa) = min(max(chem(i,k,j,p_antsoa), epsilc), soa_max)
           endif
        endif

     enddo
  enddo
  enddo

  end subroutine simple_soa_voc

  !-----------------------------------------------------------------------
  ! Reduced two-product SOA scheme,
  ! Odum's 2p model for Secondary Organic Aerosols
  ! This version is not connected to wrapper yet, so in other words 
  ! it is not working. Will be option3 in the future.
  ! 
  !
  ! Input VOC emissions must already be converted to mol m-2 s-1.
  ! Raw GRA2PES style units of mole km-2 hr-1 should be converted before
  ! this subroutine is called.
  ! Version for NEMO has not been developed yet
  !
  ! Current mapping:
  !   ALK = HC03 + HC04 + HC05 + HC06
  !   BNZ = HC38
  !   TOL = HC41
  !   XYL = HC42 + HC43 + HC44
  !   ISO = HC10
  !   TRP = HC11
  !
  ! Product formation:
  !   VOC emission -> reacted VOC mass -> semivolatile products
  !   semivolatile products partition using Mo / (Mo + Cstar_T)
  !
  ! Output:
  !   p_soa can be antsoa for now. This is initial version
  ! Minsu Choi, CIRES/NOAA GSL, May 24, 2026
  !-----------------------------------------------------------------------
  subroutine simple_soa_2p(dtstep, chem, num_chem, swdown, coszen, &
                           rho_phy, dz8w, t_phy,                  &
                           e_hc03, e_hc04, e_hc05, e_hc06,        &
                           e_hc38, e_hc41,                        &
                           e_hc42, e_hc43, e_hc44,                &
                           e_hc10, e_hc11,                        &
                           p_soa,                                 &
                           p_smoke_fine, p_unspc_fine,            &
                           ids, ide, jds, jde, kds, kde,          &
                           ims, ime, jms, jme, kms, kme,          &
                           its, ite, jts, jte, kts, kte)

    implicit none

    integer, intent(in) :: ids, ide, jds, jde, kds, kde
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    integer, intent(in) :: its, ite, jts, jte, kts, kte
    integer, intent(in) :: num_chem
    integer, intent(in) :: p_soa
    integer, intent(in), optional :: p_smoke_fine, p_unspc_fine

    real(RKIND), intent(in) :: dtstep

    real(RKIND), dimension(ims:ime,kms:kme,jms:jme,1:num_chem), intent(inout) :: chem
    real(RKIND), dimension(ims:ime,jms:jme), intent(in) :: swdown, coszen
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: rho_phy, dz8w, t_phy

    ! Anthropogenic VOC emissions, mol m-2 s-1
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc03
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc04
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc05
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc06
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc38
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc41
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc42
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc43
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc44
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc10
    real(RKIND), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: e_hc11

    integer :: i, j, k

    ! OH surrogate
    real(RKIND), parameter :: oh_ref = 1.5e6_RKIND
    real(RKIND), parameter :: sw_ref = 800.0_RKIND
    real(RKIND), parameter :: sw_min = 10.0_RKIND

    ! Constants
    real(RKIND), parameter :: rgas = 8.314_RKIND
    real(RKIND), parameter :: tref = 298.0_RKIND
    real(RKIND), parameter :: dz_min = 1.0_RKIND
    real(RKIND), parameter :: rho_min = 0.1_RKIND
    real(RKIND), parameter :: mo_floor = 0.1_RKIND
    real(RKIND), parameter :: soa_cap = 100.0_RKIND
    real(RKIND), parameter :: soa_max = 5.0e3_RKIND

    ! Approximate organic mass fraction for partitioning mass proxy
    real(RKIND), parameter :: f_smoke_oa = 0.70_RKIND
    real(RKIND), parameter :: f_unspc_oa = 0.50_RKIND

    ! Effective OH rate constants, cm3 molecule-1 s-1
    ! These are intentionally simple tunable values for this reduced scheme.
    real(RKIND), parameter :: koh_alk = 1.0e-11_RKIND
    real(RKIND), parameter :: koh_bnz = 1.2e-12_RKIND
    real(RKIND), parameter :: koh_tol = 5.6e-12_RKIND
    real(RKIND), parameter :: koh_xyl = 2.3e-11_RKIND
    real(RKIND), parameter :: koh_iso = 1.0e-10_RKIND
    real(RKIND), parameter :: koh_trp = 5.3e-11_RKIND

    ! Surrogate molecular weights, g mol-1
    real(RKIND), parameter :: mw_alk = 170.0_RKIND
    real(RKIND), parameter :: mw_bnz = 78.11_RKIND
    real(RKIND), parameter :: mw_tol = 92.14_RKIND
    real(RKIND), parameter :: mw_xyl = 106.17_RKIND
    real(RKIND), parameter :: mw_iso = 68.12_RKIND
    real(RKIND), parameter :: mw_trp = 136.24_RKIND

    ! alpha = dimensionless mass yield.
    ! Cstar = ug m-3 at 298 K.
    ! Hvap = J mol-1.
    ! Current code is based on CMAQv5.2
    ! ALK
    real(RKIND), parameter :: a_alk1 = 0.0334_RKIND
    real(RKIND), parameter :: c_alk1 = 0.1472_RKIND
    real(RKIND), parameter :: h_alk1 = 53.0e3_RKIND
    real(RKIND), parameter :: a_alk2 = 0.2164_RKIND
    real(RKIND), parameter :: c_alk2 = 51.8775_RKIND
    real(RKIND), parameter :: h_alk2 = 53.0e3_RKIND

    ! xylene
    real(RKIND), parameter :: a_xyl1 = 0.0310_RKIND
    real(RKIND), parameter :: c_xyl1 = 1.3140_RKIND
    real(RKIND), parameter :: h_xyl1 = 32.0e3_RKIND
    real(RKIND), parameter :: a_xyl2 = 0.0900_RKIND
    real(RKIND), parameter :: c_xyl2 = 34.4830_RKIND
    real(RKIND), parameter :: h_xyl2 = 32.0e3_RKIND

    ! toluene
    real(RKIND), parameter :: a_tol1 = 0.0580_RKIND
    real(RKIND), parameter :: c_tol1 = 2.3260_RKIND
    real(RKIND), parameter :: h_tol1 = 18.0e3_RKIND
    real(RKIND), parameter :: a_tol2 = 0.1130_RKIND
    real(RKIND), parameter :: c_tol2 = 21.2770_RKIND
    real(RKIND), parameter :: h_tol2 = 18.0e3_RKIND

    ! benzene
    real(RKIND), parameter :: a_bnz1 = 0.0720_RKIND
    real(RKIND), parameter :: c_bnz1 = 0.3020_RKIND
    real(RKIND), parameter :: h_bnz1 = 18.0e3_RKIND
    real(RKIND), parameter :: a_bnz2 = 0.8880_RKIND
    real(RKIND), parameter :: c_bnz2 = 111.1100_RKIND
    real(RKIND), parameter :: h_bnz2 = 18.0e3_RKIND

    ! terpene
    real(RKIND), parameter :: a_trp1 = 0.1393_RKIND
    real(RKIND), parameter :: c_trp1 = 14.7920_RKIND
    real(RKIND), parameter :: h_trp1 = 40.0e3_RKIND
    real(RKIND), parameter :: a_trp2 = 0.4542_RKIND
    real(RKIND), parameter :: c_trp2 = 133.7297_RKIND
    real(RKIND), parameter :: h_trp2 = 40.0e3_RKIND

    ! isoprene
    real(RKIND), parameter :: a_iso1 = 0.2320_RKIND
    real(RKIND), parameter :: c_iso1 = 116.0100_RKIND
    real(RKIND), parameter :: h_iso1 = 40.0e3_RKIND
    real(RKIND), parameter :: a_iso2 = 0.0288_RKIND
    real(RKIND), parameter :: c_iso2 = 0.6170_RKIND
    real(RKIND), parameter :: h_iso2 = 40.0e3_RKIND

    real(RKIND) :: oh_eff
    real(RKIND) :: rho
    real(RKIND) :: dz
    real(RKIND) :: temp
    real(RKIND) :: mo
    real(RKIND) :: soa_add_ugm3
    real(RKIND) :: soa_add_ugkg

    real(RKIND) :: alk_flux, bnz_flux, tol_flux, xyl_flux, iso_flux, trp_flux
    real(RKIND) :: alk_ugm3, bnz_ugm3, tol_ugm3, xyl_ugm3, iso_ugm3, trp_ugm3

    if (p_soa < 1 .or. p_soa > num_chem) return

    do j = jts, jte
    do i = its, ite

       if (swdown(i,j) > sw_min .and. coszen(i,j) > 0.0_RKIND) then
          oh_eff = oh_ref * min(1.0_RKIND, swdown(i,j) / sw_ref)
       else
          oh_eff = 0.0_RKIND
       endif

       if (oh_eff <= 0.0_RKIND) cycle

       do k = kts, kte

          rho  = max(rho_phy(i,k,j), rho_min)
          dz   = max(dz8w(i,k,j), dz_min)
          temp = max(t_phy(i,k,j), 180.0_RKIND)

          ! Organic aerosol mass available for partitioning, ug m-3.
          mo = max(chem(i,k,j,p_soa) * rho, mo_floor)

          if (present(p_smoke_fine)) then
             if (p_smoke_fine >= 1 .and. p_smoke_fine <= num_chem) then
                mo = mo + f_smoke_oa * max(chem(i,k,j,p_smoke_fine), 0.0_RKIND) * rho
             endif
          endif

          if (present(p_unspc_fine)) then
             if (p_unspc_fine >= 1 .and. p_unspc_fine <= num_chem) then
                mo = mo + f_unspc_oa * max(chem(i,k,j,p_unspc_fine), 0.0_RKIND) * rho
             endif
          endif

          ! Emission, mol m-2 s-1.
          ! Again, we don't have this yet!!!!!!
          ! This will only work when cheMPAS-Fire read VOCs from GRA2PES
          alk_flux = max(e_hc03(i,k,j),0.0_RKIND) + max(e_hc04(i,k,j),0.0_RKIND) + &
                     max(e_hc05(i,k,j),0.0_RKIND) + max(e_hc06(i,k,j),0.0_RKIND)

          bnz_flux = max(e_hc38(i,k,j),0.0_RKIND)
          tol_flux = max(e_hc41(i,k,j),0.0_RKIND)

          xyl_flux = max(e_hc42(i,k,j),0.0_RKIND) + max(e_hc43(i,k,j),0.0_RKIND) + &
                     max(e_hc44(i,k,j),0.0_RKIND)

          iso_flux = max(e_hc10(i,k,j),0.0_RKIND)
          trp_flux = max(e_hc11(i,k,j),0.0_RKIND)

          ! Reacted precursor mass over this timestep, ug m-3.
          alk_ugm3 = reacted_mass_ugm3(alk_flux, mw_alk, koh_alk, oh_eff, dtstep, dz)
          bnz_ugm3 = reacted_mass_ugm3(bnz_flux, mw_bnz, koh_bnz, oh_eff, dtstep, dz)
          tol_ugm3 = reacted_mass_ugm3(tol_flux, mw_tol, koh_tol, oh_eff, dtstep, dz)
          xyl_ugm3 = reacted_mass_ugm3(xyl_flux, mw_xyl, koh_xyl, oh_eff, dtstep, dz)
          iso_ugm3 = reacted_mass_ugm3(iso_flux, mw_iso, koh_iso, oh_eff, dtstep, dz)
          trp_ugm3 = reacted_mass_ugm3(trp_flux, mw_trp, koh_trp, oh_eff, dtstep, dz)

          soa_add_ugm3 = 0.0_RKIND

          soa_add_ugm3 = soa_add_ugm3 + soa_2p(alk_ugm3, mo, temp, &
                                               a_alk1, c_alk1, h_alk1, &
                                               a_alk2, c_alk2, h_alk2)

          soa_add_ugm3 = soa_add_ugm3 + soa_2p(bnz_ugm3, mo, temp, &
                                               a_bnz1, c_bnz1, h_bnz1, &
                                               a_bnz2, c_bnz2, h_bnz2)

          soa_add_ugm3 = soa_add_ugm3 + soa_2p(tol_ugm3, mo, temp, &
                                               a_tol1, c_tol1, h_tol1, &
                                               a_tol2, c_tol2, h_tol2)

          soa_add_ugm3 = soa_add_ugm3 + soa_2p(xyl_ugm3, mo, temp, &
                                               a_xyl1, c_xyl1, h_xyl1, &
                                               a_xyl2, c_xyl2, h_xyl2)

          soa_add_ugm3 = soa_add_ugm3 + soa_2p(iso_ugm3, mo, temp, &
                                               a_iso1, c_iso1, h_iso1, &
                                               a_iso2, c_iso2, h_iso2)

          soa_add_ugm3 = soa_add_ugm3 + soa_2p(trp_ugm3, mo, temp, &
                                               a_trp1, c_trp1, h_trp1, &
                                               a_trp2, c_trp2, h_trp2)

          if (soa_add_ugm3 <= 0.0_RKIND) cycle

          ! ug m-3 -> ug kg-1
          soa_add_ugkg = soa_add_ugm3 / rho
          soa_add_ugkg = min(soa_add_ugkg, soa_cap)

          chem(i,k,j,p_soa) = chem(i,k,j,p_soa) + soa_add_ugkg
          chem(i,k,j,p_soa) = min(max(chem(i,k,j,p_soa), epsilc), soa_max)

       enddo
    enddo
    enddo

  contains
    ! These functions are directly imported from CMAQv5.2
    ! Code has been slightly re-organized by Minsu Choi, CIRES/NOAA GSL
    real(RKIND) function reacted_mass_ugm3(e_flux, mw, koh, oh, dt, dzloc)
      implicit none

      real(RKIND), intent(in) :: e_flux
      real(RKIND), intent(in) :: mw
      real(RKIND), intent(in) :: koh
      real(RKIND), intent(in) :: oh
      real(RKIND), intent(in) :: dt
      real(RKIND), intent(in) :: dzloc

      real(RKIND) :: reacted_frac

      if (e_flux <= 0.0_RKIND) then
         reacted_mass_ugm3 = 0.0_RKIND
         return
      endif

      reacted_frac = 1.0_RKIND - exp(-koh * oh * dt)

      ! mol m-2 s-1 * s / m = mol m-3
      ! mol m-3 * g mol-1 * 1e6 = ug m-3
      reacted_mass_ugm3 = reacted_frac * e_flux * dt * mw * 1.0e6_RKIND / dzloc

    end function reacted_mass_ugm3

    real(RKIND) function soa_2p(parent_ugm3, mo, temp, &
                                alpha1, cstar1, hvap1, &
                                alpha2, cstar2, hvap2)
      implicit none

      real(RKIND), intent(in) :: parent_ugm3
      real(RKIND), intent(in) :: mo
      real(RKIND), intent(in) :: temp
      real(RKIND), intent(in) :: alpha1, cstar1, hvap1
      real(RKIND), intent(in) :: alpha2, cstar2, hvap2

      real(RKIND) :: c1_t
      real(RKIND) :: c2_t
      real(RKIND) :: fp1
      real(RKIND) :: fp2

      if (parent_ugm3 <= 0.0_RKIND) then
         soa_2p = 0.0_RKIND
         return
      endif

      c1_t = cstar_temp(cstar1, hvap1, temp)
      c2_t = cstar_temp(cstar2, hvap2, temp)

      fp1 = mo / max(mo + c1_t, 1.0e-12_RKIND)
      fp2 = mo / max(mo + c2_t, 1.0e-12_RKIND)

      soa_2p = parent_ugm3 * (alpha1 * fp1 + alpha2 * fp2)

    end function soa_2p

    real(RKIND) function cstar_temp(cstar298, hvap, temp)
      implicit none

      real(RKIND), intent(in) :: cstar298
      real(RKIND), intent(in) :: hvap
      real(RKIND), intent(in) :: temp

      cstar_temp = cstar298 * exp((hvap / rgas) * (1.0_RKIND / tref - 1.0_RKIND / temp))
      cstar_temp = max(cstar_temp, 1.0e-12_RKIND)

    end function cstar_temp

  end subroutine simple_soa_2p

end module module_simple_soa     
