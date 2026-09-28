!
!  This module developed by Jordan Schnell (CIRES/NOAA GSL) following
!  Druge et al. (2019)
!  https://doi.org/10.5194/acp-19-3707-2019
!
!  For serious questions contact jordan.schnell@noaa.gov
!  Code further updated by Minsu Choi, CIRES/NOAA GSL

MODULE module_tactic_sna

  USE mpas_kind_types, ONLY : RKIND
  USE mpas_smoke_init, ONLY : p_unspc_fine, p_smoke_fine, p_nox, &
                              p_nh3, p_nh4_a_fine, p_so2,        &
                              p_so4_a_fine, p_no3_a_fine

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpas_smoke_tactic_sna_driver

  ! Unspecified (anthropogenic) fine PM left after the previous call.
  ! Used to speciate only the mass added since that call.
  REAL(RKIND), ALLOCATABLE, SAVE :: unspc_tmp(:,:,:)

CONTAINS

  SUBROUTINE mpas_smoke_tactic_sna_driver(              &
                         ktau, dt, chem, num_chem,       &
                         relhum, t_phy, dz8w, rho_phy,   &
                         nifa, nwfa, hno3_bkgd,          &
                         swdown, coszen,                 &
                         ids, ide, jds, jde, kds, kde,   &
                         ims, ime, jms, jme, kms, kme,   &
                         its, ite, jts, jte, kts, kte)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ktau, num_chem,               &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte

    REAL(RKIND), INTENT(IN) :: dt

    REAL(RKIND), DIMENSION(ims:ime, kms:kme, jms:jme), &
                           INTENT(IN) ::               &
                           relhum, t_phy, dz8w,        &
                           rho_phy, nifa, nwfa

    REAL(RKIND), DIMENSION(ims:ime, jms:jme), INTENT(IN) :: &
                           swdown, coszen

    REAL(RKIND), DIMENSION(ims:ime, kms:kme, jms:jme), &
                           INTENT(INOUT) ::               &
                           hno3_bkgd

    REAL(RKIND), DIMENSION(ims:ime, kms:kme, jms:jme, &
                           1:num_chem), INTENT(INOUT) :: chem

    INTEGER :: i, j, k

    ! Term is from Mozurkewich 1993 equilibrium work (M. Mozurkewich et al 1993).
    ! https://www.sciencedirect.com/science/article/abs/pii/0960168693903564
    REAL(RKIND) :: ta, ts, tn
    REAL(RKIND) :: drh, kp, rh, rh1
    REAL(RKIND) :: p1, p2, p3, pc
    REAL(RKIND) :: ta_star, gamma
    REAL(RKIND) :: nh4_so4, nh3_so4
    REAL(RKIND) :: nh4_no3
    REAL(RKIND) :: nh3_u, nh3_eq, hno3_eq
    REAL(RKIND) :: nh3_m, nh4_m, so4_m
    REAL(RKIND) :: hno3_m, no3_m
    REAL(RKIND) :: nh3_new, hno3_new
    REAL(RKIND) :: relax, discriminant
    REAL(RKIND) :: unspc_bulk, smoke_bulk
    REAL(RKIND) :: light_frac, oh_eff, nox_aging_frac
    REAL(RKIND) :: nox_old, nox_aged
    REAL(RKIND) :: so2_aging_frac, so2_old, so2_aged

    REAL(RKIND), PARAMETER :: mw_nh3 = 17.031_RKIND
    REAL(RKIND), PARAMETER :: mw_nh4 = 18.04_RKIND
    REAL(RKIND), PARAMETER :: mw_so4 = 96.06_RKIND
    REAL(RKIND), PARAMETER :: mw_so2 = 64.066_RKIND
    REAL(RKIND), PARAMETER :: mw_no3 = 62.0049_RKIND
    REAL(RKIND), PARAMETER :: mw_no2 = 46.0055_RKIND
    REAL(RKIND), PARAMETER :: mw_air = 28.97_RKIND

    ! Gas constant in nbar m3 umol-1 K-1 (8.3145 J mol-1 K-1).
    ! Used to convert the Mozurkewich kp (nbar**2) to (umol m-3)**2.
    REAL(RKIND), PARAMETER :: r_nbar = 0.083145_RKIND

    ! HNO3 background in ppbv.
    REAL(RKIND), PARAMETER :: hno3_bkgd_ppbv = 0.001_RKIND

    ! Simple OH-dependent conversion of NOx to HNO3.
    ! Same approach implemented in SOA module
    ! Please visit module_simple_soa.F90 for more detail
    REAL(RKIND), PARAMETER :: oh_ref = 1.5e6_RKIND
    REAL(RKIND), PARAMETER :: oh_night = 1.3e4_RKIND
    REAL(RKIND), PARAMETER :: k_no2_oh = 1.25e-11_RKIND
    REAL(RKIND), PARAMETER :: sw_ref = 800._RKIND

    ! SO2 + OH -> (eventually) H2SO4, gas-phase pathway only.
    ! JPL/IUPAC recommended effective bimolecular rate near the
    ! surface (~9e-13 cm3 molecule-1 s-1); the true SO2+OH+M
    ! reaction is pressure-dependent (three-body, low-pressure
    ! falloff), so this fixed value is a near-surface approximation.
    ! This captures only the gas-phase route: aqueous/in-cloud
    ! oxidation (often the dominant SO4 source when cloud water is
    ! present, e.g. GOCART's SO2+H2O2 term) is not represented here.
    REAL(RKIND), PARAMETER :: k_so2_oh = 9.0e-13_RKIND

    REAL(RKIND), PARAMETER :: tau = 300._RKIND
    REAL(RKIND), PARAMETER :: tiny_conc = 1.0e-30_RKIND

    ! This is the part that we speciate primary inorgnics from
    ! biomass burning, and anthropogenic emission using fractional
    ! contribution from NEMO, and field campaign
    REAL(RKIND), PARAMETER :: frac_unspc_no3 = 0.01_RKIND
    REAL(RKIND), PARAMETER :: frac_unspc_so4 = 0.01_RKIND
    REAL(RKIND), PARAMETER :: frac_unspc_nh4 = 0.01_RKIND
    REAL(RKIND), PARAMETER :: frac_unspc =                  &
                              1._RKIND - frac_unspc_no3 -   &
                              frac_unspc_so4 - frac_unspc_nh4

    REAL(RKIND), PARAMETER :: frac_bb_no3 = 0.03_RKIND
    REAL(RKIND), PARAMETER :: frac_bb_so4 = 0.01_RKIND
    REAL(RKIND), PARAMETER :: frac_bb_nh4 = 0.03_RKIND
    REAL(RKIND), PARAMETER :: frac_bb_smoke =               &
                              1._RKIND - frac_bb_no3 -      &
                              frac_bb_so4 - frac_bb_nh4

    if (p_nox        < 1 .or. p_nox        > num_chem) return
    if (p_nh3        < 1 .or. p_nh3        > num_chem) return
    if (p_nh4_a_fine < 1 .or. p_nh4_a_fine > num_chem) return
    if (p_so4_a_fine < 1 .or. p_so4_a_fine > num_chem) return
    if (p_no3_a_fine < 1 .or. p_no3_a_fine > num_chem) return

    ! Bookkeeping array for the unspecified PM (see unspc_tmp above).
    ! On (re)allocation it starts from the current field, so no mass is
    ! speciated retroactively; at ktau==1 it is zeroed cell by cell below.
    IF ( p_unspc_fine>=1 .AND. p_unspc_fine<=num_chem ) THEN
      IF ( ALLOCATED(unspc_tmp) ) THEN
        IF ( ANY(LBOUND(unspc_tmp)/=(/ims,kms,jms/)) .OR. &
             ANY(UBOUND(unspc_tmp)/=(/ime,kme,jme/)) ) DEALLOCATE(unspc_tmp)
      ENDIF
      IF ( .NOT.ALLOCATED(unspc_tmp) ) THEN
        ALLOCATE(unspc_tmp(ims:ime,kms:kme,jms:jme))
        unspc_tmp = max(0._RKIND,chem(:,:,:,p_unspc_fine))
      ENDIF
    ENDIF

    ! dz8w, nifa, and nwfa are retained for a future
    ! aerosol-number-dependent relaxation timescale.
    DO j = Jts , Jte
       DO k = Kts , Kte
          DO i = Its , Ite

             IF ( Rho_phy(i,k,j)<=TINY_CONC .OR. T_phy(i,k,j)<=TINY_CONC ) CYCLE
             IF ( i==Its .AND. k==Kts .AND. j==Jts ) WRITE (*,*) 'Debug:TACTIC, ktau = ' , Ktau
        ! Initialize
        ! Now, this version only has primary species from 0 tstep
        ! Otherwise, SNA concentrations accumulated unrealistically from smoke_fine over tstep
        ! This will be resolved when we put inline emission module for biomass burning
             IF ( Ktau==1 ) Hno3_bkgd(i,k,j) = HNO3_BKGD_PPBV*1.0E-9_RKIND

        ! Anthropogenic (unspecified) PM: speciate every step, but only the
        ! mass added since the previous call (unspc_tmp), so SNA does not
        ! accumulate from the remaining bulk. At ktau==1 this is the whole field.
             IF ( p_unspc_fine>=1 .AND. p_unspc_fine<=Num_chem ) THEN
               IF ( Ktau==1 ) Unspc_tmp(i,k,j) = 0._RKIND
               unspc_bulk = max(0._RKIND,Chem(i,k,j,p_unspc_fine)-Unspc_tmp(i,k,j))
               Chem(i,k,j,p_unspc_fine) = max(0._RKIND,Chem(i,k,j,p_unspc_fine)) - (1._RKIND-FRAC_UNSPC)*unspc_bulk
               Chem(i,k,j,p_no3_a_fine) = max(0._RKIND,Chem(i,k,j,p_no3_a_fine)) + FRAC_UNSPC_NO3*unspc_bulk
               Chem(i,k,j,p_so4_a_fine) = max(0._RKIND,Chem(i,k,j,p_so4_a_fine)) + FRAC_UNSPC_SO4*unspc_bulk
               Chem(i,k,j,p_nh4_a_fine) = max(0._RKIND,Chem(i,k,j,p_nh4_a_fine)) + FRAC_UNSPC_NH4*unspc_bulk
               Unspc_tmp(i,k,j) = Chem(i,k,j,p_unspc_fine)
             ENDIF

             IF ( Ktau==1 .AND. p_smoke_fine>=1 .AND. p_smoke_fine<=Num_chem ) THEN
               smoke_bulk = max(0._RKIND,Chem(i,k,j,p_smoke_fine))
               Chem(i,k,j,p_smoke_fine) = FRAC_BB_SMOKE*smoke_bulk
               Chem(i,k,j,p_no3_a_fine) = max(0._RKIND,Chem(i,k,j,p_no3_a_fine)) + FRAC_BB_NO3*smoke_bulk
               Chem(i,k,j,p_so4_a_fine) = max(0._RKIND,Chem(i,k,j,p_so4_a_fine)) + FRAC_BB_SO4*smoke_bulk
               Chem(i,k,j,p_nh4_a_fine) = max(0._RKIND,Chem(i,k,j,p_nh4_a_fine)) + FRAC_BB_NH4*smoke_bulk
             ENDIF

        ! Use the same simple diagnostic OH field as the SOA module.
             light_frac = max(0._RKIND,min(1._RKIND,Swdown(i,j)/SW_REF))

             IF ( Coszen(i,j)>0._RKIND ) THEN
               oh_eff = max(OH_NIGHT,OH_REF*light_frac)
             ELSE
               oh_eff = OH_NIGHT
             ENDIF

        ! Treat the lumped NOx tracer as NO2-equivalent for this
        ! initial parameterization. NOx is stored in ug kg-1.
        ! Transfer oxidized NOx to HNO3 on a one-to-one molar basis.
        ! This need to be discussed further for NO,NO2 speciation
             nox_aging_frac = 1._RKIND - exp(-K_NO2_OH*oh_eff*max(0._RKIND,Dt))

             nox_old = max(0._RKIND,Chem(i,k,j,p_nox))
             nox_aged = min(nox_old,nox_aging_frac*nox_old)

             Chem(i,k,j,p_nox) = max(0._RKIND,nox_old-nox_aged)
        ! hno3_bkgd is stored as mol mol-1. Convert the aged
        ! NO2-equivalent mass mixing ratio to a molar mixing ratio.
             Hno3_bkgd(i,k,j) = Hno3_bkgd(i,k,j) + nox_aged*1.0E-9_RKIND*MW_AIR/MW_NO2

        ! SO2 + OH -> H2SO4, condensing directly onto fine SO4.
        ! Same diagnostic OH proxy as the NOx aging term above.
        ! H2SO4 has negligible vapor pressure, so unlike NH4NO3 no
        ! equilibrium partitioning is needed: the oxidized mass is
        ! simply moved from SO2 to fine SO4, done before so4_m is
        ! read below so the added sulfate correctly competes for
        ! NH3 in the same step's equilibrium.
        ! Now keep in mind, gas-phase so2 chem not a significant contributor for SO4
             IF ( p_so2>=1 .AND. p_so2<=Num_chem ) THEN
               so2_aging_frac = 1._RKIND - exp(-K_SO2_OH*oh_eff*max(0._RKIND,Dt))

               so2_old = max(0._RKIND,Chem(i,k,j,p_so2))
               so2_aged = min(so2_old,so2_aging_frac*so2_old)

               Chem(i,k,j,p_so2) = max(0._RKIND,so2_old-so2_aged)
               Chem(i,k,j,p_so4_a_fine) = max(0._RKIND,Chem(i,k,j,p_so4_a_fine)) + (MW_SO4/MW_SO2)*so2_aged
             ENDIF

        ! NH3 is stored in ug kg-1. Convert it to numerical
        ! umol m-3 using the same convention as the working
        ! reference implementation.
             nh3_m = max(0._RKIND,Chem(i,k,j,p_nh3))*Rho_phy(i,k,j)/MW_NH3

        ! Persistent HNO3 is stored as mol mol-1. Convert it to
        ! numerical umol m-3.
             hno3_m = max(0._RKIND,Hno3_bkgd(i,k,j))*Rho_phy(i,k,j)*1.0E9_RKIND/MW_AIR

        ! Aerosols are in ug kg-1. The resulting numerical
        ! concentrations are umol m-3.
             nh4_m = max(0._RKIND,Chem(i,k,j,p_nh4_a_fine))*Rho_phy(i,k,j)/MW_NH4

             so4_m = max(0._RKIND,Chem(i,k,j,p_so4_a_fine))*Rho_phy(i,k,j)/MW_SO4

             no3_m = max(0._RKIND,Chem(i,k,j,p_no3_a_fine))*Rho_phy(i,k,j)/MW_NO3

        ! Total molar ammonia, sulfate, and nitrate.
             ta = nh3_m + nh4_m
             ts = so4_m
             tn = hno3_m + no3_m

             IF ( ts>ta ) THEN
               gamma = 1._RKIND
             ELSEIF ( 2._RKIND*ts>ta ) THEN
               gamma = 1.5_RKIND
             ELSE
               gamma = 2._RKIND
             ENDIF

             ta_star = max(0._RKIND,ta-gamma*ts)

             nh4_so4 = min(ta,gamma*ts)

             nh3_so4 = max(0._RKIND,nh4_so4-nh4_m)
             nh3_u = max(0._RKIND,nh3_m-nh3_so4)

             rh = max(0._RKIND,min(1._RKIND,Relhum(i,k,j)))

        ! Druge et al. give DRH in percent.
             drh = 0.01_RKIND*exp(723.7_RKIND/T_phy(i,k,j)+1.6954_RKIND)

             rh1 = 1._RKIND - rh

        ! Dry equilibrium constant.
             kp = exp(118.87_RKIND-24084._RKIND/T_phy(i,k,j)-6.025_RKIND*log(T_phy(i,k,j)))

        ! Humidity correction above the deliquescence RH.
             IF ( rh>=drh ) THEN
               p1 = exp(-135.94_RKIND+8763._RKIND/T_phy(i,k,j)+19.12_RKIND*log(T_phy(i,k,j)))

               p2 = exp(-122.65_RKIND+9969._RKIND/T_phy(i,k,j)+16.22_RKIND*log(T_phy(i,k,j)))

               p3 = exp(-182.61_RKIND+13875._RKIND/T_phy(i,k,j)+24.46_RKIND*log(T_phy(i,k,j)))

               pc = max(0._RKIND,(p1-p2*rh1+p3*rh1**2)*rh1**1.75_RKIND)

               kp = kp*pc
             ENDIF

        ! kp is a partial-pressure product in nbar**2 (Mozurkewich 1993).
        ! Convert it to (umol m-3)**2 so it matches tn*ta_star:
        ! p [nbar] = C [umol m-3] * R_NBAR * T  ->  kp_conc = kp / (R_NBAR*T)**2
             kp = kp/(R_NBAR*T_phy(i,k,j))**2

        ! Equilibrium ammonium nitrate concentration.
             IF ( tn*ta_star<=kp ) THEN
               nh4_no3 = 0._RKIND
             ELSE
               discriminant = max(0._RKIND,(ta_star+tn)**2-4._RKIND*(tn*ta_star-kp))

               nh4_no3 = 0.5_RKIND*(ta_star+tn-sqrt(discriminant))

               nh4_no3 = min(ta_star,min(tn,max(0._RKIND,nh4_no3)))
             ENDIF

             nh3_eq = max(0._RKIND,ta_star-nh4_no3)
             hno3_eq = max(0._RKIND,tn-nh4_no3)

        ! First-order relaxation toward equilibrium.
        ! This is a placeholder for future research
        ! Purpose of this relaxation time is calcualte it from aerosol num. conc.
        ! and also size, HNO3 uptake coefficient
        ! thus we need aerosol dynamics.
        ! With this  relaxation time, we should have more rapid equilibrium in highly polluted region
        ! while slower equilibrium in less polluted region
        ! Minsu Choi, CIRES/NOAA GSL
             relax = 1.0_RKIND
!             relax = 1._RKIND - exp(-max(0._RKIND,Dt)/TAU)

             nh3_new = nh3_u + relax*(nh3_eq-nh3_u)
             hno3_new = hno3_m + relax*(hno3_eq-hno3_m)

        ! Convert NH3 from umol m-3 to ug kg-1.
             Chem(i,k,j,p_nh3) = MW_NH3/Rho_phy(i,k,j)*nh3_new

        ! Convert persistent HNO3 from numerical umol m-3
        ! back to mol mol-1.
             Hno3_bkgd(i,k,j) = hno3_new*MW_AIR/(Rho_phy(i,k,j)*1.0E9_RKIND)

        ! Convert aerosol concentrations from umol m-3 to ug kg-1.
             Chem(i,k,j,p_nh4_a_fine) = MW_NH4/Rho_phy(i,k,j)*max(0._RKIND,ta-nh3_new)

             Chem(i,k,j,p_no3_a_fine) = MW_NO3/Rho_phy(i,k,j)*max(0._RKIND,tn-hno3_new)

        ! Sulfate is conserved by the SNA partitioning.
             Chem(i,k,j,p_so4_a_fine) = MW_SO4/Rho_phy(i,k,j)*ts
          ENDDO
       ENDDO
    ENDDO

   END SUBROUTINE mpas_smoke_tactic_sna_driver
END MODULE module_tactic_sna
