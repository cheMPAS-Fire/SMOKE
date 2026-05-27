MODULE dep_dry_mod_wesley


  ! --- MPAS and External Modules ---
  USE mpas_kind_types
  USE dep_data_mod
  USE mpas_smoke_config, ONLY : n_dbg_lines
  USE mpas_smoke_init
  USE mpas_timer, ONLY : mpas_timer_start, mpas_timer_stop

  ! ======================================================================
  ! MODULE: module_gas_deposition
  ! DESCRIPTION: Calculates dry deposition velocities for trace gases
  !              using the Wesely resistance analogy, SWDOWN for diurnal
  !              stomatal response, and MODIS lookup tables.
  !              Aerodynamic resistance (Ra) is calculated internally.
  ! ======================================================================

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: calc_gas_dep_vel

CONTAINS

  SUBROUTINE calc_gas_dep_vel(                   &
       num_chem,                                 &
       ustar, wspd, temp, press, swdown, ivgtyp, &  
       ddvel,                                    &
       ims, ime, jms, jme, its, ite, jts, jte    )


    ! --- WRF/MPAS Standard Indices ---
    INTEGER, INTENT(IN) :: ims, ime, jms, jme  
    INTEGER, INTENT(IN) :: its, ite, jts, jte  
    INTEGER, INTENT(IN) :: num_chem            

    ! --- Meteorological & Surface Inputs (2D Arrays) ---
    REAL(KIND=RKIND), DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ustar  ! Friction velocity (m/s)
    REAL(KIND=RKIND), DIMENSION(ims:ime, jms:jme), INTENT(IN) :: wspd   ! Wind speed at lowest level (m/s)
    REAL(KIND=RKIND), DIMENSION(ims:ime, jms:jme), INTENT(IN) :: temp   ! Surface temperature (K)
    REAL(KIND=RKIND), DIMENSION(ims:ime, jms:jme), INTENT(IN) :: press  ! Surface pressure (Pa)
    REAL(KIND=RKIND), DIMENSION(ims:ime, jms:jme), INTENT(IN) :: swdown ! Shortwave down (W/m2)
    INTEGER, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: ivgtyp          ! MODIS category (1-21)

    ! --- Output Array ---
    REAL(KIND=RKIND), DIMENSION(ims:ime, jms:jme, num_chem), INTENT(OUT) :: ddvel 

    ! --- Local Variables ---
    INTEGER  :: i, j, nv, lu
    REAL(KIND=RKIND) :: ra_loc, rb, rc, r_total       
    REAL(KIND=RKIND) :: r_stomatal, r_mesophyll, r_cuticle, r_ground
    REAL(KIND=RKIND) :: par         
    
    ! Land-use base specific variables retrieved from tables
    REAL(KIND=RKIND) :: rs_min_base 
    REAL(KIND=RKIND) :: rcut_base   
    REAL(KIND=RKIND) :: rgrnd_base  
    
    REAL(KIND=RKIND), PARAMETER :: r_min = 1.0e-4_RKIND 
    REAL(KIND=RKIND), PARAMETER :: swdown_thresh = 10.0_RKIND 

    ! ====================================================================
    ! Main Calculation Loop
    ! ====================================================================
    
    DO j = jts, jte
       DO i = its, ite
          
          ! Ensure land use index is safely bounded within 1 to 21
          lu = MAX(1, MIN(21, ivgtyp(i,j)))
          
          ! Retrieve base resistances from the MODIS lookup tables
          rs_min_base = rs_min_tbl(lu)
          rcut_base   = rcut_tbl(lu)
          rgrnd_base  = rgrnd_tbl(lu)

          ! Estimate PAR (~45% of incoming shortwave radiation)
          par = MAX(0.0_RKIND, swdown(i,j) * 0.45_RKIND)

          ! --------------------------------------------------------------
          ! 1. Aerodynamic Resistance (Ra)
          ! --------------------------------------------------------------
          IF (ustar(i,j) > 1.0e-4_RKIND) THEN
             ra_loc = wspd(i,j) / (ustar(i,j)**2)
          ELSE
             ra_loc = 9999.0_RKIND 
          END IF

          ! Loop over chemical species
          DO nv = 1, num_chem

             ! -----------------------------------------------------------
             ! 2. Quasi-Laminar Boundary Layer Resistance (Rb)
             ! -----------------------------------------------------------
             IF (ustar(i,j) > 1.0e-4_RKIND) THEN
                rb = (5.0_RKIND / ustar(i,j)) * (diff_ratio(nv)**0.6667_RKIND)
             ELSE
                rb = 9999.0_RKIND 
             END IF

             ! -----------------------------------------------------------
             ! 3. Surface / Canopy Resistance (Rc)
             ! -----------------------------------------------------------
             
             ! A) Stomatal resistance 
             IF (rs_min_base < 9000.0_RKIND .AND. swdown(i,j) > swdown_thresh) THEN
                r_stomatal = rs_min_base * diff_ratio(nv) * &
                             (1.0_RKIND + (200.0_RKIND / (par + 0.1_RKIND)))
             ELSE
                r_stomatal = 9999.0_RKIND
             END IF
             
             ! B) Mesophyll resistance
             r_mesophyll = 1.0_RKIND / ( (henry_const(nv) / 3000.0_RKIND) + &
                                         (100.0_RKIND * reactivity(nv)) + 1.0e-9_RKIND )
             
             ! C) Cuticular & Ground resistances
             r_cuticle = rcut_base / ( (1.0e-5_RKIND * henry_const(nv)) + &
                                       reactivity(nv) + 1.0e-9_RKIND )

             r_ground  = rgrnd_base / ( (1.0e-5_RKIND * henry_const(nv)) + &
                                       reactivity(nv) + 1.0e-9_RKIND )

             ! Combine parallel pathways for total Rc
             rc = 1.0_RKIND / ( (1.0_RKIND / (r_stomatal + r_mesophyll)) + &
                                (1.0_RKIND / r_cuticle) + &
                                (1.0_RKIND / r_ground) )

             ! -----------------------------------------------------------
             ! 4. Total Deposition Velocity Calculation
             ! -----------------------------------------------------------
             r_total = ra_loc + rb + rc

             IF (r_total > r_min) THEN
                ddvel(i, j, nv) = 1.0_RKIND / r_total
             ELSE
                ddvel(i, j, nv) = 0.0_RKIND
             END IF

          END DO 
       END DO 
    END DO 

  END SUBROUTINE calc_gas_dep_vel

END MODULE dep_dry_mod_wesley
