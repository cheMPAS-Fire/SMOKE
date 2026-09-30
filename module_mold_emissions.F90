MODULE module_mold_emissions

  use mpas_smoke_init
  use mpas_kind_types
  use dep_data_mod 

!----------------------------------------------------------------------
! Module to calculate Fungal Spore (Mold) Emissions
! 
! References:
! Heald, C. L., & Spracklen, D. V. (2009). Atmospheric budget of 
! primary biological aerosol particles from fungal spores. 
! Geophysical Research Letters, 36(9).
!
! Hoose, C., et al. (2010). A global simulation of aerosols and their
! ice nucleating effects. Environmental Research Letters.
!----------------------------------------------------------------------

CONTAINS

SUBROUTINE calc_mold_emiss(dt, nsoil, num_chem, chem,             &
                           dz8w, t_phy, rho_phy, qv, lai, smois,  &
                           ims, ime, jms, jme, kms, kme,          &
                           its, ite, jts, jte, kts, kte           )

  IMPLICIT NONE

  !--------------------------------------------------------------------
  ! Input/Output variables
  !--------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte
  INTEGER, INTENT(IN) :: nsoil, num_chem
  REAL, INTENT(IN)    :: dt

  ! Meteorological and land-surface inputs
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: t_phy     ! Temperature (K)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: rho_phy     ! density (kg/m3)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: dz8w     ! layer height (m)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: qv    ! Water vapor mixing ratio (kg/kg)
  REAL, DIMENSION( ims:ime, jms:jme ),          INTENT(IN) :: lai       ! Leaf Area Index (m2/m2)
  REAL, DIMENSION( ims:ime, 1:nsoil, jms:jme ),          INTENT(IN) :: smois ! Soil moisture (m3/m3)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem)  :: chem

  ! Emission output array for the single mold tracer

  !--------------------------------------------------------------------
  ! Local variables and parameters
  !--------------------------------------------------------------------
  INTEGER :: i, j, k, ksoil
  REAL    :: t_sfc, qv_sfc, lai_val, sm_val
  REAL    :: mold_flux

  ! Empirical constants
  ! WARNING: C_SCALE must be recalibrated. Because soil moisture is a 
  ! fraction (typically 0.05 to 0.5), multiplying by it will reduce the 
  ! overall flux magnitude compared to the previous equation.
  REAL, PARAMETER :: C_SCALE = 2500.0  ! Adjusted empirical scaling factor
  REAL, PARAMETER :: T_MIN   = 273.15  ! Minimum temp for emission (K)

  !--------------------------------------------------------------------
  ! Execution
  !--------------------------------------------------------------------
  
  ! Emissions are injected into the lowest atmospheric model level
  k = kts
  ksoil = 1
  ! Loop over the spatial grid using WRF tile bounds for parallelization
  DO j = jts, jte
     DO i = its, ite
        

        ! Extract surface meteorological and soil variables
        t_sfc   = t_phy(i,k,j)
        qv_sfc  = MAX(0.0, qv(i,k,j))     ! Ensure non-negative atmospheric humidity
        lai_val = MAX(0.0, lai(i,j))          ! Ensure non-negative LAI
        sm_val  = MAX(0.0, smois(i,ksoil,j))    ! Ensure non-negative soil moisture

        ! Calculate flux if above the biological freezing threshold
        IF (t_sfc > T_MIN) THEN
           ! Linear parameterization driven by humidity, LAI, and soil moisture
           mold_flux = C_SCALE * lai_val * qv_sfc * sm_val
        ELSE
           mold_flux = 0.0
        END IF

        ! Add calculated flux to the chem array
        chem(i,k,j,p_mold_fine) = chem(i,k,j,p_mold_fine) + (mold_flux * dt) / (rho_phy(i,k,j) * dz8w(i,k,j))

     END DO
  END DO

END SUBROUTINE calc_mold_emiss

END MODULE module_mold_emissions
