!>\file  dep_data_mod.F90
!! This file contains data for the dry deposition modules.
module dep_data_mod

  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  public :: aero_dry_dep_init, aero_wet_dep_init, gas_dry_dep_init

  real(RKIND), parameter :: max_dep_vel = 0.005                   ! m/s (may need to set per species)
  real(RKIND), parameter :: dep_ref_hgt = 2.0                     ! Meters 
  real(RKIND), parameter :: pi = 3.1415926536
!  3*PI
  REAL(RKIND), PARAMETER :: threepi=3.0*pi
  real(RKIND), parameter :: gravity =  9.81
! mean gravitational acceleration [ m/sec**2 ]
  REAL(RKIND), PARAMETER :: grav=9.80622
  real(RKIND), parameter :: boltzmann = 1.3807e-16
! universal gas constant [ J/mol-K ]
  REAL(RKIND), PARAMETER :: rgasuniv=8.314510
! Avogadro's Constant [ 1/mol ]
  REAL, PARAMETER :: avo=6.0221367E23
  ! Boltzmann's Constant [ J / K ]i\
  REAL(RKIND), PARAMETER :: boltz=rgasuniv/avo
  real(RKIND), parameter :: Cb = 2., Cim = 0.4, alpha = 0.8, Cin = 2.5, vv = 0.8
  real(RKIND), parameter :: A_for = 0.1 ! forest
  real(RKIND), parameter :: A_grs = 0.2 ! grass
  real(RKIND), parameter :: A_wat = 100. ! water
  real(RKIND), parameter :: eps0_for = 0.8*0.01 ! forest
  real(RKIND), parameter :: eps0_grs = 0.4*0.01 ! grass
  real(RKIND), parameter :: eps0_wat = 0.6*0.01 ! water
  REAL(RKIND), PARAMETER :: one3=1.0/3.0
  REAL(RKIND), PARAMETER :: two3=2.0/3.0
!  SQRT( 2 )
  REAL(RKIND), PARAMETER :: sqrt2=1.4142135623731
!  SQRT( PI )
  REAL(RKIND), PARAMETER :: sqrtpi=1.7724539
  REAL(RKIND) :: karman = 0.4                             ! von Karman constant
  REAL(RKIND), PARAMETER :: conmin= 1.E-16
  REAL(RKIND), PARAMETER :: pirs=3.14159265358979324
  REAL(RKIND), PARAMETER :: f6dpi=6.0/pirs
  REAL(RKIND), PARAMETER :: f6dpim9=1.0E-9*f6dpi
!  starting standard surface temperature [ K ]
  REAL(RKIND), PARAMETER :: sigma1 = 1.8
  REAL(RKIND), PARAMETER :: mean_diameter1 = 4.e-8
  REAL(RKIND), PARAMETER :: fact_wfa = 1.e-9*6.0/pirs*exp(4.5*log(sigma1)**2)/mean_diameter1**3
! lowest particle diameter ( m )   
  REAL(RKIND), PARAMETER :: dgmin=1.0E-09
! lowest particle density ( Kg/m**3 )
  REAL(RKIND), PARAMETER :: densmin=1.0E03

! Arrays to hold density and diameters of aerosols
  real(RKIND), dimension(1000), save :: aero_dens, aero_diam, ls_frac
  REAL(KIND=RKIND), DIMENSION(1000), save :: henry_const  ! Effective Henry's Law Constant (M/atm)
  REAL(KIND=RKIND), DIMENSION(1000), save :: reactivity   ! Reactivity factor for Wesely Rc (f0)
  REAL(KIND=RKIND), DIMENSION(1000), save :: diff_ratio   ! Ratio of molecular diffusivity to H2O
! ====================================================================
  ! MODIS Land-Use Lookup Tables (Midsummer/Active Vegetation mapping)
  ! ====================================================================
  INTEGER, PARAMETER :: max_modis = 21
  INTEGER, PARAMETER :: max_season = 5
  INTEGER            :: s
  REAL(KIND=RKIND), DIMENSION(max_modis,max_season) :: rs_min_tbl ! Min Stomatal Resistance (s/m)
  REAL(KIND=RKIND), DIMENSION(max_modis,max_season) :: rcut_tbl   ! Base Cuticular Resistance (s/m)
  REAL(KIND=RKIND), DIMENSION(max_modis,max_season) :: rgrnd_tbl  ! Base Ground/Soil Resistance (s/m)

  contains
   subroutine aero_wet_dep_init

      implicit none

      ls_frac(:) = -999._RKIND

      if (p_smoke_ultrafine>0)   ls_frac(p_smoke_ultrafine)  = 0.5_RKIND 
      if (p_smoke_fine>0)        ls_frac(p_smoke_fine)       = 0.5_RKIND
      if (p_smoke_coarse>0)      ls_frac(p_smoke_coarse)     = 0.5_RKIND
!      
      if (p_dust_ultrafine>0)    ls_frac(p_dust_ultrafine)   = 0.5_RKIND
      if (p_dust_fine>0)         ls_frac(p_dust_fine)        = 0.5_RKIND
      if (p_dust_coarse>0)       ls_frac(p_dust_coarse)      = 0.5_RKIND
!      
      if (p_unspc_ultrafine>0)   ls_frac(p_unspc_ultrafine)  = 0.5_RKIND
      if (p_unspc_fine>0)        ls_frac(p_unspc_fine)       = 0.5_RKIND
      if (p_unspc_coarse>0)      ls_frac(p_unspc_coarse)     = 0.5_RKIND
!      
      if (p_ssalt_fine>0)        ls_frac(p_ssalt_fine)       = 0.5_RKIND
      if (p_ssalt_coarse>0)      ls_frac(p_ssalt_coarse)     = 0.5_RKIND
!
      if (p_bact_fine>0)         ls_frac(p_bact_fine)        = 0.5_RKIND
      if (p_polp_all>0)          ls_frac(p_polp_all)         = 0.5_RKIND
      if (p_polp_tree>0)         ls_frac(p_polp_tree)        = 0.5_RKIND
      if (p_polp_grass>0)        ls_frac(p_polp_grass)       = 0.5_RKIND
      if (p_polp_weed>0)         ls_frac(p_polp_weed)        = 0.5_RKIND
!
      if (p_pols_all>0)          ls_frac(p_pols_all)         = 0.5_RKIND
      if (p_pols_tree>0)         ls_frac(p_pols_tree)        = 0.5_RKIND
      if (p_pols_grass>0)        ls_frac(p_pols_grass)       = 0.5_RKIND
      if (p_pols_weed>0)         ls_frac(p_pols_weed)        = 0.5_RKIND

   end subroutine aero_wet_dep_init


   subroutine aero_dry_dep_init

      implicit none

      aero_dens(:) = -999._RKIND
      aero_diam(:) = -999._RKIND

    ! Aerosol densities (kg/m3)
      if (p_smoke_ultrafine>0)   aero_dens(p_smoke_ultrafine)  = 1.4E3_RKIND
      if (p_smoke_fine>0)        aero_dens(p_smoke_fine)       = 1.4E3_RKIND
      if (p_smoke_coarse>0)      aero_dens(p_smoke_coarse)     = 1.4E3_RKIND
!      
      if (p_dust_ultrafine>0)    aero_dens(p_dust_ultrafine)   = 2.6E3_RKIND
      if (p_dust_fine>0)         aero_dens(p_dust_fine)        = 2.6E3_RKIND 
      if (p_dust_coarse>0)       aero_dens(p_dust_coarse)      = 2.6E3_RKIND
!      
      if (p_unspc_ultrafine>0)   aero_dens(p_unspc_ultrafine)  = 2.6E3_RKIND
      if (p_unspc_fine>0)        aero_dens(p_unspc_fine)       = 2.6E3_RKIND
      if (p_unspc_coarse>0)      aero_dens(p_unspc_coarse)     = 2.6E3_RKIND
!      
      if (p_ssalt_fine>0)        aero_dens(p_ssalt_fine)       = 2.2E3_RKIND
      if (p_ssalt_coarse>0)      aero_dens(p_ssalt_coarse)     = 2.2E3_RKIND
!
      if (p_bact_fine>0)         aero_dens(p_bact_fine)        = 1.425E3_RKIND
      if (p_polp_all>0)          aero_dens(p_polp_all)         = 1.2E3_RKIND
      if (p_polp_tree>0)         aero_dens(p_polp_tree)        = 1.2E3_RKIND
      if (p_polp_grass>0)        aero_dens(p_polp_grass)       = 1.0E3_RKIND
      if (p_polp_weed>0)         aero_dens(p_polp_weed)        = 1.3E3_RKIND
!
      if (p_pols_all>0)          aero_dens(p_pols_all)         = 1.425E3_RKIND
      if (p_pols_tree>0)         aero_dens(p_pols_tree)        = 1.425E3_RKIND
      if (p_pols_grass>0)        aero_dens(p_pols_grass)       = 1.425E3_RKIND
      if (p_pols_weed>0)         aero_dens(p_pols_weed)        = 1.425E3_RKIND
    ! Aerosol diameters (m)
      if (p_smoke_ultrafine>0)   aero_diam(p_smoke_ultrafine)   = 4E-9_RKIND    ! JLS, check
      if (p_smoke_fine>0)        aero_diam(p_smoke_fine)        = 4E-8_RKIND    ! JLS, check
      if (p_smoke_coarse>0)      aero_diam(p_smoke_coarse)      = 10E-6_RKIND
!
      if (p_dust_ultrafine>0)    aero_diam(p_dust_ultrafine)    = 1E-9_RKIND    ! JLS, check
      if (p_dust_fine>0)         aero_diam(p_dust_fine)         = 1E-6_RKIND
      if (p_dust_coarse>0)       aero_diam(p_dust_coarse)       = 4.5E-6_RKIND

      if (p_unspc_ultrafine>0)   aero_diam(p_unspc_ultrafine)   = 1E-9_RKIND    ! JLS, check
      if (p_unspc_fine>0)        aero_diam(p_unspc_fine)        = 1E-7_RKIND
      if (p_unspc_coarse>0)      aero_diam(p_unspc_coarse)      = 4.5E-6_RKIND
!      
      if (p_ssalt_fine>0)        aero_diam(p_ssalt_fine)        = 6.32E-7_RKIND
      if (p_ssalt_coarse>0)      aero_diam(p_ssalt_coarse)      = 5.632E-6_RKIND
!
      if (p_bact_fine>0)         aero_diam(p_bact_fine)         = 5.E-6_RKIND

!
      if (p_polp_all>0)          aero_diam(p_polp_all)          = 30E-6_RKIND
      if (p_polp_tree>0)         aero_diam(p_polp_tree)         = 35E-6_RKIND
      if (p_polp_grass>0)        aero_diam(p_polp_grass)        = 35E-6_RKIND
      if (p_polp_weed>0)         aero_diam(p_polp_weed)         = 20E-6_RKIND
!
      if (p_pols_all>0)          aero_diam(p_pols_all)          = 1.5E-7_RKIND
      if (p_pols_tree>0)         aero_diam(p_pols_tree)         = 1.5E-7_RKIND
      if (p_pols_grass>0)        aero_diam(p_pols_grass)        = 1.5E-7_RKIND
      if (p_pols_weed>0)         aero_diam(p_pols_weed)         = 1.5E-7_RKIND

   end subroutine aero_dry_dep_init

   subroutine gas_dry_dep_init

      implicit none

    

      henry_const(:) = -999._RKIND  ! Effective Henry's Law Constant (M/atm)
      reactivity(:) = -999._RKIND   ! Reactivity factor for Wesely Rc (f0)
      diff_ratio(:) = -999._RKIND   ! Ratio of molecular diffusivity to H2O

! --- Inorganics ---
    
    ! SO2
    ! Ref: Wesely (1989) Table 1 & 2; Sander (2015) for H*
    henry_const(p_so2) = 1.0e5_RKIND
    reactivity(p_so2)  = 1.0_RKIND
    diff_ratio(p_so2)  = 0.53_RKIND

    ! NO
    ! Ref: Wesely (1989) Table 2; Sander (2015)
    !henry_const(p_no)  = 2.0e-3_RKIND
    !reactivity(p_no)   = 0.1_RKIND
    !diff_ratio(p_no)   = 0.80_RKIND

    ! NO2
    ! Ref: Wesely (1989) Table 2; Sander (2015)
    !henry_const(p_no2) = 1.0e-2_RKIND
    !reactivity(p_no2)  = 0.1_RKIND
    !diff_ratio(p_no2)  = 0.60_RKIND

    ! NOX Tracer
    ! Ref: Treating identically to NO2 for bulk modeling purposes
    henry_const(p_nox) = 1.0e-2_RKIND
    reactivity(p_nox)  = 0.1_RKIND
    diff_ratio(p_nox)  = 0.60_RKIND

    ! HNO3
    ! Ref: Wesely (1989) Table 2; Sander (2015) - effectively infinite H*
    !henry_const(p_hno3)= 1.0e14_RKIND
    !reactivity(p_hno3) = 1.0_RKIND
    !diff_ratio(p_hno3) = 0.53_RKIND

    ! H2O2
    ! Ref: Wesely (1989) Table 2; Sander (2015)
    !henry_const(p_h2o2)= 1.0e5_RKIND
    !reactivity(p_h2o2) = 1.0_RKIND
    !diff_ratio(p_h2o2) = 0.70_RKIND

    ! CO
    ! Ref: Massman (1998) for diffusivity; Sander (2015) for low H*
    henry_const(p_co)  = 1.0e-3_RKIND
    reactivity(p_co)   = 0.0_RKIND
    diff_ratio(p_co)   = 0.80_RKIND

    ! NH3
    ! Ref: Zhang et al. (2003) for NH3 specific parameterizations; Sander (2015)
    henry_const(p_nh3) = 7.4e1_RKIND
    reactivity(p_nh3)  = 1.0_RKIND
    diff_ratio(p_nh3)  = 0.90_RKIND

! --- Custom Bulk Proxies ---
    
    ! Bulk SOA (Bulk, Anthropogenic, Biomass Burning)
    ! Ref: Ahmadov et al. (2012) WRF-Chem SOA treatments. 
    ! Treated as highly sticky/condensable, slow diffusing due to MW.
! TODO - once SOA code is merged

!    henry_const(p_soa)    = 1.0e5_RKIND
!    reactivity(p_soa)     = 0.1_RKIND
!    diff_ratio(p_soa)     = 0.20_RKIND
!
!    henry_const(p_antsoa) = 1.0e5_RKIND
!    reactivity(p_antsoa)  = 0.1_RKIND
!    diff_ratio(p_antsoa)  = 0.20_RKIND
!
!    henry_const(p_bbsoa)  = 1.0e5_RKIND
!    reactivity(p_bbsoa)   = 0.1_RKIND
!    diff_ratio(p_bbsoa)   = 0.20_RKIND

! Bulk VOCs (Bulk, Anthropogenic, Biomass Burning)
    ! Ref: Generic CTM approximations for unspeciated non-methane hydrocarbons (NMHCs).
    ! Representative of a mix of oxygenated and non-oxygenated heavier species.
    !henry_const(p_voc)    = 1.0_RKIND
    !reactivity(p_voc)     = 0.5_RKIND
    !diff_ratio(p_voc)     = 0.40_RKIND

! TODO - once SOA code is merged
!    henry_const(p_antvoc) = 1.0_RKIND
!    reactivity(p_antvoc)  = 0.5_RKIND
!    diff_ratio(p_antvoc)  = 0.40_RKIND
!
!    henry_const(p_bbvoc)  = 1.0_RKIND
!    reactivity(p_bbvoc)   = 0.5_RKIND
!    diff_ratio(p_bbvoc)   = 0.40_RKIND

! ====================================================================
    ! 3. Populate MODIS Land-Use Parameters
    ! 
    ! PRIMARY REFERENCES:
    ! - Base mappings derived from Wesely (1989) Table 1 (Land-use categories 
    !   and seasonal categories) mapped to MODIS 21-category definitions 
    !   commonly used in standard WRF-Chem implementations.

! ====================================================================
    ! 2. Populate Seasonal MODIS Land-Use Parameters
    ! PRIMARY REFERENCES: Wesely (1989) Table 1.
    ! Seasons: 1=Summer, 2=Autumn, 3=Late Autumn, 4=Winter, 5=Spring
    ! ====================================================================
    
    ! --------------------------------------------------------------------
    ! SEASON 1: MIDSUMMER (Active, Lush Vegetation)
    ! --------------------------------------------------------------------
    ! Evergreen Needleleaf / Broadleaf
    rs_min_tbl(1:2, 1) =  130.0_RKIND
    rcut_tbl(1:2, 1)   = 2000.0_RKIND
    rgrnd_tbl(1:2, 1)  =  200.0_RKIND
    
    ! Deciduous Needleleaf / Broadleaf / Mixed Forest
    rs_min_tbl(3:5, 1) =   70.0_RKIND
    rcut_tbl(3:5, 1)   = 2000.0_RKIND
    rgrnd_tbl(3:5, 1)  =  200.0_RKIND
    
    ! Shrublands, Savannas, Grasslands
    rs_min_tbl(6:10, 1) =  120.0_RKIND
    rcut_tbl(6:10, 1)   = 2000.0_RKIND
    rgrnd_tbl(6:10, 1)  =  200.0_RKIND
    
    ! Wetlands
    rs_min_tbl(11, 1) =   80.0_RKIND
    rcut_tbl(11, 1)   = 2000.0_RKIND
    rgrnd_tbl(11, 1)  =  100.0_RKIND
    
    ! Croplands & Mosaics
    rs_min_tbl(12, 1) =   70.0_RKIND
    rs_min_tbl(14, 1) =   70.0_RKIND
    rcut_tbl(12, 1)   = 2000.0_RKIND
    rcut_tbl(14, 1)   = 2000.0_RKIND
    rgrnd_tbl(12, 1)  =  150.0_RKIND
    rgrnd_tbl(14, 1)  =  150.0_RKIND
    
    ! Tundra
    rs_min_tbl(18:19, 1) =  150.0_RKIND
    rcut_tbl(18:19, 1)   = 2000.0_RKIND
    rgrnd_tbl(18:19, 1)  =  300.0_RKIND


    ! --------------------------------------------------------------------
    ! SEASON 2: AUTUMN (Unharvested, Senescence)
    ! --------------------------------------------------------------------
    ! Evergreen types maintain stomata
    rs_min_tbl(1:2, 2) =  130.0_RKIND
    rcut_tbl(1:2, 2)   = 2000.0_RKIND
    rgrnd_tbl(1:2, 2)  =  200.0_RKIND
    
    ! Deciduous & Mixed Forest (Leaves dying, resistance up)
    rs_min_tbl(3:5, 2) =  120.0_RKIND
    rcut_tbl(3:5, 2)   = 3000.0_RKIND
    rgrnd_tbl(3:5, 2)  =  200.0_RKIND
    
    ! Shrublands, Savannas, Grasslands
    rs_min_tbl(6:10, 2) =  150.0_RKIND
    rcut_tbl(6:10, 2)   = 3000.0_RKIND
    rgrnd_tbl(6:10, 2)  =  300.0_RKIND
    
    ! Wetlands
    rs_min_tbl(11, 2) =  120.0_RKIND
    rcut_tbl(11, 2)   = 3000.0_RKIND
    rgrnd_tbl(11, 2)  =  100.0_RKIND
    
    ! Croplands & Mosaics
    rs_min_tbl(12, 2) =  120.0_RKIND
    rs_min_tbl(14, 2) =  120.0_RKIND
    rcut_tbl(12, 2)   = 3000.0_RKIND
    rcut_tbl(14, 2)   = 3000.0_RKIND
    rgrnd_tbl(12, 2)  =  200.0_RKIND
    rgrnd_tbl(14, 2)  =  200.0_RKIND
    
    ! Tundra
    rs_min_tbl(18:19, 2) =  200.0_RKIND
    rcut_tbl(18:19, 2)   = 3000.0_RKIND
    rgrnd_tbl(18:19, 2)  =  300.0_RKIND


    ! --------------------------------------------------------------------
    ! SEASON 3: LATE AUTUMN (Post-frost, dormant, no snow)
    ! --------------------------------------------------------------------
    ! Evergreen
    rs_min_tbl(1:2, 3) =  130.0_RKIND
    rcut_tbl(1:2, 3)   = 2000.0_RKIND
    rgrnd_tbl(1:2, 3)  =  200.0_RKIND
    
    ! Deciduous (Leaves gone, stomatal resistance infinite)
    rs_min_tbl(3:5, 3) = 9999.0_RKIND
    rcut_tbl(3:5, 3)   = 4000.0_RKIND
    rgrnd_tbl(3:5, 3)  =  200.0_RKIND
    
    ! Grass/Shrubs (Dormant)
    rs_min_tbl(6:10, 3) = 9999.0_RKIND
    rcut_tbl(6:10, 3)   = 4000.0_RKIND
    rgrnd_tbl(6:10, 3)  =  300.0_RKIND
    
    ! Wetlands
    rs_min_tbl(11, 3) = 9999.0_RKIND
    rcut_tbl(11, 3)   = 4000.0_RKIND
    rgrnd_tbl(11, 3)  =  100.0_RKIND
    
    ! Croplands & Mosaics (Harvested)
    rs_min_tbl(12, 3) = 9999.0_RKIND
    rs_min_tbl(14, 3) = 9999.0_RKIND
    rcut_tbl(12, 3)   = 4000.0_RKIND
    rcut_tbl(14, 3)   = 4000.0_RKIND
    rgrnd_tbl(12, 3)  =  300.0_RKIND
    rgrnd_tbl(14, 3)  =  300.0_RKIND
    
    ! Tundra
    rs_min_tbl(18:19, 3) = 9999.0_RKIND
    rcut_tbl(18:19, 3)   = 4000.0_RKIND
    rgrnd_tbl(18:19, 3)  =  300.0_RKIND


    ! --------------------------------------------------------------------
    ! SEASON 4: WINTER (Snow on ground, subfreezing)
    ! --------------------------------------------------------------------
    ! Evergreen (Stomata severely restricted by cold, ground covered in snow)
    rs_min_tbl(1:2, 4) =  250.0_RKIND
    rcut_tbl(1:2, 4)   = 4000.0_RKIND
    rgrnd_tbl(1:2, 4)  = 1000.0_RKIND
    
    ! Deciduous (No leaves, snow ground)
    rs_min_tbl(3:5, 4) = 9999.0_RKIND
    rcut_tbl(3:5, 4)   = 9999.0_RKIND
    rgrnd_tbl(3:5, 4)  = 1000.0_RKIND
    
    ! Grass/Shrubs (Buried in snow)
    rs_min_tbl(6:10, 4) = 9999.0_RKIND
    rcut_tbl(6:10, 4)   = 9999.0_RKIND
    rgrnd_tbl(6:10, 4)  = 1000.0_RKIND
    
    ! Wetlands (Frozen)
    rs_min_tbl(11, 4) = 9999.0_RKIND
    rcut_tbl(11, 4)   = 9999.0_RKIND
    rgrnd_tbl(11, 4)  = 1000.0_RKIND
    
    ! Croplands & Mosaics
    rs_min_tbl(12, 4) = 9999.0_RKIND
    rs_min_tbl(14, 4) = 9999.0_RKIND
    rcut_tbl(12, 4)   = 9999.0_RKIND
    rcut_tbl(14, 4)   = 9999.0_RKIND
    rgrnd_tbl(12, 4)  = 1000.0_RKIND
    rgrnd_tbl(14, 4)  = 1000.0_RKIND
    
    ! Tundra
    rs_min_tbl(18:19, 4) = 9999.0_RKIND
    rcut_tbl(18:19, 4)   = 9999.0_RKIND
    rgrnd_tbl(18:19, 4)  = 1000.0_RKIND


    ! --------------------------------------------------------------------
    ! SEASON 5: TRANSITIONAL SPRING (Emerging vegetation)
    ! --------------------------------------------------------------------
    ! Evergreen 
    rs_min_tbl(1:2, 5) =  130.0_RKIND
    rcut_tbl(1:2, 5)   = 2000.0_RKIND
    rgrnd_tbl(1:2, 5)  =  200.0_RKIND
    
    ! Deciduous (Buds opening)
    rs_min_tbl(3:5, 5) =  120.0_RKIND
    rcut_tbl(3:5, 5)   = 2000.0_RKIND
    rgrnd_tbl(3:5, 5)  =  200.0_RKIND
    
    ! Shrublands, Savannas, Grasslands
    rs_min_tbl(6:10, 5) =  150.0_RKIND
    rcut_tbl(6:10, 5)   = 2000.0_RKIND
    rgrnd_tbl(6:10, 5)  =  200.0_RKIND
    
    ! Wetlands
    rs_min_tbl(11, 5) =  120.0_RKIND
    rcut_tbl(11, 5)   = 2000.0_RKIND
    rgrnd_tbl(11, 5)  =  100.0_RKIND
    
    ! Croplands & Mosaics
    rs_min_tbl(12, 5) =  120.0_RKIND
    rs_min_tbl(14, 5) =  120.0_RKIND
    rcut_tbl(12, 5)   = 2000.0_RKIND
    rcut_tbl(14, 5)   = 2000.0_RKIND
    rgrnd_tbl(12, 5)  =  200.0_RKIND
    rgrnd_tbl(14, 5)  =  200.0_RKIND
    
    ! Tundra
    rs_min_tbl(18:19, 5) =  200.0_RKIND
    rcut_tbl(18:19, 5)   = 2000.0_RKIND
    rgrnd_tbl(18:19, 5)  =  300.0_RKIND


    ! --------------------------------------------------------------------
    ! NON-VEGETATED TYPES (Constant across all 5 Seasons)
    ! --------------------------------------------------------------------
    DO s = 1, max_season
       ! Urban and Built-Up (13)
       rs_min_tbl(13, s) = 9999.0_RKIND
       rcut_tbl(13, s)   = 9999.0_RKIND
       
       ! Snow and Ice (15)
       rs_min_tbl(15, s) = 9999.0_RKIND
       rcut_tbl(15, s)   = 9999.0_RKIND
       rgrnd_tbl(15, s)  = 1000.0_RKIND
       
       ! Barren or Sparsely Vegetated (16)
       rs_min_tbl(16, s) = 9999.0_RKIND
       rcut_tbl(16, s)   = 9999.0_RKIND
       
       ! Water (Oceans, Lakes) (17)
       rs_min_tbl(17, s) = 9999.0_RKIND
       rcut_tbl(17, s)   = 9999.0_RKIND
       rgrnd_tbl(17, s)  =   10.0_RKIND
       
       ! Bare Ground Tundra (20)
       rs_min_tbl(20, s) = 9999.0_RKIND
       rcut_tbl(20, s)   = 9999.0_RKIND
       
       ! Unclassified/Missing (21)
       rs_min_tbl(21, s) = 9999.0_RKIND
       rcut_tbl(21, s)   = 9999.0_RKIND
    END DO

    ! Ground resistance specific adjustments for non-vegetated types in winter (S4)
    rgrnd_tbl(13, 1:3) = 400.0_RKIND
    rgrnd_tbl(13, 4)   = 1000.0_RKIND
    rgrnd_tbl(13, 5)   = 400.0_RKIND
    
    rgrnd_tbl(16, 1:3) = 500.0_RKIND
    rgrnd_tbl(16, 4)   = 1000.0_RKIND
    rgrnd_tbl(16, 5)   = 500.0_RKIND
    
    rgrnd_tbl(20, 1:3) = 400.0_RKIND
    rgrnd_tbl(20, 4)   = 1000.0_RKIND
    rgrnd_tbl(20, 5)   = 400.0_RKIND
    
    rgrnd_tbl(21, 1:3) = 400.0_RKIND
    rgrnd_tbl(21, 4)   = 1000.0_RKIND
    rgrnd_tbl(21, 5)   = 400.0_RKIND

   end subroutine gas_dry_dep_init

end module dep_data_mod
