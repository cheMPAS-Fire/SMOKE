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
  REAL(KIND=RKIND), DIMENSION(max_modis) :: rs_min_tbl ! Min Stomatal Resistance (s/m)
  REAL(KIND=RKIND), DIMENSION(max_modis) :: rcut_tbl   ! Base Cuticular Resistance (s/m)
  REAL(KIND=RKIND), DIMENSION(max_modis) :: rgrnd_tbl  ! Base Ground/Soil Resistance (s/m)

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
    ! - Values represent Midsummer/Active Growing Season (Category 1 in 
    !   the original Wesely formulation).
    ! ====================================================================
    
    ! 1. Evergreen Needleleaf
    rs_min_tbl(1) =  130.0_RKIND
    rcut_tbl(1)   = 2000.0_RKIND
    rgrnd_tbl(1)  =  200.0_RKIND
    
    ! 2. Evergreen Broadleaf
    rs_min_tbl(2) =  130.0_RKIND
    rcut_tbl(2)   = 2000.0_RKIND
    rgrnd_tbl(2)  =  200.0_RKIND
    
    ! 3. Deciduous Needleleaf
    rs_min_tbl(3) =  100.0_RKIND
    rcut_tbl(3)   = 2000.0_RKIND
    rgrnd_tbl(3)  =  200.0_RKIND
    
    ! 4. Deciduous Broadleaf
    rs_min_tbl(4) =   70.0_RKIND
    rcut_tbl(4)   = 2000.0_RKIND
    rgrnd_tbl(4)  =  200.0_RKIND
    
    ! 5. Mixed Forest
    rs_min_tbl(5) =  100.0_RKIND
    rcut_tbl(5)   = 2000.0_RKIND
    rgrnd_tbl(5)  =  200.0_RKIND
    
    ! 6. Closed Shrublands
    rs_min_tbl(6) =  120.0_RKIND
    rcut_tbl(6)   = 2000.0_RKIND
    rgrnd_tbl(6)  =  200.0_RKIND
    
    ! 7. Open Shrublands
    rs_min_tbl(7) =  120.0_RKIND
    rcut_tbl(7)   = 2000.0_RKIND
    rgrnd_tbl(7)  =  300.0_RKIND
    
    ! 8. Woody Savannas
    rs_min_tbl(8) =  120.0_RKIND
    rcut_tbl(8)   = 2000.0_RKIND
    rgrnd_tbl(8)  =  200.0_RKIND
    
    ! 9. Savannas
    rs_min_tbl(9) =  120.0_RKIND
    rcut_tbl(9)   = 2000.0_RKIND
    rgrnd_tbl(9)  =  200.0_RKIND
    
    ! 10. Grasslands
    rs_min_tbl(10) =  120.0_RKIND
    rcut_tbl(10)   = 2000.0_RKIND
    rgrnd_tbl(10)  =  200.0_RKIND
    
    ! 11. Permanent Wetlands
    rs_min_tbl(11) =   80.0_RKIND
    rcut_tbl(11)   = 2000.0_RKIND
    rgrnd_tbl(11)  =  100.0_RKIND
    
    ! 12. Croplands
    rs_min_tbl(12) =   70.0_RKIND
    rcut_tbl(12)   = 2000.0_RKIND
    rgrnd_tbl(12)  =  150.0_RKIND
    
    ! 13. Urban and Built-Up
    rs_min_tbl(13) = 9999.0_RKIND
    rcut_tbl(13)   = 9999.0_RKIND
    rgrnd_tbl(13)  =  400.0_RKIND
    
    ! 14. Cropland/Natural Veg Mosaic
    rs_min_tbl(14) =   70.0_RKIND
    rcut_tbl(14)   = 2000.0_RKIND
    rgrnd_tbl(14)  =  150.0_RKIND
    
    ! 15. Snow and Ice
    rs_min_tbl(15) = 9999.0_RKIND
    rcut_tbl(15)   = 9999.0_RKIND
    rgrnd_tbl(15)  = 1000.0_RKIND
    
    ! 16. Barren or Sparsely Vegetated
    rs_min_tbl(16) = 9999.0_RKIND
    rcut_tbl(16)   = 9999.0_RKIND
    rgrnd_tbl(16)  =  500.0_RKIND
    
    ! 17. Water (Oceans, Lakes)
    rs_min_tbl(17) = 9999.0_RKIND
    rcut_tbl(17)   = 9999.0_RKIND
    rgrnd_tbl(17)  =   10.0_RKIND
    
    ! 18. Wooded Tundra
    rs_min_tbl(18) =  150.0_RKIND
    rcut_tbl(18)   = 2000.0_RKIND
    rgrnd_tbl(18)  =  300.0_RKIND
    
    ! 19. Mixed Tundra
    rs_min_tbl(19) =  150.0_RKIND
    rcut_tbl(19)   = 2000.0_RKIND
    rgrnd_tbl(19)  =  300.0_RKIND
    
    ! 20. Bare Ground Tundra
    rs_min_tbl(20) = 9999.0_RKIND
    rcut_tbl(20)   = 9999.0_RKIND
    rgrnd_tbl(20)  =  400.0_RKIND
    
    ! 21. Unclassified/Missing
    rs_min_tbl(21) = 9999.0_RKIND
    rcut_tbl(21)   = 9999.0_RKIND
    rgrnd_tbl(21)  =  400.0_RKIND

   end subroutine gas_dry_dep_init

end module dep_data_mod
