module module_hab_emissions

! Main references:
! The GOCART/WRF-Chem Implementation: Gong (2003). A parameterization of sea-salt aerosol source function for sub- and super-micron particles. (Published in Global Biogeochemical Cycles). 
! Monahan et al. (1986). A model of marine aerosol generation via whitecaps and wave disruption.
! Langley (2019), The Foundational Bioaerosol $\Delta T$ Model: Aerosolization of Cyanobacteria: A Mesocosm Study
! The Lake Spray Aerosol (LSA) Adaptation: May et al. (2022). Incorporation of freshwater lake spray aerosol emissions into the WRF-Chem model.
! Microcystin/Toxin Specifics: Plaas et al. (2023). (Specifically stemming from Haley Plaas's extensive doctoral work and subsequent papers, e.g., Understanding the Impacts of Harmful Cyanobacterial Blooms on Air Quality)
! Recent Field Corroboration: Carter (2022). Investigation of Multiple Cyanotoxins in Toxic Lake Aerosols.
!TODO
! ADD oxidation via Zorbas et al., 2023: https://pubs.acs.org/doi/10.1021/acsearthspacechem.3c00050#:~:text=In%20addition%2C%20aerosol%20parameters%20like,force%20for%20conducting%20this%20study.


  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  private

  public :: hab_bacteria_driver

    ! --- Tunable Physical Parameters ---
    REAL(RKIND), parameter :: RHO_WATER     = 1000.0  ! Density of source water (kg/m3). Use 1024.0 for marine.
    REAL(RKIND), parameter :: U10_THRESH    = 3.5     ! Minimum wind speed for whitecapping (m/s)
    REAL(RKIND), parameter :: EF_BASE       = 150.0   ! Maximum Enrichment Factor at calm conditions
    REAL(RKIND), parameter :: EF_DECAY      = 0.15    ! Decay rate of EF with increasing wind speed
    REAL(RKIND), parameter :: FLUX_SCALAR   = 1.0E-8  ! Empirical scaling factor for bulk water flux to match GOCART totals
    REAL(RKIND), parameter :: DT_SCALAR = 0.05_RKIND  ! Fractional flux enhancement per 1K of positive Delta T

contains

    !-----------------------------------------------------------------------
    ! SUBROUTINE: hab_bacteria_driver
    ! 
    ! DESCRIPTION:
    ! Wraps the WRF/MPAS-style 2D tile loops. Calculates the 10m wind speed,
    ! derives the aerosolized water flux, calculates the enrichment factor,
    ! and applies the resulting bacteria tendency to the lowest model level.
    !-----------------------------------------------------------------------
    subroutine hab_bacteria_driver(dt, rho_phy, dz8w, u10, v10, xland, &
                                   xice, tsk, t2m,                 &
                                   bact_water_conc,                &
                                   chem,num_chem,                  &
                                   index_bact_fine,                &
                                   ids, ide, jds, jde, kds, kde,   &
                                   ims, ime, jms, jme, kms, kme,   &
                                   its, ite, jts, jte, kts, kte)

        ! --- Standard WRF/MPAS Dimensions ---
        integer, intent(in) :: ids,ide, jds,jde, kds,kde
        integer, intent(in) :: ims,ime, jms,jme, kms,kme
        integer, intent(in) :: its,ite, jts,jte, kts,kte
        integer, intent(in) :: num_chem,index_bact_fine

        ! --- Meteorological Inputs ---
        REAL(RKIND), intent(in) :: dt                              ! Time step (s)
        REAL(RKIND), intent(in), dimension(ims:ime, kms:kme, jms:jme) :: rho_phy  ! air density (m3/kg)
        REAL(RKIND), intent(in), dimension(ims:ime, kms:kme, jms:jme) :: dz8w     ! Layer thickness (m)
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: u10      ! 10m U-wind (m/s)
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: v10      ! 10m V-wind (m/s)
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: xland    ! Land mask (1=land, 2=water)
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: xice     ! Ice mask
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: tsk      ! Skin temperature (K)
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: t2m      ! 2-meter air temperature (K) <-- CHANGED

        ! --- HAB Specific Arrays ---
        ! Input: Surface concentration map (e.g., cells/m3 or ug/m3)
        REAL(RKIND), intent(in), dimension(ims:ime, jms:jme)          :: bact_water_conc 
        
        ! In/Out: chem (tracer_units / kg_air / s)
        REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT)     :: chem

        ! --- Local Variables ---
        integer :: i, j, k
        REAL(RKIND) :: wspd10       ! Wind speed magnitude
        REAL(RKIND) :: ef           ! Enrichment factor
        REAL(RKIND) :: w_flux       ! Volumetric water spray flux (m3_water m-2 s-1)
        REAL(RKIND) :: b_flux       ! Bacteria surface flux (tracer_units m-2 s-1)
        REAL(RKIND) :: delta_t            ! Difference between water and air temp
        REAL(RKIND) :: convective_scalar  ! Multiplier for the final flux

        ! --- Execution ---
        k = kts ! Emissions only applied to the lowest model layer

        ! Loop over tiles (Outer loop J, Inner loop I for Fortran memory contiguity)
        do j = jts, jte
            do i = its, ite
                
                ! Check if grid cell is water (xland typically > 1.5 for water) 
                ! and if there is a non-zero bacteria concentration.
                if (xland(i,j) > 1.5 .and. xice(i,j) < 0.1_RKIND .and. bact_water_conc(i,j) > 0.0) then

                    ! 1. Calculate wind speed magnitude
                    wspd10 = sqrt(u10(i,j)**2 + v10(i,j)**2)
                    ! Enforce upper limit
                    if (wspd10 > 35.0_RKIND) wspd10 = 35.0_RKIND

                    ! 2. Calculate Enrichment Factor for this grid cell
                    call calc_enrichment_factor(wspd10, ef)

                    ! 3. Calculate bulk volumetric water flux (m3/m2/s)
                    call calc_bulk_spray_flux(wspd10, w_flux)

                    ! 4. Calculate Convective Enhancement Scalar
                    ! Using lowest model level (k=kts) for air temperature
                    delta_t = tsk(i,j) - t2m(i,j)
                    if (delta_t > 0.0_RKIND) then
                        convective_scalar = 1.0_RKIND + (DT_SCALAR * delta_t)
                    else
                        ! Stable or neutral conditions: no enhancement
                        convective_scalar = 1.0_RKIND
                    end if

                    ! 5. Calculate total emitted bioaerosol flux
                    ! F = (Water_Flux) * (Bulk_Concentration) * EF
                    b_flux = w_flux * bact_water_conc(i,j) * ef * convective_scalar       ! * 1000 if map is in ug/L==mg/m3

                    ! Tendency = Flux / (Air_Density * Layer_Depth)
                    ! Add to existing tendency
                    ! bact_tend(i,k,j) = bact_tend(i,k,j) + (b_flux / (air_rho * dz8w(i,k,j)))
                    chem(i,k,j,index_bact_fine) = chem(i,k,j,index_bact_fine) + ((b_flux * dt) / (rho_phy(i,k,j) * dz8w(i,k,j)))

                end if

            end do
        end do

    end subroutine hab_bacteria_driver


    !-----------------------------------------------------------------------
    ! SUBROUTINE: calc_enrichment_factor
    ! 
    ! DESCRIPTION:
    ! Calculates the enrichment factor based on wind speed. Models the 
    ! dilution of the sea surface microlayer (SML) due to turbulent 
    ! mixing at higher wind speeds.
    !-----------------------------------------------------------------------
    subroutine calc_enrichment_factor(wspd, ef_out)
        REAL(RKIND), intent(in)  :: wspd
        REAL(RKIND), intent(out) :: ef_out
        REAL(RKIND) :: excess_wind
        
        if (wspd > U10_THRESH) then
            excess_wind = wspd - U10_THRESH
            ! Exponential decay of enrichment due to wave mixing
            ef_out = EF_BASE * exp(-1.0 * EF_DECAY * excess_wind)
            
            ! Ensure EF doesn't drop below 1.0 (bulk water concentration)
            if (ef_out < 1.0) ef_out = 1.0
        else
            ! Below whitecap threshold, microlayer is fully intact/enriched
            ef_out = EF_BASE
        end if

    end subroutine calc_enrichment_factor


    !-----------------------------------------------------------------------
    ! SUBROUTINE: calc_bulk_spray_flux
    ! 
    ! DESCRIPTION:
    ! Computes the total volumetric flux of water emitted as aerosol.
    ! Uses a simplified Monahan whitecap dependence (U10^3.41).
    !-----------------------------------------------------------------------
    subroutine calc_bulk_spray_flux(wspd, flux_out)
        REAL(RKIND), intent(in)  :: wspd
        REAL(RKIND), intent(out) :: flux_out
        
        if (wspd > U10_THRESH) then
            ! Monahan whitecap coverage scales with U10^3.41
            ! FLUX_SCALAR is used to convert this relationship into a bulk volume
            flux_out = FLUX_SCALAR * (wspd ** 3.41)
        else
            flux_out = 0.0
        end if

    end subroutine calc_bulk_spray_flux

end module module_hab_emissions
