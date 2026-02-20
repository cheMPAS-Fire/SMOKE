!>\file dep_dry_mod.F90
!! This file is for the dry depostion driver.
!-------------REVISION HISTORY---------------!
! XX/XX/XXXX : original implementation (Ravan Ahmadov)
! 08/17/2023 : modified to follow Emerson et al., (2020) (Jordan Schnell)
! 08/17/2023 : gravitational settling folowing the coarse pm settling driver (Jordan Schnell)

module dep_dry_mod_emerson

  use mpas_kind_types
  use dep_data_mod     
  use mpas_smoke_config, only : n_dbg_lines
  use mpas_smoke_init
  use mpas_timer, only : mpas_timer_start, mpas_timer_stop

  implicit none

  private

  public :: dry_dep_driver_emerson, particle_settling_wrapper, &
            particle_settling_aersett

contains
    subroutine dry_dep_driver_emerson(rmol,ustar,znt,num_chem,ddvel,         &
               chem,delz,snowh,t_phy,p_phy,rho_phy,ivgtyp,gravity,dt,        &
               drydep_flux,tend_chem_settle,dbg_opt,settling_opt,vg,         &
               ids,ide, jds,jde, kds,kde,                                    &
               ims,ime, jms,jme, kms,kme,                                    &
               its,ite, jts,jte, kts,kte, curr_secs                          )
!
! compute dry deposition velocity for aerosol particles
! Based on Emerson et al. (2020), PNAS,
! www.pnas.org/cgi/doi/10.1073/pnas.2014761117
! Code adapted from Hee-Ryu and Min, (2022):
! https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2021MS002792
!----------------------------------------------------------------------
  IMPLICIT NONE
 
       INTEGER,  INTENT(IN   ) ::  num_chem,                               &
                                   ids,ide, jds,jde, kds,kde,              &
                                   ims,ime, jms,jme, kms,kme,              &
                                   its,ite, jts,jte, kts,kte

       REAL(RKIND) :: curr_secs

       REAL(RKIND), DIMENSION( ims:ime , jms:jme )        ,            &
           INTENT(IN) :: ustar, rmol, znt, snowh
       REAL(RKIND), DIMENSION( ims:ime , kms:kme , jms:jme ),          &
           INTENT(IN   ) :: t_phy, rho_phy, p_phy, delz                    
       INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  ivgtyp          
       REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ),  &
                                             INTENT(IN ) :: chem
       REAL(RKIND), INTENT(IN) :: gravity,dt
       LOGICAL, INTENT(IN) :: dbg_opt
       INTEGER, INTENT(IN) :: settling_opt
 !
 ! Output arrays 
       REAL(RKIND), DIMENSION( ims:ime, jms:jme, 1:num_chem ), INTENT(INOUT)       :: drydep_flux
       REAL(RKIND), DIMENSION( ims:ime, jms:jme, 1:num_chem ), INTENT(OUT)         :: ddvel
       REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem), INTENT(OUT) :: tend_chem_settle, vg
 ! Local
       real(RKIND), DIMENSION( its:ite, kts:kte, jts:jte)     :: airkinvisc        ! Air kinetic viscosity [cm2/s]
       real(RKIND), DIMENSION( its:ite, kts:kte, jts:jte)     :: freepath          ! Air molecular freepath [cm]
       real(RKIND), DIMENSION( its:ite, jts:jte)              :: aer_res           ! Aerodynmaic resistance
       real(RKIND), DIMENSION( its:ite, jts:jte)              :: A, eps0           ! land surface params [-]
       real(RKIND) :: Cc                      ! Cunningham/slip correction factor [-]
       real(RKIND) :: DDp, Eb                 ! Brownian diffusion []
       real(RKIND) :: Eim                     ! Impaction []
       real(RKIND) :: Ein                     ! Interception []
       real(RKIND) :: Sc                      ! Schmit number []
       real(RKIND) :: St                      ! Stokes number []
       real(RKIND) :: dp                      ! aerosol diameter [cm]
       real(RKIND) :: aerodens                ! aerosol density [g/cm3] 
       real(RKIND) :: Rs                      ! Surface resistance
       real(RKIND) :: growth_fac,vsettl,dtmax,conver,converi,dzmin
       real(RKIND) :: rmol_local
       integer :: i, j, k, ntdt, nv
       integer :: icall=0
!> -- Gas constant
       real(RKIND), parameter :: RSI = 8.314510_RKIND
       real(RKIND) :: four_ninths, g100

       logical, parameter :: do_timing = .true.

       integer, parameter :: max_iter_settle = 5

       real(RKIND) :: ust

       tend_chem_settle(:,:,:,:) = 0._RKIND
       ddvel(:,:,:)              = 0._RKIND
       vg(:,:,:,:)               = 0._RKIND
 !      growth_fac        = 1._RKIND
       conver            = 1.e-9_RKIND
       converi           = 1.e9_RKIND
 !      four_ninths       = 4._RKIND / 9._RKIND
       g100             = gravity * 1.e2_RKIND
       ntdt=INT(dt)
 
       if (mod(int(curr_secs),1800) .eq. 0) then
           icall = 0
       endif
       
! Set up 2D vars
       ! MODIS type lu, large roughness lengths (e.g., urban or forest)
       ! -----------------------------------------------------------------------
       ! *** TO DO -- set A and eps0 for all land surface types *** !!!
       ! -----------------------------------------------------------------------
       ! Set if snow greater than 1 cm
       do j = jts, jte
       do i = its, ite
          if ( ivgtyp(i,j) .eq. 13 .or. ivgtyp(i,j) .le. 5 ) then ! Forest
             A(i,j) = A_for
             eps0(i,j) = eps0_for
          else if ( ivgtyp(i,j) .eq. 17 ) then ! water
             A(i,j) = A_wat
             eps0(i,j) = eps0_wat
          else ! otherwise
             A(i,j) = A_grs
             eps0(i,j) = eps0_grs
          end if
          rmol_local = rmol(i,j)
          call depvel( rmol_local, dep_ref_hgt, znt(i,j), ustar(i,j), aer_res(i,j) )
       enddo
       enddo
!
! Set up 3D vars
!       if  (do_timing) call mpas_timer_start('dep_prep')
       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          airkinvisc(i,k,j)  = ( 1.8325e-4_RKIND * ( 416.16_RKIND / ( t_phy(i,k,j) + 120._RKIND) ) *   &
                       ( ( t_phy(i,k,j) / 296.16_RKIND )**1.5_RKIND ) ) / (rho_phy(i,k,j) * 1.e-3_RKIND) ! Convert density to mol/cm^3
          ! Air molecular freepath (cm)  
          freepath(i,k,j)    = 7.39758e-4_RKIND * airkinvisc(i,k,j) / sqrt( t_phy(i,k,j) )
!          delz_flip(i,k,j)   = delz(i,kte-kts+k,j)
       enddo
       enddo
       enddo
!      
! 3D + chem --> vg
!       if  (do_timing) call mpas_timer_start('vg_and_ddvel_calc')

       if ( settling_opt .eq. 1 ) then
          do nv = 1, num_chem
             if (aero_diam(nv) .lt. 0) cycle  ! At some point we'll do something different for gasses
             ! Convert diameter to cm and aerodens to g/cm3
             dp       = aero_diam(nv) * 100._RKIND
             aerodens = aero_dens(nv) * 1.e-3_RKIND
             do j = jts, jte
             do k = kts, kte
             do i = its, ite
                ! Cunningham correction factor
                Cc = 1._RKIND + 2._RKIND * freepath(i,k,j) / dp * ( 1.257_RKIND + 0.4_RKIND*exp( -0.55_RKIND * dp / freepath(i,k,j) ) )
                ! Gravitational Settling
                vg(i,k,j,nv) = aerodens * dp * dp * g100 * Cc / &       ! Convert gravity to cm/s^2
                       ( 18._RKIND * airkinvisc(i,k,j) * (rho_phy(i,k,j)*1.e-3_RKIND) ) ! Convert density to mol/cm^
             enddo
             enddo
             enddo
          enddo
       elseif (settling_opt .eq. 2 ) then
        ! Calculate setting velocities based on AERSETT, accounting for non-sphericity
           call particle_settling_aersett(t_phy,rho_phy,p_phy,vg,num_chem,  &
                           ids,ide, jds,jde, kds,kde,                       &
                           ims,ime, jms,jme, kms,kme,                       &
                           its,ite, jts,jte, kts,kte                        )
      
       endif

! 2D + chem (surface dep)
       k=kts 
       do nv = 1, num_chem
          if (aero_diam(nv) .lt. 0) cycle  ! At some point we'll do something different for gasses
          ! Convert diameter to cm and aerodens to g/cm3
          dp       = aero_diam(nv) * 100._RKIND
          aerodens = aero_dens(nv) * 1.e-3_RKIND
          do j = jts, jte
          do i = its, ite
                ust        = ustar(i,j)*1.e2_RKIND
                ! Cunningham correction factor
                Cc = 1._RKIND + 2._RKIND * freepath(i,k,j) / dp * ( 1.257_RKIND + 0.4_RKIND*exp( -0.55_RKIND * dp / freepath(i,k,j) ) )
                ! Brownian Diffusion
                DDp = ( boltzmann * t_phy(i,k,j) ) * Cc / (3._RKIND * pi * airkinvisc(i,k,j) * (rho_phy(i,k,j)*1.e-3_RKIND)  * dp) ! Convert density to mol/cm^3
                ! Schmit number
                Sc = airkinvisc(i,k,j) / DDp
                ! Brownian Diffusion
                Eb = Cb * Sc**(-0.666666667_RKIND)
                ! Stokes number
                St = ust*ust * vg(i,k,j,nv) / airkinvisc(i,k,j) / g100 ! Convert ustar to cm/s, gravity to cm/s^2
                ! Impaction 
                Eim = Cim * ( St / ( alpha + St ) )**1.7_RKIND
                ! Interception
                Ein = Cin * ( dp / A(i,j) )**vv
                ! Surface resistance
                Rs = 1._RKIND / ( ust * ( Eb + Eim + Ein) * eps0(i,j) ) ! Convert ustar to cm/s
                ! Compute final ddvel = aer_res + RS, set max at max_dep_vel in dep_data_mod.F[ m/s]
                ! The /100. term converts from cm/s to m/s, required for MYNN.
                if ( settling_opt .gt. 0 ) then
                      ddvel(i,j,nv) = max(min( ( vg(i,k,j,nv) + 1._RKIND / (aer_res(i,j)+Rs) )*1.e-2_RKIND, max_dep_vel),0._RKIND)
                else
                      ddvel(i,j,nv) = max(min( ( 1._RKIND / (aer_res(i,j)+Rs) )*1.e-2_RKIND, max_dep_vel),0._RKIND)
                endif
                if ( dbg_opt .and. (icall .le. n_dbg_lines) ) then
                   icall = icall + 1
                endif
                drydep_flux(i,j,nv) = drydep_flux(i,j,nv) + chem(i,k,j,nv)*rho_phy(i,k,j)*ddvel(i,j,nv)*dt*10._RKIND
          enddo
          enddo
       enddo
        
end subroutine dry_dep_driver_emerson
!
!--------------------------------------------------------------------------------
!
subroutine depvel( rmol, zr, z0, ustar, aer_res )
!--------------------------------------------------
!     THIS FUNCTION HAS BEEN DESIGNED TO EVALUATE AN UPPER LIMIT
!     FOR THE POLLUTANT DEPOSITION VELOCITY AS A FUNCTION OF THE
!     SURFACE ROUGHNESS AND METEOROLOGICAL CONDITIONS.
!     PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner  (Feb. 1991)
!                  by Winfried Seidl     (Aug. 1997)
!.....PROGRAM VARIABLES...
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     Z0       - SURFACE ROUGHNESS HEIGHT
!     USTAR    - FRICTION VELOCITY U*
!     AER_RES  - AERODYNAMIC RESISTANCE
!.....REFERENCES...
!     MCRAE, G.J. ET AL. (1983) MATHEMATICAL MODELING OF PHOTOCHEMICAL
!       AIR POLLUTION, ENVIRONMENTAL QUALITY LABORATORY REPORT 18,
!       CALIFORNIA INSTITUTE OF TECHNOLOGY, PASADENA, CALIFORNIA.
!.....RESTRICTIONS...
!     1. THE MODEL EDDY DIFFUSIVITIES ARE BASED ON MONIN-OBUKHOV
!        SIMILARITY THEORY AND SO ARE ONLY APPLICABLE IN THE
!        SURFACE LAYER, A HEIGHT OF O(30M).
!     2. ALL INPUT UNITS MUST BE CONSISTENT
!     3. THE PHI FUNCTIONS USED TO CALCULATE THE FRICTION
!        VELOCITY U* AND THE POLLUTANT INTEGRALS ARE BASED
!        ON THE WORK OF BUSINGER ET AL.(1971).
!     4. THE MOMENTUM AND POLLUTANT DIFFUSIVITIES ARE NOT
!        THE SAME FOR THE CASES L<0 AND L>0.
!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        REAL(RKIND), intent(in)    :: ustar, z0, zr
        REAL(RKIND), intent(out)   :: aer_res
        REAL(RKIND), intent(inout) :: rmol
!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: l
        REAL(RKIND)    :: ao, ar, polint, vk
!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC alog
!--------------------------------------------------
!     Set the von Karman constant
!--------------------------------------------------
        vk = karman
!--------------------------------------------------
!     DETERMINE THE STABILITY BASED ON THE CONDITIONS
!             1/L < 0 UNSTABLE
!             1/L = 0 NEUTRAL
!             1/L > 0 STABLE
!--------------------------------------------------
        if(abs(rmol) < 1.E-6 ) rmol = 0._RKIND
        IF ( rmol < 0._RKIND ) THEN
          ar = ((1._RKIND-9._RKIND*zr*rmol)**(0.25)+0.001_RKIND)**2.
          ao = ((1._RKIND-9._RKIND*z0*rmol)**(0.25)+0.001_RKIND)**2.
          polint = 0.74_RKIND*(alog((ar-1._RKIND)/(ar+1._RKIND))-alog((ao-1._RKIND)/(ao+1._RKIND)))
        ELSEIF ( rmol==0._RKIND ) THEN
          polint = 0.74_RKIND*alog(zr/z0)
        ELSE
          polint = 0.74_RKIND*alog(zr/z0) + 4.7_RKIND*rmol*(zr-z0)
        END IF
        !vgpart = ustar*vk/polint
        aer_res = polint/(karman*max(ustar,1.0e-4_RKIND)) * 1.e-2_RKIND ! convert to s/cm
end subroutine depvel
!
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
!
subroutine particle_settling_wrapper(tend_chem_settle,chem,rho_phy,delz_flip,vg,  &
                             dt, kts,kte,its,ite,jts,jte,num_chem,ims,ime, jms,jme, kms,kme       )
     IMPLICIT NONE
     
     INTEGER, INTENT(IN ) :: kts, kte,its,ite,jts,jte,num_chem,ims,ime, jms,jme, kms,kme 
     REAL(RKIND), INTENT(IN) :: dt
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT (IN)  :: rho_phy
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT (IN)  :: delz_flip
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(IN) :: chem
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(IN) :: vg
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT) :: tend_chem_settle
     
     REAL(RKIND) :: dt_settl, growth_fac, four_ninths, dzmin, vsettl, dtmax
     INTEGER     :: ntdt, ndt_settl
!
!--- Local------
     integer, parameter :: max_iter_settle = 10
     real(RKIND), parameter :: one_over_dyn_visc = 1.e5_RKIND ! 5.5248E5_RKIND
     INTEGER :: k,n,l2,i,j,nv
     REAL(RKIND) :: temp_tc, transfer_to_below_level, vd_wrk1
  
     growth_fac        = 1._RKIND
     four_ninths       = 4._RKIND / 9._RKIND
     ntdt = INT(dt)

     do nv = 1,num_chem
     if (aero_diam(nv) .lt. 0) cycle  ! At some point we'll do something different for gasses
     ! -- NOTE, diameters and densities are NOT converted to cm and g/cm3 like in Emerson
     vsettl = four_ninths * gravity * aero_dens(nv) * ( growth_fac * ( 0.5_RKIND * aero_diam(nv) ))**2.0_RKIND * one_over_dyn_visc

     do j = jts,jte
     do i = its,ite
     dzmin = delz_flip(i,kts,j)
     do k = kts,kte
        dzmin = min(dzmin,delz_flip(i,k,j))
     enddo
     ! -- Get necessary info for settling
     ! -- Determine the maximum time-step satisying the CFL condition:
     ! -- dt_settl calculations (from original coarsepm_settling)
     ! 1.5E-5 = dyn_visc --> dust_data_mod.F90

     dtmax = dzmin / vsettl
 
     ndt_settl = MAX( 1, INT( ntdt /dtmax) )
     ! Limit maximum number of iterations
     IF (ndt_settl > max_iter_settle) ndt_settl = max_iter_settle
     dt_settl = REAL(ntdt,kind=RKIND) /REAL(ndt_settl,kind=RKIND)

     do n = 1,ndt_settl
     transfer_to_below_level = 0._RKIND

     do k = kte,kts,-1

        l2 = kte - k + 1
     
        temp_tc = chem(i,k,j,nv)
     
        vd_wrk1 = dt_settl * vg(i,k,j,nv)*1.e-2_RKIND / delz_flip(i,l2,j)               ! convert vg to m/s
     
        temp_tc = temp_tc * (1._RKIND - vd_wrk1) + transfer_to_below_level
        if ( k .gt. kts ) then
        transfer_to_below_level =(temp_tc*vd_wrk1)*((delz_flip(i,l2,j) &
                        *rho_phy(i,k,j))/(delz_flip(i,l2+1,j)*rho_phy(i,k-1,j)))          ! [ug/kg]
        endif
        tend_chem_settle(i,k,j,nv) = tend_chem_settle(i,k,j,nv) + (temp_tc - chem(i,k,j,nv))

     enddo ! k
     enddo ! n
     enddo ! i
     enddo ! j
     enddo ! nv



end subroutine particle_settling_wrapper


subroutine particle_settling_aersett(t_phy,rho_phy,p_phy,vg, num_chem,          &
                                     ids,ide, jds,jde, kds,kde,                 &
                                     ims,ime, jms,jme, kms,kme,                 &
                                     its,ite, jts,jte, kts,kte                  )

   IMPLICIT NONE

   INTEGER, INTENT(IN ) :: kts, kte,its,ite,jts,jte,ims,ime, jms,jme, kms,kme, ids, ide, jds, jde, kds, kde
   INTEGER, INTENT(IN ) :: num_chem
   REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT (IN)  :: rho_phy, t_phy, p_phy
   REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT) :: vg

  real(RKIND), parameter  :: beta=1.458e-6_RKIND
     ! The beta constant in kg/(s.m.K^.5) in the expression for dynamic viscosity.
     ! Value from the US Standard Atmosphere, 1976
  real(RKIND), parameter  :: S=110.4_RKIND
       ! The Sutherland constant S=110.4 K in the expression for dynamic viscosity
     ! Value from the US Standard Atmosphere, 1976 (Eq.51)

   LOGICAL, PARAMETER :: use_LUT = .false.

   INTEGER :: i,j,k,nv

   REAL(RKIND) :: D, rho_p, rho_a, aspect_p, mu_a, mfp_a, temp_a

   do nv = 1,num_chem
    ! Particle diameter, density 
      D        = aero_diam(nv)
      rho_p    = aero_dens(nv)
      aspect_p = 1.05_RKIND !aero_aspect(nv)
    !  
      do j = jts,jte
      do k = kts,kte
      do i = its,ite
        ! Air density 
         rho_a  = rho_phy(i,k,j)
         mu_a  =  beta * t_phy(i,k,j)**1.5 / (t_phy(i,k,j)+S)
         mfp_a = sqrt(3.13159_RKIND/8._RKIND) * mu_a/0.4987445_RKIND * 1._RKIND/sqrt(p_phy(i,k,j)*rho_a) ! Jennings 1988
         temp_a = t_phy(i,k,j)
        !
         vg(i,k,j,nv) = vinfty_aersett_spheroids_avg(D, rho_p, rho_a, mu_a, mfp_a, aspect_p, temp_a, use_LUT)
        !
      end do
      end do
      end do
   end do


end subroutine particle_settling_aersett

  function vinfty_aersett_spheroids_avg(D, rho_p, rho_a, mu_a, mfp_a, aspect_p, temp_a, use_LUT)
    ! Settling speed for a spheroid following MaMa25.
    ! The probability distribution of orientation is taken into account as described in MaMa25.
    ! if use_LUT is set to true, exectution is faster to the price of a (vary) small numerical error.
    real(RKIND), intent(IN) :: D, rho_p, rho_a, mu_a, mfp_a, temp_a
    real(RKIND), intent(IN) :: aspect_p
    logical, intent(IN) :: use_LUT

    real(RKIND) :: vinfty_aersett_spheroids_avg

    real(RKIND) :: a, alpha, w0, v0, v90

    a = D * aspect_p**(2.0_RKIND / 3.0_RKIND)
    alpha = a * (1._RKIND - exp(3.0_RKIND * (1._RKIND-aspect_p))) * 9.8_RKIND * (3.14159 * D**3 / 6.0_RKIND) * (rho_p - rho_a) &
         & / (16 * boltzmann * temp_a)
 !
    w0 = ialpha_LUT(alpha)
 !
    v90 = vinfty_aersett_spheroids(D, rho_p, rho_a, mu_a, mfp_a, aspect_p, 90, use_LUT)
    if(w0<0.01) then
       vinfty_aersett_spheroids_avg = v90
    else
       v0 = vinfty_aersett_spheroids(D, rho_p, rho_a, mu_a, mfp_a, aspect_p, 0, use_LUT)
       vinfty_aersett_spheroids_avg = w0 * v0 + (1.0 - w0) * v90
    endif
  end function vinfty_aersett_spheroids_avg

  function ialpha_LUT(alpha)
    ! Calculate I(alpha) from the LUT initialized by aersett_init
    real(RKIND) :: ialpha_LUT
    real(RKIND), intent(IN) ::  alpha
    integer :: index
    if (alpha .le. alpha_min) then
       ialpha_LUT = 1. / 3.
    else if (alpha .ge. alpha_max) then
       ialpha_LUT = 0.
    else
       index = 1 + floor((alpha - alpha_min) / alpha_step)
       ialpha_LUT = ialpha_ar(index)
    endif
  end function ialpha_LUT

  function vinfty_aersett_spheroids(D, rho_p, rho_a, mu_a, mfp_a, aspect_p, angle_p, use_LUT)
    ! Settling speed for a spheroid following MaMa23. the spheroid has to be either horizontally (angle_p=90)
    ! or vertically (angle_p=0) oriented so.
    ! if use_LUT is set to true, exectution is faster to the price of a (vary) small numerical error. 
    real(RKIND), intent(IN) :: D, rho_p, rho_a, mu_a, mfp_a
    real(RKIND), intent(IN) :: aspect_p
    integer, intent(IN) :: angle_p
    logical, intent(IN) :: use_LUT

    real(RKIND) :: vinfty_aersett_spheroids, a_lam_phi, Cc
    real(RKIND) :: Kn, R_tilde, phi

    call phi_and_a(angle_p, aspect_p, phi, a_lam_phi)


    Kn = 2. * mfp_a/( D*aspect_p**(2./3.) * phi )
    Cc = 1. + Kn * ( 1.257 + 0.4 * exp( -1.1/Kn ) )
    vinfty_aersett_spheroids = 4. * Cc * D**2. * (rho_p - rho_a) * 9.8_RKIND / (3. * mu_a * a_lam_phi)
    R_tilde = Cc * D**3 * rho_a * (rho_p - rho_a) * 9.8_RKIND / (18. * mu_a**2)
    if(R_tilde>Re_0)then
      vinfty_aersett_spheroids = vinfty_aersett_spheroids * sfunc(R_tilde)
    endif

  end function vinfty_aersett_spheroids

  function sfunc(R)
    ! the S function described in MaMa23
    real(RKIND), intent(IN) :: R
    real(RKIND)             :: sfunc
    sfunc = 1. - (1. + (R /4.880)**(-0.4335))**(-1.905)
  end function sfunc


  subroutine phi_and_a(angle_p, aspect_p, phi, a)
    ! Calculate Phi and A from the formulae in MaMa23
    integer, intent(IN) :: angle_p
    real(RKIND), intent(IN) :: aspect_p
    real(RKIND), intent(OUT) :: phi, a
    real(RKIND) :: part1, part2, Ep, Gp, eccen, aspect_round
    real(RKIND), parameter  :: f = 0.9113

    eccen = (1. - aspect_p**(-2._RKIND))**(0.5)
    Ep = asin(eccen)/eccen
    Gp = (1./aspect_p) - Ep
    aspect_round = aspect_p**2. - 1.

    if ( angle_p==0 ) then

        part1 = ( 2. * aspect_p**2. - 1. )/( ( aspect_round )**0.5 ) * log( aspect_p + ( aspect_round )**0.5 ) - aspect_p
        part2 = 2. * Ep * f + (Gp/eccen**2.) * (eccen**2. * ( 4. - 2.*f) - 4. + ( 3. - ( Pi/( 2. * aspect_p**2. ) )) * f)
        phi = ( 1.657 / ( 8. * ( aspect_round ) ) ) * part1 * part2
        a   = 64 * aspect_p**(2./3.) * eccen**3. / ( -2. * eccen + ( 1 + eccen**2) * log((1+eccen)/(1-eccen)) )
    else if ( angle_p==90 ) then

        part1 = ( 2. * aspect_p**2. - 3. )/( ( aspect_round )**0.5 ) * log( aspect_p + ( aspect_round )**0.5 ) + aspect_p
        part2 = Ep * ( 4. + ( Pi/2 - 1. ) * f ) + ( Gp/eccen**2. ) * ( 2. + ( ( 4. * eccen**2 + Pi - 6. )/4. ) * f )
        phi   = ( 1.657 / ( 16. * ( aspect_round ) ) ) * part1 * part2
        a     = 64*aspect_p**(2./3.) * 2. * eccen**3. / (2. * eccen + ( 3. * eccen**2. - 1. ) * log((1+eccen)/(1-eccen)) )
    else
       print*, "Angle is either 0 or 90."
    end if

  end subroutine phi_and_a


end module dep_dry_mod_emerson
