!WRF:MEDIATION_LAYER:PHYSICS
!
! Contains initialization subroutine lightning_init and driver subroutine
! lightning_driver.
!
! History:
!   3.5?  - rewritten and added init, separate out flash rate
!           parameterization from emission
!   3.4.0 - Added cpm option
!   3.3.x - lightning_driver written by M. Barth called in
!           emission_driver in chem
!   4.0   - implementation in MPAS
!
! Contact: Jordan Schnell <jordan.schnell@noaa.gov>
!
!**********************************************************************

 MODULE module_lightning_driver
  use mpas_kind_types
  use mpas_smoke_init
  use mpas_pool_routines
  use module_peg_util

  implicit none

  private

  public :: lightning_init, lightning_driver
 CONTAINS

!**********************************************************************
!
! SUBROUTINE lightning_init
!
! Performs compatibility checks and zero out flash arrays at first timestep.
!
!**********************************************************************


 SUBROUTINE lightning_init ( dt,           &
                            ! Namelist control options
                             config_convection_scheme,config_microp_scheme,  &
                             config_lightning_option, lightning_dt,           &
                             lightning_start_seconds,                  &
                             lightning_activation_time,                              &
                             iccg_prescribed_num, iccg_prescribed_den, &
                             lightning_cellcount_method,                         &
                            ! Order dependent args for domain, mem, and tile dims
                             ids, ide, jds, jde, kds, kde,             &
                             ims, ime, jms, jme, kms, kme,             &
                             its, ite, jts, jte, kts, kte,             &
                            ! IC and CG flash rates and accumulated flash count
                             ic_flashcount, ic_flashrate,              &
                             cg_flashcount, cg_flashrate               &
                            )
!-----------------------------------------------------------------
 USE module_peg_util, only : peg_error_fatal
 IMPLICIT NONE
!-----------------------------------------------------------------

 REAL,     INTENT(IN)        :: dt
 CHARACTER(LEN=*),  INTENT(IN)  :: config_convection_scheme,config_microp_scheme,config_lightning_option
 REAL,     INTENT(IN)        :: lightning_dt, lightning_start_seconds
 REAL,     INTENT(INOUT)     :: lightning_activation_time
 REAL,     INTENT(IN)        :: iccg_prescribed_num, iccg_prescribed_den
 INTEGER,  INTENT(IN)        :: lightning_cellcount_method
 INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte

! Making these optional just in case qualitative lightning indices get implemented
 REAL, OPTIONAL, DIMENSION( ims:ime,jms:jme ), &
                 INTENT(OUT) :: ic_flashcount, ic_flashrate, &
                                cg_flashcount, cg_flashrate

 CHARACTER (LEN=80) :: message

!-----------------------------------------------------------------

!-- do not reset unless it is the first timestep or config_lightning_option is on
 IF (trim(config_lightning_option) .eq. 'off') THEN
   return
 ENDIF
!-- check to see if lightning_dt is less than zero
 IF ( lightning_dt < dt ) then
    call peg_error_fatal(6,' lightning_init: lightning_dt < model_timetep')
 ENDIF
 !
 ! Ligtning activation time, for when lightning_dt .ne. dt
 lightning_activation_time = lightning_start_seconds
!--  check to see if the prescribed IC:CG ratio is valid (0/0 and -1 are not allowed)
 IF (iccg_prescribed_den .eq. 0. .and. iccg_prescribed_num .eq. 0.) THEN
    call peg_error_fatal(6,' lightning_init: iccg_prescribed cannot be 0.0/0.0')
 ENDIF
 IF (iccg_prescribed_den .ne. 0.) THEN
    IF (iccg_prescribed_num/iccg_prescribed_den .eq. -1.) THEN
        call peg_error_fatal(6, ' lightning_init: iccg_prescribed cannot be -1')
    ENDIF
 ENDIF
 !-- check to see if config_lightning_option is valid
 !
 !   Add new schemes here so it is recognized and proper checks are performed
 !
 ! Convective resolved/permitted
 if ( trim(config_lightning_option) .eq. 'ltng_crm_PR92w' .or. trim(config_lightning_option) .eq. 'ltng_crm_PR92z' ) then
     IF ( trim(config_microp_scheme) .eq. 'off') THEN
       CALL peg_error_fatal(6,' lightning_init: Selected lightning option requires microphysics scheme' )
     ENDIF
     call mpas_log_write(' lightning_init: CRM lightning option used: ltng_crm_PR92w or ltng_crm_PR92z')
 ! Convective parameterized
 else if ( trim(config_lightning_option) .eq. 'ltng_cpm_PR92z' ) then
     IF ( trim(config_convection_scheme) .ne. 'cu_grell_freitas' .and. trim(config_convection_scheme) .ne. 'cu_grell_freitas_li' ) THEN
       call peg_error_fatal( 6, ' lightning_init: Selected lightning option requires GF or GFL convective parameterization' )
     ENDIF
     call mpas_log_write( ' lightning_init: CPM lightning option selected: ltng_cpm_PR92z')
 else if ( trim(config_lightning_option) .eq. 'ltng_lpi' ) then
     call mpas_log_write( ' lightning_init: LPIM lightning option selected: ltng_lpi')
 ! Non-existing options
 else
     CALL peg_error_fatal(6,'lightning_init: invalid config_lightning_option')
 endif 

!-- zero out arrays
 IF ( PRESENT( ic_flashcount ) .and. PRESENT( ic_flashrate ) .and. &
      PRESENT( cg_flashcount ) .and. PRESENT( cg_flashrate ) ) THEN
    ic_flashrate(:,:)  = 0.
    ic_flashcount(:,:) = 0.
    cg_flashrate(:,:)  = 0.
    cg_flashcount(:,:) = 0.
 ELSE
    CALL peg_error_fatal (6, ' lightning_init: flash arrays not present' )
 ENDIF

!-- Resolve auto-cellcount method option (lightning_cellcount_method=0)
 IF ( ( lightning_cellcount_method .eq. 0 ) .and. (trim(config_lightning_option) .eq. 'ltng_crm_PR92w' )) THEN
   CALL peg_error_fatal(6,'lightning_init:lightning_cellcount_method = 0 and needs to be set, 1 for dx > 10km, 2 for dx<12km')
 ENDIF

 END SUBROUTINE lightning_init


!**********************************************************************
!
! SUBROUTINE lightning_driver
!
! Redirect to the appropriate lightning subroutine.
!
!**********************************************************************

 SUBROUTINE lightning_driver ( &
                          ! Frequently used prognostics
                            curr_secs, dt, area,                  &
                            xlat, xlon, xland, ht,                &
                            t_phy, p_phy, rho, u, v, w,           &
                            th_phy, pi_phy,dz8w,                  &  
                            z,                                    &
                            qv,qc,qr,qi,qs,qg,                    &
                          ! Scheme specific prognostics
                            ktop_deep,                            &
                            refl,                                 &
                          ! Mandatory namelist inputs
                            config_lightning_option,                     &
                            lightning_dt,                         &
                            lightning_start_seconds,              &
                            lightning_activation_time,                          &
                            flashrate_factor,                     &
                          ! IC:CG namelist settings
                            iccg_method,                          &
                            iccg_prescribed_num,                  &
                            iccg_prescribed_den,                  &
                          ! Scheme specific namelist inputs
                            lightning_cellcount_method,                     &
                            cldtop_adjustment,                    &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
                          ! Mandatory outputs for all quantitative schemes
                            ic_flashcount, ic_flashrate,          &
                            cg_flashcount, cg_flashrate,           &
                            lpi                                   &
                            )
!-----------------------------------------------------------------
! Framework
! Parameterization options
 USE module_ltng_crmpr92       ! config_lightning_option == 1,   ltng_crm_PR92w
                               ! config_lightning_option == 2,   ltng_crm_PR92z
 USE module_ltng_cpmpr92z      ! config_lightning_option == 11,  ltng_cpm_PR92z

! IC:CG methods
 USE module_ltng_iccg

! LPI 
  USE module_ltng_lpi

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL(RKIND), INTENT(IN   )    ::       curr_secs
 REAL(RKIND),    INTENT(IN   )    ::       dt

 REAL(RKIND),    DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: xlat, xlon, xland, ht, area
 REAL(RKIND),    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: t_phy, p_phy, rho
 REAL(RKIND),    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: th_phy, pi_phy, dz8w  
 REAL(RKIND),    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: u, v, w, z
 REAL(RKIND),    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: qv,qc,qr,qi,qs,qg

! Scheme specific prognostics
 INTEGER, DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: ktop_deep     ! model LNB from config_convection_scheme
 REAL(RKIND),    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: refl          ! reflectivity from config_microp_scheme

 REAL(RKIND), DIMENSION( ims:ime, jms:jme ), INTENT(INOUT), OPTIONAL :: LPI
! Mandatory namelist inputs
 CHARACTER(LEN=*), INTENT(IN   )    ::       config_lightning_option
 REAL(RKIND),    INTENT(IN   )    ::       lightning_dt, lightning_start_seconds, flashrate_factor
 REAL(RKIND),    INTENT(INOUT)    ::       lightning_activation_time

! IC:CG namelist settings
 INTEGER, INTENT(IN   )    ::       iccg_method
 REAL(RKIND),    INTENT(IN   )    ::       iccg_prescribed_num, iccg_prescribed_den

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       lightning_cellcount_method                    ! used in CRM
 REAL(RKIND),    INTENT(IN   )    ::       cldtop_adjustment                   ! used in CPM

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       its,ite, jts,jte, kts,kte

! Mandatory outputs for all quantitative schemes
 REAL(RKIND), OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: ic_flashcount , cg_flashcount
 REAL(RKIND), OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(  OUT) :: ic_flashrate  , cg_flashrate


! Local variables
 REAL(RKIND) :: LtngActivationTime
 REAL(RKIND) :: nextTime
 REAL(RKIND), DIMENSION( ims:ime, jms:jme ) :: total_flashrate
 CHARACTER (LEN=80) :: message
 LOGICAL :: do_ltng

 REAL(RKIND),    DIMENSION( ims:ime, kms:kme, jms:jme) :: qh

 REAL(RKIND), PARAMETER            :: reflthreshold = 20. ! reflectivity threshold for CRM schemes
 REAL(RKIND), DIMENSION( kms:kme ) :: cellcount

!-----------------------------------------------------------------

 qh(:,:,:) = 0._RKIND


 IF ( trim(config_lightning_option) .eq. 'off' ) RETURN

 nextTime = curr_secs + REAL(dt,8)
 LtngActivationTime = REAL(lightning_activation_time,8)
 do_ltng = LtngActivationTime >= curr_secs .and. LtngActivationTime <= nextTime

 IF( .not. do_ltng ) THEN
   RETURN
 ENDIF

!-----------------------------------------------------------------
! This driver performs several steps in order to produce lightning
! flash rate and flash count diagnostics:
!
! 1. Determine cloud extents for specific CRM schemes
! 2. Total flash rate assignment to 2D array
! 3. Partitioning of total lightning into IC & CG
! 4. Scale flash rate by flashrate_factor and lightning_dt
!
!-----------------------------------------------------------------

 IF ( trim(config_lightning_option) .eq. 'ltng_crm_PR92w' .or. &
      trim(config_lightning_option) .eq. 'ltng_crm_PR92z' ) THEN
   CALL countCells( &
          ! Inputs
            refl, reflthreshold, lightning_cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            its, ite, jts, jte, kts, kte,              &
          ! Outputs
            cellcount )
 ENDIF

!-----------------------------------------------------------------

    ! CRM lightning options
    if ( trim(config_lightning_option) .eq. 'ltng_crm_PR92w' ) then

        CALL ltng_crmpr92w ( &
                  ! Frequently used prognostics
                    area, xland, ht, z, t_phy,          &
                  ! Scheme specific prognostics
                    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    lightning_cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate                       &
                  )
    elseif ( trim(config_lightning_option) .eq. 'ltng_crm_PR92z' ) then
        CALL ltng_crmpr92z ( &
                  ! Frequently used prognostics
                    area, xland, ht, z, t_phy,          &
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    lightning_cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate                       &
                  )

    ! CPM lightning options
    elseif ( trim(config_lightning_option) .eq. 'ltng_cpm_PR92z' ) then

        CALL ltng_cpmpr92z ( &
                  ! Frequently used prognostics
                    area, xland, ht, z, t_phy,      &
                    ktop_deep, cldtop_adjustment,     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,     &
                    ims, ime, jms, jme, kms, kme,     &
                    its, ite, jts, jte, kts, kte,     &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate                   &
                  )

    ! LPI lightning options
    elseif ( trim(config_lightning_option) .eq. 'ltng_lpi' ) then
        CALL   calclpi(W=w,                              &
                     Z=z,                              &
                     PI_PHY=pi_phy, RHO_PHY=rho,             &
                     TH_PHY=TH_PHY,P_PHY=p_phy,                  &
                     DZ8w=dz8w,                          &
!qc_vis,qr_vis,qi_vis,qs_vis,qg_vis,   &
                     QV=qv,         &   !Qv=qv_curr,                         &
                     QC=qc,         &   !Qc=qc_curr,                         &
                     QR=qr,         &   !QR=qr_curr,                         &
                     QI=qi,         &   !QI=qi_curr,                         &
                     QS=qs,         &   !qs_curr,                         &
                     QG=qg,         &   !qg_curr,                         &
                     QH=qh,       &   !qh_curr,                         &
                  lpi=lpi &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
                 ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte)
    else
  
        CALL peg_error_fatal (6,' lightning_driver: The lightning option does not exist')
    endif

!-----------------------------------------------------------------

    IF ( iccg_method == 0 ) THEN
        ! Flash rate option defaults
        
        IF ( TRIM(config_lightning_option) == 'ltng_crm_PR92w' .OR. &
             TRIM(config_lightning_option) == 'ltng_crm_PR92z' .OR. &
             TRIM(config_lightning_option) == 'ltng_cpm_PR92z' ) THEN
             
            CALL iccg_boccippio( &
                        xlat, xlon,                                &
                        iccg_prescribed_num, iccg_prescribed_den,  &
                      ! Order dependent args for domain, mem, and tile dims
                        ids, ide, jds, jde, kds, kde,              &
                        ims, ime, jms, jme, kms, kme,              &
                        its, ite, jts, jte, kts, kte,              &
                      ! Input
                        total_flashrate,                           &
                      ! Output
                        ic_flashrate, cg_flashrate                 &
                      )
    
        ELSE
            CALL iccg_user_prescribed( &
                        iccg_prescribed_num,                  &
                        iccg_prescribed_den,                  &
                      ! Order dependent args for domain, mem, and tile dims
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte,         &
                      ! Input
                        total_flashrate,                      &
                      ! Output
                        ic_flashrate, cg_flashrate            &
                      )
        END IF
    
    ELSE IF ( iccg_method == 1 ) THEN
        ! Used-prescribed constant
        CALL iccg_user_prescribed( &
                    iccg_prescribed_num,                  &
                    iccg_prescribed_den,                  &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Input
                    total_flashrate,                      &
                  ! Output
                    ic_flashrate, cg_flashrate            &
                  )
    
    ELSE IF ( iccg_method == 2 ) THEN
        ! Boccippio et al, 2001
        CALL iccg_boccippio( &
                    xlat, xlon,                                &
                    iccg_prescribed_num, iccg_prescribed_den,  &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                  )
    
    ELSE IF ( iccg_method == 3 ) THEN
        ! Price and Rind, 1993
        IF ( TRIM(config_lightning_option) == 'ltng_crm_PR92w' .OR. &
             TRIM(config_lightning_option) == 'ltng_crm_PR92z' ) THEN
             
            CALL iccg_crm_pr93( &
                        refl, reflthreshold, t_phy, z,             &
                      ! Order dependent args for domain, mem, and tile dims
                        ids, ide, jds, jde, kds, kde,              &
                        ims, ime, jms, jme, kms, kme,              &
                        its, ite, jts, jte, kts, kte,              &
                      ! Input
                        total_flashrate,                           &
                      ! Output
                        ic_flashrate, cg_flashrate                 &
                    )
        ELSE
            CALL iccg_pr93( &
                        ktop_deep, cldtop_adjustment, t_phy, z,    &
                      ! Order dependent args for domain, mem, and tile dims
                        ids, ide, jds, jde, kds, kde,              &
                        ims, ime, jms, jme, kms, kme,              &
                        its, ite, jts, jte, kts, kte,              &
                      ! Input
                        total_flashrate,                           &
                      ! Output
                        ic_flashrate, cg_flashrate                 &
                    )
        END IF
    
    ELSE
        ! Invalid IC:CG method
        CALL peg_error_fatal (6,' lightning_driver: Invalid IC:CG method (iccg_method) = ')
    
    END IF
!-----------------------------------------------------------------

 ic_flashrate(its:ite,jts:jte) = ic_flashrate(its:ite,jts:jte) * flashrate_factor
 cg_flashrate(its:ite,jts:jte) = cg_flashrate(its:ite,jts:jte) * flashrate_factor

 ic_flashcount(its:ite,jts:jte) = ic_flashcount(its:ite,jts:jte) + ic_flashrate(its:ite,jts:jte) * lightning_dt
 cg_flashcount(its:ite,jts:jte) = cg_flashcount(its:ite,jts:jte) + cg_flashrate(its:ite,jts:jte) * lightning_dt
 
 do
   if( REAL(lightning_activation_time,8) <= nextTime ) then
     lightning_activation_time = lightning_activation_time + lightning_dt
   else
     exit
   endif
 enddo

!-----------------------------------------------------------------

 END SUBROUTINE lightning_driver


!**********************************************************************
!
! SUBROUTINE countCells
!
! For counting number of cells where reflectivity exceeds a certain
! threshold. Typically used by CRM schemes to redistribute lightning
! within convective cores.
!
! Departure from original implementation:
! Output includes domain-wide cellcounts if lightning_cellcount_method = 2
!
!**********************************************************************

 SUBROUTINE countCells( &
          ! Inputs
            refl, reflthreshold, lightning_cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            its, ite, jts, jte, kts, kte,              &
          ! Outputs
            cellcount )


 IMPLICIT NONE
!-----------------------------------------------------------------

! Inputs
 REAL,    DIMENSION( ims:ime,kms:kme,jms:jme ), INTENT(IN   ) :: refl
 REAL,    INTENT(IN   ) :: reflthreshold
 INTEGER, INTENT(IN   ) :: lightning_cellcount_method

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       its,ite, jts,jte, kts,kte


! Outputs
 REAL,    DIMENSION( kms:kme ), INTENT(  OUT) :: cellcount

! Local vars
 INTEGER :: i,k,j

!-----------------------------------------------------------------

 cellcount(kts:kte) = 0.
 DO j=jts,jte
   DO k=kts,kte
     DO i=its,ite
       IF ( refl(i,k,j) .gt. reflthreshold ) THEN
         cellcount(k) = cellcount(k) + 1
       ENDIF
     ENDDO
   ENDDO
 ENDDO


 END SUBROUTINE

 END MODULE module_lightning_driver
