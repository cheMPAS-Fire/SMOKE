MODULE module_smoke_diagnostics
!
!  Jordan Schnell (NOAA GSL)
!  For serious questions contact jordan.schnell@noaa.gov
!

   USE mpas_kind_types
   USE mpas_smoke_init
   USE rad_data_mod
   USE mpas_mie_module , only:optical_averaging
   USE module_data_rrtmgaeropt
   IMPLICIT NONE

   PRIVATE

   PUBLIC :: mpas_aod_diag , mpas_visibility_diag

CONTAINS

   SUBROUTINE mpas_aod_diag(Id,Curr_secs,Dtstep,Nbin_o,Chem,Aod3d,Aod3d_Simple,Rho_phy,Relhum,Dz8w,Num_chem,Tauaersw,Extaersw,Gaersw,Waersw,    &
                          & Bscoefsw,L2aer,L3aer,L4aer,L5aer,L6aer,L7aer,Tauaerlw,Extaerlw,Ids,Ide,Jds,Jde,Kds,Kde,Ims,Ime,Jms,Jme,&
                          & Kms,Kme,Its,Ite,Jts,Jte,Kts,Kte)

      INTEGER , INTENT(IN) :: Num_chem
      INTEGER , INTENT(IN) :: Ims
      INTEGER , INTENT(IN) :: Ime
      INTEGER , INTENT(IN) :: Jms
      INTEGER , INTENT(IN) :: Jme
      INTEGER , INTENT(IN) :: Kms
      INTEGER , INTENT(IN) :: Kme
      INTEGER , INTENT(IN) :: Id
      REAL(rkind) , INTENT(IN) :: Curr_secs
      REAL(rkind) , INTENT(IN) :: Dtstep
      INTEGER , INTENT(IN) :: Nbin_o
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:Num_chem) :: Chem
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Aod3d_Simple
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Aod3d
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Rho_phy
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Relhum
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Dz8w
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: Tauaersw
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: Extaersw
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: Gaersw
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: Waersw
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: Bscoefsw
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: L2aer
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: L3aer
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: L4aer
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: L5aer
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: L6aer
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) :: L7aer
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:NLWBANDS) :: Tauaerlw
      REAL(rkind) , INTENT(INOUT) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:NLWBANDS) :: Extaerlw
      INTEGER , INTENT(IN) :: Ids
      INTEGER , INTENT(IN) :: Ide
      INTEGER , INTENT(IN) :: Jds
      INTEGER , INTENT(IN) :: Jde
      INTEGER , INTENT(IN) :: Kds
      INTEGER , INTENT(IN) :: Kde
      INTEGER , INTENT(IN) :: Its
      INTEGER , INTENT(IN) :: Ite
      INTEGER , INTENT(IN) :: Jts
      INTEGER , INTENT(IN) :: Jte
      INTEGER , INTENT(IN) :: Kts
      INTEGER , INTENT(IN) :: Kte
!
! Local variable declarations
!
      REAL(rkind) :: alpha , ext , tau400 , tau550 , tau600, y_factor
      INTEGER :: i , j , k , nv

      CALL optical_averaging(Id,Curr_secs,Dtstep,Nbin_o,Chem,Num_chem,Dz8w,Rho_phy,Relhum,Tauaersw,Extaersw,Gaersw,Waersw,Bscoefsw,&
                           & L2aer,L3aer,L4aer,L5aer,L6aer,L7aer,Tauaerlw,Extaerlw,Ids,Ide,Jds,Jde,Kds,Kde,Ims,Ime,Jms,Jme,Kms,Kme,&
                           & Its,Ite,Jts,Jte,Kts,Kte)

! New AOD implemented here, Minsu Choi CIRES/NOAA GSL
      Aod3d(:,:,:) = 0._RKIND
      Aod3d_Simple(:,:,:) = 0._RKIND
      DO nv = 1 , Num_chem
         IF ( nv==p_ch4 ) CYCLE
         ext = sc_eff(nv) + ab_eff(nv)
         DO j = Jts , Jte
            DO k = Kts , Kte
               DO i = Its , Ite
                 ! ---------------------------------------------------------
                 ! Method 1: Aod3d_Simple (Accumulated per species)
                 ! Implements RH method, still reference missing  
                 ! Unit here is very confusing, need to confirm during PR
                 ! ---------------------------------------------------------
                 IF ( Relhum(i,k,j) > 0.3_RKIND ) THEN
                    y_factor = ((1.0_RKIND - 0.3_RKIND) / &
                               (1.0_RKIND - MIN(0.99_RKIND, Relhum(i,k,j))))**0.18_RKIND
                 ELSE
                    y_factor = 1.0_RKIND
                 ENDIF
                 Aod3d_Simple(i,k,j) = Aod3d_Simple(i,k,j) + &
                                       1.e-6_RKIND * (ext * y_factor) * &
                                       chem(i,k,j,nv) * rho_phy(i,k,j) * dz8w(i,k,j)
                 ! ---------------------------------------------------------
                 ! Method 2: Aod3d (Mie Calculation)
                 ! Implements Angstrom Exponent interpolation (400nm -> 550nm)
                 ! ---------------------------------------------------------
                 tau400 = Tauaersw(i,k,j,2)
                 tau600 = Tauaersw(i,k,j,3)
                 IF ( tau400 > 1.E-23_RKIND .AND. tau600 > 1.E-23_RKIND ) THEN
                    alpha = LOG(tau400 / tau600) / LOG(600.0_RKIND / 400.0_RKIND)
                    tau550 = tau400 * (400.0_RKIND / 550.0_RKIND)**alpha
                 ELSE
                    tau550 = 0.0_RKIND
                 ENDIF
                 Aod3d(i,k,j) = MAX(tau550, 0.0_RKIND)
               ENDDO
            ENDDO
         ENDDO
      ENDDO


   END SUBROUTINE mpas_aod_diag

   SUBROUTINE mpas_visibility_diag(Qcloud,Qrain,Qice,Qsnow,Qgrpl,Blcldw,Blcldi,Rho_phy,Wind10m,Wind,Rh2m,Rh,Qv,T2m,T,Coszen,       &
                                 & Extcoef55,Vis,Ids,Ide,Jds,Jde,Kds,Kde,Ims,Ime,Jms,Jme,Kms,Kme,Its,Ite,Jts,Jte,Kts,Kte)

      REAL(rkind) , PARAMETER :: VISFACTOR = 3.912_RKIND , RH_THRESHOLD = 0.3_RKIND
      INTEGER , INTENT(IN) :: Ims
      INTEGER , INTENT(IN) :: Ime
      INTEGER , INTENT(IN) :: Jms
      INTEGER , INTENT(IN) :: Jme
      INTEGER , INTENT(IN) :: Kms
      INTEGER , INTENT(IN) :: Kme
      INTEGER , INTENT(IN) :: Its
      INTEGER , INTENT(IN) :: Ite
      INTEGER , INTENT(IN) :: Jts
      INTEGER , INTENT(IN) :: Jte
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Qcloud
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Qrain
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Qice
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Qsnow
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Qgrpl
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Blcldw
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Blcldi
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Rho_phy
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Jms:Jme) :: Wind10m
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Wind
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Jms:Jme) :: Rh2m
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Rh
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Qv
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Jms:Jme) :: T2m
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: T
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Jms:Jme) :: Coszen
      REAL(rkind) , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Extcoef55
      REAL(rkind) , INTENT(OUT) , DIMENSION(Ims:Ime,Jms:Jme) :: Vis
      INTEGER , INTENT(IN) :: Ids
      INTEGER , INTENT(IN) :: Ide
      INTEGER , INTENT(IN) :: Jds
      INTEGER , INTENT(IN) :: Jde
      INTEGER , INTENT(IN) :: Kds
      INTEGER , INTENT(IN) :: Kde
      INTEGER , INTENT(IN) :: Kts
      INTEGER , INTENT(IN) :: Kte
      REAL(rkind) :: alpha_haze , bc , bg , bi , blcldi2 , blcldw2 , br , bs , extcoeff , extcoeff552 , ext_rh , haze_ext_coeff ,  &
                   & hydro_extcoeff , prob_ext_coeff_gt_p29 , qcloud2 , qgrpl2 , qice2 , qrain2 , qrh , qsnow2 , rhmax , rh_tmp ,  &
                   & sm_mse , tiny_number , tvd , vis_haze , vis_hydlith , vis_night , y , zen_fac
      INTEGER :: d , i , j , k
      REAL(rkind) , DIMENSION(Its:Ite,Jts:Jte) :: vis_alpha

  ! local
! Variables for simple RH calculation

!   sm_mse = 4.5_RKIND
      tiny_number = 1.E-12_RKIND

      DO j = Jts , Jte
         DO i = Its , Ite
      !Initialize
            qcloud2 = 0._RKIND
            blcldw2 = 0._RKIND
            qrain2 = 0._RKIND
            qice2 = 0._RKIND
            blcldi2 = 0._RKIND
            qsnow2 = 0._RKIND
            qgrpl2 = 0._RKIND
            extcoeff552 = 0._RKIND
      ! Follwowing UPP: CALVIS_GSD.f, take max of hydrometeors in lowest 3 levels
      ! - in UPP, only bottom rho_phy is used, shouldn't we use rho_phy from that level (as below)?
            k = Kts
!      do k = 1,3
            qcloud2 = Qcloud(i,k,j)*Rho_phy(i,k,j)*1000._RKIND !max(qcloud2,qcloud(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            blcldw2 = Blcldw(i,k,j)*Rho_phy(i,k,j)*1000._RKIND !max(blcldw2,blcldw(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            qrain2 = Qrain(i,k,j)*Rho_phy(i,k,j)*1000._RKIND  !max(qrain2,qrain(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            qice2 = Qice(i,k,j)*Rho_phy(i,k,j)*1000._RKIND   ! max(qice2,qice(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            blcldi2 = Blcldi(i,k,j)*Rho_phy(i,k,j)*1000._RKIND !max(blcldi2,blcldi(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            qsnow2 = Qsnow(i,k,j)*Rho_phy(i,k,j)*1000._RKIND  !max(qsnow2,qsnow(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            qgrpl2 = Qgrpl(i,k,j)*Rho_phy(i,k,j)*1000._RKIND  !max(qgrpl2,qgrpl(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
            extcoeff552 = Extcoef55(i,k,j)*1.E-3_RKIND
                                                   !max(extcoeff552,extcoeff55(i,k,j)) ! JLS - EXT55 is in units = 1/km, covert to 1/m
!      enddo

            bc = 144.7_RKIND*(qcloud2+blcldw2)**0.88_RKIND
            br = 2.24_RKIND*qrain2**0.75_RKIND
            bi = 327.8_RKIND*(qice2+blcldi2)**1.00_RKIND
            bs = 10.36_RKIND*qsnow2**0.7776_RKIND
            bg = 8.0_RKIND*qgrpl2**0.7500_RKIND

            hydro_extcoeff = (bc+br+bi+bs+bg)*1.E-3_RKIND
                                                  ! m^-1

! TODO, JLS, bring in QV 2M
            vis_haze = 999999._RKIND
            IF ( Qv(i,1,j)>0._RKIND ) vis_haze = 1500._RKIND*(105._RKIND-Rh2m(i,j))*(5._RKIND/min(1000._RKIND*Qv(i,1,j),5._RKIND))

      ! Follow UPP/CALVIS_GSD, First compute max RH of lowest 2 layers
            rhmax = max(Rh(i,1,j),Rh(i,2,j))
            qrh = max(0._RKIND,min(0.8_RKIND,(rhmax*0.01_RKIND-0.15_RKIND)))
            vis_haze = 90000._RKIND*exp(-2.5_RKIND*qrh)

      ! Calculate a Weibull function "alpha" term.  This can be
      ! used later with visibility (which acts as the "beta" term
      ! in the Weibull function) to create a probability distribution
      ! for visibility. Alpha can be thought of as the "level of
      ! certainty" that beta (model visibility) is correct. Fog is
      ! notoriously difficult to model. In the below algorithm,
      ! the alpha value (certainty) decreases as PWAT, mixing ratio,
      ! or winds decrease (possibly foggy conditions), but increases
      ! if RH decreases (more certainly not foggy).  If PWAT is lower
      ! then there is a higher chance of radiation fog because there
      ! is less insulating cloud above.
      ! -------------------------------------------------------------

!      alpha_haze=3.6
!      IF (q2m(i,j) .gt. 0.) THEN
!        alpha_haze=0.1 + pwater(i,j)/25.     + wind125m(i,j)/3. + &
!                        (100.-rh2m(i,j))/10. + 1./(1000.*q2m(i,j))
!        alpha_haze=min(alpha_haze,3.6)
!      ENDIF

      ! Calculate visibility from hydro/lithometeors
      ! Maximum visibility -> 999999 meters
      ! --------------------------------------------
            extcoeff = hydro_extcoeff + extcoeff552
            IF ( extcoeff>tiny_number ) THEN
               vis_hydlith = min(VISFACTOR/extcoeff,999999._RKIND)
            ELSE
               vis_hydlith = 999999._RKIND
            ENDIF

      ! -- Dec 2003 - Roy Rasmussen (NCAR) expression for night vs. day vis
      !   1.609 factor is number of km in mile.
            vis_night = 1.69*((vis_hydlith/1.609)**0.86)*1.609
            zen_fac = min(0.1,max(Coszen(i,j),0.))*10._RKIND
            vis_hydlith = zen_fac*vis_hydlith + (1.-zen_fac)*vis_night

      ! Calculate total visibility
      ! Take minimum visibility from hydro/lithometeors and haze
      ! Set alpha to be alpha_haze if haze dominates, or 3.6
      ! (a Guassian distribution) when hydro/lithometeors dominate
      ! ----------------------------------------------------------
            IF ( vis_hydlith<vis_haze ) THEN
               Vis(i,j) = vis_hydlith
         !vis_alpha(i,j)=3.6
            ELSE
               Vis(i,j) = vis_haze
         !vis_alpha(i,j)=alpha_haze
            ENDIF

         ENDDO
         ! i
      ENDDO
         ! j

   END SUBROUTINE mpas_visibility_diag

END MODULE module_smoke_diagnostics
