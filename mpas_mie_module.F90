MODULE mpas_mie_module
   USE module_data_rrtmgaeropt
   USE mpas_kind_types
   USE mpas_smoke_init
   IMPLICIT NONE
   INTEGER , PARAMETER , PRIVATE :: LUNERR = -1
CONTAINS
!---------------------------------------------------------------------------------------
! Original Source code was part of WRF-Chem (optical_averaging.F)
! Code has been originally contributed by PNNL and other WRF-Chem developers
!
! Original Source Code for cheMPAS-Fire has been written by Minsu Choi CIRES/NOAA GSL
! December, 2025, Minsu.Choi@noaa.gov, Minsu.Choi@colorado.edu
! Current version of the source code = simplified version
! For 2 size bins with 4 tracers
! Code has been further modernized/re-engineered by SPAG
!---------------------------------------------------------------------------------------
   SUBROUTINE optical_averaging(Id,Curr_secs,Dtstep,Nbin_o,Chem,Num_chem,Dz8w,Rho_phy,Relhum,Tauaersw,Extaersw,Gaersw,Waersw,      &
                              & Bscoefsw,L2aer,L3aer,L4aer,L5aer,L6aer,L7aer,Tauaerlw,Extaerlw,Ids,Ide,Jds,Jde,Kds,Kde,Ims,Ime,Jms,&
                              & Jme,Kms,Kme,Its,Ite,Jts,Jte,Kts,Kte)
 
      USE module_data_rrtmgaeropt
      USE mpas_kind_types
      USE mpas_smoke_init
 
      INTEGER , INTENT(IN) :: Id , Nbin_o
      INTEGER , INTENT(IN) :: Num_chem
      INTEGER , INTENT(IN) :: Ids , Ide , Jds , Jde , Kds , Kde
      INTEGER , INTENT(IN) :: Ims , Ime , Jms , Jme , Kms , Kme
      INTEGER , INTENT(IN) :: Its , Ite , Jts , Jte , Kts , Kte
      REAL(KIND=rkind) , INTENT(IN) :: Curr_secs
      REAL , INTENT(IN) :: Dtstep
      REAL , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,Num_chem) , INTENT(IN) :: Chem
      REAL , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) , INTENT(IN) :: Dz8w , Rho_phy , Relhum
 
      REAL , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) , INTENT(INOUT) :: Tauaersw , Extaersw , Gaersw , Waersw , Bscoefsw
      REAL , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:4) , INTENT(INOUT) :: L2aer , L3aer , L4aer , L5aer , L6aer , L7aer
      REAL , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,1:NLWBANDS) , INTENT(INOUT) :: Tauaerlw , Extaerlw
 
      REAL , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o) :: radius_wet , number_bin , radius_core
      COMPLEX , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o,1:4) :: swrefindx
      COMPLEX , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o,1:NLWBANDS) :: lwrefindx
 
      REAL , DIMENSION(1:Nbin_o,Kts:Kte) :: radius_wet_col , number_bin_col
      COMPLEX , DIMENSION(1:Nbin_o,Kts:Kte,1:4) :: swrefindx_col
      COMPLEX , DIMENSION(1:Nbin_o,Kts:Kte,1:NLWBANDS) :: lwrefindx_col
      REAL , DIMENSION(Kts:Kte) :: dz
 
      REAL , DIMENSION(4,Kts:Kte) :: swsizeaer , swextaer , swwaer , swgaer , swtauaer , swbscoef
      REAL , DIMENSION(4,Kts:Kte) :: l2 , l3 , l4 , l5 , l6 , l7
      REAL , DIMENSION(NLWBANDS,Kts:Kte) :: lwextaer , lwtauaer
 
      INTEGER :: i , j , k , isize , ns
      INTEGER :: ii , jj , kk
 
      ii = Its
      jj = Jts
      kk = Kts
 
      CALL optical_prep_simple(Nbin_o,Chem,Num_chem,Rho_phy,Relhum,radius_core,radius_wet,number_bin,swrefindx,lwrefindx,Ids,Ide,  &
                             & Jds,Jde,Kds,Kde,Ims,Ime,Jms,Jme,Kms,Kme,Its,Ite,Jts,Jte,Kts,Kte)
 
      DO j = Jts , Jte
         DO i = Its , Ite
 
            DO k = Kts , Kte
               dz(k) = Dz8w(i,k,j)
            ENDDO
 
            DO k = Kts , Kte
               DO isize = 1 , Nbin_o
                  number_bin_col(isize,k) = number_bin(i,k,j,isize)
                  radius_wet_col(isize,k) = radius_wet(i,k,j,isize)
                  swrefindx_col(isize,k,:) = swrefindx(i,k,j,isize,:)
                  lwrefindx_col(isize,k,:) = lwrefindx(i,k,j,isize,:)
               ENDDO
            ENDDO
 
            CALL mieaer(Id,i,j,Nbin_o,number_bin_col,radius_wet_col,swrefindx_col,lwrefindx_col,dz,Curr_secs,Kts,Kte,swsizeaer,    &
                      & swextaer,swwaer,swgaer,swtauaer,lwextaer,lwtauaer,l2,l3,l4,l5,l6,l7,swbscoef)
 
 
            DO k = Kts , Kte
               DO ns = 1 , 4
                  Tauaersw(i,k,j,ns) = max(swtauaer(ns,k),1.E-20)
                  Extaersw(i,k,j,ns) = max(swextaer(ns,k),1.E-20)
                  Gaersw(i,k,j,ns) = max(min(swgaer(ns,k),1.-1.E-8),1.E-20)
                  Waersw(i,k,j,ns) = max(min(swwaer(ns,k),1.-1.E-8),1.E-20)
                  Bscoefsw(i,k,j,ns) = max(swbscoef(ns,k),1.E-20)
               ENDDO
               L2aer(i,k,j,:) = l2(:,k)
               L3aer(i,k,j,:) = l3(:,k)
               L4aer(i,k,j,:) = l4(:,k)
               L5aer(i,k,j,:) = l5(:,k)
               L6aer(i,k,j,:) = l6(:,k)
               L7aer(i,k,j,:) = l7(:,k)
               DO ns = 1 , NLWBANDS
                  Tauaerlw(i,k,j,ns) = max(lwtauaer(ns,k),1.E-20)
                  Extaerlw(i,k,j,ns) = max(lwextaer(ns,k),1.E-20)
               ENDDO
            ENDDO
 
 
         ENDDO
      ENDDO
   END SUBROUTINE optical_averaging
 
   SUBROUTINE optical_prep_simple(Nbin_o,Chem,Num_chem,Rho_phy,Relhum,Radius_core,Radius_wet,Number_bin,Swrefindx,Lwrefindx,Ids,   &
                                & Ide,Jds,Jde,Kds,Kde,Ims,Ime,Jms,Jme,Kms,Kme,Its,Ite,Jts,Jte,Kts,Kte)
!---------------------------------------------------------------------------------------
! Original Source Code for cheMPAS-Fire has been written by Minsu Choi CIRES/NOAA GSL
! December, 2025, Minsu.Choi@noaa.gov, Minsu.Choi@colorado.edu
! Current version of the source code = simplified version
! For 2 size bins with 4 tracers
! Code has been further modernized/re-engineered by SPAG
!---------------------------------------------------------------------------------------
 
      USE module_data_rrtmgaeropt
      REAL , PARAMETER :: PI = 3.14159265358979324 , TINY = 1.0E-30 , KAPPA_SMOKE = 0.14 , KAPPA_DUST = 0.1

      INTEGER , INTENT(IN) :: Nbin_o
      INTEGER , INTENT(IN) :: Num_chem
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
      INTEGER , INTENT(IN) :: Kts
      INTEGER , INTENT(IN) :: Kte
      REAL , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme,Num_chem) :: Chem
      REAL , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Rho_phy
      REAL , INTENT(IN) , DIMENSION(Ims:Ime,Kms:Kme,Jms:Jme) :: Relhum
      REAL , INTENT(OUT) , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o) :: Radius_core
      REAL , INTENT(OUT) , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o) :: Radius_wet
      REAL , INTENT(OUT) , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o) :: Number_bin
      COMPLEX , INTENT(OUT) , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o,1:NSWBANDS) :: Swrefindx
      COMPLEX , INTENT(OUT) , DIMENSION(Its:Ite,Kts:Kte,Jts:Jte,1:Nbin_o,1:NLWBANDS) :: Lwrefindx
      INTEGER , INTENT(IN) :: Ids
      INTEGER , INTENT(IN) :: Ide
      INTEGER , INTENT(IN) :: Jds
      INTEGER , INTENT(IN) :: Jde
      INTEGER , INTENT(IN) :: Kds
      INTEGER , INTENT(IN) :: Kde
!
! Local variable declarations rewritten by SPAG
!
      REAL :: conv1a , dens_dust , dens_smoke , dgmin_cm , dgnum_um , dg_c_cm , dg_f_cm , dhi_um , dlo_tmp , dlo_um , dp_dry_a ,   &
            & dp_wet_a , duma , kappa_avg , lnsg , mass_dust , mass_dust_c , mass_dust_f , mass_smoke , mass_smoke_c ,             &
            & mass_smoke_f , mp_g , num_a , num_ac , num_aj , rho_eff , rh_clamped , sginia , sginic , sginin , sixpi , ss1 , ss2 ,&
            & ss3 , vbar_cm3 , vol_ac , vol_aj , vol_dry_a , vol_dust , vol_h2o , vol_oc , vol_smoke , vol_wet_a
      INTEGER :: i , iflag , isize , j , k , ns
      COMPLEX , DIMENSION(1:NLWBANDS) :: lwref_index_dust , lwref_index_h2o , lwref_index_smoke
      COMPLEX :: ri_ave_a , ri_dum
      COMPLEX , DIMENSION(1:NSWBANDS) :: swref_index_dust , swref_index_h2o , swref_index_smoke
      REAL , DIMENSION(1:Nbin_o) :: xdia_cm , xdia_um , xmas_sectc , xmas_secti , xmas_sectj , xnum_sectc , xnum_secti , xnum_sectj
!
! End of declarations rewritten by SPAG
!
  ! Aerosol hygroscopicitiy (dimensionless)
  ! Simple aerosol hygroscopic growth
 
      sixpi = 6.0/PI                    ! This is for OC
      dlo_um = 0.0390625
      dhi_um = 10.0
      dgmin_cm = 1.0E-07
      iflag = 2
      duma = 1.0
 
      dens_smoke = 1.6
                   ! was 1.8 in WRF-Chem, but 1.6 used here for fresh smoke
      dens_dust = 2.6
 
      sginin = 2.0
      sginia = 2.0
      sginic = 2.2
 
      dg_f_cm = 0.2E-4
      dg_c_cm = 2.0E-4
 
      dlo_tmp = dlo_um
 
      DO isize = 1 , Nbin_o
         xdia_um(isize) = (dlo_tmp+dlo_tmp*2.0)/2.0
         dlo_tmp = dlo_tmp*2.0
      ENDDO
 
      DO isize = 1 , Nbin_o
         xdia_cm(isize) = xdia_um(isize)*1.0E-4
      ENDDO
 
 
!  do ns = 1, nswbands
!    swref_index_smoke(ns) = cmplx(refrsw_smoke(ns), refisw_smoke(ns))
!    swref_index_dust(ns)  = cmplx(refrsw_dust(ns),  refisw_dust(ns))
!    swref_index_h2o(ns)  = cmplx(refrwsw(ns),  refiwsw(ns))
!  end do
!  do ns = 1, nlwbands
!    lwref_index_smoke(ns) = cmplx(refrlw_smoke(ns), refilw_smoke(ns))
!    lwref_index_dust(ns)  = cmplx(refrlw_dust(ns),  refilw_dust(ns))
!    lwref_index_h2o(ns)  = cmplx(refrwlw(ns),  refiwlw(ns))
!  end do
 
      DO ns = 1 , NSWBANDS
 
         swref_index_smoke(ns) = cmplx(Aer_optics(ID_SMOKE)%sw_real(ns),Aer_optics(ID_SMOKE)%sw_imag(ns))
         swref_index_dust(ns) = cmplx(Aer_optics(ID_DUST)%sw_real(ns),Aer_optics(ID_DUST)%sw_imag(ns))
 
         swref_index_h2o(ns) = cmplx(Aer_optics(ID_WATER)%sw_real(ns),Aer_optics(ID_WATER)%sw_imag(ns))
      ENDDO
 
      DO ns = 1 , NLWBANDS
         lwref_index_smoke(ns) = cmplx(Aer_optics(ID_SMOKE)%lw_real(ns),Aer_optics(ID_SMOKE)%lw_imag(ns))
 
         lwref_index_dust(ns) = cmplx(Aer_optics(ID_DUST)%lw_real(ns),Aer_optics(ID_DUST)%lw_imag(ns))
 
         lwref_index_h2o(ns) = cmplx(Aer_optics(ID_WATER)%lw_real(ns),Aer_optics(ID_WATER)%lw_imag(ns))
      ENDDO
 
! initialize - ms
      Swrefindx = cmplx(0.0,0.0)
      Lwrefindx = cmplx(0.0,0.0)
      Radius_wet = 0.0
      Number_bin = 0.0
      Radius_core = 0.0
 
      DO j = Jts , Jte
         DO k = Kts , Kte
            DO i = Its , Ite
 
! initialize - ms
               mass_smoke_f = 0.0
               mass_smoke_c = 0.0
               mass_dust_f = 0.0
               mass_dust_c = 0.0
               vol_aj = 0.0
               vol_ac = 0.0
               num_aj = 0.0
               num_ac = 0.0
 
        ! alt was used in WRF-Chem.
        ! here in mpas, we use rho
! Check, Minsu Choi CIRES/NOAA GSL
               conv1a = Rho_phy(i,k,j)*1.0E-12
                                          ! unit is currently not clear of the rho_phy
 
               mass_smoke_f = Chem(i,k,j,p_smoke_fine)*conv1a
               mass_dust_f = Chem(i,k,j,p_dust_fine)*conv1a
!        mass_smoke_c = chem(i,k,j,p_smoke_fine) * conv1a    ! for initial testing purpose
               mass_smoke_c = 0.0_RKIND
               mass_dust_c = Chem(i,k,j,p_dust_coarse)*conv1a
 
 
               vol_aj = (mass_smoke_f/dens_smoke) + (mass_dust_f/dens_dust)
               vol_ac = (mass_smoke_c/dens_smoke) + (mass_dust_c/dens_dust)
 
               rho_eff = (mass_smoke_f+mass_dust_f)/max(TINY,vol_aj)
               lnsg = log(max(1.01,sginia))
               vbar_cm3 = (PI/6.0)*(max(1.0E-12,dg_f_cm)**3)*exp(4.5*lnsg*lnsg)
               mp_g = rho_eff*vbar_cm3
 
               IF ( mp_g>TINY ) THEN
                  num_aj = (mass_smoke_f+mass_dust_f)/mp_g
               ELSE
                  num_aj = 0.0
               ENDIF
 
               rho_eff = (mass_smoke_c+mass_dust_c)/max(TINY,vol_ac)
               lnsg = log(max(1.01,sginic))
               vbar_cm3 = (PI/6.0)*(max(1.0E-12,dg_c_cm)**3)*exp(4.5*lnsg*lnsg)
               mp_g = rho_eff*vbar_cm3
               IF ( mp_g>TINY ) THEN
                  num_ac = (mass_smoke_c+mass_dust_c)/mp_g
               ELSE
                  num_ac = 0.0
               ENDIF
 
               xnum_secti = 0.0
               xmas_secti = 0.0
 
               IF ( num_aj>1.0E-20 .AND. vol_aj>1.0E-20 ) THEN
                  ss1 = log(sginia)
                  ss2 = exp(ss1*ss1*36.0/8.0)
                  ss3 = (sixpi*vol_aj/(num_aj*ss2))**0.3333333
                  dgnum_um = max(dgmin_cm,ss3)*1.0E+04
                  CALL sect02(dgnum_um,sginia,dens_smoke,iflag,duma,Nbin_o,dlo_um,dhi_um,xnum_sectj,xmas_sectj)
               ELSE
                  xnum_sectj = 0.0
                  xmas_sectj = 0.0
               ENDIF
 
               IF ( num_ac>1.0E-20 .AND. vol_ac>1.0E-20 ) THEN
                  ss1 = log(sginic)
                  ss2 = exp(ss1*ss1*36.0/8.0)
                  ss3 = (sixpi*vol_ac/(num_ac*ss2))**0.3333333
                  dgnum_um = max(dgmin_cm,ss3)*1.0E+04
                  CALL sect02(dgnum_um,sginic,dens_dust,iflag,duma,Nbin_o,dlo_um,dhi_um,xnum_sectc,xmas_sectc)
               ELSE
                  xnum_sectc = 0.0
                  xmas_sectc = 0.0
               ENDIF
 
               DO isize = 1 , Nbin_o
 
                  mass_smoke = mass_smoke_f*xmas_sectj(isize) + mass_smoke_c*xmas_sectc(isize)
                  mass_dust = mass_dust_f*xmas_sectj(isize) + mass_dust_c*xmas_sectc(isize)
 
                  vol_smoke = mass_smoke/dens_smoke
                  vol_oc = 0.9*mass_smoke/dens_smoke
                  vol_dust = mass_dust/dens_dust
 
                  vol_dry_a = vol_smoke + vol_dust
! Initial version of Mie code.
! Now we turn off this
! Treat hygroscopic growth
!          vol_wet_a = vol_dry_a
                  IF ( vol_dry_a>1.0E-20 ) THEN
                     kappa_avg = (KAPPA_SMOKE*vol_oc+KAPPA_DUST*vol_dust)/vol_dry_a
                  ELSE
                     kappa_avg = 0.0
                  ENDIF
                  rh_clamped = min(0.95,max(0.0,Relhum(i,k,j)))
                  vol_h2o = vol_dry_a*kappa_avg*(rh_clamped/(1.0-rh_clamped))
                  vol_wet_a = vol_dry_a + vol_h2o
 
 
                  num_a = num_aj*xnum_sectj(isize) + num_ac*xnum_sectc(isize)
 
                  IF ( num_a>1.0E-20 .AND. vol_dry_a>1.0E-20 ) THEN
                     dp_dry_a = (sixpi*vol_dry_a/num_a)**0.3333333
                     dp_wet_a = (sixpi*vol_wet_a/num_a)**0.3333333
                  ELSE
                     dp_dry_a = xdia_cm(isize)
                     dp_wet_a = xdia_cm(isize)
                  ENDIF
 
                  DO ns = 1 , NSWBANDS
                     IF ( vol_dry_a>1.0E-20 ) THEN
                        ri_dum = swref_index_smoke(ns)*vol_smoke + swref_index_dust(ns)*vol_dust + swref_index_h2o(ns)*vol_h2o
 
                        ri_ave_a = ri_dum/vol_wet_a
                     ELSE
                        ri_ave_a = cmplx(1.5,0.0)
                     ENDIF
                     Swrefindx(i,k,j,isize,ns) = ri_ave_a
                  ENDDO
 
                  DO ns = 1 , NLWBANDS
                     IF ( vol_dry_a>1.0E-20 ) THEN
                        ri_dum = lwref_index_smoke(ns)*vol_smoke + lwref_index_dust(ns)*vol_dust + lwref_index_h2o(ns)*vol_h2o
                        ri_ave_a = ri_dum/vol_wet_a
                     ELSE
                        ri_ave_a = cmplx(1.5,0.0)
                     ENDIF
                     Lwrefindx(i,k,j,isize,ns) = ri_ave_a
                  ENDDO
 
                  IF ( dp_wet_a/2.0<(dlo_um*1.0E-4)/2.0 ) THEN
                     Radius_wet(i,k,j,isize) = (dlo_um*1.0E-4)/2.0
                  ELSE
                     Radius_wet(i,k,j,isize) = dp_wet_a/2.0
                  ENDIF
                  Number_bin(i,k,j,isize) = num_a
                  Radius_core(i,k,j,isize) = 0.0
 
               ENDDO
 
            ENDDO
         ENDDO
      ENDDO
 
   END SUBROUTINE optical_prep_simple
 
 
   SUBROUTINE mieaer(Id,Iclm,Jclm,Nbin_a,Number_bin_col,Radius_wet_col,Swrefindx_col,Lwrefindx_col,Dz,Curr_secs,Kts,Kte,Swsizeaer, &
                   & Swextaer,Swwaer,Swgaer,Swtauaer,Lwextaer,Lwtauaer,L2,L3,L4,L5,L6,L7,Swbscoef)
                                           ! added bscoef JCB 2007/02/01
 
 
      USE module_peg_util , only:peg_error_fatal , peg_message
      USE module_data_rrtmgaeropt , only:Aer_optics , Num_aer_species
 
      INTEGER , PARAMETER :: NSPINT = 4 ! Num of spectral for FAST-J
      INTEGER , INTENT(IN) :: Kts , Kte
      INTEGER , INTENT(IN) :: Id , Iclm , Jclm , Nbin_a
      REAL(KIND=rkind) , INTENT(IN) :: Curr_secs
 
      REAL , DIMENSION(1:NSPINT,Kts:Kte) , INTENT(OUT) :: Swsizeaer , Swextaer , Swwaer , Swgaer , Swtauaer
      REAL , DIMENSION(1:NLWBANDS,Kts:Kte) , INTENT(OUT) :: Lwextaer , Lwtauaer
      REAL , DIMENSION(1:NSPINT,Kts:Kte) , INTENT(OUT) :: L2 , L3 , L4 , L5 , L6 , L7
      REAL , DIMENSION(1:NSPINT,Kts:Kte) , INTENT(OUT) :: Swbscoef  !JCB 2007/02/01
      REAL , INTENT(IN) , DIMENSION(1:Nbin_a,Kts:Kte) :: Number_bin_col
      REAL , INTENT(INOUT) , DIMENSION(1:Nbin_a,Kts:Kte) :: Radius_wet_col
      COMPLEX , INTENT(IN) , DIMENSION(1:Nbin_a,Kts:Kte,NSPINT) :: Swrefindx_col
      COMPLEX , INTENT(IN) , DIMENSION(1:Nbin_a,Kts:Kte,NLWBANDS) :: Lwrefindx_col
      REAL , INTENT(IN) , DIMENSION(Kts:Kte) :: Dz
 
        !fitting variables
      INTEGER LTYPE   ! total number of indicies of refraction
      PARAMETER (LTYPE=1)      ! bracket refractive indices based on information from Rahul, 2002/11/07
      INTEGER nrefr , nrefi , nr , ni
      SAVE nrefr , nrefi
      COMPLEX*16 sforw , sback , tforw(2) , tback(2)
      REAL*8 pmom(0:7,1)
      LOGICAL , SAVE :: ini_fit   ! initial mie fit only for the first time step
      DATA ini_fit/.TRUE./
        ! nsiz = number of wet particle sizes
      INTEGER , PARAMETER :: NSIZ = 200 , NLOG = 30
                                                !,ncoef=5
      REAL p2(NSIZ) , p3(NSIZ) , p4(NSIZ) , p5(NSIZ)
      REAL p6(NSIZ) , p7(NSIZ)
      LOGICAL perfct , anyang , prnt(2)
      REAL*8 xmu(1)
      DATA xmu/1./ , anyang/.FALSE./
      DATA prnt/.FALSE. , .FALSE./
      INTEGER numang , nmom , ipolzn , momdim
      DATA numang/0/
      COMPLEX*16 s1(1) , s2(1)
      REAL*8 mimcut
      DATA perfct/.FALSE./ , mimcut/0.0/
      DATA nmom/7/ , ipolzn/0/ , momdim/7/
      INTEGER n
      REAL*8 thesize      ! 2 pi radpart / waveleng = size parameter
      REAL*8 qext(NSIZ)   ! array of extinction efficiencies
      REAL*8 qsca(NSIZ)   ! array of scattering efficiencies
      REAL*8 gqsc(NSIZ)   ! array of asymmetry factor * scattering efficiency
      REAL qext4(NSIZ)            !  extinction, real*4
      REAL qsca4(NSIZ)            !  extinction, real*4
      REAL qabs4(NSIZ)            !  extinction, real*4
      REAL asymm(NSIZ)    ! array of asymmetry factor
      REAL sb2(NSIZ)       ! JCB 2007/02/01 - 4*abs(sback)^2/(size parameter)^2 backscattering efficiency
      COMPLEX*16 crefin , crefd , crefw
      SAVE crefw
      REAL , SAVE :: rmin , rmax  ! min, max aerosol size bin
      REAL bma , bpa
      REAL refr       ! real part of refractive index
      REAL refi       ! imaginary part of refractive index
      REAL refrmin   ! minimum of real part of refractive index
      REAL refrmax   ! maximum of real part of refractive index
      REAL refimin   ! minimum of imag part of refractive index
      REAL refimax   ! maximum of imag part of refractive index
      REAL drefr   ! increment in real part of refractive index
      REAL drefi   ! increment in imag part of refractive index
      COMPLEX specrefndx(LTYPE)   ! refractivr indices
      INTEGER , PARAMETER :: NAEROSOLS = 5
 
        !parameterization variables
      REAL weighte , weights , weighta
      REAL x
      REAL thesum   ! for normalizing things
      REAL sizem   ! size in microns
      INTEGER m , j , nc , klevel
      REAL pext             ! parameterized specific extinction (cm2/g)
      REAL pscat        !scattering cross section
      REAL pabs             ! parameterized specific extinction (cm2/g)
      REAL pasm         ! parameterized asymmetry factor
      REAL ppmom2       ! 2 Lengendre expansion coefficient (numbered 0,1,2,...)
      REAL ppmom3       ! 3     ...
      REAL ppmom4       ! 4     ...
      REAL ppmom5       ! 5     ...
      REAL ppmom6       ! 6     ...
      REAL ppmom7       ! 7     ...
      REAL sback2       ! JCB 2007/02/01 sback*conjg(sback)
      REAL cext(NCOEF) , casm(NCOEF) , cpmom2(NCOEF) , cabs(NCOEF)
      REAL cscat(NCOEF)    ! JCB 2004/02/09
      REAL cpmom3(NCOEF) , cpmom4(NCOEF) , cpmom5(NCOEF)
      REAL cpmom6(NCOEF) , cpmom7(NCOEF)
      REAL cpsback2p(NCOEF)   ! JCB 2007/02/09  - backscatter
      INTEGER itab , jtab
      REAL ttab , utab
      REAL , SAVE :: xrmin , xrmax , xr
      REAL rs(NSIZ)   ! surface mode radius (cm)
      REAL xrad   ! normalized aerosol radius
      REAL ch(NCOEF)   ! chebychev polynomial
      REAL , PARAMETER :: RMMIN = 0.005E-4 , RMMAX = 50.E-4
      INTEGER , PARAMETER :: LUNERR = -1
 
 
        !others
      INTEGER i , k , l , ns
      REAL pie , third
      INTEGER ibin
      CHARACTER*150 msg
      INTEGER kcallmieaer , kcallmieaer2
 
 
      pie = 4.*atan(1.)
      third = 1./3.
      rmin = RMMIN
      rmax = RMMAX
 
!######################################################################
!initial fitting to mie calculation based on Ghan et al. 2002 and 2007
!#####################################################################
      IF ( ini_fit ) THEN
         ini_fit = .FALSE.
 
  !----------------------------------------------------------------------
  !shortwave
  !---------------------------------------------------------------------
   ! wavelength loop
         DO ns = 1 , NSPINT
            Crefwsw(ns) = cmplx(Refrwsw(ns),Refiwsw(ns))
            refrmin = real(Crefwsw(ns))
            refrmax = real(Crefwsw(ns))
            refimin = -imag(Crefwsw(ns))
            refimax = -imag(Crefwsw(ns))
            DO l = 1 , Num_aer_species
               refr = Aer_optics(l)%sw_real(ns)
               refi = -Aer_optics(l)%sw_imag(ns)
 
               refrmin = min(refrmin,refr)
               refrmax = max(refrmax,refr)
               refimin = min(refimin,refi)
               refimax = max(refimax,refi)
            ENDDO
 
            drefr = (refrmax-refrmin)
            IF ( drefr>1.E-4 ) THEN
               nrefr = PREFR
               drefr = drefr/(nrefr-1)
            ELSE
               nrefr = 1
            ENDIF
 
            drefi = (refimax-refimin)
            IF ( drefi>1.E-4 ) THEN
               nrefi = PREFI
               drefi = drefi/(nrefi-1)
            ELSE
               nrefi = 1
            ENDIF
 
            bma = 0.5*log(rmax/rmin)
                                ! JCB
            bpa = 0.5*log(rmax*rmin)
                                ! JCB
 
            DO nr = 1 , nrefr
               DO ni = 1 , nrefi
 
                  Refrtabsw(nr,ns) = refrmin + (nr-1)*drefr
                  Refitabsw(ni,ns) = refimin/0.2*(0.2**real(ni))
                                                             !slightly different from Ghan and Zaveri
                  IF ( ni==nrefi ) Refitabsw(ni,ns) = -1.0E-20
                                                          ! JCB change
                  crefd = cmplx(Refrtabsw(nr,ns),Refitabsw(ni,ns))
 
!              mie calculations of optical efficiencies
                  DO n = 1 , NSIZ
                     xr = cos(pie*(float(n)-0.5)/float(NSIZ))
                     rs(n) = exp(xr*bma+bpa)
 
!                size parameter and weighted refractive index
                     thesize = 2.*pie*rs(n)/Wavmidsw(ns)
                     thesize = min(thesize,10000.D0)
 
                     CALL miev0(thesize,crefd,perfct,mimcut,anyang,numang,xmu,nmom,ipolzn,momdim,prnt,qext(n),qsca(n),gqsc(n),pmom,&
                              & sforw,sback,s1,s2,tforw,tback)
                     qext4(n) = qext(n)
                     qsca4(n) = min(qsca(n),qext(n))
                     qabs4(n) = qext4(n) - qsca4(n)
                     qabs4(n) = max(qabs4(n),1.E-20)
                                                ! avoid 0
                     asymm(n) = gqsc(n)/qsca4(n)
                                            ! assume always greater than zero
! coefficients of phase function expansion; note modification by JCB of miev0 coefficients
                     p2(n) = pmom(2,1)/pmom(0,1)*5.0
                     p3(n) = pmom(3,1)/pmom(0,1)*7.0
                     p4(n) = pmom(4,1)/pmom(0,1)*9.0
                     p5(n) = pmom(5,1)/pmom(0,1)*11.0
                     p6(n) = pmom(6,1)/pmom(0,1)*13.0
                     p7(n) = pmom(7,1)/pmom(0,1)*15.0
! backscattering efficiency, Bohren and Huffman, page 122
! as stated by Bohren and Huffman, this is 4*pie times what is should be
! may need to be smoothed - a very rough function - for the time being we won't apply smoothing
! and let the integration over the size distribution be the smoothing
                     sb2(n) = 4.0*sback*dconjg(sback)/(thesize*thesize)
                                                                   ! JCB 2007/02/01
!
!PMA makes sure it is greater than zero
!
                     qext4(n) = max(qext4(n),1.E-20)
                     qabs4(n) = max(qabs4(n),1.E-20)
                     qsca4(n) = max(qsca4(n),1.E-20)
                     asymm(n) = max(asymm(n),1.E-20)
                     sb2(n) = max(sb2(n),1.E-20)
 
                  ENDDO
 
                  CALL fitcurv(rs,qext4,Extpsw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv(rs,qabs4,Abspsw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv(rs,qsca4,Ascatpsw(1,nr,ni,ns),NCOEF,NSIZ)
                                                                      ! scattering efficiency
                  CALL fitcurv(rs,asymm,Asmpsw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv(rs,sb2,Sbackpsw(1,nr,ni,ns),NCOEF,NSIZ)
                                                                    ! backscattering efficiency
                  CALL fitcurv_nolog(rs,p2,Pmom2psw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv_nolog(rs,p3,Pmom3psw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv_nolog(rs,p4,Pmom4psw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv_nolog(rs,p5,Pmom5psw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv_nolog(rs,p6,Pmom6psw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv_nolog(rs,p7,Pmom7psw(1,nr,ni,ns),NCOEF,NSIZ)
 
               ENDDO
            ENDDO
         ENDDO        ! ns for shortwave
 
 
  !----------------------------------------------------------------------
  !longwave
  !---------------------------------------------------------------------
   ! wavelength loop
         DO ns = 1 , NLWBANDS
        !wavelength for longwave 1/cm --> cm
            Wavmidlw(ns) = 0.5*(1./Wavenumber1_longwave(ns)+1./Wavenumber2_longwave(ns))
 
            Crefwlw(ns) = cmplx(Refrwlw(ns),Refiwlw(ns))
            refrmin = real(Crefwlw(ns))
            refrmax = real(Crefwlw(ns))
            refimin = -imag(Crefwlw(ns))
            refimax = -imag(Crefwlw(ns))
 
            DO l = 1 , Num_aer_species
               refr = Aer_optics(l)%lw_real(ns)
               refi = -Aer_optics(l)%lw_imag(ns)
               refrmin = min(refrmin,refr)
               refrmax = max(refrmax,refr)
               refimin = min(refimin,refi)
               refimax = max(refimax,refi)
            ENDDO
 
            drefr = (refrmax-refrmin)
            IF ( drefr>1.E-4 ) THEN
               nrefr = PREFR
               drefr = drefr/(nrefr-1)
            ELSE
               nrefr = 1
            ENDIF
 
            drefi = (refimax-refimin)
            IF ( drefi>1.E-4 ) THEN
               nrefi = PREFI
               drefi = drefi/(nrefi-1)
            ELSE
               nrefi = 1
            ENDIF
 
            bma = 0.5*log(rmax/rmin)
                                ! JCB
            bpa = 0.5*log(rmax*rmin)
                                ! JCB
 
            DO nr = 1 , nrefr
               DO ni = 1 , nrefi
 
                  Refrtablw(nr,ns) = refrmin + (nr-1)*drefr
                  Refitablw(ni,ns) = refimin/0.2*(0.2**real(ni))
                                                             !slightly different from Ghan and Zaveri
                  IF ( ni==nrefi ) Refitablw(nrefi,ns) = -1.0E-21
                                                             ! JCB change
                  crefd = cmplx(Refrtablw(nr,ns),Refitablw(ni,ns))
 
!              mie calculations of optical efficiencies
                  DO n = 1 , NSIZ
                     xr = cos(pie*(float(n)-0.5)/float(NSIZ))
                     rs(n) = exp(xr*bma+bpa)
 
!                size parameter and weighted refractive index
                     thesize = 2.*pie*rs(n)/Wavmidlw(ns)
                     thesize = min(thesize,10000.D0)
 
                     CALL miev0(thesize,crefd,perfct,mimcut,anyang,numang,xmu,nmom,ipolzn,momdim,prnt,qext(n),qsca(n),gqsc(n),pmom,&
                              & sforw,sback,s1,s2,tforw,tback)
                     qext4(n) = qext(n)
                     qext4(n) = max(qext4(n),1.E-20)
                                                ! avoid 0
                     qsca4(n) = min(qsca(n),qext(n))
                     qsca4(n) = max(qsca4(n),1.E-20)
                                                ! avoid 0
                     qabs4(n) = qext4(n) - qsca4(n)
                     qabs4(n) = max(qabs4(n),1.E-20)
                                                ! avoid 0
                     asymm(n) = gqsc(n)/qsca4(n)
                                            ! assume always greater than zero
                     asymm(n) = max(asymm(n),1.E-20)
                                                !PMA makes sure it is greater than zero
                  ENDDO
!
               !if (nr==1.and.ni==1) then
               !endif
                  CALL fitcurv(rs,qext4,Extplw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv(rs,qabs4,Absplw(1,nr,ni,ns),NCOEF,NSIZ)
                  CALL fitcurv(rs,qsca4,Ascatplw(1,nr,ni,ns),NCOEF,NSIZ)
                                                                      ! scattering efficiency
                  CALL fitcurv(rs,asymm,Asmplw(1,nr,ni,ns),NCOEF,NSIZ)
               ENDDO
            ENDDO
         ENDDO        ! ns for longwave
 
 
      ENDIF !ini_fit
 
 
      xrmin = log(rmin)
      xrmax = log(rmax)
 
!######################################################################
!parameterization of mie calculation for shortwave
!#####################################################################
 
! begin level loop
      DO klevel = 1 , Kte
! sum densities for normalization
         thesum = 0.0
         DO m = 1 , Nbin_a
            thesum = thesum + Number_bin_col(m,klevel)
         ENDDO
! Begin shortwave spectral loop
         DO ns = 1 , NSWBANDS
 
!        aerosol optical properties
            Swtauaer(ns,klevel) = 0.
            Swwaer(ns,klevel) = 0.
            Swgaer(ns,klevel) = 0.
            Swsizeaer(ns,klevel) = 0.0
            Swextaer(ns,klevel) = 0.0
            L2(ns,klevel) = 0.0
            L3(ns,klevel) = 0.0
            L4(ns,klevel) = 0.0
            L5(ns,klevel) = 0.0
            L6(ns,klevel) = 0.0
            L7(ns,klevel) = 0.0
            Swbscoef(ns,klevel) = 0.0 ! JCB 2007/02/01 - backscattering coefficient
            IF ( thesum>1E-21 ) THEN      ! set everything = 0 if no aerosol !wig changed 0.0 to 1e-21, 31-Oct-2005
 
! loop over the bins
               DO m = 1 , Nbin_a
                             ! nbin_a is number of bins
! here's the size
                  sizem = Radius_wet_col(m,klevel)
                                               ! radius in cm
 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !ec  diagnostics
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ! check limits of particle size
          ! rce 2004-dec-07 - use klevel in write statements
                  IF ( Radius_wet_col(m,klevel)<=rmin ) THEN
                     Radius_wet_col(m,klevel) = rmin
                     WRITE (msg,'(a, 5i4,1x, e11.4)') 'mieaer: radius_wet set to rmin,'//'id,i,j,k,m,rm(m,k)' , Id , Iclm , Jclm , &
                          & klevel , m , Radius_wet_col(m,klevel)
                     CALL peg_message(LUNERR,msg)
                  ENDIF
                  IF ( Radius_wet_col(m,klevel)>rmax ) THEN
                     Radius_wet_col(m,klevel) = rmax
                 !only print when the number is significant
                     IF ( Number_bin_col(m,klevel)>=1.E-10 ) THEN
                        WRITE (msg,'(a, 5i4,1x, 2e11.4)') 'mieaer: radius_wet set to rmax,'//'id,i,j,k,m,rm(m,k),number' , Id ,    &
                             & Iclm , Jclm , klevel , m , Radius_wet_col(m,klevel) , Number_bin_col(m,klevel)
                        CALL peg_message(LUNERR,msg)
                     ENDIF
                  ENDIF
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
                  x = log(Radius_wet_col(m,klevel))
                                                ! radius in cm
                  crefin = Swrefindx_col(m,klevel,ns)
                  refr = real(crefin)
                  refi = -imag(crefin)
                  xrad = x
                  thesize = 2.0*pie*exp(x)/Wavmidsw(ns)
                ! normalize size parameter
                  xrad = (2*xrad-xrmax-xrmin)/(xrmax-xrmin)
 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !ec  diagnostics
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ! retain this diagnostic code
                  IF ( abs(refr)>10.0 .OR. abs(refr)<=0.001 ) THEN
                     WRITE (msg,'(a,1x, e14.5)') 'mieaer /refr/ outside range 1e-3 - 10 '//'refr= ' , refr
                     CALL peg_error_fatal(LUNERR,msg)
                  ENDIF
                  IF ( abs(refi)>10. ) THEN
                     WRITE (msg,'(a,1x, e14.5)') 'mieaer /refi/ >10 '//'refi' , refi
                     CALL peg_error_fatal(LUNERR,msg)
                  ENDIF
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
! interpolate coefficients linear in refractive index
! first call calcs itab,jtab,ttab,utab
                  itab = 0
                  CALL binterp(Extpsw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab,   &
                             & cext)
 
! JCB 2004/02/09  -- new code for scattering cross section
                  CALL binterp(Ascatpsw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cscat)
                  CALL binterp(Asmpsw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab,   &
                             & casm)
                  CALL binterp(Pmom2psw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpmom2)
                  CALL binterp(Pmom3psw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpmom3)
                  CALL binterp(Pmom4psw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpmom4)
                  CALL binterp(Pmom5psw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpmom5)
                  CALL binterp(Pmom6psw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpmom6)
                  CALL binterp(Pmom7psw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpmom7)
                  CALL binterp(Sbackpsw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtabsw(1,ns),Refitabsw(1,ns),itab,jtab,ttab,utab, &
                             & cpsback2p)
 
!                 chebyshev polynomials
                  ch(1) = 1.
                  ch(2) = xrad
                  DO nc = 3 , NCOEF
                     ch(nc) = 2.*xrad*ch(nc-1) - ch(nc-2)
                  ENDDO
!                 parameterized optical properties
 
                  pext = 0.5*cext(1)
                  DO nc = 2 , NCOEF
                     pext = pext + ch(nc)*cext(nc)
                  ENDDO
                  pext = exp(pext)
 
! JCB 2004/02/09 -- for scattering efficiency
                  pscat = 0.5*cscat(1)
                  DO nc = 2 , NCOEF
                     pscat = pscat + ch(nc)*cscat(nc)
                  ENDDO
                  pscat = exp(pscat)
!
                  pasm = 0.5*casm(1)
                  DO nc = 2 , NCOEF
                     pasm = pasm + ch(nc)*casm(nc)
                  ENDDO
                  pasm = exp(pasm)
!
                  ppmom2 = 0.5*cpmom2(1)
                  DO nc = 2 , NCOEF
                     ppmom2 = ppmom2 + ch(nc)*cpmom2(nc)
                  ENDDO
                  IF ( ppmom2<=0.0 ) ppmom2 = 0.0
!
                  ppmom3 = 0.5*cpmom3(1)
                  DO nc = 2 , NCOEF
                     ppmom3 = ppmom3 + ch(nc)*cpmom3(nc)
                  ENDDO
                  IF ( ppmom3<=0.0 ) ppmom3 = 0.0
!
                  ppmom4 = 0.5*cpmom4(1)
                  DO nc = 2 , NCOEF
                     ppmom4 = ppmom4 + ch(nc)*cpmom4(nc)
                  ENDDO
                  IF ( ppmom4<=0.0 .OR. sizem<=0.03E-04 ) ppmom4 = 0.0
!
                  ppmom5 = 0.5*cpmom5(1)
                  DO nc = 2 , NCOEF
                     ppmom5 = ppmom5 + ch(nc)*cpmom5(nc)
                  ENDDO
                  IF ( ppmom5<=0.0 .OR. sizem<=0.03E-04 ) ppmom5 = 0.0
!
                  ppmom6 = 0.5*cpmom6(1)
                  DO nc = 2 , NCOEF
                     ppmom6 = ppmom6 + ch(nc)*cpmom6(nc)
                  ENDDO
                  IF ( ppmom6<=0.0 .OR. sizem<=0.03E-04 ) ppmom6 = 0.0
!
                  ppmom7 = 0.5*cpmom7(1)
                  DO nc = 2 , NCOEF
                     ppmom7 = ppmom7 + ch(nc)*cpmom7(nc)
                  ENDDO
                  IF ( ppmom7<=0.0 .OR. sizem<=0.03E-04 ) ppmom7 = 0.0
!
                  sback2 = 0.5*cpsback2p(1)
                                          ! JCB 2007/02/01 - backscattering efficiency
                  DO nc = 2 , NCOEF
                     sback2 = sback2 + ch(nc)*cpsback2p(nc)
                  ENDDO
                  sback2 = exp(sback2)
                  IF ( sback2<=0.0 ) sback2 = 0.0
!
!
! weights:
                  pscat = min(pscat,pext)
                                !czhao
                  weighte = pext*pie*exp(x)**2
                                   ! JCB, extinction cross section
                  weights = pscat*pie*exp(x)**2
                                    ! JCB, scattering cross section
                  Swtauaer(ns,klevel) = Swtauaer(ns,klevel) + weighte*Number_bin_col(m,klevel)
                                                                                  ! must be multiplied by deltaZ
!      if (iclm==30.and.jclm==49.and.klevel==2.and.m==5) then
!      write(0,*) 'czhao check swtauaer calculation in MIE',ns,m,weighte,number_bin_col(m,klevel),swtauaer(ns,klevel)*dz(klevel)*100
!      print*, 'czhao check swtauaer calculation in MIE',ns,m,weighte,number_bin_col(m,klevel),swtauaer(ns,klevel)*dz(klevel)*100
!      endif
                  Swsizeaer(ns,klevel) = Swsizeaer(ns,klevel) + exp(x)*10000.0*Number_bin_col(m,klevel)
                  Swwaer(ns,klevel) = Swwaer(ns,klevel) + weights*Number_bin_col(m,klevel)
                                                                             !JCB
                  Swgaer(ns,klevel) = Swgaer(ns,klevel) + pasm*weights*Number_bin_col(m,klevel)
                                                                                  !JCB
! need weighting by scattering cross section ?  JCB 2004/02/09
                  L2(ns,klevel) = L2(ns,klevel) + weights*ppmom2*Number_bin_col(m,klevel)
                  L3(ns,klevel) = L3(ns,klevel) + weights*ppmom3*Number_bin_col(m,klevel)
                  L4(ns,klevel) = L4(ns,klevel) + weights*ppmom4*Number_bin_col(m,klevel)
                  L5(ns,klevel) = L5(ns,klevel) + weights*ppmom5*Number_bin_col(m,klevel)
                  L6(ns,klevel) = L6(ns,klevel) + weights*ppmom6*Number_bin_col(m,klevel)
                  L7(ns,klevel) = L7(ns,klevel) + weights*ppmom7*Number_bin_col(m,klevel)
! convert backscattering efficiency to backscattering coefficient, units (cm)^-1
                  Swbscoef(ns,klevel) = Swbscoef(ns,klevel) + pie*exp(x)**2*sback2*Number_bin_col(m,klevel)
                                                                                             ! backscatter
 
               ENDDO
               ! end of nbin_a loop
 
! take averages - weighted by cross section - new code JCB 2004/02/09
               Swsizeaer(ns,klevel) = Swsizeaer(ns,klevel)/thesum
               Swgaer(ns,klevel) = Swgaer(ns,klevel)/Swwaer(ns,klevel)
                                                              ! JCB removed *3 factor 2/9/2004
! because factor is applied in subroutine opmie, file zz01fastj_mod.f
               L2(ns,klevel) = L2(ns,klevel)/Swwaer(ns,klevel)
               L3(ns,klevel) = L3(ns,klevel)/Swwaer(ns,klevel)
               L4(ns,klevel) = L4(ns,klevel)/Swwaer(ns,klevel)
               L5(ns,klevel) = L5(ns,klevel)/Swwaer(ns,klevel)
               L6(ns,klevel) = L6(ns,klevel)/Swwaer(ns,klevel)
               L7(ns,klevel) = L7(ns,klevel)/Swwaer(ns,klevel)
! backscatter coef, divide by 4*Pie to get units of (km*ster)^-1 JCB 2007/02/01
               Swbscoef(ns,klevel) = Swbscoef(ns,klevel)*1.0E5
                                                       ! units are now (km)^-1
               Swextaer(ns,klevel) = Swtauaer(ns,klevel)*1.0E5
                                                       ! now true extincion, units (km)^-1
! this must be last!!
               Swwaer(ns,klevel) = Swwaer(ns,klevel)/Swtauaer(ns,klevel)
                                                                ! JCB
            ENDIF
 
!70   continue ! bail out if no aerosol;go on to next wavelength bin
 
         ENDDO   ! end of wavelength loop
 
      ENDDO      ! end of klevel loop
!
! before returning, multiply tauaer by depth of individual cells.
! tauaer is in cm-1, dz in m; multiply dz by 100 to convert from m to cm.
      DO ns = 1 , NSWBANDS
         DO klevel = 1 , Kte
            Swtauaer(ns,klevel) = Swtauaer(ns,klevel)*Dz(klevel)*100.
         ENDDO
      ENDDO
 
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ec  fastj diagnostics
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Iclm==chem_dbg_i ) THEN
         IF ( Jclm==chem_dbg_j ) THEN
!   initial entry
            IF ( kcallmieaer==0 ) THEN
               WRITE (*,99001) chem_dbg_i , chem_dbg_j
99001          FORMAT (' for cell i = ',i3,' j = ',i3)
               WRITE (*,99002)
99002          FORMAT ('curr_secs',3x,'i',3x,'j',3x,'k',3x,'dzmfastj',8x,'tauaer(1,k)',1x,'tauaer(2,k)',1x,'tauaer(3,k)',3x,       &
                      &'tauaer(4,k)',5x,'waer(1,k)',7x,'waer(2,k)',7x,'waer(3,k)',7x,'waer(4,k)',7x,'gaer(1,k)',7x,'gaer(2,k)',7x, &
                      &'gaer(3,k)',7x,'gaer(4,k)',7x,'extaer(1,k)',5x,'extaer(2,k)',5x,'extaer(3,k)',5x,'extaer(4,k)',5x,          &
                      &'sizeaer(1,k)',4x,'sizeaer(2,k)',4x,'sizeaer(3,k)',4x,'sizeaer(4,k)')
            ENDIF
!ec output for run_out.30
            DO k = 1 , Kte
               WRITE (*,99003) Curr_secs , Iclm , Jclm , k , Dz(k) , (Swtauaer(n,k),n=1,4) , (Swwaer(n,k),n=1,4) ,                 &
                             & (Swgaer(n,k),n=1,4) , (Swextaer(n,k),n=1,4) , (Swsizeaer(n,k),n=1,4)
99003          FORMAT (i7,3(2x,i4),2x,21(e14.6,2x))
            ENDDO
            kcallmieaer = kcallmieaer + 1
         ENDIF
      ENDIF
!ec end print of fastj diagnostics
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif
 
 
!######################################################################
!parameterization of mie calculation for longwave
!#####################################################################
 
! begin level loop
      DO klevel = 1 , Kte
! sum densities for normalization
         thesum = 0.0
         DO m = 1 , Nbin_a
            thesum = thesum + Number_bin_col(m,klevel)
         ENDDO
! Begin longwave spectral loop
         DO ns = 1 , NLWBANDS
 
!        aerosol optical properties
            Lwtauaer(ns,klevel) = 0.
            Lwextaer(ns,klevel) = 0.0
            IF ( thesum>1E-21 ) THEN      ! set everything = 0 if no aerosol !wig changed 0.0 to 1e-21, 31-Oct-2005
 
! loop over the bins
               DO m = 1 , Nbin_a
                             ! nbin_a is number of bins
! here's the size
                  sizem = Radius_wet_col(m,klevel)
                                               ! radius in cm
                  x = log(Radius_wet_col(m,klevel))
                                                ! radius in cm
                  crefin = Lwrefindx_col(m,klevel,ns)
                  refr = real(crefin)
                  refi = -imag(crefin)
                  xrad = x
                  thesize = 2.0*pie*exp(x)/Wavmidlw(ns)
                ! normalize size parameter
                  xrad = (2*xrad-xrmax-xrmin)/(xrmax-xrmin)
 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !ec  diagnostics
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ! retain this diagnostic code
                  IF ( abs(refr)>10.0 .OR. abs(refr)<=0.001 ) THEN
                     WRITE (msg,'(a,1x, e14.5)') 'mieaer /refr/ outside range 1e-3 - 10 '//'refr= ' , refr
                     CALL peg_error_fatal(LUNERR,msg)
                  ENDIF
                  IF ( abs(refi)>10. ) THEN
                     WRITE (msg,'(a,1x, e14.5)') 'mieaer /refi/ >10 '//'refi' , refi
                     CALL peg_error_fatal(LUNERR,msg)
                  ENDIF
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
! interpolate coefficients linear in refractive index
! first call calcs itab,jtab,ttab,utab
                  itab = 0
                  CALL binterp(Absplw(1,1,1,ns),NCOEF,nrefr,nrefi,refr,refi,Refrtablw(1,ns),Refitablw(1,ns),itab,jtab,ttab,utab,   &
                             & cabs)
 
!                 chebyshev polynomials
                  ch(1) = 1.
                  ch(2) = xrad
                  DO nc = 3 , NCOEF
                     ch(nc) = 2.*xrad*ch(nc-1) - ch(nc-2)
                  ENDDO
!                 parameterized optical properties
                  pabs = 0.5*cabs(1)
                  DO nc = 2 , NCOEF
                     pabs = pabs + ch(nc)*cabs(nc)
                  ENDDO
                  pabs = exp(pabs)
 
!
! weights:
                  weighta = pabs*pie*exp(x)**2
                                   ! JCB, extinction cross section
        !weighta: cm2 and number_bin_col #/cm3 -> /cm
                  Lwtauaer(ns,klevel) = Lwtauaer(ns,klevel) + weighta*Number_bin_col(m,klevel)
                                                                                 ! must be multiplied by deltaZ
 
               ENDDO
               ! end of nbin_a loop
 
! take averages - weighted by cross section - new code JCB 2004/02/09
               Lwextaer(ns,klevel) = Lwtauaer(ns,klevel)*1.0E5
                                                       ! now true extincion, units (km)^-1
            ENDIF
 
         ENDDO   ! end of wavelength loop
 
      ENDDO      ! end of klevel loop
!
! before returning, multiply tauaer by depth of individual cells.
! tauaer is in cm-1, dz in m; multiply dz by 100 to convert from m to cm.
      DO ns = 1 , NLWBANDS
         DO klevel = 1 , Kte
            Lwtauaer(ns,klevel) = Lwtauaer(ns,klevel)*Dz(klevel)*100.
         ENDDO
      ENDDO
 
   END SUBROUTINE mieaer
!****************************************************************
 
!****************************************************************
 
   SUBROUTINE fitcurv(Rs,Yin,Coef,Ncoef,Maxm)
 
!     fit y(x) using Chebychev polynomials
!     wig 7-Sep-2004: Removed dependency on pre-determined maximum
!                     array size and replaced with f90 array info.
 
 
      USE module_peg_util , only:peg_message
!      integer nmodes, nrows, maxm, ncoef
!      parameter (nmodes=500,nrows=8)
      INTEGER , INTENT(IN) :: Maxm , Ncoef
 
!      real rs(nmodes),yin(nmodes),coef(ncoef)
!      real x(nmodes),y(nmodes)
      REAL , DIMENSION(Ncoef) :: Coef
      REAL , DIMENSION(:) :: Rs , Yin
      REAL x(size(Rs)) , y(size(Yin))
 
      INTEGER m
      REAL xmin , xmax
      CHARACTER*80 msg
 
!!$      if(maxm.gt.nmodes)then
!!$        write ( msg, '(a, 1x,i6)' )  &
!!$           'FASTJ mie nmodes too small in fitcurv, '  //  &
!!$           'maxm ', maxm
!!$!        write(*,*)'nmodes too small in fitcurv',maxm
!!$        call peg_error_fatal( lunerr, msg )
!!$      endif
 
      DO m = 1 , Maxm
! To prevent the log of 0 or negative values, as the code was blowing up when compile with intel
! Added by Manish Shrivastava
! Need to be checked
         x(m) = log(max(Rs(m),1D-20))
         y(m) = log(max(Yin(m),1D-20))
      ENDDO
 
      xmin = x(1)
      xmax = x(Maxm)
      DO m = 1 , Maxm
         x(m) = (2*x(m)-xmax-xmin)/(xmax-xmin)
      ENDDO
 
      CALL chebft(Coef,Ncoef,Maxm,y)
 
   END SUBROUTINE fitcurv
!**************************************************************
   SUBROUTINE fitcurv_nolog(Rs,Yin,Coef,Ncoef,Maxm)
 
!     fit y(x) using Chebychev polynomials
!     wig 7-Sep-2004: Removed dependency on pre-determined maximum
!                     array size and replaced with f90 array info.
 
      USE module_peg_util , only:peg_message
 
!      integer nmodes, nrows, maxm, ncoef
!      parameter (nmodes=500,nrows=8)
      INTEGER , INTENT(IN) :: Maxm , Ncoef
 
!      real rs(nmodes),yin(nmodes),coef(ncoef)
      REAL , DIMENSION(:) :: Rs , Yin
      REAL , DIMENSION(Ncoef) :: Coef(Ncoef)
      REAL x(size(Rs)) , y(size(Yin))
 
      INTEGER m
      REAL xmin , xmax
      CHARACTER*80 msg
 
!!$      if(maxm.gt.nmodes)then
!!$        write ( msg, '(a,1x, i6)' )  &
!!$           'FASTJ mie nmodes too small in fitcurv '  //  &
!!$           'maxm ', maxm
!!$!        write(*,*)'nmodes too small in fitcurv',maxm
!!$        call peg_error_fatal( lunerr, msg )
!!$      endif
 
      DO m = 1 , Maxm
         x(m) = log(Rs(m))
         y(m) = Yin(m)
                  ! note, no "log" here
      ENDDO
 
      xmin = x(1)
      xmax = x(Maxm)
      DO m = 1 , Maxm
         x(m) = (2*x(m)-xmax-xmin)/(xmax-xmin)
      ENDDO
 
      CALL chebft(Coef,Ncoef,Maxm,y)
 
   END SUBROUTINE fitcurv_nolog
!************************************************************************
   SUBROUTINE chebft(C,Ncoef,N,F)
!     given a function f with values at zeroes x_k of Chebychef polynomial
!     T_n(x), calculate coefficients c_j such that
!     f(x)=sum(k=1,n) c_k t_(k-1)(y) - 0.5*c_1
!     where y=(x-0.5*(xmax+xmin))/(0.5*(xmax-xmin))
!     See Numerical Recipes, pp. 148-150.
 
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER :: PI = 3.14159265
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Ncoef
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(OUT) , DIMENSION(Ncoef) :: C
      REAL , INTENT(IN) , DIMENSION(N) :: F
!
! Local variable declarations rewritten by SPAG
!
      REAL :: fac , thesum
      INTEGER :: j , k
!
! End of declarations rewritten by SPAG
!
 
! local variables
 
      fac = 2./N
      DO j = 1 , Ncoef
         thesum = 0
         DO k = 1 , N
            thesum = thesum + F(k)*cos((PI*(j-1))*((k-0.5)/N))
         ENDDO
         C(j) = fac*thesum
      ENDDO
   END SUBROUTINE chebft
!*************************************************************************
   SUBROUTINE binterp(Table,Km,Im,Jm,X,Y,Xtab,Ytab,Ix,Jy,T,U,Out)
 
!     bilinear interpolation of table
!
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Km
      INTEGER , INTENT(IN) :: Im
      INTEGER , INTENT(IN) :: Jm
      REAL , INTENT(IN) , DIMENSION(Km,Im,Jm) :: Table
      REAL , INTENT(IN) :: X
      REAL , INTENT(IN) :: Y
      REAL , INTENT(IN) , DIMENSION(Im) :: Xtab
      REAL , INTENT(IN) , DIMENSION(Jm) :: Ytab
      INTEGER , INTENT(INOUT) :: Ix
      INTEGER , INTENT(INOUT) :: Jy
      REAL , INTENT(INOUT) :: T
      REAL , INTENT(INOUT) :: U
      REAL , INTENT(OUT) , DIMENSION(Km) :: Out
!
! Local variable declarations rewritten by SPAG
!
      REAL :: dx , dy , tcu , tcuc , tu , tuc
      INTEGER :: i , ip1 , j , jp1 , k
!
! End of declarations rewritten by SPAG
!
 
      IF ( Ix<=0 ) THEN
         IF ( Im>1 ) THEN
            SPAG_Loop_1_1: DO i = 1 , Im
               IF ( X<Xtab(i) ) EXIT SPAG_Loop_1_1
            ENDDO SPAG_Loop_1_1
            Ix = max0(i-1,1)
            ip1 = min0(Ix+1,Im)
            dx = (Xtab(ip1)-Xtab(Ix))
            IF ( abs(dx)>1.E-20 ) THEN
               T = (X-Xtab(Ix))/(Xtab(Ix+1)-Xtab(Ix))
            ELSE
               T = 0
            ENDIF
         ELSE
            Ix = 1
            ip1 = 1
            T = 0
         ENDIF
         IF ( Jm>1 ) THEN
            SPAG_Loop_1_2: DO j = 1 , Jm
               IF ( Y<Ytab(j) ) EXIT SPAG_Loop_1_2
            ENDDO SPAG_Loop_1_2
            Jy = max0(j-1,1)
            jp1 = min0(Jy+1,Jm)
            dy = (Ytab(jp1)-Ytab(Jy))
            IF ( abs(dy)>1.E-20 ) THEN
               U = (Y-Ytab(Jy))/dy
            ELSE
               U = 0
            ENDIF
         ELSE
            Jy = 1
            jp1 = 1
            U = 0
         ENDIF
      ENDIF
      jp1 = min(Jy+1,Jm)
      ip1 = min(Ix+1,Im)
      tu = T*U
      tuc = T - tu
      tcuc = 1 - tuc - U
      tcu = U - tu
      DO k = 1 , Km
         Out(k) = tcuc*Table(k,Ix,Jy) + tuc*Table(k,ip1,Jy) + tu*Table(k,ip1,jp1) + tcu*Table(k,Ix,jp1)
      ENDDO
   END SUBROUTINE binterp
!***************************************************************
   SUBROUTINE miev0(Xx,Crefin,Perfct,Mimcut,Anyang,Numang,Xmu,Nmom,Ipolzn,Momdim,Prnt,Qext,Qsca,Gqsc,Pmom,Sforw,Sback,S1,S2,Tforw, &
                  & Tback)
!
!    computes mie scattering and extinction efficiencies; asymmetry
!    factor;  forward- and backscatter amplitude;  scattering
!    amplitudes for incident polarization parallel and perpendicular
!    to the plane of scattering, as functions of scattering angle;
!    coefficients in the legendre polynomial expansions of either the
!    unpolarized phase function or the polarized phase matrix;
!    and some quantities needed in polarized radiative transfer.
!
!      calls :  biga, ckinmi, small1, small2, testmi, miprnt,
!               lpcoef, errmsg
!
!      i n t e r n a l   v a r i a b l e s
!      -----------------------------------
!
!  an,bn           mie coefficients  little-a-sub-n, little-b-sub-n
!                     ( ref. 1, eq. 16 )
!  anm1,bnm1       mie coefficients  little-a-sub-(n-1),
!                     little-b-sub-(n-1);  used in -gqsc- sum
!  anp             coeffs. in s+ expansion ( ref. 2, p. 1507 )
!  bnp             coeffs. in s- expansion ( ref. 2, p. 1507 )
!  anpm            coeffs. in s+ expansion ( ref. 2, p. 1507 )
!                     when  mu  is replaced by  - mu
!  bnpm            coeffs. in s- expansion ( ref. 2, p. 1507 )
!                     when  mu  is replaced by  - mu
!  calcmo(k)       true, calculate moments for k-th phase quantity
!                     (derived from -ipolzn-; used only in 'lpcoef')
!  cbiga(n)        bessel function ratio capital-a-sub-n (ref. 2, eq. 2)
!                     ( complex version )
!  cior            complex index of refraction with negative
!                     imaginary part (van de hulst convention)
!  cioriv          1 / cior
!  coeff           ( 2n + 1 ) / ( n ( n + 1 ) )
!  fn              floating point version of index in loop performing
!                     mie series summation
!  lita,litb(n)    mie coefficients -an-, -bn-, saved in arrays for
!                     use in calculating legendre moments *pmom*
!  maxtrm          max. possible no. of terms in mie series
!  mm              + 1 and  - 1,  alternately.
!  mim             magnitude of imaginary refractive index
!  mre             real part of refractive index
!  maxang          max. possible value of input variable -numang-
!  nangd2          (numang+1)/2 ( no. of angles in 0-90 deg; anyang=f )
!  noabs           true, sphere non-absorbing (determined by -mimcut-)
!  np1dn           ( n + 1 ) / n
!  npquan          highest-numbered phase quantity for which moments are
!                     to be calculated (the largest digit in -ipolzn-
!                     if  ipolzn .ne. 0)
!  ntrm            no. of terms in mie series
!  pass1           true on first entry, false thereafter; for self-test
!  pin(j)          angular function little-pi-sub-n ( ref. 2, eq. 3 )
!                     at j-th angle
!  pinm1(j)        little-pi-sub-(n-1) ( see -pin- ) at j-th angle
!  psinm1          ricatti-bessel function psi-sub-(n-1), argument -xx-
!  psin            ricatti-bessel function psi-sub-n of argument -xx-
!                     ( ref. 1, p. 11 ff. )
!  rbiga(n)        bessel function ratio capital-a-sub-n (ref. 2, eq. 2)
!                     ( real version, for when imag refrac index = 0 )
!  rioriv          1 / mre
!  rn              1 / n
!  rtmp            (real) temporary variable
!  sp(j)           s+  for j-th angle  ( ref. 2, p. 1507 )
!  sm(j)           s-  for j-th angle  ( ref. 2, p. 1507 )
!  sps(j)          s+  for (numang+1-j)-th angle ( anyang=false )
!  sms(j)          s-  for (numang+1-j)-th angle ( anyang=false )
!  taun            angular function little-tau-sub-n ( ref. 2, eq. 4 )
!                     at j-th angle
!  tcoef           n ( n+1 ) ( 2n+1 ) (for summing tforw,tback series)
!  twonp1          2n + 1
!  yesang          true if scattering amplitudes are to be calculated
!  zetnm1          ricatti-bessel function  zeta-sub-(n-1) of argument
!                     -xx-  ( ref. 2, eq. 17 )
!  zetn            ricatti-bessel function  zeta-sub-n of argument -xx-
!
! ----------------------------------------------------------------------
! --------  i / o specifications for subroutine miev0  -----------------
! ----------------------------------------------------------------------
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER :: MAXANG = 501 , MXANG2 = MAXANG/2 + 1 , MAXTRM = 1100
      REAL*8 , PARAMETER :: ONETHR = 1./3.
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Momdim
      REAL*8 , INTENT(INOUT) :: Xx
      COMPLEX*16 , INTENT(INOUT) :: Crefin
      LOGICAL , INTENT(INOUT) :: Perfct
      REAL*8 , INTENT(INOUT) :: Mimcut
      LOGICAL , INTENT(INOUT) :: Anyang
      INTEGER , INTENT(INOUT) :: Numang
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xmu
      INTEGER , INTENT(INOUT) :: Nmom
      INTEGER , INTENT(INOUT) :: Ipolzn
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Prnt
      REAL*8 , INTENT(INOUT) :: Qext
      REAL*8 , INTENT(INOUT) :: Qsca
      REAL*8 , INTENT(INOUT) :: Gqsc
      REAL*8 , DIMENSION(0:Momdim,*) :: Pmom
      COMPLEX*16 , INTENT(INOUT) :: Sforw
      COMPLEX*16 , INTENT(INOUT) :: Sback
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(*) :: S1
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(*) :: S2
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(*) :: Tforw
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(*) :: Tback
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: an , anm1 , anp , anpm , bn , bnm1 , bnp , bnpm , cior , cioriv , cresav , ctmp , zet , zetn , zetnm1
      LOGICAL :: anysav , noabs , ok , persav , yesang
      LOGICAL , DIMENSION(4) :: calcmo
      COMPLEX*16 , DIMENSION(MAXTRM) :: cbiga , lita , litb
      REAL*8 :: chin , chinm1 , coeff , fn , mim , mimsav , mm , mre , np1dn , psin , psinm1 , rioriv , rn , rtmp , taun , tcoef , &
              & twonp1 , xinv , xmusav , xxsav
      INTEGER :: i , iposav , j , n , nangd2 , nmosav , npquan , ntrm , numsav
      LOGICAL , SAVE :: pass1
      REAL*8 , DIMENSION(MAXANG) :: pin , pinm1
      REAL*8 , DIMENSION(MAXTRM) :: rbiga
      COMPLEX*16 , DIMENSION(MAXANG) :: sm , sp
      COMPLEX*16 , DIMENSION(MXANG2) :: sms , sps
      REAL*8 :: sq
!
! End of declarations rewritten by SPAG
!
! ----------------------------------------------------------------------
!
!
!                                  ** note --  maxtrm = 10100  is neces-
!                                  ** sary to do some of the test probs,
!                                  ** but 1100 is sufficient for most
!                                  ** conceivable applications
!
      EQUIVALENCE (cbiga,rbiga)
      DATA pass1/.TRUE./
      sq(ctmp) = dble(ctmp)**2 + dimag(ctmp)**2
      INTEGER :: spag_nextblock_1
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
!
!
            IF ( pass1 ) THEN
!                                   ** save certain user input values
               xxsav = Xx
               cresav = Crefin
               mimsav = Mimcut
               persav = Perfct
               anysav = Anyang
               nmosav = Nmom
               iposav = Ipolzn
               numsav = Numang
               xmusav = Xmu(1)
!                              ** reset input values for test case
               Xx = 10.0
               Crefin = (1.5,-0.1)
               Perfct = .FALSE.
               Mimcut = 0.0
               Anyang = .TRUE.
               Numang = 1
               Xmu(1) = -0.7660444
               Nmom = 1
               Ipolzn = -1
!
            ENDIF
            spag_nextblock_1 = 2
         CASE (2)
!                                        ** check input and calculate
!                                        ** certain variables from input
!
            CALL ckinmi(Numang,MAXANG,Xx,Perfct,Crefin,Momdim,Nmom,Ipolzn,Anyang,Xmu,calcmo,npquan)
!
            IF ( Perfct .AND. Xx<=0.1 ) THEN
!                                            ** use totally-reflecting
!                                            ** small-particle limit
!
               CALL small1(Xx,Numang,Xmu,Qext,Qsca,Gqsc,Sforw,Sback,S1,S2,Tforw,Tback,lita,litb)
               ntrm = 2
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!
            IF ( .NOT.Perfct ) THEN
!
               cior = Crefin
               IF ( dimag(cior)>0.0 ) cior = dconjg(cior)
               mre = dble(cior)
               mim = -dimag(cior)
               noabs = mim<=Mimcut
               cioriv = 1.0/cior
               rioriv = 1.0/mre
!
               IF ( Xx*dmax1(1.D0,cdabs(cior))<=0.D1 ) THEN
!
!                                    ** use general-refractive-index
!                                    ** small-particle limit
!                                    ** ( ref. 2, p. 1508 )
!
                  CALL small2(Xx,cior,.NOT.noabs,Numang,Xmu,Qext,Qsca,Gqsc,Sforw,Sback,S1,S2,Tforw,Tback,lita,litb)
                  ntrm = 2
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!
            ENDIF
!
            nangd2 = (Numang+1)/2
            yesang = Numang>0
!                              ** estimate number of terms in mie series
!                              ** ( ref. 2, p. 1508 )
            IF ( Xx<=8.0 ) THEN
               ntrm = Xx + 4.*Xx**ONETHR + 1.
            ELSEIF ( Xx<4200. ) THEN
               ntrm = Xx + 4.05*Xx**ONETHR + 2.
            ELSE
               ntrm = Xx + 4.*Xx**ONETHR + 2.
            ENDIF
            IF ( ntrm+1>MAXTRM ) CALL errmsg('miev0--parameter maxtrm too small',.TRUE.)
!
!                            ** calculate logarithmic derivatives of
!                            ** j-bessel-fcn., big-a-sub-(1 to ntrm)
            IF ( .NOT.Perfct ) CALL biga(cior,Xx,ntrm,noabs,yesang,rbiga,cbiga)
!
!                            ** initialize ricatti-bessel functions
!                            ** (psi,chi,zeta)-sub-(0,1) for upward
!                            ** recurrence ( ref. 1, eq. 19 )
            xinv = 1.0/Xx
            psinm1 = dsin(Xx)
            chinm1 = dcos(Xx)
            psin = psinm1*xinv - chinm1
            chin = chinm1*xinv + psinm1
            zetnm1 = dcmplx(psinm1,chinm1)
            zetn = dcmplx(psin,chin)
!                                     ** initialize previous coeffi-
!                                     ** cients for -gqsc- series
            anm1 = (0.0,0.0)
            bnm1 = (0.0,0.0)
!                             ** initialize angular function little-pi
!                             ** and sums for s+, s- ( ref. 2, p. 1507 )
            IF ( Anyang ) THEN
               DO j = 1 , Numang
                  pinm1(j) = 0.0
                  pin(j) = 1.0
                  sp(j) = (0.0,0.0)
                  sm(j) = (0.0,0.0)
               ENDDO
            ELSE
               DO j = 1 , nangd2
                  pinm1(j) = 0.0
                  pin(j) = 1.0
                  sp(j) = (0.0,0.0)
                  sm(j) = (0.0,0.0)
                  sps(j) = (0.0,0.0)
                  sms(j) = (0.0,0.0)
               ENDDO
            ENDIF
!                         ** initialize mie sums for efficiencies, etc.
            Qsca = 0.0
            Gqsc = 0.0
            Sforw = (0.,0.)
            Sback = (0.,0.)
            Tforw(1) = (0.,0.)
            Tback(1) = (0.,0.)
!
!
! ---------  loop to sum mie series  -----------------------------------
!
            mm = +1.0
            DO n = 1 , ntrm
!                           ** compute various numerical coefficients
               fn = n
               rn = 1.0/fn
               np1dn = 1.0 + rn
               twonp1 = 2*n + 1
               coeff = twonp1/(fn*(n+1))
               tcoef = twonp1*(fn*(n+1))
!
!                              ** calculate mie series coefficients
               IF ( Perfct ) THEN
!                                   ** totally-reflecting case
!
                  an = ((fn*xinv)*psin-psinm1)/((fn*xinv)*zetn-zetnm1)
                  bn = psin/zetn
!
               ELSEIF ( noabs ) THEN
!                                      ** no-absorption case
!
                  an = ((rioriv*rbiga(n)+(fn*xinv))*psin-psinm1)/((rioriv*rbiga(n)+(fn*xinv))*zetn-zetnm1)
                  bn = ((mre*rbiga(n)+(fn*xinv))*psin-psinm1)/((mre*rbiga(n)+(fn*xinv))*zetn-zetnm1)
               ELSE
!                                       ** absorptive case
!
                  an = ((cioriv*cbiga(n)+(fn*xinv))*psin-psinm1)/((cioriv*cbiga(n)+(fn*xinv))*zetn-zetnm1)
                  bn = ((cior*cbiga(n)+(fn*xinv))*psin-psinm1)/((cior*cbiga(n)+(fn*xinv))*zetn-zetnm1)
                  Qsca = Qsca + twonp1*(sq(an)+sq(bn))
!
               ENDIF
!                       ** save mie coefficients for *pmom* calculation
               lita(n) = an
               litb(n) = bn
!                            ** increment mie sums for non-angle-
!                            ** dependent quantities
!
               Sforw = Sforw + twonp1*(an+bn)
               Tforw(1) = Tforw(1) + tcoef*(an-bn)
               Sback = Sback + (mm*twonp1)*(an-bn)
               Tback(1) = Tback(1) + (mm*tcoef)*(an+bn)
               Gqsc = Gqsc + (fn-rn)*dble(anm1*dconjg(an)+bnm1*dconjg(bn)) + coeff*dble(an*dconjg(bn))
!
               IF ( yesang ) THEN
!                                      ** put mie coefficients in form
!                                      ** needed for computing s+, s-
!                                      ** ( ref. 2, p. 1507 )
                  anp = coeff*(an+bn)
                  bnp = coeff*(an-bn)
!                                      ** increment mie sums for s+, s-
!                                      ** while upward recursing
!                                      ** angular functions little pi
!                                      ** and little tau
                  IF ( Anyang ) THEN
!                                         ** arbitrary angles
!
!                                              ** vectorizable loop
                     DO j = 1 , Numang
                        rtmp = (Xmu(j)*pin(j)) - pinm1(j)
                        taun = fn*rtmp - pinm1(j)
                        sp(j) = sp(j) + anp*(pin(j)+taun)
                        sm(j) = sm(j) + bnp*(pin(j)-taun)
                        pinm1(j) = pin(j)
                        pin(j) = (Xmu(j)*pin(j)) + np1dn*rtmp
                     ENDDO
!
                  ELSE
!                                  ** angles symmetric about 90 degrees
                     anpm = mm*anp
                     bnpm = mm*bnp
!                                          ** vectorizable loop
                     DO j = 1 , nangd2
                        rtmp = (Xmu(j)*pin(j)) - pinm1(j)
                        taun = fn*rtmp - pinm1(j)
                        sp(j) = sp(j) + anp*(pin(j)+taun)
                        sms(j) = sms(j) + bnpm*(pin(j)+taun)
                        sm(j) = sm(j) + bnp*(pin(j)-taun)
                        sps(j) = sps(j) + anpm*(pin(j)-taun)
                        pinm1(j) = pin(j)
                        pin(j) = (Xmu(j)*pin(j)) + np1dn*rtmp
                     ENDDO
!
                  ENDIF
               ENDIF
!                          ** update relevant quantities for next
!                          ** pass through loop
               mm = -mm
               anm1 = an
               bnm1 = bn
!                           ** upward recurrence for ricatti-bessel
!                           ** functions ( ref. 1, eq. 17 )
!
               zet = (twonp1*xinv)*zetn - zetnm1
               zetnm1 = zetn
               zetn = zet
               psinm1 = psin
               psin = dble(zetn)
            ENDDO
!
! ---------- end loop to sum mie series --------------------------------
!
!
            Qext = 2./Xx**2*dble(Sforw)
            IF ( Perfct .OR. noabs ) THEN
               Qsca = Qext
            ELSE
               Qsca = 2./Xx**2*Qsca
            ENDIF
!
            Gqsc = 4./Xx**2*Gqsc
            Sforw = 0.5*Sforw
            Sback = 0.5*Sback
            Tforw(2) = 0.5*(Sforw+0.25*Tforw(1))
            Tforw(1) = 0.5*(Sforw-0.25*Tforw(1))
            Tback(2) = 0.5*(Sback+0.25*Tback(1))
            Tback(1) = 0.5*(-Sback+0.25*Tback(1))
!
            IF ( yesang ) THEN
!                                ** recover scattering amplitudes
!                                ** from s+, s- ( ref. 1, eq. 11 )
               IF ( Anyang ) THEN
!                                         ** vectorizable loop
                  DO j = 1 , Numang
                     S1(j) = 0.5*(sp(j)+sm(j))
                     S2(j) = 0.5*(sp(j)-sm(j))
                  ENDDO
!
               ELSE
!                                         ** vectorizable loop
                  DO j = 1 , nangd2
                     S1(j) = 0.5*(sp(j)+sm(j))
                     S2(j) = 0.5*(sp(j)-sm(j))
                  ENDDO
!                                         ** vectorizable loop
                  DO j = 1 , nangd2
                     S1(Numang+1-j) = 0.5*(sps(j)+sms(j))
                     S2(Numang+1-j) = 0.5*(sps(j)-sms(j))
                  ENDDO
               ENDIF
!
            ENDIF
            spag_nextblock_1 = 3
         CASE (3)
!                                         ** calculate legendre moments
            IF ( Nmom>0 ) CALL lpcoef(ntrm,Nmom,Ipolzn,Momdim,calcmo,npquan,lita,litb,Pmom)
!
            IF ( dimag(Crefin)>0.0 ) THEN
!                                         ** take complex conjugates
!                                         ** of scattering amplitudes
               Sforw = dconjg(Sforw)
               Sback = dconjg(Sback)
               DO i = 1 , 2
                  Tforw(i) = dconjg(Tforw(i))
                  Tback(i) = dconjg(Tback(i))
               ENDDO
!
               DO j = 1 , Numang
                  S1(j) = dconjg(S1(j))
                  S2(j) = dconjg(S2(j))
               ENDDO
!
            ENDIF
!
            IF ( pass1 ) THEN
!                             ** compare test case results with
!                             ** correct answers and abort if bad
!
               CALL testmi(Qext,Qsca,Gqsc,Sforw,Sback,S1,S2,Tforw,Tback,Pmom,Momdim,ok)
               IF ( .NOT.ok ) THEN
                  Prnt(1) = .FALSE.
                  Prnt(2) = .FALSE.
                  CALL miprnt(Prnt,Xx,Perfct,Crefin,Numang,Xmu,Qext,Qsca,Gqsc,Nmom,Ipolzn,Momdim,calcmo,Pmom,Sforw,Sback,Tforw,    &
                            & Tback,S1,S2)
                  CALL errmsg('miev0 -- self-test failed',.TRUE.)
               ENDIF
!                                       ** restore user input values
               Xx = xxsav
               Crefin = cresav
               Mimcut = mimsav
               Perfct = persav
               Anyang = anysav
               Nmom = nmosav
               Ipolzn = iposav
               Numang = numsav
               Xmu(1) = xmusav
               pass1 = .FALSE.
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
!
            ENDIF
!
            IF ( Prnt(1) .OR. Prnt(2) ) CALL miprnt(Prnt,Xx,Perfct,Crefin,Numang,Xmu,Qext,Qsca,Gqsc,Nmom,Ipolzn,Momdim,calcmo,Pmom,&
               & Sforw,Sback,Tforw,Tback,S1,S2)
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
!
!
   END SUBROUTINE miev0
!****************************************************************************
   SUBROUTINE ckinmi(Numang,Maxang,Xx,Perfct,Crefin,Momdim,Nmom,Ipolzn,Anyang,Xmu,Calcmo,Npquan)
!
!        check for bad input to 'miev0' and calculate -calcmo,npquan-
!
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Numang
      INTEGER , INTENT(IN) :: Maxang
      REAL*8 , INTENT(IN) :: Xx
      LOGICAL , INTENT(IN) :: Perfct
      COMPLEX*16 , INTENT(IN) :: Crefin
      INTEGER , INTENT(IN) :: Momdim
      INTEGER , INTENT(IN) :: Nmom
      INTEGER , INTENT(IN) :: Ipolzn
      LOGICAL , INTENT(IN) :: Anyang
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xmu
      LOGICAL , INTENT(OUT) , DIMENSION(*) :: Calcmo
      INTEGER , INTENT(INOUT) :: Npquan
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ip , j , l
      LOGICAL :: inperr
      CHARACTER(4) :: string
!
! End of declarations rewritten by SPAG
!
!
!
      inperr = .FALSE.
!
      IF ( Numang>Maxang ) THEN
         CALL errmsg('miev0--parameter maxang too small',.TRUE.)
         inperr = .TRUE.
      ENDIF
      IF ( Numang<0 ) CALL wrtbad('numang',inperr)
      IF ( Xx<0. ) CALL wrtbad('xx',inperr)
      IF ( .NOT.Perfct .AND. dble(Crefin)<=0. ) CALL wrtbad('crefin',inperr)
      IF ( Momdim<1 ) CALL wrtbad('momdim',inperr)
!
      IF ( Nmom/=0 ) THEN
         IF ( Nmom<0 .OR. Nmom>Momdim ) CALL wrtbad('nmom',inperr)
         IF ( iabs(Ipolzn)>4444 ) CALL wrtbad('ipolzn',inperr)
         Npquan = 0
         DO l = 1 , 4
            Calcmo(l) = .FALSE.
         ENDDO
         IF ( Ipolzn/=0 ) THEN
!                                 ** parse out -ipolzn- into its digits
!                                 ** to find which phase quantities are
!                                 ** to have their moments calculated
!
            WRITE (string,'(i4)') iabs(Ipolzn)
            DO j = 1 , 4
               ip = ichar(string(j:j)) - ichar('0')
               IF ( ip>=1 .AND. ip<=4 ) Calcmo(ip) = .TRUE.
               IF ( ip==0 .OR. (ip>=5 .AND. ip<=9) ) CALL wrtbad('ipolzn',inperr)
               Npquan = max0(Npquan,ip)
            ENDDO
         ENDIF
      ENDIF
!
      IF ( Anyang ) THEN
!                                ** allow for slight imperfections in
!                                ** computation of cosine
         DO i = 1 , Numang
            IF ( Xmu(i)<-1.00001 .OR. Xmu(i)>1.00001 ) CALL wrtbad('xmu',inperr)
         ENDDO
      ELSE
         DO i = 1 , (Numang+1)/2
            IF ( Xmu(i)<-0.00001 .OR. Xmu(i)>1.00001 ) CALL wrtbad('xmu',inperr)
         ENDDO
      ENDIF
!
      IF ( inperr ) CALL errmsg('miev0--input error(s).  aborting...',.TRUE.)
!
      IF ( Xx>20000.0 .OR. dble(Crefin)>10.0 .OR. dabs(dimag(Crefin))>10.0 )                                                       &
         & CALL errmsg('miev0--xx or crefin outside tested range',.FALSE.)
!
   END SUBROUTINE ckinmi
!***********************************************************************
   SUBROUTINE lpcoef(Ntrm,Nmom,Ipolzn,Momdim,Calcmo,Npquan,A,B,Pmom)
!
!         calculate legendre polynomial expansion coefficients (also
!         called moments) for phase quantities ( ref. 5 formulation )
!
!     input:  ntrm                    number terms in mie series
!             nmom, ipolzn, momdim    'miev0' arguments
!             calcmo                  flags calculated from -ipolzn-
!             npquan                  defined in 'miev0'
!             a, b                    mie series coefficients
!
!     output: pmom                   legendre moments ('miev0' argument)
!
!     *** notes ***
!
!         (1)  eqs. 2-5 are in error in dave, appl. opt. 9,
!         1888 (1970).  eq. 2 refers to m1, not m2;  eq. 3 refers to
!         m2, not m1.  in eqs. 4 and 5, the subscripts on the second
!         term in square brackets should be interchanged.
!
!         (2)  the general-case logic in this subroutine works correctly
!         in the two-term mie series case, but subroutine  'lpco2t'
!         is called instead, for speed.
!
!         (3)  subroutine  'lpco1t', to do the one-term case, is never
!         called within the context of 'miev0', but is included for
!         complete generality.
!
!         (4)  some improvement in speed is obtainable by combining the
!         310- and 410-loops, if moments for both the third and fourth
!         phase quantities are desired, because the third phase quantity
!         is the real part of a complex series, while the fourth phase
!         quantity is the imaginary part of that very same series.  but
!         most users are not interested in the fourth phase quantity,
!         which is related to circular polarization, so the present
!         scheme is usually more efficient.
!
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER :: MAXTRM = 1102 , MAXMOM = 2*MAXTRM , MXMOM2 = MAXMOM/2 , MAXRCP = 4*MAXTRM + 2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Momdim
      INTEGER , INTENT(IN) :: Ntrm
      INTEGER :: Nmom
      INTEGER :: Ipolzn
      LOGICAL , DIMENSION(*) :: Calcmo
      INTEGER , INTENT(IN) :: Npquan
      COMPLEX*16 , DIMENSION(*) :: A
      COMPLEX*16 , DIMENSION(*) :: B
      REAL*8 , INTENT(INOUT) , DIMENSION(0:Momdim,*) :: Pmom
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 , DIMENSION(0:MAXTRM) :: am
      REAL*8 , DIMENSION(0:MXMOM2) :: bi , bidel
      COMPLEX*16 , DIMENSION(MAXTRM) :: c , cm , cs , d , dm , ds
      LOGICAL :: evenl
      INTEGER :: i , idel , imax , j , k , l , ld2 , m , mmax , nummom
      LOGICAL , SAVE :: pass1
      REAL*8 , DIMENSION(MAXRCP) , SAVE :: recip
      REAL*8 :: thesum
!
! End of declarations rewritten by SPAG
!
!
!           ** specification of local variables
!
!      am(m)       numerical coefficients  a-sub-m-super-l
!                     in dave, eqs. 1-15, as simplified in ref. 5.
!
!      bi(i)       numerical coefficients  b-sub-i-super-l
!                     in dave, eqs. 1-15, as simplified in ref. 5.
!
!      bidel(i)    1/2 bi(i) times factor capital-del in dave
!
!      cm,dm()     arrays c and d in dave, eqs. 16-17 (mueller form),
!                     calculated using recurrence derived in ref. 5
!
!      cs,ds()     arrays c and d in ref. 4, eqs. a5-a6 (sekera form),
!                     calculated using recurrence derived in ref. 5
!
!      c,d()       either -cm,dm- or -cs,ds-, depending on -ipolzn-
!
!      evenl       true for even-numbered moments;  false otherwise
!
!      idel        1 + little-del  in dave
!
!      maxtrm      max. no. of terms in mie series
!
!      maxmom      max. no. of non-zero moments
!
!      nummom      number of non-zero moments
!
!      recip(k)    1 / k
!
      EQUIVALENCE (c,cm) , (d,dm)
      DATA pass1/.TRUE./
!
!
      IF ( pass1 ) THEN
!
         DO k = 1 , MAXRCP
            recip(k) = 1.0/k
         ENDDO
         pass1 = .FALSE.
!
      ENDIF
!
      DO j = 1 , max0(1,Npquan)
         DO l = 0 , Nmom
            Pmom(l,j) = 0.0
         ENDDO
      ENDDO
!
      IF ( Ntrm==1 ) THEN
         CALL lpco1t(Nmom,Ipolzn,Momdim,Calcmo,A,B,Pmom)
         RETURN
      ELSEIF ( Ntrm==2 ) THEN
         CALL lpco2t(Nmom,Ipolzn,Momdim,Calcmo,A,B,Pmom)
         RETURN
      ENDIF
!
      IF ( Ntrm+2>MAXTRM ) CALL errmsg('lpcoef--parameter maxtrm too small',.TRUE.)
!
!                                     ** calculate mueller c, d arrays
      cm(Ntrm+2) = (0.,0.)
      dm(Ntrm+2) = (0.,0.)
      cm(Ntrm+1) = (1.-recip(Ntrm+1))*B(Ntrm)
      dm(Ntrm+1) = (1.-recip(Ntrm+1))*A(Ntrm)
      cm(Ntrm) = (recip(Ntrm)+recip(Ntrm+1))*A(Ntrm) + (1.-recip(Ntrm))*B(Ntrm-1)
      dm(Ntrm) = (recip(Ntrm)+recip(Ntrm+1))*B(Ntrm) + (1.-recip(Ntrm))*A(Ntrm-1)
!
      DO k = Ntrm - 1 , 2 , -1
         cm(k) = cm(k+2) - (1.+recip(k+1))*B(k+1) + (recip(k)+recip(k+1))*A(k) + (1.-recip(k))*B(k-1)
         dm(k) = dm(k+2) - (1.+recip(k+1))*A(k+1) + (recip(k)+recip(k+1))*B(k) + (1.-recip(k))*A(k-1)
      ENDDO
      cm(1) = cm(3) + 1.5*(A(1)-B(2))
      dm(1) = dm(3) + 1.5*(B(1)-A(2))
!
      IF ( Ipolzn>=0 ) THEN
!
         DO k = 1 , Ntrm + 2
            c(k) = (2*k-1)*cm(k)
            d(k) = (2*k-1)*dm(k)
         ENDDO
!
      ELSE
!                                    ** compute sekera c and d arrays
         cs(Ntrm+2) = (0.,0.)
         ds(Ntrm+2) = (0.,0.)
         cs(Ntrm+1) = (0.,0.)
         ds(Ntrm+1) = (0.,0.)
!
         DO k = Ntrm , 1 , -1
            cs(k) = cs(k+2) + (2*k+1)*(cm(k+1)-B(k))
            ds(k) = ds(k+2) + (2*k+1)*(dm(k+1)-A(k))
         ENDDO
!
         DO k = 1 , Ntrm + 2
            c(k) = (2*k-1)*cs(k)
            d(k) = (2*k-1)*ds(k)
         ENDDO
!
      ENDIF
!
!
      IF ( Ipolzn<0 ) nummom = min0(Nmom,2*Ntrm-2)
      IF ( Ipolzn>=0 ) nummom = min0(Nmom,2*Ntrm)
      IF ( nummom>MAXMOM ) CALL errmsg('lpcoef--parameter maxtrm too small',.TRUE.)
!
!                               ** loop over moments
      SPAG_Loop_1_1: DO l = 0 , nummom
         ld2 = l/2
         evenl = mod(l,2)==0
!                                    ** calculate numerical coefficients
!                                    ** a-sub-m and b-sub-i in dave
!                                    ** double-sums for moments
         IF ( l==0 ) THEN
!
            idel = 1
            DO m = 0 , Ntrm
               am(m) = 2.0*recip(2*m+1)
            ENDDO
            bi(0) = 1.0
!
         ELSEIF ( evenl ) THEN
!
            idel = 1
            DO m = ld2 , Ntrm
               am(m) = (1.+recip(2*m-l+1))*am(m)
            ENDDO
            DO i = 0 , ld2 - 1
               bi(i) = (1.-recip(l-2*i))*bi(i)
            ENDDO
            bi(ld2) = (2.-recip(l))*bi(ld2-1)
!
         ELSE
!
            idel = 2
            DO m = ld2 , Ntrm
               am(m) = (1.-recip(2*m+l+2))*am(m)
            ENDDO
            DO i = 0 , ld2
               bi(i) = (1.-recip(l+2*i+1))*bi(i)
            ENDDO
!
         ENDIF
!                                     ** establish upper limits for sums
!                                     ** and incorporate factor capital-
!                                     ** del into b-sub-i
         mmax = Ntrm - idel
         IF ( Ipolzn>=0 ) mmax = mmax + 1
         imax = min0(ld2,mmax-ld2)
         IF ( imax<0 ) EXIT SPAG_Loop_1_1
         DO i = 0 , imax
            bidel(i) = bi(i)
         ENDDO
         IF ( evenl ) bidel(0) = 0.5*bidel(0)
!
!                                    ** perform double sums just for
!                                    ** phase quantities desired by user
         IF ( Ipolzn==0 ) THEN
!
            DO i = 0 , imax
!                                           ** vectorizable loop (cray)
               thesum = 0.0
               DO m = ld2 , mmax - i
                  thesum = thesum + am(m)*(dble(c(m-i+1)*dconjg(c(m+i+idel)))+dble(d(m-i+1)*dconjg(d(m+i+idel))))
               ENDDO
               Pmom(l,1) = Pmom(l,1) + bidel(i)*thesum
            ENDDO
            Pmom(l,1) = 0.5*Pmom(l,1)
            CYCLE
!
         ENDIF
!
         IF ( Calcmo(1) ) THEN
            DO i = 0 , imax
!                                           ** vectorizable loop (cray)
               thesum = 0.0
               DO m = ld2 , mmax - i
                  thesum = thesum + am(m)*dble(c(m-i+1)*dconjg(c(m+i+idel)))
               ENDDO
               Pmom(l,1) = Pmom(l,1) + bidel(i)*thesum
            ENDDO
         ENDIF
!
!
         IF ( Calcmo(2) ) THEN
            DO i = 0 , imax
!                                           ** vectorizable loop (cray)
               thesum = 0.0
               DO m = ld2 , mmax - i
                  thesum = thesum + am(m)*dble(d(m-i+1)*dconjg(d(m+i+idel)))
               ENDDO
               Pmom(l,2) = Pmom(l,2) + bidel(i)*thesum
            ENDDO
         ENDIF
!
!
         IF ( Calcmo(3) ) THEN
            DO i = 0 , imax
!                                           ** vectorizable loop (cray)
               thesum = 0.0
               DO m = ld2 , mmax - i
                  thesum = thesum + am(m)*(dble(c(m-i+1)*dconjg(d(m+i+idel)))+dble(c(m+i+idel)*dconjg(d(m-i+1))))
               ENDDO
               Pmom(l,3) = Pmom(l,3) + bidel(i)*thesum
            ENDDO
            Pmom(l,3) = 0.5*Pmom(l,3)
         ENDIF
!
!
         IF ( Calcmo(4) ) THEN
            DO i = 0 , imax
!                                           ** vectorizable loop (cray)
               thesum = 0.0
               DO m = ld2 , mmax - i
                  thesum = thesum + am(m)*(dimag(c(m-i+1)*dconjg(d(m+i+idel)))+dimag(c(m+i+idel)*dconjg(d(m-i+1))))
               ENDDO
               Pmom(l,4) = Pmom(l,4) + bidel(i)*thesum
            ENDDO
            Pmom(l,4) = -0.5*Pmom(l,4)
         ENDIF
!
      ENDDO SPAG_Loop_1_1
!
!
   END SUBROUTINE lpcoef
!*********************************************************************
   SUBROUTINE lpco1t(Nmom,Ipolzn,Momdim,Calcmo,A,B,Pmom)
!
!         calculate legendre polynomial expansion coefficients (also
!         called moments) for phase quantities in special case where
!         no. terms in mie series = 1
!
!        input:  nmom, ipolzn, momdim     'miev0' arguments
!                calcmo                   flags calculated from -ipolzn-
!                a(1), b(1)               mie series coefficients
!
!        output: pmom                     legendre moments
!
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Momdim
      INTEGER , INTENT(IN) :: Nmom
      INTEGER , INTENT(IN) :: Ipolzn
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Calcmo
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: A
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: B
      REAL*8 , INTENT(OUT) , DIMENSION(0:Momdim,*) :: Pmom
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: a1b1c , ctmp
      REAL*8 :: a1sq , b1sq
      INTEGER :: l , nummom
      REAL*8 :: sq
!
! End of declarations rewritten by SPAG
!
      sq(ctmp) = dble(ctmp)**2 + dimag(ctmp)**2
!
!
      a1sq = sq(A(1))
      b1sq = sq(B(1))
      a1b1c = A(1)*dconjg(B(1))
!
      IF ( Ipolzn<0 ) THEN
!
         IF ( Calcmo(1) ) Pmom(0,1) = 2.25*b1sq
         IF ( Calcmo(2) ) Pmom(0,2) = 2.25*a1sq
         IF ( Calcmo(3) ) Pmom(0,3) = 2.25*dble(a1b1c)
         IF ( Calcmo(4) ) Pmom(0,4) = 2.25*dimag(a1b1c)
!
      ELSE
!
         nummom = min0(Nmom,2)
!                                   ** loop over moments
         DO l = 0 , nummom
!
            IF ( Ipolzn==0 ) THEN
               IF ( l==0 ) Pmom(l,1) = 1.5*(a1sq+b1sq)
               IF ( l==1 ) Pmom(l,1) = 1.5*dble(a1b1c)
               IF ( l==2 ) Pmom(l,1) = 0.15*(a1sq+b1sq)
               CYCLE
            ENDIF
!
            IF ( Calcmo(1) ) THEN
               IF ( l==0 ) Pmom(l,1) = 2.25*(a1sq+b1sq/3.)
               IF ( l==1 ) Pmom(l,1) = 1.5*dble(a1b1c)
               IF ( l==2 ) Pmom(l,1) = 0.3*b1sq
            ENDIF
!
            IF ( Calcmo(2) ) THEN
               IF ( l==0 ) Pmom(l,2) = 2.25*(b1sq+a1sq/3.)
               IF ( l==1 ) Pmom(l,2) = 1.5*dble(a1b1c)
               IF ( l==2 ) Pmom(l,2) = 0.3*a1sq
            ENDIF
!
            IF ( Calcmo(3) ) THEN
               IF ( l==0 ) Pmom(l,3) = 3.0*dble(a1b1c)
               IF ( l==1 ) Pmom(l,3) = 0.75*(a1sq+b1sq)
               IF ( l==2 ) Pmom(l,3) = 0.3*dble(a1b1c)
            ENDIF
!
            IF ( Calcmo(4) ) THEN
               IF ( l==0 ) Pmom(l,4) = -1.5*dimag(a1b1c)
               IF ( l==1 ) Pmom(l,4) = 0.0
               IF ( l==2 ) Pmom(l,4) = 0.3*dimag(a1b1c)
            ENDIF
!
         ENDDO
!
      ENDIF
!
   END SUBROUTINE lpco1t
!********************************************************************
   SUBROUTINE lpco2t(Nmom,Ipolzn,Momdim,Calcmo,A,B,Pmom)
!
!         calculate legendre polynomial expansion coefficients (also
!         called moments) for phase quantities in special case where
!         no. terms in mie series = 2
!
!        input:  nmom, ipolzn, momdim     'miev0' arguments
!                calcmo                   flags calculated from -ipolzn-
!                a(1-2), b(1-2)           mie series coefficients
!
!        output: pmom                     legendre moments
!
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Momdim
      INTEGER , INTENT(IN) :: Nmom
      INTEGER , INTENT(IN) :: Ipolzn
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Calcmo
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: A
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: B
      REAL*8 , INTENT(OUT) , DIMENSION(0:Momdim,*) :: Pmom
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: a2c , b2c , ca , cac , cat , cb , cbc , cbt , cg , ch , ctmp
      REAL*8 :: a2sq , b2sq , pm1 , pm2
      INTEGER :: l , nummom
      REAL*8 :: sq
!
! End of declarations rewritten by SPAG
!
      sq(ctmp) = dble(ctmp)**2 + dimag(ctmp)**2
!
!
      ca = 3.*A(1) - 5.*B(2)
      cat = 3.*B(1) - 5.*A(2)
      cac = dconjg(ca)
      a2sq = sq(A(2))
      b2sq = sq(B(2))
      a2c = dconjg(A(2))
      b2c = dconjg(B(2))
!
      IF ( Ipolzn<0 ) THEN
!                                   ** loop over sekera moments
         nummom = min0(Nmom,2)
         DO l = 0 , nummom
!
            IF ( Calcmo(1) ) THEN
               IF ( l==0 ) Pmom(l,1) = 0.25*(sq(cat)+(100./3.)*b2sq)
               IF ( l==1 ) Pmom(l,1) = (5./3.)*dble(cat*b2c)
               IF ( l==2 ) Pmom(l,1) = (10./3.)*b2sq
            ENDIF
!
            IF ( Calcmo(2) ) THEN
               IF ( l==0 ) Pmom(l,2) = 0.25*(sq(ca)+(100./3.)*a2sq)
               IF ( l==1 ) Pmom(l,2) = (5./3.)*dble(ca*a2c)
               IF ( l==2 ) Pmom(l,2) = (10./3.)*a2sq
            ENDIF
!
            IF ( Calcmo(3) ) THEN
               IF ( l==0 ) Pmom(l,3) = 0.25*dble(cat*cac+(100./3.)*B(2)*a2c)
               IF ( l==1 ) Pmom(l,3) = 5./6.*dble(B(2)*cac+cat*a2c)
               IF ( l==2 ) Pmom(l,3) = 10./3.*dble(B(2)*a2c)
            ENDIF
!
            IF ( Calcmo(4) ) THEN
               IF ( l==0 ) Pmom(l,4) = -0.25*dimag(cat*cac+(100./3.)*B(2)*a2c)
               IF ( l==1 ) Pmom(l,4) = -5./6.*dimag(B(2)*cac+cat*a2c)
               IF ( l==2 ) Pmom(l,4) = -10./3.*dimag(B(2)*a2c)
            ENDIF
!
         ENDDO
!
      ELSE
!
         cb = 3.*B(1) + 5.*A(2)
         cbt = 3.*A(1) + 5.*B(2)
         cbc = dconjg(cb)
         cg = (cbc*cbt+10.*(cac*A(2)+b2c*cat))/3.
         ch = 2.*(cbc*A(2)+b2c*cbt)
!
!                                   ** loop over mueller moments
         nummom = min0(Nmom,4)
         DO l = 0 , nummom
!
            IF ( Ipolzn==0 .OR. Calcmo(1) ) THEN
               IF ( l==0 ) pm1 = 0.25*sq(ca) + sq(cb)/12. + (5./3.)*dble(ca*b2c) + 5.*b2sq
               IF ( l==1 ) pm1 = dble(cb*(cac/6.+b2c))
               IF ( l==2 ) pm1 = sq(cb)/30. + (20./7.)*b2sq + (2./3.)*dble(ca*b2c)
               IF ( l==3 ) pm1 = (2./7.)*dble(cb*b2c)
               IF ( l==4 ) pm1 = (40./63.)*b2sq
               IF ( Calcmo(1) ) Pmom(l,1) = pm1
            ENDIF
!
            IF ( Ipolzn==0 .OR. Calcmo(2) ) THEN
               IF ( l==0 ) pm2 = 0.25*sq(cat) + sq(cbt)/12. + (5./3.)*dble(cat*a2c) + 5.*a2sq
               IF ( l==1 ) pm2 = dble(cbt*(dconjg(cat)/6.+a2c))
               IF ( l==2 ) pm2 = sq(cbt)/30. + (20./7.)*a2sq + (2./3.)*dble(cat*a2c)
               IF ( l==3 ) pm2 = (2./7.)*dble(cbt*a2c)
               IF ( l==4 ) pm2 = (40./63.)*a2sq
               IF ( Calcmo(2) ) Pmom(l,2) = pm2
            ENDIF
!
            IF ( Ipolzn==0 ) THEN
               Pmom(l,1) = 0.5*(pm1+pm2)
               CYCLE
            ENDIF
!
            IF ( Calcmo(3) ) THEN
               IF ( l==0 ) Pmom(l,3) = 0.25*dble(cac*cat+cg+20.*b2c*A(2))
               IF ( l==1 ) Pmom(l,3) = dble(cac*cbt+cbc*cat+3.*ch)/12.
               IF ( l==2 ) Pmom(l,3) = 0.1*dble(cg+(200./7.)*b2c*A(2))
               IF ( l==3 ) Pmom(l,3) = dble(ch)/14.
               IF ( l==4 ) Pmom(l,3) = 40./63.*dble(b2c*A(2))
            ENDIF
!
            IF ( Calcmo(4) ) THEN
               IF ( l==0 ) Pmom(l,4) = 0.25*dimag(cac*cat+cg+20.*b2c*A(2))
               IF ( l==1 ) Pmom(l,4) = dimag(cac*cbt+cbc*cat+3.*ch)/12.
               IF ( l==2 ) Pmom(l,4) = 0.1*dimag(cg+(200./7.)*b2c*A(2))
               IF ( l==3 ) Pmom(l,4) = dimag(ch)/14.
               IF ( l==4 ) Pmom(l,4) = 40./63.*dimag(b2c*A(2))
            ENDIF
!
         ENDDO
!
      ENDIF
!
   END SUBROUTINE lpco2t
!*********************************************************************
   SUBROUTINE biga(Cior,Xx,Ntrm,Noabs,Yesang,Rbiga,Cbiga)
!
!        calculate logarithmic derivatives of j-bessel-function
!
!     input :  cior, xx, ntrm, noabs, yesang  (defined in 'miev0')
!
!    output :  rbiga or cbiga  (defined in 'miev0')
!
!    internal variables :
!
!       confra     value of lentz continued fraction for -cbiga(ntrm)-,
!                     used to initialize downward recurrence.
!       down       = true, use down-recurrence.  false, do not.
!       f1,f2,f3   arithmetic statement functions used in determining
!                     whether to use up-  or down-recurrence
!                     ( ref. 2, eqs. 6-8 )
!       mre        real refractive index
!       mim        imaginary refractive index
!       rezinv     1 / ( mre * xx ); temporary variable for recurrence
!       zinv       1 / ( cior * xx ); temporary variable for recurrence
!
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX*16 , INTENT(IN) :: Cior
      REAL*8 :: Xx
      INTEGER :: Ntrm
      LOGICAL , INTENT(IN) :: Noabs
      LOGICAL , INTENT(IN) :: Yesang
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Rbiga
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(*) :: Cbiga
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: ctmp , zinv
      LOGICAL :: down
      REAL*8 :: f1 , f2 , f3
      REAL*8 :: mim , mre , rezinv , rtmp
      INTEGER :: n
!
! End of declarations rewritten by SPAG
!
!      complex*16  cior, ctmp, confra, cbiga(*), zinv
      f1(mre) = -8.0 + mre**2*(26.22+mre*(-0.4474+mre**3*(0.00204-0.000175*mre)))
      f2(mre) = 3.9 + mre*(-10.8+13.78*mre)
      f3(mre) = -15.04 + mre*(8.42+16.35*mre)
!
!                                  ** decide whether 'biga' can be
!                                  ** calculated by up-recurrence
      mre = dble(Cior)
      mim = dabs(dimag(Cior))
      IF ( mre<1.0 .OR. mre>10.0 .OR. mim>10.0 ) THEN
         down = .TRUE.
      ELSEIF ( Yesang ) THEN
         down = .TRUE.
         IF ( mim*Xx<f2(mre) ) down = .FALSE.
      ELSE
         down = .TRUE.
         IF ( mim*Xx<f1(mre) ) down = .FALSE.
      ENDIF
!
      zinv = 1.0/(Cior*Xx)
      rezinv = 1.0/(mre*Xx)
!
      IF ( down ) THEN
!                          ** compute initial high-order 'biga' using
!                          ** lentz method ( ref. 1, pp. 17-20 )
!
         ctmp = confra(Ntrm,zinv,Xx)
!
!                                   *** downward recurrence for 'biga'
!                                   *** ( ref. 1, eq. 22 )
         IF ( Noabs ) THEN
!                                            ** no-absorption case
            Rbiga(Ntrm) = dble(ctmp)
            DO n = Ntrm , 2 , -1
               Rbiga(n-1) = (n*rezinv) - 1.0/((n*rezinv)+Rbiga(n))
            ENDDO
!
         ELSE
!                                            ** absorptive case
            Cbiga(Ntrm) = ctmp
            DO n = Ntrm , 2 , -1
               Cbiga(n-1) = (n*zinv) - 1.0/((n*zinv)+Cbiga(n))
            ENDDO
!
         ENDIF
!
!                              *** upward recurrence for 'biga'
!                              *** ( ref. 1, eqs. 20-21 )
      ELSEIF ( Noabs ) THEN
!                                            ** no-absorption case
         rtmp = dsin(mre*Xx)
         Rbiga(1) = -rezinv + rtmp/(rtmp*rezinv-dcos(mre*Xx))
         DO n = 2 , Ntrm
            Rbiga(n) = -(n*rezinv) + 1.0/((n*rezinv)-Rbiga(n-1))
         ENDDO
!
      ELSE
!                                                ** absorptive case
         ctmp = cdexp(-dcmplx(0.D0,2.D0)*Cior*Xx)
         Cbiga(1) = -zinv + (1.-ctmp)/(zinv*(1.-ctmp)-dcmplx(0.D0,1.D0)*(1.+ctmp))
         DO n = 2 , Ntrm
            Cbiga(n) = -(n*zinv) + 1.0/((n*zinv)-Cbiga(n-1))
         ENDDO
!
      ENDIF
!
   END SUBROUTINE biga
!**********************************************************************
   FUNCTION confra(N,Zinv,Xx)
!
!         compute bessel function ratio capital-a-sub-n from its
!         continued fraction using lentz method ( ref. 1, pp. 17-20 )
!
!         zinv = reciprocal of argument of capital-a
!
!    i n t e r n a l    v a r i a b l e s
!    ------------------------------------
!
!    cak      term in continued fraction expansion of capital-a
!                ( ref. 1, eq. 25 )
!    capt     factor used in lentz iteration for capital-a
!                ( ref. 1, eq. 27 )
!    cdenom   denominator in -capt-  ( ref. 1, eq. 28b )
!    cnumer   numerator   in -capt-  ( ref. 1, eq. 28a )
!    cdtd     product of two successive denominators of -capt-
!                factors  ( ref. 1, eq. 34c )
!    cntn     product of two successive numerators of -capt-
!                factors  ( ref. 1, eq. 34b )
!    eps1     ill-conditioning criterion
!    eps2     convergence criterion
!    kk       subscript k of -cak-  ( ref. 1, eq. 25b )
!    kount    iteration counter ( used only to prevent runaway )
!    maxit    max. allowed no. of iterations
!    mm       + 1  and - 1, alternately
!
!
! Function and Dummy argument declarations rewritten by SPAG
!
      COMPLEX*16 :: confra
      INTEGER , INTENT(IN) :: N
      COMPLEX*16 , INTENT(IN) :: Zinv
      REAL*8 :: Xx
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: cak , capt , cdenom , cdtd , cntn , cnumer
      REAL*8 , SAVE :: eps1 , eps2
      INTEGER :: kk , kount , mm
      INTEGER , SAVE :: maxit
!
! End of declarations rewritten by SPAG
!
      DATA maxit/10000/
      DATA eps1/1.D-2/ , eps2/1.D-8/
!
!                                      *** ref. 1, eqs. 25a, 27
      confra = (N+1)*Zinv
      mm = -1
      kk = 2*N + 3
      cak = (mm*kk)*Zinv
      cdenom = cak
      cnumer = cdenom + 1.0/confra
      kount = 1
      SPAG_Loop_1_1: DO
!
         kount = kount + 1
         IF ( kount>maxit ) CALL errmsg('confra--iteration failed to converge$',.TRUE.)
!
!                                         *** ref. 2, eq. 25b
         mm = -mm
         kk = kk + 2
         cak = (mm*kk)*Zinv
!                                         *** ref. 2, eq. 32
         IF ( cdabs(cnumer/cak)<=eps1 .OR. cdabs(cdenom/cak)<=eps1 ) THEN
!
!                                  ** ill-conditioned case -- stride
!                                  ** two terms instead of one
!
!                                         *** ref. 2, eqs. 34
            cntn = cak*cnumer + 1.0
            cdtd = cak*cdenom + 1.0
            confra = (cntn/cdtd)*confra
!                                             *** ref. 2, eq. 25b
            mm = -mm
            kk = kk + 2
            cak = (mm*kk)*Zinv
!                                         *** ref. 2, eqs. 35
            cnumer = cak + cnumer/cntn
            cdenom = cak + cdenom/cdtd
            kount = kount + 1
            CYCLE
!
         ELSE
!                                ** well-conditioned case
!
!                                        *** ref. 2, eqs. 26, 27
            capt = cnumer/cdenom
            confra = capt*confra
!                                    ** check for convergence
!                                    ** ( ref. 2, eq. 31 )
!
            IF ( dabs(dble(capt)-1.0)>=eps2 .OR. dabs(dimag(capt))>=eps2 ) THEN
!
!                                        *** ref. 2, eqs. 30a-b
               cnumer = cak + 1.0/cnumer
               cdenom = cak + 1.0/cdenom
               CYCLE
            ENDIF
         ENDIF
         EXIT SPAG_Loop_1_1
      ENDDO SPAG_Loop_1_1
!
!
   END FUNCTION confra
!********************************************************************
   SUBROUTINE miprnt(Prnt,Xx,Perfct,Crefin,Numang,Xmu,Qext,Qsca,Gqsc,Nmom,Ipolzn,Momdim,Calcmo,Pmom,Sforw,Sback,Tforw,Tback,S1,S2)
!
!         print scattering quantities of a single particle
!
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Momdim
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Prnt
      REAL*8 , INTENT(IN) :: Xx
      LOGICAL , INTENT(IN) :: Perfct
      COMPLEX*16 , INTENT(IN) :: Crefin
      INTEGER , INTENT(IN) :: Numang
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xmu
      REAL*8 , INTENT(IN) :: Qext
      REAL*8 , INTENT(IN) :: Qsca
      REAL*8 , INTENT(IN) :: Gqsc
      INTEGER , INTENT(IN) :: Nmom
      INTEGER , INTENT(IN) :: Ipolzn
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Calcmo
      REAL*8 , INTENT(IN) , DIMENSION(0:Momdim,*) :: Pmom
      COMPLEX*16 , INTENT(IN) :: Sforw
      COMPLEX*16 , INTENT(IN) :: Sback
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: Tforw
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: Tback
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: S1
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: S2
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 :: fi1 , fi2 , fnorm
      CHARACTER(22) :: fmt
      INTEGER :: i , j , m
!
! End of declarations rewritten by SPAG
!
!
!
      IF ( Perfct ) WRITE (*,'(''1'',10x,a,1p,e11.4)') 'perfectly conducting case, size parameter =' , Xx
      IF ( .NOT.Perfct ) WRITE (*,'(''1'',10x,3(a,1p,e11.4))') 'refractive index:  real ' , dble(Crefin) , '  imag ' ,             &
                              & dimag(Crefin) , ',   size parameter =' , Xx
!
      IF ( Prnt(1) .AND. Numang>0 ) THEN
!
         WRITE (*,'(/,a)') '    cos(angle)  ------- s1 ---------  ------- s2 ---------'//                                          &
                          &'  --- s1*conjg(s2) ---   i1=s1**2   i2=s2**2  (i1+i2)/2'//'  deg polzn'
         DO i = 1 , Numang
            fi1 = dble(S1(i))**2 + dimag(S1(i))**2
            fi2 = dble(S2(i))**2 + dimag(S2(i))**2
            WRITE (*,'( i4, f10.6, 1p,10e11.3 )') i , Xmu(i) , S1(i) , S2(i) , S1(i)*dconjg(S2(i)) , fi1 , fi2 , 0.5*(fi1+fi2) ,   &
                 & (fi2-fi1)/(fi2+fi1)
         ENDDO
!
      ENDIF
!
!
      IF ( Prnt(2) ) THEN
!
         WRITE (*,'(/,a,9x,a,17x,a,17x,a,/,(0p,f7.2, 1p,6e12.3) )') '  angle' , 's-sub-1' , 't-sub-1' , 't-sub-2' , 0.0 , Sforw ,  &
              & Tforw(1) , Tforw(2) , 180. , Sback , Tback(1) , Tback(2)
         WRITE (*,'(/,4(a,1p,e11.4))') ' efficiency factors,  extinction:' , Qext , '   scattering:' , Qsca , '   absorption:' ,   &
                                     & Qext - Qsca , '   rad. pressure:' , Qext - Gqsc
!
         IF ( Nmom>0 ) THEN
!
            WRITE (*,'(/,a)') ' normalized moments of :'
            IF ( Ipolzn==0 ) WRITE (*,'(''+'',27x,a)') 'phase fcn'
            IF ( Ipolzn>0 ) WRITE (*,'(''+'',33x,a)') 'm1           m2          s21          d21'
            IF ( Ipolzn<0 ) WRITE (*,'(''+'',33x,a)') 'r1           r2           r3           r4'
!
            fnorm = 4./(Xx**2*Qsca)
            DO m = 0 , Nmom
               WRITE (*,'(a,i4)') '      moment no.' , m
               DO j = 1 , 4
                  IF ( Calcmo(j) ) THEN
                     WRITE (fmt,99001) 24 + (j-1)*13
!
!
99001                FORMAT ('( ''+'', t',i2,', 1p,e13.4 )')
                     WRITE (*,fmt) fnorm*Pmom(m,j)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
      ENDIF
   END SUBROUTINE miprnt
!**************************************************************************
   SUBROUTINE small1(Xx,Numang,Xmu,Qext,Qsca,Gqsc,Sforw,Sback,S1,S2,Tforw,Tback,A,B)
!
!       small-particle limit of mie quantities in totally reflecting
!       limit ( mie series truncated after 2 terms )
!
!        a,b       first two mie coefficients, with numerator and
!                  denominator expanded in powers of -xx- ( a factor
!                  of xx**3 is missing but is restored before return
!                  to calling program )  ( ref. 2, p. 1508 )
!
!
! PARAMETER definitions rewritten by SPAG
!
      REAL*8 , PARAMETER :: TWOTHR = 2./3. , FIVTHR = 5./3. , FIVNIN = 5./9.
!
! Dummy argument declarations rewritten by SPAG
!
      REAL*8 , INTENT(IN) :: Xx
      INTEGER , INTENT(IN) :: Numang
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xmu
      REAL*8 , INTENT(OUT) :: Qext
      REAL*8 , INTENT(INOUT) :: Qsca
      REAL*8 , INTENT(OUT) :: Gqsc
      COMPLEX*16 , INTENT(OUT) :: Sforw
      COMPLEX*16 , INTENT(OUT) :: Sback
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: S1
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: S2
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: Tforw
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: Tback
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(2) :: A
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(2) :: B
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: ctmp
      INTEGER :: j
      REAL*8 :: rtmp
      REAL*8 :: sq
!
! End of declarations rewritten by SPAG
!
!
      sq(ctmp) = dble(ctmp)**2 + dimag(ctmp)**2
!
!
      A(1) = dcmplx(0.D0,TWOTHR*(1.-0.2*Xx**2))/dcmplx(1.D0-0.5*Xx**2,TWOTHR*Xx**3)
!
      B(1) = dcmplx(0.D0,-(1.-0.1*Xx**2)/3.)/dcmplx(1.D0+0.5*Xx**2,-Xx**3/3.)
!
      A(2) = dcmplx(0.D0,Xx**2/30.)
      B(2) = dcmplx(0.D0,-Xx**2/45.)
!
      Qsca = 6.*Xx**4*(sq(A(1))+sq(B(1))+FIVTHR*(sq(A(2))+sq(B(2))))
      Qext = Qsca
      Gqsc = 6.*Xx**4*dble(A(1)*dconjg(A(2)+B(1))+(B(1)+FIVNIN*A(2))*dconjg(B(2)))
!
      rtmp = 1.5*Xx**3
      Sforw = rtmp*(A(1)+B(1)+FIVTHR*(A(2)+B(2)))
      Sback = rtmp*(A(1)-B(1)-FIVTHR*(A(2)-B(2)))
      Tforw(1) = rtmp*(B(1)+FIVTHR*(2.*B(2)-A(2)))
      Tforw(2) = rtmp*(A(1)+FIVTHR*(2.*A(2)-B(2)))
      Tback(1) = rtmp*(B(1)-FIVTHR*(2.*B(2)+A(2)))
      Tback(2) = rtmp*(A(1)-FIVTHR*(2.*A(2)+B(2)))
!
      DO j = 1 , Numang
         S1(j) = rtmp*(A(1)+B(1)*Xmu(j)+FIVTHR*(A(2)*Xmu(j)+B(2)*(2.*Xmu(j)**2-1.)))
         S2(j) = rtmp*(B(1)+A(1)*Xmu(j)+FIVTHR*(B(2)*Xmu(j)+A(2)*(2.*Xmu(j)**2-1.)))
      ENDDO
!                                     ** recover actual mie coefficients
      A(1) = Xx**3*A(1)
      A(2) = Xx**3*A(2)
      B(1) = Xx**3*B(1)
      B(2) = Xx**3*B(2)
!
   END SUBROUTINE small1
!*************************************************************************
   SUBROUTINE small2(Xx,Cior,Calcqe,Numang,Xmu,Qext,Qsca,Gqsc,Sforw,Sback,S1,S2,Tforw,Tback,A,B)
!
!       small-particle limit of mie quantities for general refractive
!       index ( mie series truncated after 2 terms )
!
!        a,b       first two mie coefficients, with numerator and
!                  denominator expanded in powers of -xx- ( a factor
!                  of xx**3 is missing but is restored before return
!                  to calling program )  ( ref. 2, p. 1508 )
!
!        ciorsq    square of refractive index
!
!
! PARAMETER definitions rewritten by SPAG
!
      REAL*8 , PARAMETER :: TWOTHR = 2./3. , FIVTHR = 5./3.
!
! Dummy argument declarations rewritten by SPAG
!
      REAL*8 , INTENT(IN) :: Xx
      COMPLEX*16 , INTENT(IN) :: Cior
      LOGICAL , INTENT(IN) :: Calcqe
      INTEGER , INTENT(IN) :: Numang
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xmu
      REAL*8 , INTENT(OUT) :: Qext
      REAL*8 , INTENT(INOUT) :: Qsca
      REAL*8 , INTENT(OUT) :: Gqsc
      COMPLEX*16 , INTENT(OUT) :: Sforw
      COMPLEX*16 , INTENT(OUT) :: Sback
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: S1
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: S2
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(*) :: Tforw
      COMPLEX*16 , INTENT(OUT) , DIMENSION(*) :: Tback
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(2) :: A
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(2) :: B
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 :: ciorsq , ctmp
      INTEGER :: j
      REAL*8 :: rtmp
      REAL*8 :: sq
!
! End of declarations rewritten by SPAG
!
!
      sq(ctmp) = dble(ctmp)**2 + dimag(ctmp)**2
!
!
      ciorsq = Cior**2
      ctmp = dcmplx(0.D0,TWOTHR)*(ciorsq-1.0)
      A(1) = ctmp*(1.0-0.1*Xx**2+(ciorsq/350.+1./280.)*Xx**4)/(ciorsq+2.0+(1.0-0.7*ciorsq)*Xx**2-(ciorsq**2/175.-0.275*ciorsq+0.25)&
           & *Xx**4+Xx**3*ctmp*(1.0-0.1*Xx**2))
!
      B(1) = (Xx**2/30.)*ctmp*(1.0+(ciorsq/35.-1./14.)*Xx**2)/(1.0-(ciorsq/15.-1./6.)*Xx**2)
!
      A(2) = (0.1*Xx**2)*ctmp*(1.0-Xx**2/14.)/(2.*ciorsq+3.-(ciorsq/7.-0.5)*Xx**2)
!
      Qsca = 6.*Xx**4*(sq(A(1))+sq(B(1))+FIVTHR*sq(A(2)))
      Gqsc = 6.*Xx**4*dble(A(1)*dconjg(A(2)+B(1)))
      Qext = Qsca
      IF ( Calcqe ) Qext = 6.*Xx*dble(A(1)+B(1)+FIVTHR*A(2))
!
      rtmp = 1.5*Xx**3
      Sforw = rtmp*(A(1)+B(1)+FIVTHR*A(2))
      Sback = rtmp*(A(1)-B(1)-FIVTHR*A(2))
      Tforw(1) = rtmp*(B(1)-FIVTHR*A(2))
      Tforw(2) = rtmp*(A(1)+2.*FIVTHR*A(2))
      Tback(1) = Tforw(1)
      Tback(2) = rtmp*(A(1)-2.*FIVTHR*A(2))
!
      DO j = 1 , Numang
         S1(j) = rtmp*(A(1)+(B(1)+FIVTHR*A(2))*Xmu(j))
         S2(j) = rtmp*(B(1)+A(1)*Xmu(j)+FIVTHR*A(2)*(2.*Xmu(j)**2-1.))
      ENDDO
!                                     ** recover actual mie coefficients
      A(1) = Xx**3*A(1)
      A(2) = Xx**3*A(2)
      B(1) = Xx**3*B(1)
      B(2) = (0.,0.)
!
   END SUBROUTINE small2
!***********************************************************************
   SUBROUTINE testmi(Qext,Qsca,Gqsc,Sforw,Sback,S1,S2,Tforw,Tback,Pmom,Momdim,Ok)
!
!         compare mie code test case results with correct answers
!         and return  ok=false  if even one result is inaccurate.
!
!         the test case is :  mie size parameter = 10
!                             refractive index   = 1.5 - 0.1 i
!                             scattering angle = 140 degrees
!                             1 sekera moment
!
!         results for this case may be found among the test cases
!         at the end of reference (1).
!
!         *** note *** when running on some computers, esp. in single
!         precision, the 'accur' criterion below may have to be relaxed.
!         however, if 'accur' must be set larger than 10**-3 for some
!         size parameters, your computer is probably not accurate
!         enough to do mie computations for those size parameters.
!
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Momdim
      REAL*8 , INTENT(IN) :: Qext
      REAL*8 , INTENT(IN) :: Qsca
      REAL*8 , INTENT(IN) :: Gqsc
      COMPLEX*16 , INTENT(IN) :: Sforw
      COMPLEX*16 , INTENT(IN) :: Sback
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: S1
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: S2
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: Tforw
      COMPLEX*16 , INTENT(IN) , DIMENSION(*) :: Tback
      REAL*8 , INTENT(IN) , DIMENSION(0:Momdim,*) :: Pmom
      LOGICAL :: Ok
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 , SAVE :: accur , testgq , testqe , testqs
      REAL*8 :: calc , exact
      INTEGER :: m , n
      REAL*8 , DIMENSION(0:1) , SAVE :: testpm
      COMPLEX*16 , SAVE :: tests1 , tests2 , testsb , testsf
      COMPLEX*16 , DIMENSION(2) , SAVE :: testtb , testtf
      LOGICAL :: wrong
!
! End of declarations rewritten by SPAG
!
!
      DATA testqe/2.459791/ , testqs/1.235144/ , testgq/1.139235/ , testsf/(61.49476,-3.177994)/ , testsb/(1.493434,0.2963657)/ ,  &
         & tests1/(-0.1548380,-1.128972)/ , tests2/(0.05669755,0.5425681)/ , testtf/(12.95238,-136.6436) , (48.54238,133.4656)/ ,  &
         & testtb/(41.88414,-15.57833) , (43.37758,-15.28196)/ , testpm/227.1975 , 183.6898/
!      data   accur / 1.e-5 /
      DATA accur/1.E-4/
      wrong(calc,exact) = dabs((calc-exact)/exact)>accur
!
!
      Ok = .TRUE.
      IF ( wrong(Qext,testqe) ) CALL tstbad('qext',abs((Qext-testqe)/testqe),Ok)
      IF ( wrong(Qsca,testqs) ) CALL tstbad('qsca',abs((Qsca-testqs)/testqs),Ok)
      IF ( wrong(Gqsc,testgq) ) CALL tstbad('gqsc',abs((Gqsc-testgq)/testgq),Ok)
!
      IF ( wrong(dble(Sforw),dble(testsf)) .OR. wrong(dimag(Sforw),dimag(testsf)) )                                                &
         & CALL tstbad('sforw',cdabs((Sforw-testsf)/testsf),Ok)
!
      IF ( wrong(dble(Sback),dble(testsb)) .OR. wrong(dimag(Sback),dimag(testsb)) )                                                &
         & CALL tstbad('sback',cdabs((Sback-testsb)/testsb),Ok)
!
      IF ( wrong(dble(S1(1)),dble(tests1)) .OR. wrong(dimag(S1(1)),dimag(tests1)) ) CALL tstbad('s1',cdabs((S1(1)-tests1)/tests1), &
         & Ok)
!
      IF ( wrong(dble(S2(1)),dble(tests2)) .OR. wrong(dimag(S2(1)),dimag(tests2)) ) CALL tstbad('s2',cdabs((S2(1)-tests2)/tests2), &
         & Ok)
!
      DO n = 1 , 2
         IF ( wrong(dble(Tforw(n)),dble(testtf(n))) .OR. wrong(dimag(Tforw(n)),dimag(testtf(n))) )                                 &
            & CALL tstbad('tforw',cdabs((Tforw(n)-testtf(n))/testtf(n)),Ok)
         IF ( wrong(dble(Tback(n)),dble(testtb(n))) .OR. wrong(dimag(Tback(n)),dimag(testtb(n))) )                                 &
            & CALL tstbad('tback',cdabs((Tback(n)-testtb(n))/testtb(n)),Ok)
      ENDDO
!
      DO m = 0 , 1
         IF ( wrong(Pmom(m,1),testpm(m)) ) CALL tstbad('pmom',dabs((Pmom(m,1)-testpm(m))/testpm(m)),Ok)
      ENDDO
!
!
   END SUBROUTINE testmi
!**************************************************************************
   SUBROUTINE errmsg(Messag,Fatal)
!
!        print out a warning or error message;  abort if error
!
 
      USE module_peg_util , only:peg_message , peg_error_fatal
      LOGICAL Fatal
      LOGICAL , SAVE :: once
      DATA once/.FALSE./
      CHARACTER*80 msg
      CHARACTER*(*) Messag
      INTEGER , SAVE :: maxmsg , nummsg
      DATA nummsg/0/ , maxmsg/100/
      INTEGER , PARAMETER :: LUNERR = -1
!
!
      IF ( Fatal ) THEN
         WRITE (msg,'(a)') 'optical averaging mie fatal error '//Messag
         CALL peg_message(LUNERR,msg)
         CALL peg_error_fatal(LUNERR,msg)
      ENDIF
!
      nummsg = nummsg + 1
      IF ( nummsg>maxmsg ) THEN
!         if ( .not.once )  write ( *,99 )
         IF ( .NOT.once ) THEN
            WRITE (msg,'(a)') 'optical averaging mie: too many warning messages -- no longer printing '
            CALL peg_message(LUNERR,msg)
         ENDIF
         once = .TRUE.
      ELSE
         msg = 'optical averaging mie warning '//Messag
         CALL peg_message(LUNERR,msg)
!         write ( *, '(2a)' )  ' ******* warning >>>>>>  ', messag
      ENDIF
!
!
!   99 format( ///,' >>>>>>  too many warning messages --  ',   &
!         'they will no longer be printed  <<<<<<<', /// )
   END SUBROUTINE errmsg
!********************************************************************
   SUBROUTINE wrtbad(Varnam,Erflag)
!
!          write names of erroneous variables
!
!      input :   varnam = name of erroneous variable to be written
!                         ( character, any length )
!
!      output :  erflag = logical flag, set true by this routine
! ----------------------------------------------------------------------
 
      USE module_peg_util , only:peg_message
      CHARACTER*(*) Varnam
      LOGICAL Erflag
      CHARACTER*80 msg
      INTEGER , SAVE :: maxmsg , nummsg
      DATA nummsg/0/ , maxmsg/50/
      INTEGER , PARAMETER :: LUNERR = -1
!
!
      nummsg = nummsg + 1
!      write ( *, '(3a)' )  ' ****  input variable  ', varnam,   &
!                           '  in error  ****'
      msg = 'optical averaging mie input variable in error '//Varnam
      CALL peg_message(LUNERR,msg)
      Erflag = .TRUE.
      IF ( nummsg==maxmsg ) CALL errmsg('too many input variable errors.  aborting...$',.TRUE.)
!
   END SUBROUTINE wrtbad
!******************************************************************
   SUBROUTINE tstbad(Varnam,Relerr,Ok)
!
!       write name (-varnam-) of variable failing self-test and its
!       percent error from the correct value.  return  ok = false.
!
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER(*) , INTENT(IN) :: Varnam
      REAL*8 , INTENT(IN) :: Relerr
      LOGICAL , INTENT(OUT) :: Ok
!
! End of declarations rewritten by SPAG
!
!
!
      Ok = .FALSE.
      WRITE (*,'(/,3a,1p,e11.2,a)') ' output variable  ' , Varnam , '  differed by' , 100.*Relerr ,                                &
                                   &'  per cent from correct value.  self-test failed.'
!
   END SUBROUTINE tstbad
!******************************************************************
!
   SUBROUTINE sect02(Dgnum_um,Sigmag,Drydens,Iflag,Duma,Nbin,Dlo_um,Dhi_um,Xnum_sect,Xmas_sect)
!
!   user specifies a single log-normal mode and a set of section boundaries
!   prog calculates mass and number for each section
!
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER :: PI = 3.1415926536
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Nbin
      REAL , INTENT(IN) :: Dgnum_um
      REAL , INTENT(IN) :: Sigmag
      REAL , INTENT(IN) :: Drydens
      INTEGER , INTENT(IN) :: Iflag
      REAL , INTENT(IN) :: Duma
      REAL , INTENT(IN) :: Dlo_um
      REAL , INTENT(IN) :: Dhi_um
      REAL , INTENT(OUT) , DIMENSION(Nbin) :: Xnum_sect
      REAL , INTENT(OUT) , DIMENSION(Nbin) :: Xmas_sect
!
! Local variable declarations rewritten by SPAG
!
      REAL :: dgnum , dhi , dlo , dstar , dumfrac , dx , summas , sumnum , sx , sxroot2 , thi , tlo , x0 , x3 , xhi , xlo , xmtot ,&
            & xntot , xvtot
      REAL , DIMENSION(Nbin) :: dhi_sect , dlo_sect
      INTEGER :: n
!
! End of declarations rewritten by SPAG
!
!       real erfc_num_recipes
!
      IF ( Iflag<=1 ) THEN
         xntot = Duma
      ELSE
         xmtot = Duma
         xntot = Duma      !czhao
      ENDIF
!   compute total volume and number for mode
!       dgnum = dgnum_um*1.0e-4
!       sx = log( sigmag )
!       x0 = log( dgnum )
!       x3 = x0 + 3.*sx*sx
!       dstar = dgnum * exp(1.5*sx*sx)
!       if (iflag .le. 1) then
!           xvtot = xntot*(pi/6.0)*dstar*dstar*dstar
!           xmtot = xvtot*drydens*1.0e12
!       else
!           xvtot = xmtot/(drydens*1.0e12)
!           xntot = xvtot/((pi/6.0)*dstar*dstar*dstar)
!       end if
!   compute section boundaries
      dlo = Dlo_um*1.0E-4
      dhi = Dhi_um*1.0E-4
      xlo = log(dlo)
      xhi = log(dhi)
      dx = (xhi-xlo)/Nbin
      DO n = 1 , Nbin
         dlo_sect(n) = exp(xlo+dx*(n-1))
         dhi_sect(n) = exp(xlo+dx*n)
      ENDDO
!   compute modal "working" parameters including total num/vol/mass
      dgnum = Dgnum_um*1.0E-4
      sx = log(Sigmag)
      x0 = log(dgnum)
      x3 = x0 + 3.*sx*sx
      dstar = dgnum*exp(1.5*sx*sx)
      IF ( Iflag<=1 ) THEN
         xvtot = xntot*(PI/6.0)*dstar*dstar*dstar
         xmtot = xvtot*Drydens*1.0E12
!czhao      xvtot = xmtot/(drydens*1.0e12)
!czhao      xntot = xvtot/((pi/6.0)*dstar*dstar*dstar)
      ENDIF
!   compute number and mass for each section
      sxroot2 = sx*sqrt(2.0)
      sumnum = 0.
      summas = 0.
!       write(22,*)
!       write(22,*) 'dgnum_um, sigmag = ', dgnum_um, sigmag
!       write(22,*) 'drydens =', drydens
!       write(22,*) 'ntot (#/cm3), mtot (ug/m3) = ', xntot, xmtot
!        write(22,9220)
!9220    format( /   &
!        '  n   dlo(um)   dhi(um)       number         mass' / )
!9225    format(   i3, 2f10.6, 2(1pe13.4) )
!9230    format( / 'sum over all sections  ', 2(1pe13.4) )
!9231    format(   'modal totals           ', 2(1pe13.4) )
      DO n = 1 , Nbin
         xlo = log(dlo_sect(n))
         xhi = log(dhi_sect(n))
         tlo = (xlo-x0)/sxroot2
         thi = (xhi-x0)/sxroot2
         IF ( tlo<=0. ) THEN
            dumfrac = 0.5*(erfc_num_recipes(-thi)-erfc_num_recipes(-tlo))
         ELSE
            dumfrac = 0.5*(erfc_num_recipes(tlo)-erfc_num_recipes(thi))
         ENDIF
         Xnum_sect(n) = xntot*dumfrac
         tlo = (xlo-x3)/sxroot2
         thi = (xhi-x3)/sxroot2
         IF ( tlo<=0. ) THEN
            dumfrac = 0.5*(erfc_num_recipes(-thi)-erfc_num_recipes(-tlo))
         ELSE
            dumfrac = 0.5*(erfc_num_recipes(tlo)-erfc_num_recipes(thi))
         ENDIF
         Xmas_sect(n) = xmtot*dumfrac
         sumnum = sumnum + Xnum_sect(n)
         summas = summas + Xmas_sect(n)
!           write(22,9225) n, 1.e4*dlo_sect(n), 1.e4*dhi_sect(n),   &
!               xnum_sect(n), xmas_sect(n)
      ENDDO
!       write(22,9230) sumnum, summas
!       write(22,9231) xntot, xmtot
 
   END SUBROUTINE sect02
!-----------------------------------------------------------------------
   FUNCTION erfc_num_recipes(X)
!
!   from press et al, numerical recipes, 1990, page 164
!
!
! Function and Dummy argument declarations rewritten by SPAG
!
      REAL :: erfc_num_recipes
      REAL , INTENT(IN) :: X
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 :: dum , erfc_dbl , t , z
!
! End of declarations rewritten by SPAG
!
      z = abs(X)
      t = 1.0/(1.0+0.5*z)
!       erfc_num_recipes =
!     &   t*exp( -z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
!     &   t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
!     &                                    t*(-1.13520398 +
!     &   t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))
      dum = (-z*z-1.26551223+                                                                                                      &
          & t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+    &
          & t*0.17087277)))))))))
      erfc_dbl = t*exp(dum)
      IF ( X<0.0 ) erfc_dbl = 2.0D0 - erfc_dbl
      erfc_num_recipes = erfc_dbl
 
   END FUNCTION erfc_num_recipes
!-----------------------------------------------------------------------
 
   SUBROUTINE miedriver(Dp_wet_a,Dp_core_a,Ri_shell_a,Ri_core_a,Vlambc,Qextc,Qscatc,Gscac,Extc,Scatc,Qbackc,Backc,Pmom)
!
! Dummy argument declarations rewritten by SPAG
!
      REAL*8 , INTENT(IN) :: Dp_wet_a
      REAL*8 , INTENT(IN) :: Dp_core_a
      COMPLEX*16 , INTENT(IN) :: Ri_shell_a
      COMPLEX*16 , INTENT(IN) :: Ri_core_a
      REAL*8 :: Vlambc
      REAL*8 :: Qextc
      REAL*8 :: Qscatc
      REAL*8 :: Gscac
      REAL*8 :: Extc
      REAL*8 :: Scatc
      REAL*8 :: Qbackc
      REAL*8 :: Backc
      REAL*8 , DIMENSION(0:7,1) :: Pmom
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 , DIMENSION(200) :: anglesc , s11 , s12 , s1c , s1r , s2c , s2r , s33 , s34 , sp , spol
      REAL*8 :: coreic , corerc , rgc , rgcmax , rgcmin , rinc , s11n , shelic , shelrc , sigmagc
      INTEGER*4 :: nang , nrgflagc
!
! End of declarations rewritten by SPAG
!
! MOSAIC INPUTS
!	dp_wet_a = diameter (cm) of aerosol
!	dp_core_a = diameter (cm) of the aerosol's core
!	ri_shell_a = refractive index (complex) of shell
!	ri_core_a = refractve index (complex ) of core (usually assumed to be LAC)
!	vlambc = wavelength of calculation (um, convert to cm)
! MOSAIC outputs
! 	qextc = scattering efficiency
!	qscac = scattering efficiency
!	gscac = asymmetry parameter
!	extc = extinction cross section  (cm^2)
!	scac = scattering cross section (cm^2)
! drives concentric sphere program
! /*---------------------------------------------------------------*/
! /* INPUTS:                                                       */
! /*---------------------------------------------------------------*/
 
!   VLAMBc: Wavelength of the radiation
!   NRGFLAGc: Flag to indicate a number density of volume radius
!   RGc: Number (RGN = Rm) or volume (RGV) weighted mean radius of
!        the particle size distribution
!   SIGMAGc: Geometric standard deviation of the distribution
!   SHELRc: Real part of the index of refraction for the shell
!   SHELIc: Imaginary part of the index of refraction for the shell
!   RINc: Inner core radius as a fraction of outer shell radius
!   CORERc: Real part of the index of refraction for the core
!   COREIc: Imaginary part of the index of refraction for the core
!   NANG: Number of scattering angles between 0 and 90 degrees,
!         inclusive
 
! /*---------------------------------------------------------------*/
! /* OUTPUTS:                                                      */
! /*---------------------------------------------------------------*/
 
!   QEXTc: Extinction efficiency of the particle
!   QSCAc: Scattering efficiency of the particle
!   QBACKc: Backscatter efficiency of the particle
!   EXTc: Extinction cross section of the particle
!   SCAc: Scattering cross section of the particle
!   BACKc: Backscatter cross section of the particle
!   GSCA: Asymmetry parameter of the particles phase function
!   ANGLES(NAN): Scattering angles in degrees
!   S1R(NAN): Real part of the amplitude scattering matrix
!   S1C(NAN): Complex part of the amplitude scattering matrix
!   S2R(NAN): Real part of the amplitude scattering matrix
!   S2C(NAN): Complex part of the amplitude scattering matrix
!   S11N: Normalization coefficient of the scattering matrix
!   S11(NAN): S11 scattering coefficients
!   S12(NAN): S12 scattering coefficients
!   S33(NAN): S33 scattering coefficients
!   S34(NAN): S34 scattering coefficients
!   SPOL(NAN): Degree of polarization of unpolarized, incident light
!   SP(NAN): Phase function
!
! NOTE: NAN=2*NANG-1 is the number of scattering angles between
!       0 and 180 degrees, inclusive.
! /*---------------------------------------------------------------*/
!
      nang = 2 ! only one angle
      nrgflagc = 0 ! size distribution
!
      rgc = Dp_wet_a/2.0 ! radius of particle
      rinc = Dp_core_a/Dp_wet_a ! fraction of radius that is the core
      rgcmin = 0.001
      rgcmax = 5.0
      sigmagc = 1.0  ! no particle size dispersion
      shelrc = real(Ri_shell_a)
      shelic = aimag(Ri_shell_a)
      corerc = real(Ri_core_a)
      coreic = aimag(Ri_core_a)
      CALL ackmieparticle(Vlambc,nrgflagc,rgcmin,rgcmax,rgc,sigmagc,shelrc,shelic,rinc,corerc,coreic,nang,Qextc,Qscatc,Qbackc,Extc,&
                        & Scatc,Backc,Gscac,anglesc,s1r,s1c,s2r,s2c,s11n,s11,s12,s33,s34,spol,sp,Pmom)
                                                                           ! jcb
!	write(6,1010)rgc,qextc,qscatc,qscatc/qextc,gscac
99001 FORMAT (5F20.12)
99002 FORMAT (2F12.6)
   END SUBROUTINE miedriver
!
!     /*--------------------------------------------------------*/
!     /* The Toon-Ackerman SUBROUTINE DMIESS for calculating the*/
!     /* scatter off of a coated sphere of some sort.           */
!     /* Toon and Ackerman, Applied Optics, Vol. 20, Pg. 3657   */
!     /*--------------------------------------------------------*/
 
!**********************************************
 
   SUBROUTINE dmiess(Ro,Rfr,Rfi,Thetd,Jx,Qext,Qscat,Ctbrqs,Eltrmx,Pie,Tau,Cstht,Si2tht,Acap,Qbs,It,Ll,R,Re2,Tmag2,Wvno,An,Bn,Ntrm)
!
! COMMON variable declarations rewritten by SPAG
!
      COMPLEX*16 , DIMENSION(3,9000) :: W
      COMMON /warray/ W
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER*4 , INTENT(IN) :: It
      INTEGER*4 , INTENT(IN) :: Ll
      REAL*8 , INTENT(IN) :: Ro
      REAL*8 , INTENT(IN) :: Rfr
      REAL*8 , INTENT(IN) :: Rfi
      REAL*8 , INTENT(INOUT) , DIMENSION(It) :: Thetd
      INTEGER*4 , INTENT(IN) :: Jx
      REAL*8 , INTENT(INOUT) :: Qext
      REAL*8 , INTENT(INOUT) :: Qscat
      REAL*8 , INTENT(INOUT) :: Ctbrqs
      REAL*8 , DIMENSION(4,It,2) :: Eltrmx
      REAL*8 , INTENT(INOUT) , DIMENSION(3,It) :: Pie
      REAL*8 , INTENT(INOUT) , DIMENSION(3,It) :: Tau
      REAL*8 , INTENT(INOUT) , DIMENSION(It) :: Cstht
      REAL*8 , INTENT(INOUT) , DIMENSION(It) :: Si2tht
      COMPLEX*16 , INTENT(INOUT) , DIMENSION(Ll) :: Acap
      REAL*8 , INTENT(OUT) :: Qbs
      REAL*8 , INTENT(IN) :: R
      REAL*8 , INTENT(IN) :: Re2
      REAL*8 , INTENT(IN) :: Tmag2
      REAL*8 , INTENT(IN) :: Wvno
      COMPLEX*16 , INTENT(OUT) , DIMENSION(500) :: An
      COMPLEX*16 , INTENT(OUT) , DIMENSION(500) :: Bn
      INTEGER*4 , INTENT(INOUT) :: Ntrm
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 :: aa , amagp , bb , cc , cosx1 , cosx4 , dd , denom , e2y1 , ey1 , ey1my4 , ey1py4 , ey4 , pig , qbsi , qbsr ,       &
              & realp , rmm , rx , rxp4 , sinx1 , sinx4 , x , x1 , x4 , y1 , y4
      COMPLEX*16 :: dh1 , dh2 , dh4 , dummy , dumsq , fn1 , fn2 , fna , fnap , fnb , fnbp , hstore , k1 , k2 , k3 , p24h21 ,       &
                  & p24h24 , pstore , rc , rf , rrf , rrfx , tc1 , tc2 , wm1
      INTEGER*4 :: iflag , j , k , m , n , nmx1 , nmx2 , nn
      REAL*8 , DIMENSION(5) :: t
      REAL*8 , DIMENSION(4) :: ta
      REAL*8 , DIMENSION(2) :: tb , tc , td , te
      COMPLEX*16 , DIMENSION(8) :: u
      COMPLEX*16 , DIMENSION(2) :: wfn
      COMPLEX*16 , DIMENSION(4) :: z
!
! End of declarations rewritten by SPAG
!
!
! **********************************************************************
!    THIS SUBROUTINE COMPUTES MIE SCATTERING BY A STRATIFIED SPHERE,
!    I.E. A PARTICLE CONSISTING OF A SPHERICAL CORE SURROUNDED BY A
!    SPHERICAL SHELL.  THE BASIC CODE USED WAS THAT DESCRIBED IN THE
!    REPORT: SUBROUTINES FOR COMPUTING THE PARAMETERS OF THE
!    ELECTROMAGNETIC RADIATION SCATTERED BY A SPHERE J.V. DAVE,
!    I B M SCIENTIFIC CENTER, PALO ALTO , CALIFORNIA.
!    REPORT NO. 320 - 3236 .. MAY 1968 .
!
!    THE MODIFICATIONS FOR STRATIFIED SPHERES ARE DESCRIBED IN
!        TOON AND ACKERMAN, APPL. OPTICS, IN PRESS, 1981
!
!    THE PARAMETERS IN THE CALLING STATEMENT ARE DEFINED AS FOLLOWS :
!      RO IS THE OUTER (SHELL) RADIUS;
!      R  IS THE CORE RADIUS;
!      RFR, RFI  ARE THE REAL AND IMAGINARY PARTS OF THE SHELL INDEX
!          OF REFRACTION IN THE FORM (RFR - I* RFI);
!      RE2, TMAG2  ARE THE INDEX PARTS FOR THE CORE;
!          ( WE ASSUME SPACE HAS UNIT INDEX. )
!      THETD(J): ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE INCIDENT
!          AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
!          IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
!          SUPPLEMENTARY VALUE, SEE COMMENTS BELOW ON ELTRMX;
!      JX: TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
!          REQUIRED.  JX SHOULD NOT EXCEED IT UNLESS THE DIMENSIONS
!          STATEMENTS ARE APPROPRIATEDLY MODIFIED;
!
!      THE DEFINITIONS FOR THE FOLLOWING SYMBOLS CAN BE FOUND IN LIGHT
!          SCATTERING BY SMALL PARTICLES, H.C.VAN DE HULST, JOHN WILEY
!          SONS, INC., NEW YORK, 1957.
!      QEXT: EFFIECIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 127.
!      QSCAT: EFFIECINCY FACTOR FOR SCATTERING,V.D. HULST,P.14 127.
!      CTBRQS: AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128
!      ELTRMX(I,J,K): ELEMENTS OF THE TRANSFORMATION MATRIX F,V.D.HULST
!          ,P.34,45 125. I=1: ELEMENT M SUB 2..I=2: ELEMENT M SUB 1..
!          I = 3: ELEMENT S SUB 21.. I = 4: ELEMENT D SUB 21..
!      ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
!          THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
!          OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
!      QBS IS THE BACK SCATTER CROSS SECTION.
!
!      IT: IS THE DIMENSION OF THETD, ELTRMX, CSTHT, PIE, TAU, SI2THT,
!          IT MUST CORRESPOND EXACTLY TO THE SECOND DIMENSION OF ELTRMX.
!      LL: IS THE DIMENSION OF ACAP
!          IN THE ORIGINAL PROGRAM THE DIMENSION OF ACAP WAS 7000.
!          FOR CONSERVING SPACE THIS SHOULD BE NOT MUCH HIGHER THAN
!          THE VALUE, N=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1
!      WVNO: 2*PIE / WAVELENGTH
!
!    ALSO THE SUBROUTINE COMPUTES THE CAPITAL A FUNCTION BY MAKING USE O
!    DOWNWARD RECURRENCE RELATIONSHIP.
!
!      TA(1): REAL PART OF WFN(1).  TA(2): IMAGINARY PART OF WFN(1).
!      TA(3): REAL PART OF WFN(2).  TA(4): IMAGINARY PART OF WFN(2).
!      TB(1): REAL PART OF FNA.     TB(2): IMAGINARY PART OF FNA.
!      TC(1): REAL PART OF FNB.     TC(2): IMAGINARY PART OF FNB.
!      TD(1): REAL PART OF FNAP.    TD(2): IMAGINARY PART OF FNAP.
!      TE(1): REAL PART OF FNBP.    TE(2): IMAGINARY PART OF FNBP.
!      FNAP, FNBP  ARE THE PRECEDING VALUES OF FNA, FNB RESPECTIVELY.
! **********************************************************************
 
!     /*--------------------------------------------------------------*/
!     /* Initially, make all types undefined.                         */
!     /*--------------------------------------------------------------*/
 
!      IMPLICIT UNDEFINED(A-Z)
 
!     /*--------------------------------------------------------*/
!     /* Dimension statements.                                  */
!     /*--------------------------------------------------------*/
 
 
 
 
!     /*--------------------------------------------------------*/
!     /* Variables used in the calculations below.              */
!     /*--------------------------------------------------------*/
 
 
!
! jcb
!
!     /*--------------------------------------------------------*/
!     /* Define the common block.                               */
!     /*--------------------------------------------------------*/
 
      INTEGER :: spag_nextblock_1
                                 ! a,b Mie coefficients, jcb  Hansen and Travis, eqn 2.44
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
 
 
 
!
!      EQUIVALENCE   (FNA,TB(1)),(FNB,TC(1)),(FNAP,TD(1)),(FNBP,TE(1))
!
!   IF THE CORE IS SMALL SCATTERING IS COMPUTED FOR THE SHELL ONLY
!
 
!     /*--------------------------------------------------------*/
!     /* Begin the Mie calculations.                            */
!     /*--------------------------------------------------------*/
            iflag = 1
            Ntrm = 0
             ! jcb
            IF ( R/Ro<1.0D-06 ) iflag = 2
            IF ( Jx>It ) THEN
               WRITE (6,99004)
99004          FORMAT (//10X,'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
               WRITE (6,99003)
               CALL errmsg('DMIESS: 30',.TRUE.)
            ENDIF
            rf = cmplx(Rfr,-Rfi)
            rc = cmplx(Re2,-Tmag2)
            x = Ro*Wvno
            k1 = rc*Wvno
            k2 = rf*Wvno
            k3 = cmplx(Wvno,0.0D0)
            z(1) = k2*Ro
            z(2) = k3*Ro
            z(3) = k1*R
            z(4) = k2*R
            x1 = dreal(z(1))
            y1 = dimag(z(1))
            x4 = dreal(z(4))
            y4 = dimag(z(4))
            rrf = 1.0D0/rf
            rx = 1.0D0/x
            rrfx = rrf*rx
            t(1) = (x**2)*(Rfr**2+Rfi**2)
            t(1) = dsqrt(t(1))
            nmx1 = 1.30D0*t(1)
!
            IF ( nmx1>Ll-1 ) THEN
               WRITE (6,99005)
               CALL errmsg('DMIESS: 32',.TRUE.)
            ENDIF
            nmx2 = t(1)*1.2
            nmx1 = min(nmx1+5,150)
                              ! jcb
            nmx2 = min(nmx2+5,135)
                              ! jcb
!   	write(6,*)x,nmx1,nmx2,ll  ! jcb
!	stop
            IF ( nmx1>150 ) THEN
            ENDIF
!        NMX1 = 150
!        NMX2 = 135
!
            Acap(nmx1+1) = (0.0D0,0.0D0)
            IF ( iflag/=2 ) THEN
               DO n = 1 , 3
                  W(n,nmx1+1) = (0.0D0,0.0D0)
               ENDDO
            ENDIF
            DO n = 1 , nmx1
               nn = nmx1 - n + 1
               Acap(nn) = (nn+1)*rrfx - 1.0D0/((nn+1)*rrfx+Acap(nn+1))
               IF ( iflag/=2 ) THEN
                  DO m = 1 , 3
                     W(m,nn) = (nn+1)/z(m+1) - 1.0D0/((nn+1)/z(m+1)+W(m,nn+1))
                  ENDDO
               ENDIF
            ENDDO
!
            DO j = 1 , Jx
               IF ( Thetd(j)<0.0D0 ) Thetd(j) = dabs(Thetd(j))
               IF ( Thetd(j)<=0.0D0 ) THEN
                  Cstht(j) = 1.0D0
                  Si2tht(j) = 0.0D0
               ELSEIF ( Thetd(j)<90.0D0 ) THEN
                  t(1) = (3.14159265359*Thetd(j))/180.0D0
                  Cstht(j) = dcos(t(1))
                  Si2tht(j) = 1.0D0 - Cstht(j)**2
               ELSEIF ( Thetd(j)>90.0 ) THEN
                  WRITE (6,99002) Thetd(j)
!
99002             FORMAT (10X,' THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN 90.0 DEGREES. IT IS ',E15.4)
                  WRITE (6,99003)
                  CALL errmsg('DMIESS: 34',.TRUE.)
               ELSE
                  Cstht(j) = 0.0D0
                  Si2tht(j) = 1.0D0
               ENDIF
            ENDDO
!
            DO j = 1 , Jx
               Pie(1,j) = 0.0D0
               Pie(2,j) = 1.0D0
               Tau(1,j) = 0.0D0
               Tau(2,j) = Cstht(j)
            ENDDO
!
! INITIALIZATION OF HOMOGENEOUS SPHERE
!
            t(1) = dcos(x)
            t(2) = dsin(x)
            wm1 = cmplx(t(1),-t(2))
            wfn(1) = cmplx(t(2),t(1))
            ta(1) = t(2)
            ta(2) = t(1)
            wfn(2) = rx*wfn(1) - wm1
            ta(3) = dreal(wfn(2))
            ta(4) = dimag(wfn(2))
!
            n = 1
            ! jcb, bug???
            IF ( iflag==2 ) THEN
               tc1 = Acap(1)*rrf + rx
               tc2 = Acap(1)*rf + rx
               fna = (tc1*ta(3)-ta(1))/(tc1*wfn(2)-wfn(1))
               fnb = (tc2*ta(3)-ta(1))/(tc2*wfn(2)-wfn(1))
               tb(1) = dreal(fna)
               tb(2) = dimag(fna)
               tc(1) = dreal(fnb)
               tc(2) = dimag(fnb)
            ELSE
               n = 1
!
! INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
!
               sinx1 = dsin(x1)
               sinx4 = dsin(x4)
               cosx1 = dcos(x1)
               cosx4 = dcos(x4)
               ey1 = dexp(y1)
               e2y1 = ey1*ey1
               ey4 = dexp(y4)
               ey1my4 = dexp(y1-y4)
               ey1py4 = ey1*ey4
               ey1my4 = dexp(y1-y4)
               aa = sinx4*(ey1py4+ey1my4)
               bb = cosx4*(ey1py4-ey1my4)
               cc = sinx1*(e2y1+1.0D0)
               dd = cosx1*(e2y1-1.0D0)
               denom = 1.0D0 + e2y1*(4.0D0*sinx1*sinx1-2.0D0+e2y1)
               realp = (aa*cc+bb*dd)/denom
               amagp = (bb*cc-aa*dd)/denom
               dummy = cmplx(realp,amagp)
               aa = sinx4*sinx4 - 0.5D0
               bb = cosx4*sinx4
               p24h24 = 0.5D0 + cmplx(aa,bb)*ey4*ey4
               aa = sinx1*sinx4 - cosx1*cosx4
               bb = sinx1*cosx4 + cosx1*sinx4
               cc = sinx1*sinx4 + cosx1*cosx4
               dd = -sinx1*cosx4 + cosx1*sinx4
               p24h21 = 0.5D0*cmplx(aa,bb)*ey1*ey4 + 0.5D0*cmplx(cc,dd)*ey1my4
               dh4 = z(4)/(1.0D0+(0.0D0,1.0D0)*z(4)) - 1.0D0/z(4)
               dh1 = z(1)/(1.0D0+(0.0D0,1.0D0)*z(1)) - 1.0D0/z(1)
               dh2 = z(2)/(1.0D0+(0.0D0,1.0D0)*z(2)) - 1.0D0/z(2)
               pstore = (dh4+n/z(4))*(W(3,n)+n/z(4))
               p24h24 = p24h24/pstore
               hstore = (dh1+n/z(1))*(W(3,n)+n/z(4))
               p24h21 = p24h21/hstore
               pstore = (Acap(n)+n/z(1))/(W(3,n)+n/z(4))
               dummy = dummy*pstore
               dumsq = dummy*dummy
!
! NOTE:  THE DEFINITIONS OF U(I) IN THIS PROGRAM ARE NOT THE SAME AS
!        THE USUBI DEFINED IN THE ARTICLE BY TOON AND ACKERMAN.  THE
!        CORRESPONDING TERMS ARE:
!          USUB1 = U(1)                       USUB2 = U(5)
!          USUB3 = U(7)                       USUB4 = DUMSQ
!          USUB5 = U(2)                       USUB6 = U(3)
!          USUB7 = U(6)                       USUB8 = U(4)
!          RATIO OF SPHERICAL BESSEL FTN TO SPHERICAL HENKAL FTN = U(8)
!
               u(1) = k3*Acap(n) - k2*W(1,n)
               u(2) = k3*Acap(n) - k2*dh2
               u(3) = k2*Acap(n) - k3*W(1,n)
               u(4) = k2*Acap(n) - k3*dh2
               u(5) = k1*W(3,n) - k2*W(2,n)
               u(6) = k2*W(3,n) - k1*W(2,n)
               u(7) = (0.0D0,-1.0D0)*(dummy*p24h21-p24h24)
               u(8) = ta(3)/wfn(2)
!
               fna = u(8)*(u(1)*u(5)*u(7)+k1*u(1)-dumsq*k3*u(5))/(u(2)*u(5)*u(7)+k1*u(2)-dumsq*k3*u(5))
               fnb = u(8)*(u(3)*u(6)*u(7)+k2*u(3)-dumsq*k2*u(6))/(u(4)*u(6)*u(7)+k2*u(4)-dumsq*k2*u(6))
!
!  Explicit equivalences added by J. Francis
 
               tb(1) = dreal(fna)
               tb(2) = dimag(fna)
               tc(1) = dreal(fnb)
               tc(2) = dimag(fnb)
            ENDIF
!
! jcb
            Ntrm = Ntrm + 1
            An(n) = fna
            Bn(n) = fnb
!	write(6,1010)ntrm,n,an(n),bn(n)
99001       FORMAT (2I5,4E15.6)
! jcb
            fnap = fna
            fnbp = fnb
            td(1) = dreal(fnap)
            td(2) = dimag(fnap)
            te(1) = dreal(fnbp)
            te(2) = dimag(fnbp)
            t(1) = 1.50D0
!
!    FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
!    FOLLOWING MEANING:
!    ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
!    ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
!    ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
!    ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
!    K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
!    DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
!
            tb(1) = t(1)*tb(1)
            tb(2) = t(1)*tb(2)
            tc(1) = t(1)*tc(1)
            tc(2) = t(1)*tc(2)
!      DO 60 J = 1,JX
!          ELTRMX(1,J,1) = TB(1) * PIE(2,J) + TC(1) * TAU(2,J)
!          ELTRMX(2,J,1) = TB(2) * PIE(2,J) + TC(2) * TAU(2,J)
!          ELTRMX(3,J,1) = TC(1) * PIE(2,J) + TB(1) * TAU(2,J)
!          ELTRMX(4,J,1) = TC(2) * PIE(2,J) + TB(2) * TAU(2,J)
!          ELTRMX(1,J,2) = TB(1) * PIE(2,J) - TC(1) * TAU(2,J)
!          ELTRMX(2,J,2) = TB(2) * PIE(2,J) - TC(2) * TAU(2,J)
!          ELTRMX(3,J,2) = TC(1) * PIE(2,J) - TB(1) * TAU(2,J)
!          ELTRMX(4,J,2) = TC(2) * PIE(2,J) - TB(2) * TAU(2,J)
!
            Qext = 2.0D0*(tb(1)+tc(1))
            Qscat = (tb(1)**2+tb(2)**2+tc(1)**2+tc(2)**2)/0.75D0
            Ctbrqs = 0.0D0
            qbsr = -2.0D0*(tc(1)-tb(1))
            qbsi = -2.0D0*(tc(2)-tb(2))
            rmm = -1.0D0
            n = 2
            spag_nextblock_1 = 2
         CASE (2)
            t(1) = 2*n - 1
                       ! start of loop, JCB
            t(2) = n - 1
            t(3) = 2*n + 1
            DO j = 1 , Jx
               Pie(3,j) = (t(1)*Pie(2,j)*Cstht(j)-n*Pie(1,j))/t(2)
               Tau(3,j) = Cstht(j)*(Pie(3,j)-Pie(1,j)) - t(1)*Si2tht(j)*Pie(2,j) + Tau(1,j)
            ENDDO
!
! HERE SET UP HOMOGENEOUS SPHERE
!
            wm1 = wfn(1)
            wfn(1) = wfn(2)
            ta(1) = dreal(wfn(1))
            ta(2) = dimag(wfn(1))
            wfn(2) = t(1)*rx*wfn(1) - wm1
            ta(3) = dreal(wfn(2))
            ta(4) = dimag(wfn(2))
!
            IF ( iflag/=2 ) THEN
!
! HERE SET UP STRATIFIED SPHERE
!
               dh2 = -n/z(2) + 1.0D0/(n/z(2)-dh2)
               dh4 = -n/z(4) + 1.0D0/(n/z(4)-dh4)
               dh1 = -n/z(1) + 1.0D0/(n/z(1)-dh1)
               pstore = (dh4+n/z(4))*(W(3,n)+n/z(4))
               p24h24 = p24h24/pstore
               hstore = (dh1+n/z(1))*(W(3,n)+n/z(4))
               p24h21 = p24h21/hstore
               pstore = (Acap(n)+n/z(1))/(W(3,n)+n/z(4))
               dummy = dummy*pstore
               dumsq = dummy*dummy
!
               u(1) = k3*Acap(n) - k2*W(1,n)
               u(2) = k3*Acap(n) - k2*dh2
               u(3) = k2*Acap(n) - k3*W(1,n)
               u(4) = k2*Acap(n) - k3*dh2
               u(5) = k1*W(3,n) - k2*W(2,n)
               u(6) = k2*W(3,n) - k1*W(2,n)
               u(7) = (0.0D0,-1.0D0)*(dummy*p24h21-p24h24)
               u(8) = ta(3)/wfn(2)
!
               fna = u(8)*(u(1)*u(5)*u(7)+k1*u(1)-dumsq*k3*u(5))/(u(2)*u(5)*u(7)+k1*u(2)-dumsq*k3*u(5))
               fnb = u(8)*(u(3)*u(6)*u(7)+k2*u(3)-dumsq*k2*u(6))/(u(4)*u(6)*u(7)+k2*u(4)-dumsq*k2*u(6))
               tb(1) = dreal(fna)
               tb(2) = dimag(fna)
               tc(1) = dreal(fnb)
               tc(2) = dimag(fnb)
            ENDIF
!
            tc1 = Acap(n)*rrf + n*rx
            tc2 = Acap(n)*rf + n*rx
            fn1 = (tc1*ta(3)-ta(1))/(tc1*wfn(2)-wfn(1))
            fn2 = (tc2*ta(3)-ta(1))/(tc2*wfn(2)-wfn(1))
            m = Wvno*R
            IF ( n>=m ) THEN
               IF ( iflag/=2 ) THEN
                  IF ( abs((fn1-fna)/fn1)<1.0D-09 .AND. abs((fn2-fnb)/fn2)<1.0D-09 ) iflag = 2
                  IF ( iflag==1 ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ENDIF
               fna = fn1
               fnb = fn2
               tb(1) = dreal(fna)
               tb(2) = dimag(fna)
               tc(1) = dreal(fnb)
               tc(2) = dimag(fnb)
            ENDIF
            spag_nextblock_1 = 3
         CASE (3)
!
! jcb
            Ntrm = Ntrm + 1
            An(n) = fna
            Bn(n) = fnb
!	write(6,1010)ntrm,n,an(n),bn(n)
! jcb
            t(5) = n
            t(4) = t(1)/(t(5)*t(2))
            t(2) = (t(2)*(t(5)+1.0D0))/t(5)
!
            Ctbrqs = Ctbrqs + t(2)*(td(1)*tb(1)+td(2)*tb(2)+te(1)*tc(1)+te(2)*tc(2)) + t(4)*(td(1)*te(1)+td(2)*te(2))
            Qext = Qext + t(3)*(tb(1)+tc(1))
            t(4) = tb(1)**2 + tb(2)**2 + tc(1)**2 + tc(2)**2
            Qscat = Qscat + t(3)*t(4)
            rmm = -rmm
            qbsr = qbsr + t(3)*rmm*(tc(1)-tb(1))
            qbsi = qbsi + t(3)*rmm*(tc(2)-tb(2))
!
            t(2) = n*(n+1)
            t(1) = t(3)/t(2)
            k = (n/2)*2
!      DO 80 J = 1,JX
!       ELTRMX(1,J,1)=ELTRMX(1,J,1)+T(1)*(TB(1)*PIE(3,J)+TC(1)*TAU(3,J))
!       ELTRMX(2,J,1)=ELTRMX(2,J,1)+T(1)*(TB(2)*PIE(3,J)+TC(2)*TAU(3,J))
!       ELTRMX(3,J,1)=ELTRMX(3,J,1)+T(1)*(TC(1)*PIE(3,J)+TB(1)*TAU(3,J))
!       ELTRMX(4,J,1)=ELTRMX(4,J,1)+T(1)*(TC(2)*PIE(3,J)+TB(2)*TAU(3,J))
!      IF ( K .EQ. N )  THEN
!       ELTRMX(1,J,2)=ELTRMX(1,J,2)+T(1)*(-TB(1)*PIE(3,J)+TC(1)*TAU(3,J))
!       ELTRMX(2,J,2)=ELTRMX(2,J,2)+T(1)*(-TB(2)*PIE(3,J)+TC(2)*TAU(3,J))
!       ELTRMX(3,J,2)=ELTRMX(3,J,2)+T(1)*(-TC(1)*PIE(3,J)+TB(1)*TAU(3,J))
!       ELTRMX(4,J,2)=ELTRMX(4,J,2)+T(1)*(-TC(2)*PIE(3,J)+TB(2)*TAU(3,J))
!      ELSE
!       ELTRMX(1,J,2)=ELTRMX(1,J,2)+T(1)*(TB(1)*PIE(3,J)-TC(1)*TAU(3,J))
!       ELTRMX(2,J,2)=ELTRMX(2,J,2)+T(1)*(TB(2)*PIE(3,J)-TC(2)*TAU(3,J))
!       ELTRMX(3,J,2)=ELTRMX(3,J,2)+T(1)*(TC(1)*PIE(3,J)-TB(1)*TAU(3,J))
!       ELTRMX(4,J,2)=ELTRMX(4,J,2)+T(1)*(TC(2)*PIE(3,J)-TB(2)*TAU(3,J))
!      END IF
!   80 CONTINUE
!
!      IF ( T(4) .LT. 1.0D-14 )   GO TO 100         ! bail out of loop
            IF ( t(4)>=1.0D-10 .AND. n<nmx2 ) THEN                 ! bail out of loop, JCB
               n = n + 1
!      DO 90 J = 1,JX
!         PIE(1,J)  =   PIE(2,J)
!         PIE(2,J)  =   PIE(3,J)
!         TAU(1,J)  =  TAU(2,J)
!         TAU(2,J)  =  TAU(3,J)
               fnap = fna
               fnbp = fnb
               td(1) = dreal(fnap)
               td(2) = dimag(fnap)
               te(1) = dreal(fnbp)
               te(2) = dimag(fnbp)
               IF ( n<=nmx2 ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               WRITE (6,99005)
               CALL errmsg('DMIESS: 36',.TRUE.)
            ENDIF
!      DO 120 J = 1,JX
!      DO 120 K = 1,2
!         DO  115  I= 1,4
!         T(I)  =  ELTRMX(I,J,K)
!  115    CONTINUE
!         ELTRMX(2,J,K)  =      T(1)**2  +  T(2)**2
!         ELTRMX(1,J,K)  =      T(3)**2  +  T(4)**2
!         ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
!         ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
!  120 CONTINUE
            t(1) = 2.0D0*rx**2
            Qext = Qext*t(1)
            Qscat = Qscat*t(1)
            Ctbrqs = 2.0D0*Ctbrqs*t(1)
!
! QBS IS THE BACK SCATTER CROSS SECTION
!
            pig = dacos(-1.0D0)
            rxp4 = rx*rx/(4.0D0*pig)
            Qbs = rxp4*(qbsr**2+qbsi**2)
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
99003 FORMAT (//10X,'PLEASE READ COMMENTS.'//)
99005 FORMAT (//10X,'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST GET DETAILED OUTPUT AND MODIFY SUBROUTINE'//)
!
   END SUBROUTINE dmiess
!
! /*****************************************************************/
! /* SUBROUTINE ACKMIEPARICLE                                      */
! /*****************************************************************/
 
! THIS PROGRAM COMPUTES SCATTERING PROPERTIES FOR DISTRIBUTIONS OF
! PARTICLES COMPOSED OF A CORE OF ONE MATERIAL AND A SHELL OF ANOTHER.
 
! /*---------------------------------------------------------------*/
! /* INPUTS:                                                       */
! /*---------------------------------------------------------------*/
 
!   VLAMBc: Wavelength of the radiation
!   NRGFLAGc: Flag to indicate a number density of volume radius
!   RGc: Geometric mean radius of the particle distribution
!   SIGMAGc: Geometric standard deviation of the distribution
!   SHELRc: Real part of the index of refraction for the shell
!   SHELIc: Imaginary part of the index of refraction for the shell
!   RINc: Inner core radius as a fraction of outer shell radius
!   CORERc: Real part of the index of refraction for the core
!   COREIc: Imaginary part of the index of refraction for the core
!   NANG: Number of scattering angles between 0 and 90 degrees,
!         inclusive
 
! /*---------------------------------------------------------------*/
! /* OUTPUTS:                                                      */
! /*---------------------------------------------------------------*/
 
!   QEXTc: Extinction efficiency of the particle
!   QSCAc: Scattering efficiency of the particle
!   QBACKc: Backscatter efficiency of the particle
!   EXTc: Extinction cross section of the particle
!   SCAc: Scattering cross section of the particle
!   BACKc: Backscatter cross section of the particle
!   GSCA: Asymmetry parameter of the particle's phase function
!   ANGLES(NAN): Scattering angles in degrees
!   S1R(NAN): Real part of the amplitude scattering matrix
!   S1C(NAN): Complex part of the amplitude scattering matrix
!   S2R(NAN): Real part of the amplitude scattering matrix
!   S2C(NAN): Complex part of the amplitude scattering matrix
!   S11N: Normalization coefficient of the scattering matrix
!   S11(NAN): S11 scattering coefficients
!   S12(NAN): S12 scattering coefficients
!   S33(NAN): S33 scattering coefficients
!   S34(NAN): S34 scattering coefficients
!   SPOL(NAN): Degree of polarization of unpolarized, incident light
!   SP(NAN): Phase function
!
! NOTE: NAN=2*NANG-1 is the number of scattering angles between
!       0 and 180 degrees, inclusive.
! /*---------------------------------------------------------------*/
 
   SUBROUTINE ackmieparticle(Vlambc,Nrgflagc,Rgcmin,Rgcmax,Rgc,Sigmagc,Shelrc,Shelic,Rinc,Corerc,Coreic,Nang,Qextc,Qscatc,Qbackc,  &
                           & Extc,Scatc,Backc,Gscac,Anglesc,S1r,S1c,S2r,S2c,S11n,S11,S12,S33,S34,Spol,Sp,Pmom)
                                                                           ! jcb
 
!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/
 
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER*4 , PARAMETER :: MXNANG = 501
!
! Dummy argument declarations rewritten by SPAG
!
      REAL*8 , INTENT(IN) :: Vlambc
      INTEGER*4 :: Nrgflagc
      REAL*8 , INTENT(IN) :: Rgcmin
      REAL*8 , INTENT(IN) :: Rgcmax
      REAL*8 , INTENT(IN) :: Rgc
      REAL*8 , INTENT(IN) :: Sigmagc
      REAL*8 , INTENT(IN) :: Shelrc
      REAL*8 , INTENT(IN) :: Shelic
      REAL*8 , INTENT(IN) :: Rinc
      REAL*8 , INTENT(IN) :: Corerc
      REAL*8 , INTENT(IN) :: Coreic
      INTEGER*4 , INTENT(IN) :: Nang
      REAL*8 , INTENT(OUT) :: Qextc
      REAL*8 , INTENT(OUT) :: Qscatc
      REAL*8 , INTENT(OUT) :: Qbackc
      REAL*8 , INTENT(OUT) :: Extc
      REAL*8 , INTENT(OUT) :: Scatc
      REAL*8 , INTENT(OUT) :: Backc
      REAL*8 , INTENT(OUT) :: Gscac
      REAL*8 , DIMENSION(*) :: Anglesc
      REAL*8 , DIMENSION(*) :: S1r
      REAL*8 , DIMENSION(*) :: S1c
      REAL*8 , DIMENSION(*) :: S2r
      REAL*8 , DIMENSION(*) :: S2c
      REAL*8 :: S11n
      REAL*8 , DIMENSION(*) :: S11
      REAL*8 , DIMENSION(*) :: S12
      REAL*8 , DIMENSION(*) :: S33
      REAL*8 , DIMENSION(*) :: S34
      REAL*8 , DIMENSION(*) :: Spol
      REAL*8 , DIMENSION(*) :: Sp
      REAL*8 , DIMENSION(0:7,1) :: Pmom
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 :: alamb , rfic , rfis , rfrc , rfrs , rgcfrac , rgmax , rgmin , rgv , sigmag
      COMPLEX*16 , DIMENSION(500) :: an , bn
      REAL*8 :: asy , bscat , ext , qbs , qext , qscat , scat
      REAL*8 , DIMENSION(2*MXNANG-1) :: cosphi , sctphs
      INTEGER*4 :: iphasemie , ipolzn , momdim , nmom , nscath , ntrm
!
! End of declarations rewritten by SPAG
!
 
!     /*--------------------------------------------------------*/
!     /* Parameter statements.                                  */
!     /*--------------------------------------------------------*/
 
 
 
!     /*--------------------------------------------------------*/
!     /* Dimension statements.                                  */
!     /*--------------------------------------------------------*/
 
 
!     /*--------------------------------------------------------*/
!     /* Define the types of the common block.                  */
!     /*--------------------------------------------------------*/
 
 
 
! for calculating the Legendre coefficient, jcb
 
!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/
 
!      IMPLICIT REAL*8 (A-H, O-Z)
 
!     /*--------------------------------------------------------*/
!     /* Input common block for scattering calculations.        */
!     /*--------------------------------------------------------*/
 
!jdf  COMMON / PHASE  / IPHASEmie
 
!jdf  COMMON / INPUTS / ALAMB, RGmin, RGmax, RGV, SIGMAG, &
!jdf                    RGCFRAC, RFRS,  RFIS,  RFRC,  RFIC
 
!     /*--------------------------------------------------------*/
!     /* Output common block for scattering calculations.       */
!     /*--------------------------------------------------------*/
 
!jdf  COMMON / OUTPUTS / QEXT, QSCAT, QBS, EXT, SCAT, BSCAT, ASY
 
!     /*--------------------------------------------------------*/
!     /* Arrays to hold the results of the scattering and       */
!     /* moment calculations.                                   */
!     /*--------------------------------------------------------*/
 
 
 
!     /*--------------------------------------------------------*/
!     /* Copy the input parameters into the common block INPUTS */
!     /*--------------------------------------------------------*/
 
      nscath = Nang                ! a,b Mie coefficients, jcb  Hansen and Travis, eqn 2.44
 
      iphasemie = 0
      alamb = Vlambc
      rgmin = Rgcmin
      rgmax = Rgcmax
      rgv = Rgc
      sigmag = Sigmagc
      rgcfrac = Rinc
      rfrs = Shelrc
      rfis = Shelic
      rfrc = Corerc
      rfic = Coreic
 
!     /*--------------------------------------------------------*/
!     /* Calculate the particle scattering properties for the   */
!     /* given wavelength, particle distribution and indices of */
!     /* refraction of inner and outer material.                */
!     /*--------------------------------------------------------*/
 
                                                                                 ! ! jcb
                                                                                ! ! jdf
                                                                                ! ! jdf
      CALL pfcnparticle(nscath,cosphi,sctphs,Anglesc,S1r,S1c,S2r,S2c,S11n,S11,S12,S33,S34,Spol,Sp,an,bn,ntrm,alamb,rgmin,rgmax,rgv,&
                      & sigmag,rgcfrac,rfrs,rfis,rfrc,rfic,qext,qscat,qbs,ext,scat,bscat,asy,iphasemie)
                                                                                  ! jdf
 
!     /*--------------------------------------------------------*/
!     /* If IPHASE = 1, then the full phase function was        */
!     /* calculated; now, go calculate its moments.             */
!     /*--------------------------------------------------------*/
 
!        IF (IPHASE .gt. 0) CALL DISMOM (NSCATA,COSPHI,SCTPHS,RMOMS)
 
!     /*--------------------------------------------------------*/
!     /* Copy the variables in the common block OUTPUTS to the  */
!     /* variable addresses passed into this routine.           */
!     /*--------------------------------------------------------*/
 
      Qextc = qext
      Qscatc = qscat
      Qbackc = qbs
 
      Extc = ext
      Scatc = scat
      Backc = bscat
 
      Gscac = asy
! jcb
!	ntrmj = number of terms in Mie series, jcb
      nmom = 7  ! largest Legendre coefficient to calculate 0:7  (8 total), jcb
      ipolzn = 0 ! unpolarized light, jcb
      momdim = 7 ! dimension of pmom, pmom(0:7), jcb
! a, b = Mie coefficients
! pmom = output of Legendre coefficients, pmom(0:7)
!	write(6,*)ntrm
!	do ii=1,ntrm
!	write(6,1030)ii,an(ii),bn(ii)
99001 FORMAT (i5,4E15.6)
!	enddo
 
      CALL lpcoefjcb(ntrm,nmom,ipolzn,momdim,an,bn,Pmom)
!	do ii=0,7
!	write(6,1040)ii,pmom(ii,1),pmom(ii,1)/pmom(0,1)
!1040	format(i5,2e15.6)
!	enddo
!     /*--------------------------------------------------------*/
!     /* FORMAT statements.                                     */
!     /*--------------------------------------------------------*/
 
99002 FORMAT (///,1X,I6,' IS AN INVALID MEAN RADIUS FLAG')
 
!     /*--------------------------------------------------------*/
!     /* DONE with this subroutine so exit.                     */
!     /*--------------------------------------------------------*/
 
   END SUBROUTINE ackmieparticle
 
! /*****************************************************************/
! /* SUBROUTINE PFCNPARTICLE                                       */
! /*****************************************************************/
!
! THIS SUBROUTINE COMPUTES THE PHASE FUNCTION SCTPHS(I) AT NSCATA
! ANGLES BETWEEN 0.0 AND 180.0 DEGREES SPECIFIED BY COSPHI(I) WHICH
! CONTAINS THE ANGLE COSINES. THESE VALUES ARE RETURNED TO FUNCTION
! PHASFN FOR AZIMUTHAL AVERAGING.
! INPUT DATA FOR THIS ROUTINE IS PASSED THROUGH COMMON /SIZDIS/
! AND CONSISTS OF THE FOLLOWING PARAMETERS
! NEWSD  = 1 IF SIZE DIS VALUES HAVE NOT PREVIOUSLY BEEN USED IN
!          THIS ROUTINE,  = 0 OTHERWISE.
! RGV    = GEOMETRIC MEAN RADIUS FOR THE VOLUME DISTRIBUTION OF THE
!          SPECIFIED PARTICLES
! SIGMAG = GEOMETRIC STANDARD DEVIATION
! RFR,I  = REAL AND IMAGINARY INDEX OF REFRACTION OF PARTICLES
! ALAMB  = WAVELENGTH AT WHICH CALCULATIONS ARE TO BE PERFORMED
!
! /*---------------------------------------------------------------*/
 
                                                                           ! ! jcb
                                                                           ! ! jdf
                                                                           ! ! jdf
   SUBROUTINE pfcnparticle(Nscath,Cosphi,Sctphs,Anglesc,S1r,S1c,S2r,S2c,S11n,S11,S12,S33,S34,Spol,Sp,An,Bn,Ntrm,Alamb,Rgmin,Rgmax, &
                         & Rgv,Sigmag,Rgcfrac,Rfrs,Rfis,Rfrc,Rfic,Qext,Qscat,Qbs,Ext,Scat,Bscat,Asy,Iphasemie)
                                                                             ! jdf
 
!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/
 
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER*4 , PARAMETER :: MXNANG = 501 , MXNWORK = 500000
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER*4 , INTENT(INOUT) :: Nscath
      REAL*8 , INTENT(IN) , DIMENSION(2*MXNANG-1) :: Cosphi
      REAL*8 , DIMENSION(2*MXNANG-1) :: Sctphs
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: Anglesc
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S1r
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S1c
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S2r
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S2c
      REAL*8 :: S11n
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S11
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S12
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: S33
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: S34
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: Spol
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: Sp
      COMPLEX*16 , DIMENSION(500) :: An
      COMPLEX*16 , DIMENSION(500) :: Bn
      INTEGER*4 :: Ntrm
      REAL*8 , INTENT(IN) :: Alamb
      REAL*8 :: Rgmin
      REAL*8 :: Rgmax
      REAL*8 , INTENT(IN) :: Rgv
      REAL*8 :: Sigmag
      REAL*8 , INTENT(IN) :: Rgcfrac
      REAL*8 , INTENT(IN) :: Rfrs
      REAL*8 , INTENT(IN) :: Rfis
      REAL*8 , INTENT(IN) :: Rfrc
      REAL*8 , INTENT(IN) :: Rfic
      REAL*8 , INTENT(OUT) :: Qext
      REAL*8 , INTENT(OUT) :: Qscat
      REAL*8 , INTENT(OUT) :: Qbs
      REAL*8 , INTENT(OUT) :: Ext
      REAL*8 , INTENT(INOUT) :: Scat
      REAL*8 , INTENT(OUT) :: Bscat
      REAL*8 , INTENT(OUT) :: Asy
      INTEGER*4 , INTENT(IN) :: Iphasemie
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX*16 , DIMENSION(MXNWORK) :: acap
      REAL*8 , DIMENSION(MXNANG) :: cstht , si2tht , theta
      REAL*8 :: ctbrqs , dqbs , dqext , dqscat , rfii , rfio , rfri , rfro , rin , rout , wnum
      REAL*8 , DIMENSION(4,MXNANG,2) :: eltrmx
      INTEGER*4 :: it , it2 , jx , ll
      INTEGER :: j , jj , nindex , nscata
      REAL*8 :: pie , x
      REAL*8 , DIMENSION(3,MXNANG) :: pii , tau
!
! End of declarations rewritten by SPAG
!
 
!     /*--------------------------------------------------------*/
!     /* Parameter statements.                                  */
!     /*--------------------------------------------------------*/
 
 
 
!     /*--------------------------------------------------------*/
!     /* Dimension statements.                                  */
!     /*--------------------------------------------------------*/
 
 
!     /*--------------------------------------------------------*/
!     /* Define the types of the common block.                  */
!     /*--------------------------------------------------------*/
 
 
 
 
!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/
 
!      IMPLICIT REAL*8 (A-H, O-Z)
 
!     /*--------------------------------------------------------*/
!     /* Input common block for scattering calculations.        */
!     /*--------------------------------------------------------*/
 
!jdf  COMMON / PHASE  / IPHASEmie
 
!jdf  COMMON / INPUTS / ALAMB, RGmin, RGmax, RGV, SIGMAG, &
!jdf                   RGCFRAC, RFRS,  RFIS,  RFRC,  RFIC
 
!     /*--------------------------------------------------------*/
!     /* Output common block for scattering calculations.       */
!     /*--------------------------------------------------------*/
 
!jdf  COMMON / OUTPUTS / QEXT, QSCAT, QBS, EXT, SCAT, BSCAT, ASY
 
!     /*--------------------------------------------------------*/
!     /* Arrays to perform the scattering calculations and to   */
!     /* hold the subsequent results.                           */
!     /*--------------------------------------------------------*/
 
 
 
 
 
 
!     /*--------------------------------------------------------*/
!     /* Obvious variable initializations.                      */
!     /*--------------------------------------------------------*/
 
      pie = dacos(-1.0D0)
 
!     /*--------------------------------------------------------*/
!     /* Maximum number of scattering angles between 0 and 90   */
!     /* degrees, inclusive.                                    */
!     /*--------------------------------------------------------*/
 
      it = MXNANG
 
!     /*--------------------------------------------------------*/
!     /* Maximum number of scattering angles between 0 and 180  */
!     /* degrees, inclusive.                                    */
!     /*--------------------------------------------------------*/
 
      it2 = 2*it - 1
 
!     /*--------------------------------------------------------*/
!     /* Dimension of the work array ACAP.                      */
!     /*--------------------------------------------------------*/
 
      ll = MXNWORK
 
!     /*--------------------------------------------------------*/
!     /* NSCATA is the actual user-requested number of          */
!     /* scattering angles between 0 and 90 degrees, inclusive. */
!     /*--------------------------------------------------------*/
 
      nscata = 2*Nscath - 1
 
!     /*--------------------------------------------------------*/
!     /* If the user did not request a phase function, then we  */
!     /* can set NSCATA and NSCATH to 0.                        */
!     /*--------------------------------------------------------*/
 
      IF ( Iphasemie<=0 ) THEN
         Nscath = 0
         nscata = 0
      ENDIF
 
!     /*--------------------------------------------------------*/
!     /* Check to make sure that the user-requested number of   */
!     /* scattering angles does not excede the current maximum  */
!     /* limit.                                                 */
!     /*--------------------------------------------------------*/
 
      IF ( nscata>it2 .OR. Nscath>it ) THEN
         WRITE (6,99002) nscata , Nscath , it2 , it
99002    FORMAT (///,1X,'NUMBER OF ANGLES SPECIFIED =',2I6,/10X,'EXCEEDS ARRAY DIMENSIONS =',2I6)
         CALL errmsg('PFCNPARTICLE: 11',.TRUE.)
      ENDIF
 
!     /*--------------------------------------------------------*/
!     /* Subroutine SCATANGLES was added by EEC[0495] in order  */
!     /* to facilitate changing the scattering angle locations  */
!     /* output by the Ackerman and Toon Mie code.              */
!     /*--------------------------------------------------------*/
 
!         CALL SCATANGLES(NSCATH,THETA,COSPHI)
 
!     /*--------------------------------------------------------*/
!     /* COMPUTE SCATTERING PROPERTIES OF THE PARTICLE.         */
!     /*--------------------------------------------------------*/
 
!     /*--------------------------------------------------------*/
!     /* DMIESS expects a wavenumber.                           */
!     /*--------------------------------------------------------*/
 
      wnum = (2.D0*pie)/Alamb
 
!     /*--------------------------------------------------------*/
!     /* DMIESS assignments of the indices of refraction of the */
!     /* core and shell materials.                              */
!     /*--------------------------------------------------------*/
 
      rfro = Rfrs
      rfio = Rfis
      rfri = Rfrc
      rfii = Rfic
 
!     /*--------------------------------------------------------*/
!     /* DMIESS core and shell radii.                           */
!     /*--------------------------------------------------------*/
 
      rout = Rgv
      rin = Rgcfrac*rout
 
!     /*--------------------------------------------------------*/
!     /* Scattering angles are symmetric about 90 degrees.      */
!     /*--------------------------------------------------------*/
 
      IF ( Nscath==0.0 ) THEN
         jx = 1
      ELSE
         jx = Nscath
      ENDIF
 
!     /*--------------------------------------------------------*/
!     /* Compute the scattering properties for this particle.   */
!     /*--------------------------------------------------------*/
 
      CALL dmiess(rout,rfro,rfio,theta,jx,dqext,dqscat,ctbrqs,eltrmx,pii,tau,cstht,si2tht,acap,dqbs,it,ll,rin,rfri,rfii,wnum,An,Bn,&
                & Ntrm)                                                        ! jcb
 
!     /*--------------------------------------------------------*/
!     /* Compute total cross-sectional area of the particle.    */
!     /*--------------------------------------------------------*/
 
      x = pie*Rgv*Rgv
 
!     /*--------------------------------------------------------*/
!     /* Assign the final extinction efficiency.                */
!     /*--------------------------------------------------------*/
 
      Qext = dqext
 
!     /*--------------------------------------------------------*/
!     /* Compute total extinction cross-section due to particle.*/
!     /*--------------------------------------------------------*/
 
      Ext = dqext*x
 
!     /*--------------------------------------------------------*/
!     /* Assign the final scattering efficiency.                */
!     /*--------------------------------------------------------*/
 
      Qscat = dqscat
 
!     /*--------------------------------------------------------*/
!     /* Compute total scattering cross-section due to particle.*/
!     /*--------------------------------------------------------*/
 
      Scat = dqscat*x
 
!     /*--------------------------------------------------------*/
!     /* Assign the final backscatter efficiency.               */
!     /*--------------------------------------------------------*/
 
      Qbs = dqbs
 
!     /*--------------------------------------------------------*/
!     /* Compute backscatter due to particle.                   */
!     /*--------------------------------------------------------*/
 
      Bscat = dqbs*x
 
!     /*--------------------------------------------------------*/
!     /* Compute asymmetry parameter due to particle.           */
!     /*--------------------------------------------------------*/
 
      Asy = (ctbrqs*x)/Scat
 
!     /*--------------------------------------------------------*/
!     /* If IPHASE is 1, compute the phase function.            */
!     /* S33 and S34 matrix elements are normalized by S11. S11 */
!     /* is normalized to 1.0 in the forward direction.  The    */
!     /* variable SPOL is the degree of polarization for        */
!     /* incident unpolarized light.                            */
!     /*--------------------------------------------------------*/
 
      IF ( Iphasemie>0 ) THEN
 
         DO j = 1 , nscata
 
            IF ( j<=jx ) THEN
               jj = j
               nindex = 1
            ELSE
               jj = nscata - j + 1
               nindex = 2
            ENDIF
 
            Anglesc(j) = Cosphi(j)
 
            S1r(j) = eltrmx(1,jj,nindex)
            S1c(j) = eltrmx(2,jj,nindex)
            S2r(j) = eltrmx(3,jj,nindex)
            S2c(j) = eltrmx(4,jj,nindex)
 
            S11(j) = 0.5D0*(S1r(j)**2+S1c(j)**2+S2r(j)**2+S2c(j)**2)
            S12(j) = 0.5D0*(S2r(j)**2+S2c(j)**2-S1r(j)**2-S1c(j)**2)
            S33(j) = S2r(j)*S1r(j) + S2c(j)*S1c(j)
            S34(j) = S2r(j)*S1c(j) - S1r(j)*S2c(j)
 
            Spol(j) = -S12(j)/S11(j)
 
            Sp(j) = (4.D0*pie)*(S11(j)/(Scat*wnum**2))
 
         ENDDO
 
!        /*-----------------------------------------------------*/
!        /* DONE with the phase function so exit the IF.        */
!        /*-----------------------------------------------------*/
 
      ENDIF
 
!     /*--------------------------------------------------------*/
!     /* END of the computations so exit the routine.           */
!     /*--------------------------------------------------------*/
 
 
 
!     /*--------------------------------------------------------*/
!     /* FORMAT statements.                                     */
!     /*--------------------------------------------------------*/
 
99001 FORMAT (7X,I3)
 
99003 FORMAT (/10X,'INTEGRATED VOLUME',T40,'=',1PE14.5,/15X,'PERCENT VOLUME IN CORE',T40,'=',0PF10.5,/15X,                         &
            & 'PERCENT VOLUME IN SHELL',T40,'=',0PF10.5,/10X,'INTEGRATED SURFACE AREA',T40,'=',1PE14.5,/10X,                       &
             &'INTEGRATED NUMBER DENSITY',T40,'=',1PE14.5)
99004 FORMAT (10X,'CORE RADIUS COMPUTED FROM :',/,20X,9A8,/)
 
99005 FORMAT (///,1X,'* * * WARNING * * *',/10X,'PHASE FUNCTION CALCULATION MAY NOT HAVE CONVERGED'/10X,                           &
             &'VALUES OF S1 AT NSDI-1 AND NSDI ARE :',2E14.6,/10X,'VALUE OF X AT NSDI =',E14.6)
 
!     /*--------------------------------------------------------*/
!     /* DONE with this subroutine so exit.                     */
!     /*--------------------------------------------------------*/
 
   END SUBROUTINE pfcnparticle
 
! /*****************************************************************/
! /*****************************************************************/
   SUBROUTINE lpcoefjcb(Ntrm,Nmom,Ipolzn,Momdim,A,B,Pmom)
!
!         calculate legendre polynomial expansion coefficients (also
!         called moments) for phase quantities ( ref. 5 formulation )
!
!     input:  ntrm                    number terms in mie series
!             nmom, ipolzn, momdim    'miev0' arguments
!             a, b                    mie series coefficients
!
!     output: pmom                   legendre moments ('miev0' argument)
!
!     *** notes ***
!
!         (1)  eqs. 2-5 are in error in dave, appl. opt. 9,
!         1888 (1970).  eq. 2 refers to m1, not m2;  eq. 3 refers to
!         m2, not m1.  in eqs. 4 and 5, the subscripts on the second
!         term in square brackets should be interchanged.
!
!         (2)  the general-case logic in this subroutine works correctly
!         in the two-term mie series case, but subroutine  'lpco2t'
!         is called instead, for speed.
!
!         (3)  subroutine  'lpco1t', to do the one-term case, is never
!         called within the context of 'miev0', but is included for
!         complete generality.
!
!         (4)  some improvement in speed is obtainable by combining the
!         310- and 410-loops, if moments for both the third and fourth
!         phase quantities are desired, because the third phase quantity
!         is the real part of a complex series, while the fourth phase
!         quantity is the imaginary part of that very same series.  but
!         most users are not interested in the fourth phase quantity,
!         which is related to circular polarization, so the present
!         scheme is usually more efficient.
!
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER :: MAXTRM = 1102 , MAXMOM = 2*MAXTRM , MXMOM2 = MAXMOM/2 , MAXRCP = 4*MAXTRM + 2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Momdim
      INTEGER , INTENT(IN) :: Ntrm
      INTEGER , INTENT(IN) :: Nmom
      INTEGER , INTENT(IN) :: Ipolzn
      COMPLEX*16 , INTENT(IN) , DIMENSION(500) :: A
      COMPLEX*16 , INTENT(IN) , DIMENSION(500) :: B
      REAL*8 , INTENT(INOUT) , DIMENSION(0:Momdim,1) :: Pmom
!
! Local variable declarations rewritten by SPAG
!
      REAL*8 , DIMENSION(0:MAXTRM) :: am
      REAL*8 , DIMENSION(0:MXMOM2) :: bi , bidel
      COMPLEX*16 , DIMENSION(MAXTRM) :: cm , cs , dm , ds
      LOGICAL :: evenl
      INTEGER :: i , idel , imax , k , l , ld2 , m , mmax , nummom
      LOGICAL , SAVE :: pass1
      REAL*8 , DIMENSION(MAXRCP) , SAVE :: recip
      REAL*8 :: sum
!
! End of declarations rewritten by SPAG
!
!
!           ** specification of local variables
!
!      am(m)       numerical coefficients  a-sub-m-super-l
!                     in dave, eqs. 1-15, as simplified in ref. 5.
!
!      bi(i)       numerical coefficients  b-sub-i-super-l
!                     in dave, eqs. 1-15, as simplified in ref. 5.
!
!      bidel(i)    1/2 bi(i) times factor capital-del in dave
!
!      cm,dm()     arrays c and d in dave, eqs. 16-17 (mueller form),
!                     calculated using recurrence derived in ref. 5
!
!      cs,ds()     arrays c and d in ref. 4, eqs. a5-a6 (sekera form),
!                     calculated using recurrence derived in ref. 5
!
!      c,d()       either -cm,dm- or -cs,ds-, depending on -ipolzn-
!
!      evenl       true for even-numbered moments;  false otherwise
!
!      idel        1 + little-del  in dave
!
!      maxtrm      max. no. of terms in mie series
!
!      maxmom      max. no. of non-zero moments
!
!      nummom      number of non-zero moments
!
!      recip(k)    1 / k
!
      DATA pass1/.TRUE./            ! the ",1" dimension is for historical reasons
!
!
      IF ( pass1 ) THEN
!
         DO k = 1 , MAXRCP
            recip(k) = 1.0/k
         ENDDO
         pass1 = .FALSE.
!
      ENDIF
!
      DO l = 0 , Nmom
         Pmom(l,1) = 0.0
      ENDDO
! these will never be called
!      if ( ntrm.eq.1 )  then
!         call  lpco1t ( nmom, ipolzn, momdim, calcmo, a, b, pmom )
!         return
!      else if ( ntrm.eq.2 )  then
!         call  lpco2t ( nmom, ipolzn, momdim, calcmo, a, b, pmom )
!         return
!      end if
!
      IF ( Ntrm+2>MAXTRM ) WRITE (6,99001)
99001 FORMAT (' lpcoef--parameter maxtrm too small')
!
!                                     ** calculate mueller c, d arrays
      cm(Ntrm+2) = (0.,0.)
      dm(Ntrm+2) = (0.,0.)
      cm(Ntrm+1) = (1.-recip(Ntrm+1))*B(Ntrm)
      dm(Ntrm+1) = (1.-recip(Ntrm+1))*A(Ntrm)
      cm(Ntrm) = (recip(Ntrm)+recip(Ntrm+1))*A(Ntrm) + (1.-recip(Ntrm))*B(Ntrm-1)
      dm(Ntrm) = (recip(Ntrm)+recip(Ntrm+1))*B(Ntrm) + (1.-recip(Ntrm))*A(Ntrm-1)
!
      DO k = Ntrm - 1 , 2 , -1
         cm(k) = cm(k+2) - (1.+recip(k+1))*B(k+1) + (recip(k)+recip(k+1))*A(k) + (1.-recip(k))*B(k-1)
         dm(k) = dm(k+2) - (1.+recip(k+1))*A(k+1) + (recip(k)+recip(k+1))*B(k) + (1.-recip(k))*A(k-1)
      ENDDO
      cm(1) = cm(3) + 1.5*(A(1)-B(2))
      dm(1) = dm(3) + 1.5*(B(1)-A(2))
!
      IF ( Ipolzn>=0 ) THEN
!
         DO k = 1 , Ntrm + 2
            cm(k) = (2*k-1)*cm(k)
            dm(k) = (2*k-1)*dm(k)
         ENDDO
!
      ELSE
!                                    ** compute sekera c and d arrays
         cs(Ntrm+2) = (0.,0.)
         ds(Ntrm+2) = (0.,0.)
         cs(Ntrm+1) = (0.,0.)
         ds(Ntrm+1) = (0.,0.)
!
         DO k = Ntrm , 1 , -1
            cs(k) = cs(k+2) + (2*k+1)*(cm(k+1)-B(k))
            ds(k) = ds(k+2) + (2*k+1)*(dm(k+1)-A(k))
         ENDDO
!
         DO k = 1 , Ntrm + 2
            cm(k) = (2*k-1)*cs(k)
            dm(k) = (2*k-1)*ds(k)
         ENDDO
!
      ENDIF
!
!
      IF ( Ipolzn<0 ) nummom = min0(Nmom,2*Ntrm-2)
      IF ( Ipolzn>=0 ) nummom = min0(Nmom,2*Ntrm)
      IF ( nummom>MAXMOM ) WRITE (6,99002)
99002 FORMAT (' lpcoef--parameter maxtrm too small')
!
!                               ** loop over moments
      SPAG_Loop_1_1: DO l = 0 , nummom
         ld2 = l/2
         evenl = mod(l,2)==0
!                                    ** calculate numerical coefficients
!                                    ** a-sub-m and b-sub-i in dave
!                                    ** double-sums for moments
         IF ( l==0 ) THEN
!
            idel = 1
            DO m = 0 , Ntrm
               am(m) = 2.0*recip(2*m+1)
            ENDDO
            bi(0) = 1.0
!
         ELSEIF ( evenl ) THEN
!
            idel = 1
            DO m = ld2 , Ntrm
               am(m) = (1.+recip(2*m-l+1))*am(m)
            ENDDO
            DO i = 0 , ld2 - 1
               bi(i) = (1.-recip(l-2*i))*bi(i)
            ENDDO
            bi(ld2) = (2.-recip(l))*bi(ld2-1)
!
         ELSE
!
            idel = 2
            DO m = ld2 , Ntrm
               am(m) = (1.-recip(2*m+l+2))*am(m)
            ENDDO
            DO i = 0 , ld2
               bi(i) = (1.-recip(l+2*i+1))*bi(i)
            ENDDO
!
         ENDIF
!                                     ** establish upper limits for sums
!                                     ** and incorporate factor capital-
!                                     ** del into b-sub-i
         mmax = Ntrm - idel
         IF ( Ipolzn>=0 ) mmax = mmax + 1
         imax = min0(ld2,mmax-ld2)
         IF ( imax<0 ) EXIT SPAG_Loop_1_1
         DO i = 0 , imax
            bidel(i) = bi(i)
         ENDDO
         IF ( evenl ) bidel(0) = 0.5*bidel(0)
!
!                                    ** perform double sums just for
!                                    ** phase quantities desired by user
         IF ( Ipolzn==0 ) THEN
            DO i = 0 , imax
               sum = 0.0
               DO m = ld2 , mmax - i
                  sum = sum + am(m)*(dble(cm(m-i+1)*dconjg(cm(m+i+idel)))+dble(dm(m-i+1)*dconjg(dm(m+i+idel))))
               ENDDO
               Pmom(l,1) = Pmom(l,1) + bidel(i)*sum
            ENDDO
            Pmom(l,1) = 0.5*Pmom(l,1)
         ENDIF
      ENDDO SPAG_Loop_1_1
   END SUBROUTINE lpcoefjcb
!
END MODULE mpas_mie_module
