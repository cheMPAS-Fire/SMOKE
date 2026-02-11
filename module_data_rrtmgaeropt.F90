MODULE module_data_rrtmgaeropt
 
!*************************************************************
!*************************************************************
!  Original source code copied from WRF-Chem.
!  This code has been written by Minsu Choi, CIRES and NOAA GSL.
!  December, 2025, Minsu.Choi@noaa.gov, Minsu.Choi@colorado.edu
!
!  The code now supports aerosol species implemented in
!  cheMPAS-Fire, developed by the Atmospheric Composition Branch
!  at NOAA GSL.
!
!  A total of four (4) shortwave bands and fourteen (14) longwave bands
!  are used for radiative calculations.
!
!  Aerosol species are mapped to integer identifiers
!  (1, 2, 3, ..., n).
!
!  The file AERO_OPT.TBL is required in the model run directory
!  to read refractive indices for individual species. Users may
!  also specify the number of aerosol species through this table.
!
!  The variables nswbands and nlwbands must be consistent with
!  the definitions in AERO_OPT.TBL and in the source code.
!  Inconsistencies will result in model failure.
!
! Code has been further modernized/re-engineered by SPAG
!*************************************************************
!*************************************************************
   USE mpas_kind_types
   USE module_peg_util , only:peg_error_fatal
   IMPLICIT NONE
 
   INTEGER , PARAMETER :: NSWBANDS = 4
   INTEGER , PARAMETER :: NLWBANDS = 16
 
   TYPE :: aer_opt_type
      CHARACTER(LEN=64) :: name      ! e.g., "WATER", "DUST"
      REAL , DIMENSION(NSWBANDS) :: sw_real , sw_imag
      REAL , DIMENSION(NLWBANDS) :: lw_real , lw_imag
   END TYPE aer_opt_type
 
   TYPE(aer_opt_type) , ALLOCATABLE , SAVE :: aer_optics(:)
   INTEGER , SAVE :: num_aer_species
 
   INTEGER , PARAMETER :: ID_WATER = 1
   INTEGER , PARAMETER :: ID_DUST = 2
   INTEGER , PARAMETER :: ID_SMOKE = 3
      ! INTEGER, PARAMETER :: id_sulfate = 4
 
   REAL , DIMENSION(NSWBANDS) , SAVE :: refrwsw , refiwsw
   REAL , DIMENSION(NLWBANDS) , SAVE :: refrwlw , refiwlw
 
 
!*************************************************************
      !wavelength
   REAL , SAVE :: wavmin(NSWBANDS)   ! Min wavelength (um) of interval
      !data wavmin /3.077,2.500,2.150,1.942,1.626,1.299, &
   DATA wavmin/0.25 , 0.35 , 0.55 , 0.998/
   REAL , SAVE :: wavmax(NSWBANDS)   ! Max wavelength (um) of interval
      !data wavmax/3.846,3.077,2.500,2.150,1.942,1.626, &
   DATA wavmax/0.35 , 0.45 , 0.65 , 1.000/
   REAL , SAVE :: wavenumber1_longwave(NLWBANDS)   !Longwave limits (cm-1)
   DATA wavenumber1_longwave/10. , 350. , 500. , 630. , 700. , 820. , 980. , 1080. , 1180. , 1390. , 1480. , 1800. , 2080. ,       &
      & 2250. , 2390. , 2600./
   REAL , SAVE :: wavenumber2_longwave(NLWBANDS)   !Longwave limits (cm-1)
   DATA wavenumber2_longwave/350. , 500. , 630. , 700. , 820. , 980. , 1080. , 1180. , 1390. , 1480. , 1800. , 2080. , 2250. ,     &
      & 2390. , 2600. , 3250./
 
      !mode or size bin
   INTEGER , PARAMETER :: MAXD_AMODE = 3
   INTEGER , PARAMETER :: NTOT_AMODE = 3
   INTEGER , PARAMETER :: MAXD_BIN = 8
   INTEGER , PARAMETER :: NTOT_BIN = 8
 
      !Chebychev polynomial
      !integer,parameter :: prefr=7,prefi=10
   INTEGER , PARAMETER :: PREFR = 7 , PREFI = 7
   INTEGER , PARAMETER :: NCOEF = 50
   REAL , PARAMETER :: RMMIN = 0.005E-4 , RMMAX = 50.E-4
                                                     ! cm
   REAL , SAVE :: refrtabsw(PREFR,NSWBANDS)
   REAL , SAVE :: refitabsw(PREFI,NSWBANDS)
   REAL , SAVE :: refrtablw(PREFR,NLWBANDS)
   REAL , SAVE :: refitablw(PREFI,NLWBANDS)
      !coefficients for parameterizing aerosol radiative properties
      !in terms of refractive index and wet radius
   REAL , SAVE :: extpsw(NCOEF,PREFR,PREFI,NSWBANDS) !specific extinction
   REAL , SAVE :: abspsw(NCOEF,PREFR,PREFI,NSWBANDS) !specific absorption
   REAL , SAVE :: ascatpsw(NCOEF,PREFR,PREFI,NSWBANDS) !specific scattering
   REAL , SAVE :: asmpsw(NCOEF,PREFR,PREFI,NSWBANDS) !asymmetry factor
   REAL , SAVE :: sbackpsw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: pmom2psw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: pmom3psw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: pmom4psw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: pmom5psw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: pmom6psw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: pmom7psw(NCOEF,PREFR,PREFI,NSWBANDS)
   REAL , SAVE :: extplw(NCOEF,PREFR,PREFI,NLWBANDS) !specific extinction
   REAL , SAVE :: absplw(NCOEF,PREFR,PREFI,NLWBANDS) !specific absorption
   REAL , SAVE :: ascatplw(NCOEF,PREFR,PREFI,NLWBANDS) !specific scattering
   REAL , SAVE :: asmplw(NCOEF,PREFR,PREFI,NLWBANDS) !asymmetry factor
 
   REAL , SAVE :: wavmidsw(NSWBANDS)
   DATA wavmidsw/0.30E-4 , 0.40E-4 , 0.60E-4 , 0.999E-04/
   REAL , SAVE :: wavmidlw(NLWBANDS)
   COMPLEX , SAVE :: crefwsw(NSWBANDS)   ! complex refractive index fro water
   COMPLEX , SAVE :: crefwlw(NLWBANDS)
 
CONTAINS
 
   SUBROUTINE rrtmgaeropt_init()
      !=================================================================
      ! Purpose: Reads AERO_OPT.TBL and populates optical arrays
      !=================================================================
      IMPLICIT NONE
      INTEGER , PARAMETER :: OPEN_OK = 0
      INTEGER :: istat , iunit , k , n , tbl_nlw , tbl_nspecies , tbl_nsw
      CHARACTER(256) :: message
 
 
      OPEN (NEWUNIT=iunit,FILE='AERO_OPT.TBL',STATUS='OLD',ACTION='READ',IOSTAT=istat)
 
      IF ( istat/=OPEN_OK ) CALL peg_error_fatal(6,'failure opening AERO_OPT.TBL')
 
      READ (iunit,*,IOSTAT=istat) tbl_nsw , tbl_nlw , tbl_nspecies
      IF ( istat/=0 ) CALL peg_error_fatal(6,'Error reading header from AERO_OPT.TBL')
 
      IF ( tbl_nsw/=NSWBANDS .OR. tbl_nlw/=NLWBANDS ) THEN
        WRITE (message,*) 'AERO_OPT.TBL mismatch. Code expects:' , NSWBANDS , NLWBANDS , ' File has:' , tbl_nsw , tbl_nlw
        CALL peg_error_fatal(6,trim(message))
      ENDIF
 
      Num_aer_species = tbl_nspecies
      IF ( .NOT.allocated(Aer_optics) ) ALLOCATE (Aer_optics(Num_aer_species))
 
      DO n = 1 , Num_aer_species
 
         READ (iunit,*,IOSTAT=istat) Aer_optics(n)%name
         IF ( istat/=0 ) CALL peg_error_fatal(6,'Error reading species name in TBL')
 
         DO k = 1 , NSWBANDS
            READ (iunit,*) Aer_optics(n)%sw_real(k) , Aer_optics(n)%sw_imag(k)
         ENDDO

         READ (iunit,*)
         DO k = 1 , NLWBANDS
            READ (iunit,*) Aer_optics(n)%lw_real(k) , Aer_optics(n)%lw_imag(k)
         ENDDO
 
      ENDDO
 
      CLOSE (iunit)
 
      IF ( Num_aer_species<ID_WATER ) THEN
        CALL peg_error_fatal(6,'AERO_OPT.TBL Water missing!')
      ENDIF
 
      ! This part for mieaer subroutine in mpas_mie_module
      ! explicit refractive indices for water are need to get bounds for
      ! miev0
      Refrwsw = Aer_optics(ID_WATER)%sw_real
      Refiwsw = Aer_optics(ID_WATER)%sw_imag
 
      Refrwlw = Aer_optics(ID_WATER)%lw_real
      Refiwlw = Aer_optics(ID_WATER)%lw_imag
 
   END SUBROUTINE rrtmgaeropt_init
 
END MODULE module_data_rrtmgaeropt
