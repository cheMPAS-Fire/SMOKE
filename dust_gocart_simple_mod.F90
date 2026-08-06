!>\file  dust_gocart_simple_mod.F90
!! This file contains the GOCART simple dust scheme.

module dust_gocart_simple_mod
!
!  This module developed by Kyoung-Min Kim (NOAA GSL)
!  For serious questions contact kyoungmin.kim@noaa.gov
!
!  07/16/2026 - GOCART is added to MPAS, K.-M. Kim

  use mpas_kind_types
  use dust_data_mod
  use mpas_smoke_init, only: p_dust_fine, p_dust_coarse
  use mpas_timer, only : mpas_timer_start, mpas_timer_stop

  implicit none

  private

  public :: gocart_dust_simple_driver

contains

  subroutine gocart_dust_simple_driver(dt,ktau,              &
       u_phy,v_phy,chem,rho_phy,dz8w,smois,u10,v10,          &
       erod_in,isltyp,xland,area,g,                          &
       e_dust_out, num_e_dust_out,                           &
       ch_dust,                                              &
       index_e_dust_out_dust_fine,                           & 
       index_e_dust_out_dust_coarse,                         &
       num_emis_dust,num_chem,num_soil_layers,num_soil_types,&
       ids,ide, jds,jde, kds,kde,                            &
       ims,ime, jms,jme, kms,kme,                            &
       its,ite, jts,jte, kts,kte                             )

    IMPLICIT NONE

    INTEGER,      INTENT(IN   ) ::                             &
         ids,ide, jds,jde, kds,kde,                            &
         ims,ime, jms,jme, kms,kme,                            &
         its,ite, jts,jte, kts,kte,                            &
         num_emis_dust,num_chem,num_soil_layers,num_soil_types,&
         num_e_dust_out, index_e_dust_out_dust_fine, index_e_dust_out_dust_coarse,ktau

    ! 2d input variables
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: xland     ! dominant land use type
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: area      ! area of grid cell [m2]
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: u10, v10  ! 10m wind speed [m/s]
    REAL(RKIND), INTENT(IN) :: ch_dust

    INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: isltyp  ! soil type

    ! 3d input variables
    REAL(RKIND), DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN) :: rho_phy, dz8w, u_phy, v_phy
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ), INTENT(INOUT) :: chem
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme,1:num_e_dust_out), INTENT(INOUT) :: e_dust_out
    REAL(RKIND), DIMENSION( ims:ime, 1:num_soil_layers, jms:jme ), INTENT(IN) :: smois
    REAL(RKIND), DIMENSION( ims:ime, 1:num_soil_types, jms:jme ), INTENT(IN) :: erod_in

    !0d input variables 
    REAL(RKIND), INTENT(IN) :: dt ! time step
    REAL(RKIND), INTENT(IN) :: g  ! gravity (m/s**2)

    ! Local variables
    integer :: nmx,i,j,k,imx,jmx,lmx
    integer :: ilwi
    real(RKIND), DIMENSION(num_soil_types) :: erod
    real(RKIND) :: airden ! air density
    real(RKIND) :: dz_lowest
    REAL(RKIND) :: w10m, gwet ! 
    real(RKIND) :: dxy
    real(RKIND) :: conver,converi ! conversion values 
    real(RKIND) :: ch_fix
    real(RKIND), DIMENSION (num_emis_dust) :: tc
    real(RKIND), DIMENSION (num_emis_dust) :: bems


    ! conversion values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    imx=1
    jmx=1
    lmx=1
    nmx=5

    ! maching scale from the namelist value
    ch_fix = ch_dust*1.0E-8

    k=kts
    do j=jts,jte
       do i=its,ite

          ! Don't do dust over water!!!
          ilwi=0
          if( (xland(i,j) - 1.5) .lt. 0.)then
             ilwi=1

             ! Total concentration at lowest model level. This is still hardcoded for 5 bins.

             tc(1)=0._RKIND !chem(i,kts,j,p_dust_fine)*conver
             tc(2)=0._RKIND
             tc(3)=0._RKIND
             tc(4)=0._RKIND
             tc(5)=0._RKIND !chem(i,kts,j,p_dust_coarse)*conver
             w10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))

             !
             ! don't trust the u10,v10 values, is model layers are very thin near surface
             !
             if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
             erod(1)=erod_in(i,1,j)
             erod(2)=erod_in(i,2,j)
             erod(3)=erod_in(i,3,j)
             
             !
             !  volumetric soil moisture over porosity
             !
             gwet=smois(i,1,j)/porosity(isltyp(i,j))
             airden=rho_phy(i,kts,j)
             dz_lowest = dz8w(i,1,j)

             ! Call dust emission routine.
             bems(:) = 0._RKIND
             call source_du( imx,jmx,lmx,nmx,num_soil_types, dt, tc,                 &
                            erod, ilwi, w10m, gwet, airden,        &
                            dz_lowest,bems,g,ch_fix)       
             
             chem(i,kts,j,p_dust_fine)  = chem(i,kts,j,p_dust_fine) + (tc(1)+0.38*tc(2))*converi ! 0.38 = ln(2.5/2)/ln(3.6/2)
             chem(i,kts,j,p_dust_coarse)= chem(i,kts,j,p_dust_coarse) + (0.62*tc(2)+tc(3)+0.737*tc(4))*converi ! 0.737 = ln(10/6)/ln(12/6)

             e_dust_out(i,kts,j,index_e_dust_out_dust_fine  ) = bems(1)+0.38*bems(2) !bems(1)
             e_dust_out(i,kts,j,index_e_dust_out_dust_coarse) = 0.62*bems(2)+bems(3)+0.737*bems(4)
          endif
       enddo
    enddo
    !

  end subroutine gocart_dust_simple_driver


  
  
  SUBROUTINE source_du( imx,jmx,lmx,nmx,num_soil_types, dt1, tc, &
                     erod, ilwi, w10m, gwet, airden, &
                     dz_lowest,bems,g0,ch_fix)

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         DT1       Time step                                      (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_FIX    Constant to fudge the total emission of dust  (s2/m2)
! *         dz_lowest heigth of the lowest layer                     (m)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/m2)
! *         BEMS      Source of each dust type           (kg/timestep/m2)
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

  INTEGER,   INTENT(IN)    :: nmx,imx,jmx,lmx,num_soil_types
  REAL(RKIND),    INTENT(IN)    :: erod(num_soil_types)
  INTEGER,   INTENT(IN)    :: ilwi

  REAL(RKIND),    INTENT(IN)    :: w10m, gwet
  REAL(RKIND),    INTENT(IN)    :: airden
  REAL(RKIND),    INTENT(INOUT) :: tc(nmx)
  REAL(RKIND),    INTENT(OUT)   :: bems(nmx) 
  REAL(RKIND),    INTENT(IN  )  :: dz_lowest

  REAL(RKIND)    :: den(nmx), diam(nmx)
  REAL(RKIND)    :: u_ts0, u_ts, dsrc, srce
  REAL(RKIND), intent(in)    :: g0
  REAL(RKIND), intent(in)    :: ch_fix
  REAL(RKIND)    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k

  ! executable statemenst

  DO n = 1, nmx
     ! Threshold velocity as a function of the dust density and the diameter
     ! from Bagnold (1941)
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     DO k = 1, ndsrc
        ! No flux if wet soil 
              rhoa = airden*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 
!             write(0,*)u_ts0,den(n),diam(n),rhoa,g
              
              ! Case of surface dry enough to erode
              IF (gwet < 0.5) THEN  
                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(m)   ! (kg s^2 m^-5)
              IF (ilwi == 1 ) THEN
                                    !(kg s^2 m^-5)*(m^3 s^-3)*s = (kg/m2) per dt1
                 dsrc = ch_fix*srce*w10m**2 * (w10m - u_ts)*dt1 
              ELSE 
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0
              
              ! Update dust mixing ratio at first model level.
              tc(n) = tc(n) + dsrc/dz_lowest/airden  ! (kg/kg)
              bems(n) = dsrc/dt1                   ! diagnostic (kg/m2/s) - kmk
     END DO
  END DO
  
END SUBROUTINE source_du
  
end module dust_gocart_simple_mod
