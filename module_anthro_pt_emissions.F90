!>\file  module_wildfire_smoke_emissions.F90
!! This file contains the MPAS-Aerosols/RRFS wildfire emission module

module module_anthro_pt_emissions
!
!  This module developed by Johana Romero-Alvarez and Jordan Schnell (NOAA GSL)
!  For serious questions contact johana.romero-alvarez@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init
  use mpas_smoke_config, only : pi 
  use module_rwc_emissions, only : plume_rise_briggs_rwc 

  implicit none

  private

  public :: mpas_smoke_anthro_pt_emis_driver

contains


  subroutine mpas_smoke_anthro_pt_emis_driver(dt,gmt,julday,ktau,                    &
                           xlat,xlong,xland, chem,num_chem,dz8w,t_phy,rho_phy,       &
                           z_at_w,zmid,pblh,wind10m,area,                            &
                           e_ant_pt_in,num_anthro_pt,num_e_ant_pt_in,                & 
                           e_ant_stack_groups_in, num_e_ant_stack_groups_in,         & 
                           anthro_pt_emis_scale_factor,                              &
                           index_e_ant_pt_in_unspc_fine,                             &
                           index_STKHT, index_STKDM, index_STKTK, index_STKVE,       &
                           index_STKLT, index_STKLG,                                 &
                           ant_pt_local_cell_idx,ant_pt_rank,myrank,                 &
                           ids,ide, jds,jde, kds,kde,                                &
                           ims,ime, jms,jme, kms,kme,                                &
                           its,ite, jts,jte, kts,kte                                 )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday, num_chem,ktau,             &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_e_ant_pt_in, num_anthro_pt,     &
                                  num_e_ant_stack_groups_in,         &
                                  index_e_ant_pt_in_unspc_fine,   &
                                  index_STKHT, index_STKDM,          &
                                  index_STKTK, index_STKVE,          &
                                  index_STKLT, index_STKLG

   REAL(RKIND), INTENT(IN    ) :: dt,gmt
   REAL(RKIND), INTENT(IN    ) :: anthro_pt_emis_scale_factor

   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: xlat,xlong,wind10m,pblh,xland,area
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: dz8w,rho_phy,t_phy,zmid,z_at_w
   REAL(RKIND),DIMENSION(25,1:num_anthro_pt,1:num_e_ant_pt_in),INTENT(IN)        :: e_ant_pt_in
   REAL(RKIND),DIMENSION(1:num_anthro_pt,1:num_e_ant_stack_groups_in),INTENT(IN)        :: e_ant_stack_groups_in
   INTEGER,DIMENSION(1:num_anthro_pt),INTENT(IN) :: ant_pt_local_cell_idx,ant_pt_rank
   INTEGER,INTENT(IN) :: myrank

   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT)     :: chem
                                                                               
  ! local
   INTEGER :: i,j,k,n,ii,kemit,kk,h,isearch,ilocal
   REAL(RKIND) :: conv_aer, conv_gas, emis
   REAL(RKIND) :: STACK_HT,STACK_DIA,STACK_VEL,STACK_TEMP,STACK_LAT,STACK_LON
   REAL(RKIND) :: T_1,T_2,wind_10m,EFF_H, PBL_H
   REAL(RKIND), DIMENSION(kts:kte) :: z_mid
 
! J-index is always 1
   j = 1
! Set the hour index to use
   h = FLOOR(gmt) + 1
! Loop over the EGUs
   do ii = 1,num_anthro_pt
! Skip if we are on the wrong rank for this EGU
      if (myrank /= ant_pt_rank(ii)) cycle
! Otherwise the index is just the local index
      i = ant_pt_local_cell_idx(ii)
! Final check in case the cell wasn't actually on any rank
      if ( i .le. 0 .or. i .gt. ite ) cycle
! Get the locations of this stack
      STACK_LAT  = e_ant_stack_groups_in(ii,index_STKLT)
      STACK_LON  = e_ant_stack_groups_in(ii,index_STKLG)
! Set the stack parameters
      STACK_HT   = e_ant_stack_groups_in(ii,index_STKHT)
      STACK_DIA  = e_ant_stack_groups_in(ii,index_STKDM)
      STACK_VEL  = e_ant_stack_groups_in(ii,index_STKVE)
      STACK_TEMP = e_ant_stack_groups_in(ii,index_STKTK)
! Set the Brigg's parameters
      T_1      = t_phy(i,kts,j)
      T_2      = t_phy(i,kts+1,j)
      wind_10m = wind10m(i,j)
      PBL_H    = pblh(i,j)
! Get the Vertical layers for this grid - TODO Move to wrapper
      do kk = kts,kte
         z_mid(kk)     = zmid(i,kk,j) - z_at_w(i,kts,j)
      enddo
! Call Briggs for this EGU and timestep      
      call plume_rise_briggs_rwc( wind_10m, T_1, T_2, PBL_H, z_mid, EFF_H, kemit, kts, kte, &
                                        STACK_HT,STACK_DIA,STACK_VEL,STACK_TEMP )
! Set the calculated emission layer
      k = kemit
!     Conversion factor for aerosol emissions (g/s) --> ug/kg --> Need to divide by area
      conv_aer = anthro_pt_emis_scale_factor * 1.e6_RKIND * dt / (rho_phy(i,k,j) *  dz8w(i,k,j) * area(i,j))
!     Conversion factor for gas phase emissions (moles/s) --> ppm/ppm --> Need to divide by area
      conv_gas = anthro_pt_emis_scale_factor * 60._RKIND * 1.E6_RKIND * 4.828E-4_RKIND * dt / ( rho_phy(i,k,j) * dz8w(i,k,j) * area(i,j) )
!
! Start applyting the emissions, selecting the correct hour, h
      if (p_unspc_fine .gt. 0 .and. index_e_ant_pt_in_unspc_fine .gt. 0 ) then
         emis = conv_aer*e_ant_pt_in(h,ii,index_e_ant_pt_in_unspc_fine)
         !print*,'pt,kemit,',kemit,' emis = ',emis
         chem(i,k,j,p_unspc_fine)   = chem(i,k,j,p_unspc_fine) + emis
         !e_ant_out(i,k,j,index_e_ant_out_unspc_fine) = e_ant_out(i,k,j,index_e_ant_out_unspc_fine) + emis
      endif
!     
   enddo ! ii

  end subroutine mpas_smoke_anthro_pt_emis_driver

end module module_anthro_pt_emissions
