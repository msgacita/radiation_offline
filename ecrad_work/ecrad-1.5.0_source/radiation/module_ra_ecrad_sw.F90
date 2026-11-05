 
 !------------------------------------------------------------------
MODULE module_ra_ecrad_sw

!PK #if defined(mpas)
!PK !MPAS specific (Laura D. Fowler):
!PK use mpas_atmphys_constants,only : cp,g=>gravity
!PK use module_ra_rrtmg_vinterp,only: vinterp_ozn

!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * updated the sourcecode to WRF revision 3.5, except for the implementation
!>      of time-varying trace gases: added multiple layers above the model-top
!>      following Cavalo et al. (2010). added option to use the ozone climatology
!>      from the CAM radiation codes.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-07-17.
!>    * in subroutine taugb7, corrected line number 6145 "if (specparm .lt. 0.125_rb) then"
!>      with if (specparm1 .lt. 0.125_rb) then, in accordance with the bug fix made in WRF
!>      revision 3.7.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2015-05-04.
CONTAINS


 ! ==================================================================
 ! DEPRECATED / NOT USED in the MONAN-MPAS coupling path.
 !
 ! This subroutine is a leftover WRF-physics-style wrapper (single
 ! precision "real" WRF tile arrays with ids/ide/ims/ime/jds/jde/
 ! kds/kde indices) that is no longer called from
 ! radiation/radiation_interface.F90:driver_radiation_lw_sw. That
 ! routine now calls ECRAD_Computing_shortwave_fluxes() below
 ! directly, using the double-precision (jprb) column arrays it
 ! already computed via gas_optics_*/cloud_optics/add_aerosol_optics,
 ! which avoids:
 !   (a) REAL(8)->REAL(4) argument type mismatches that used to
 !       prevent this file from compiling with the rest of ecRad,
 !   (b) the fake single-column "j=1" WRF tile loop, which does not
 !       correspond to the MPAS/MONAN unstructured mesh, and
 !   (c) several placeholder/incomplete assignments left in this
 !       routine (e.g. gas%mixing_ratio/aerosol%mixing_ratio filling
 !       loops with no right-hand side, xlat/xlon hard-coded to 0,
 !       dz8w hard-coded to -9999, and tauaer3d/ssaaer3d/asyaer3d
 !       being overwritten with od_sw/ssa_sw/g_sw -- i.e. radiation
 !       *output* -- instead of real aerosol optical input).
 !
 ! Kept here only for reference/back-compatibility with any WRF-style
 ! caller. Do not call this from new MONAN/MPAS code; use
 ! radiation_driver_monan_mpas.F90 instead, which fills the native
 ! ecRad column types (single_level_type, thermodynamics_type,
 ! gas_type, cloud_type, aerosol_type) directly from MPAS pools.
 ! ==================================================================
 subroutine ecrad_swrad(config,cloud,flux,single_level,thermodynamics ,aerosol,gas,&
                       p3d,p8w,pi3d,&
                       t3d,t8w,dz8w,&
                       qv3d,qc3d,qr3d,       &
                       qi3d,qs3d,qg3d,       &
                       cldfra3d,o33d,tsk,albedo,        &
                       xland,xice,snow,coszr,xtime,gmt,julday,radt,    &
                       degrad,declin,solcon,xlat,xlong,icloud,o3input, &
                       noznlevels,pin,o3clim,gsw,swcf,rthratensw,      &
                       has_reqc,has_reqi,has_reqs,re_cloud,            &
                       re_ice,re_snow,                                 &
                       aer_opt,tauaer3d,ssaaer3d,asyaer3d,             &
                       swupt,swuptc,swdnt,swdntc,                      &
                       swupb,swupbc,swdnb,swdnbc,                      &
                       swupflx, swupflxc, swdnflx, swdnflxc,           &
                       swddir,swddni,swddif,                           &
                       od_sw,ssa_sw,g_sw,                              &
                       od_sw_cloud, ssa_sw_cloud, g_sw_cloud,          &
                       incoming_sw,sw_albedo_direct,sw_albedo_diffuse, &
                       nbndsw,istartcol, iendcol,                      &
                       nlev,                                           &
                       ids,ide, jds,jde, kds,kde,                      & 
                       ims,ime, jms,jme, kms,kme,                      &
                       its,ite, jts,jte, kts,kte                  &
                       )


!=================================================================================================================
    use parkind1,                 only : jprb
    use radiation_io,             only : nulout
    use radiation_config,         only : config_type
    use radiation_cloud,          only : cloud_type
    use radiation_flux,           only : flux_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_gas,            only : gas_type,IO3
    use radiation_gas_constants,  only : NMaxGases,GasName
    use radiation_aerosol_optics_data, only : aerosol_optics_type
    use radiation_config,         only : config_type, &
         &   IGasModelMonochromatic, IGasModelIFSRRTMG, &
         &   ISolverMcICA, ISolverSpartacus, ISolverHomogeneous, &
         &   ISolverTripleclouds
 implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(single_level_type),  intent(in)      :: single_level
    type(thermodynamics_type),intent(inout)   :: thermodynamics

    type(config_type),        intent(in)   :: config
    type(cloud_type),         intent(inout):: cloud
    type(aerosol_type),       intent(in)   :: aerosol
    type(gas_type),           intent(in)   :: gas

    ! Output
    type(flux_type),          intent(inout):: flux

    ! Local variables 
    !                 
    type(aerosol_optics_type)    :: aerosol_optics
    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb),intent(inout), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw
    ! Layer in-cloud optical depth, single scattering albedo and
    ! asymmetry factor of hydrometeors in each shortwave band
    real(jprb),intent(inout), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud
    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb),intent(inout), dimension(config%n_g_sw,istartcol:iendcol) :: incoming_sw
    ! Direct and diffuse shortwave surface albedo in each shortwave
    ! g-point; note that these are weighted averages of the values
    ! from individual tiles
    real(jprb),intent(inout), dimension(config%n_g_sw, istartcol:iendcol) :: sw_albedo_direct
    real(jprb),intent(inout), dimension(config%n_g_sw, istartcol:iendcol) :: sw_albedo_diffuse

!--- input arguments:
 integer,intent(in):: ids,ide, jds,jde, kds,kde, &
                      ims,ime, jms,jme, kms,kme, &
                      its,ite, jts,jte, kts,kte,nbndsw

 integer,intent(in):: icloud,has_reqc,has_reqi,has_reqs
 integer,intent(in):: julday
 integer,intent(in),optional:: o3input

 real,intent(in):: radt,degrad,xtime,declin,solcon,gmt
 real,intent(in),dimension(ims:ime,jms:jme):: xlat,xlong
 real,intent(in),dimension(ims:ime,jms:jme):: albedo,tsk,snow,xice,xland
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: t3d,p3d,pi3d
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: dz8w
 real,intent(in),dimension(ims:ime,kms:kme+1,jms:jme):: p8w,t8w

 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: re_cloud,re_ice,re_snow
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    cldfra3d,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,o33d

!--- additional input arguments to use the CAM ozone climatology:
 integer,intent(in):: noznlevels
 real,intent(in),dimension(1:noznlevels),optional:: pin
 real,intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: o3clim

!--- additional input arguments of the aerosol optical depth, single scattering albedo, and asymmetry factor. to
!    date, the only kind of aerosols included in MPAS are the "water-friendly" and "ice-friendly" aerosols used
!    in the Thompson cloud microphysics scheme:
 integer,intent(in),optional:: aer_opt
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme,1:nbndsw),optional:: tauaer3d,ssaaer3d,asyaer3d

!--- inout arguments:
 real,intent(inout),dimension(ims:ime,jms:jme):: coszr,gsw,swcf
 real,intent(inout),dimension(ims:ime,jms:jme),optional:: &
    swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc

 real,intent(inout),dimension(ims:ime,kms:kme,jms:jme):: rthratensw

!--- output arguments:
  real,intent(out),dimension(ims:ime,jms:jme),optional:: &
     swddir,swddni,swddif
  real,intent(out),dimension(ims:ime,kms:kme+2,jms:jme ),optional:: &
     swupflx,swupflxc,swdnflx,swdnflxc
  integer :: i,k,j,n
   ! Fill reversed thermodynamic arrays
   DO k=1,nlev+1 
      do j = jts,jte
         DO i=istartcol,iendcol
            thermodynamics%pressure_hl   (i, k) = p8w(i,k,j)
            thermodynamics%temperature_hl(i, k) = t8w(i,k,j)
         END DO
      END DO
   END DO   
   if (config%do_clouds) then
      DO k=1,nlev 
         do j = jts,jte
            DO i=istartcol,iendcol
            ! Fill reversed cloud arrays
               cloud%q_liq(i,k) = qc3d(i,k,j)
               cloud%re_liq(i,k)= re_cloud(i,k,j)
               cloud%q_ice(i,k) = qi3d(i,k,j)
               cloud%re_ice(i,k)= re_ice(i,k,j)
               cloud%fraction(i,k)= cldfra3d(i,k,j)
               if (allocated(cloud%fractional_std)) then
               !    cloud%fractional_std(i,k)=
               else
                  cloud%fractional_std(i,k) = 0.0_jprb
               end if
               !print*,'allocated(cloud%inv_cloud_effective_size)',&
               !       allocated(cloud%inv_cloud_effective_size)
               if (allocated(cloud%inv_cloud_effective_size)) then
               !    cloud%inv_cloud_effective_size(i,k)=
               else
                !cloud%inv_cloud_effective_size(i,k) = 0.0_jprb
               end if
            END DO 
         END DO
      END DO 
      DO k=1,nlev-1
         DO i=istartcol,iendcol
            ! Fill reversed cloud arrays
            ! cloud%overlap_param(i,k) = 
         END DO
      END DO 
   end if
    ! Fill reversed gas arrays
    if (allocated(gas%mixing_ratio)) then
      DO n=1,NMaxGases ! Number of aerosol types
         DO k=1,nlev
            DO i=istartcol,iendcol
               !gas%mixing_ratio(i,k,n) =
               !gas%mixing_ratio(i,k,IO3)
               !print*,i,k,n,GasName(n), gas%mixing_ratio(i,k,n)
            END DO
         END DO
      END DO 
    end if
!    call gas%reverse(istartcol, iendcol, gas_rev)

    if (allocated(aerosol%mixing_ratio)) then
      DO n=1,aerosol_optics%ntype ! Number of aerosol types
         DO k=aerosol%iendlev,aerosol%istartlev
            DO i=istartcol,iendcol
               !aerosol%mixing_ratio(i,k,n) =
            END DO
         END DO
      END DO 
    end if

   CALL ECRAD_Computing_shortwave_fluxes(nlev, istartcol, iendcol, &
                                       config,cloud,flux,single_level,&
                                       thermodynamics,&
                                       sw_albedo_direct,sw_albedo_diffuse,incoming_sw,&
                                       od_sw, ssa_sw, g_sw,od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

end subroutine ecrad_swrad

!
!=================================================================================================================
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
!!                    'Computing shortwave fluxes'                             !!
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ECRAD_Computing_shortwave_fluxes(nlev, istartcol, iendcol,&
                                              config,cloud,flux,single_level,&
                                              thermodynamics,&
                                              sw_albedo_direct,sw_albedo_diffuse,incoming_sw,&
                                              od_sw, ssa_sw, g_sw,od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
    use parkind1,                 only : jprb
    use radiation_io,             only : nulout
    use radiation_config,         only : config_type
    use radiation_cloud,          only : cloud_type
    use radiation_flux,           only : flux_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_spartacus_sw,   only : solver_spartacus_sw
    use radiation_tripleclouds_sw,only : solver_tripleclouds_sw
    use radiation_mcica_sw,       only : solver_mcica_sw
    use radiation_cloudless_sw,   only : solver_cloudless_sw
    use radiation_homogeneous_sw, only : solver_homogeneous_sw

    use radiation_config,         only : config_type, &
         &   IGasModelMonochromatic, IGasModelIFSRRTMG, &
         &   ISolverMcICA, ISolverSpartacus, ISolverHomogeneous, &
         &   ISolverTripleclouds

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(single_level_type),  intent(in)   :: single_level
    type(thermodynamics_type),intent(in)   :: thermodynamics

    type(config_type),        intent(in)   :: config
    type(cloud_type),         intent(inout):: cloud
    ! Output
    type(flux_type),          intent(inout):: flux

    ! Local variables

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb),intent(inout), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw
    ! Layer in-cloud optical depth, single scattering albedo and
    ! asymmetry factor of hydrometeors in each shortwave band
    real(jprb),intent(inout), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud
    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb),intent(inout), dimension(config%n_g_sw,istartcol:iendcol) :: incoming_sw
    ! Direct and diffuse shortwave surface albedo in each shortwave
    ! g-point; note that these are weighted averages of the values
    ! from individual tiles
    real(jprb),intent(inout), dimension(config%n_g_sw, istartcol:iendcol) :: sw_albedo_direct
    real(jprb),intent(inout), dimension(config%n_g_sw, istartcol:iendcol) :: sw_albedo_diffuse

          write(*,'(a)') 'Computing shortwave fluxes'


      if (config%do_sw) then
        if (config%iverbose >= 2) then
          write(nulout,'(a)') 'Computing shortwave fluxes'
        end if

        if (config%i_solver_sw == ISolverMcICA) then
          ! Compute fluxes using the McICA shortwave solver
          call solver_mcica_sw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        else if (config%i_solver_sw == ISolverSPARTACUS) then
          ! Compute fluxes using the SPARTACUS shortwave solver
          call solver_spartacus_sw(nlev,istartcol,iendcol, &
               &  config, single_level, thermodynamics, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        else if (config%i_solver_sw == ISolverTripleclouds) then
          ! Compute fluxes using the Tripleclouds shortwave solver
          call solver_tripleclouds_sw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        elseif (config%i_solver_sw == ISolverHomogeneous) then
          ! Compute fluxes using the homogeneous solver
          call solver_homogeneous_sw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
               &  g_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        else
          ! Compute fluxes using the cloudless solver
          call solver_cloudless_sw(nlev,istartcol,iendcol, &
               &  config, single_level, od_sw, ssa_sw, g_sw, &
               &  sw_albedo_direct, sw_albedo_diffuse, &
               &  incoming_sw, flux)
        end if
      end if

  end subroutine ECRAD_Computing_shortwave_fluxes


!------------------------------------------------------------------

END MODULE module_ra_ecrad_sw
