 !------------------------------------------------------------------
MODULE module_ra_ecrad_lw

!PK #if defined(mpas)
!PK !MPAS specific (Laura D. Fowler):
!PK use mpas_atmphys_constants,only : cp,g=>gravity
!PK use module_ra_ecrad_vinterp,only: vinterp_ozn

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

!=================================================================================================================
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
!!                    'Computing longwave fluxes'                              !!
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ECRAD_Computing_longwave_fluxes(ncol, nlev, istartcol, iendcol, &
                                             config,cloud,flux,single_level,thermodynamics,&
                                             od_lw,ssa_lw, g_lw, lw_emission,lw_albedo,planck_hl,&
                                             od_lw_cloud,ssa_lw_cloud, g_lw_cloud)
    use parkind1,                 only : jprb
    use radiation_io,             only : nulout
    use radiation_config,         only : config_type
    use radiation_cloud,          only : cloud_type
    use radiation_flux,           only : flux_type
    use radiation_config,         only : &
         &   IGasModelMonochromatic, IGasModelIFSRRTMG, &
         &   ISolverMcICA, ISolverSpartacus, ISolverHomogeneous, &
         &   ISolverTripleclouds
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloudless_lw,   only : solver_cloudless_lw
    use radiation_mcica_lw,       only : solver_mcica_lw
    use radiation_homogeneous_lw, only : solver_homogeneous_lw
    use radiation_spartacus_lw,   only : solver_spartacus_lw
    use radiation_tripleclouds_lw,only : solver_tripleclouds_lw
    ! Inputs
    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in)   :: config
    type(single_level_type),  intent(in)   :: single_level
    type(thermodynamics_type),intent(in)   :: thermodynamics
    type(cloud_type),         intent(inout):: cloud
    ! Output
    type(flux_type),          intent(inout):: flux
    ! Local variables

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each longwave g-point, where the latter
    ! two variables are only defined if aerosol longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb),intent(inout), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_lw
    real(jprb),intent(inout), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw, g_lw

    ! The longwave emission from and albedo of the surface in each
    ! longwave g-point; note that these are weighted averages of the
    ! values from individual tiles
    real(jprb),intent(inout), dimension(config%n_g_lw, istartcol:iendcol) :: lw_emission
    real(jprb),intent(inout), dimension(config%n_g_lw, istartcol:iendcol) :: lw_albedo
    ! The Planck function (emitted flux from a black body) at half
    ! levels
    real(jprb),intent(inout), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! Layer in-cloud optical depth, single scattering albedo and
    ! asymmetry factor of hydrometeors in each longwave band, where
    ! the latter two variables are only defined if hydrometeor
    ! longwave scattering is enabled (otherwise both are treated as
    ! zero).
    real(jprb),intent(inout), dimension(config%n_bands_lw,nlev,istartcol:iendcol) :: od_lw_cloud
    real(jprb),intent(inout), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw_cloud, g_lw_cloud



      if (config%do_lw) then
        if (config%iverbose >= 2) then
          write(nulout,'(a)') 'Computing longwave fluxes'
        end if

        if (config%i_solver_lw == ISolverMcICA) then
          ! Compute fluxes using the McICA longwave solver
          call solver_mcica_lw(nlev,istartcol,iendcol, &
               &  config, single_level, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, &
               &  g_lw_cloud, planck_hl, lw_emission, lw_albedo, flux)
        else if (config%i_solver_lw == ISolverSPARTACUS) then
          ! Compute fluxes using the SPARTACUS longwave solver
          call solver_spartacus_lw(nlev,istartcol,iendcol, &
               &  config, thermodynamics, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  planck_hl, lw_emission, lw_albedo, flux)
        else if (config%i_solver_lw == ISolverTripleclouds) then
          ! Compute fluxes using the Tripleclouds longwave solver
          call solver_tripleclouds_lw(nlev,istartcol,iendcol, &
               &  config, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
               &  planck_hl, lw_emission, lw_albedo, flux)
        elseif (config%i_solver_lw == ISolverHomogeneous) then
          ! Compute fluxes using the homogeneous solver
          call solver_homogeneous_lw(nlev,istartcol,iendcol, &
               &  config, cloud, & 
               &  od_lw, ssa_lw, g_lw, od_lw_cloud, ssa_lw_cloud, &
               &  g_lw_cloud, planck_hl, lw_emission, lw_albedo, flux)
        else
          ! Compute fluxes using the cloudless solver
          call solver_cloudless_lw(nlev,istartcol,iendcol, &
               &  config, od_lw, ssa_lw, g_lw, &
               &  planck_hl, lw_emission, lw_albedo, flux)
        end if
      end if


  end subroutine ECRAD_Computing_longwave_fluxes
!------------------------------------------------------------------

END MODULE module_ra_ecrad_lw
