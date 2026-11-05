! era5_ecrad_bridge.F90
!
! Fills ecRad's native column types (single_level_type,
! thermodynamics_type, gas_type, cloud_type, flux_type) directly from
! one era5_column_type (see era5_bridge.F90) and calls
! driver_radiation_lw_sw. This is the ecRad counterpart to
! rrtmg_offline_driver.F90's RRTMG-SW/RRTMG-LW calls, so that the same
! ERA5 input can drive either radiation scheme, selected at run time
! via the &offline_control radiation_scheme namelist entry (see
! namelist_offline.F90 and rad_offline_driver.F90).
!
! Level-order note (this is SIMPLER than the RRTMG path): ERA5
! pressure levels are stored top-to-bottom (index 1 = lowest pressure/
! highest altitude), which is exactly the "increasing pressure with
! index" convention ecRad's own examples use natively -- so, unlike
! rrtmg_offline_driver.F90 (which reverses ERA5's order to match
! MPAS/WRF's bottom-up convention) or radiation_driver_monan_mpas.F90
! (which fills bottom-up and relies on ecRad's automatic
! radiation_reverse), this bridge uses ERA5's OWN order directly, with
! no reversal anywhere. Only the "valid" (above-ground) levels,
! col%p(1:col%nlev_valid), are used.
!
! Ozone note (also simpler than the RRTMG path): ecRad's gas_type%put()
! takes a full-column mass mixing ratio array directly, with no
! restriction analogous to RRTMG's o3input flag -- so ERA5's real
! ozone (col%o3) is used as-is, with none of the "deliberate
! modification" workaround documented in module_ra_rrtmg_sw.F/_lw.F.
!
! What is NOT done here (documented gaps, consistent with the RRTMG
! path so the two schemes are compared on equal footing):
!   - aerosols (config%use_aerosols stays .false.)
!   - cloud effective radius (fixed 10/40 micron placeholders, ERA5
!     has no cloud microphysics diagnostics)
!   - cloud overlap decorrelation length (left at whatever
!     overlap_scheme_name/i_overlap_scheme says in the ecRad
!     config.nam; no per-column latitude/day-of-year adjustment the
!     way radiation_driver_monan_mpas.F90 does, since that formula was
!     itself taken from MONAN's own RRTMG coupling and there is no
!     obviously-more-correct ERA5 equivalent to reach for here)
!
module era5_ecrad_bridge

  use parkind1,        only : jprb, jpim
  use mpas_kind_types,  only : RKIND
  use era5_bridge,      only : era5_column_type

  use radiation_config,         only : config_type
  use radiation_interface,      only : setup_radiation, driver_radiation_lw_sw
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, IVolumeMixingRatio
  use radiation_gas_constants,  only : IH2O, ICO2, IO3, IN2O, ICH4, IO2
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_flux,           only : flux_type

  implicit none
  private

  public :: era5_ecrad_init
  public :: era5_ecrad_run_column

  logical, save :: is_initialized = .false.

contains

  !=================================================================
  ! era5_ecrad_init
  !
  ! Call once, before any era5_ecrad_run_column call. Reads the ecRad
  ! namelist (same config.nam format as ecRad's own standalone driver)
  ! and loads the gas/cloud/aerosol optics look-up tables it points
  ! to (config%directory_name, typically "data/" -- these NetCDF files
  ! are NOT bundled with the ecRad source and must be obtained
  ! separately; see DADOS_ECRAD.md).
  !=================================================================
  subroutine era5_ecrad_init(config, namelist_path)

    type(config_type), intent(inout) :: config
    character(len=*),   intent(in)   :: namelist_path

    call config%read(file_name=trim(namelist_path))
    call setup_radiation(config)
    is_initialized = .true.

  end subroutine era5_ecrad_init


  !=================================================================
  ! era5_ecrad_run_column
  !
  ! Runs ecRad for one ERA5 column/timestep and returns SW/LW
  ! radiative heating rate profiles (K/s, uncoupled -- same convention
  ! as rrtmg_offline_driver.F90's rthratensw/rthratenlw), on the
  ! col%nlev_valid above-ground ERA5 levels, in ERA5's OWN top-to-
  ! bottom order (index 1 = highest level) so the two schemes' CSV
  ! output lines up directly.
  !=================================================================
  subroutine era5_ecrad_run_column(config, col, declin, solcon, gmt, &
       &  degrad, rthratensw, rthratenlw, flux_sw_up, flux_sw_dn, &
       &  flux_lw_up, flux_lw_dn)

    type(config_type),       intent(in)  :: config
    type(era5_column_type),  intent(in)  :: col
    real(kind=RKIND), intent(in) :: declin  ! radians
    real(kind=RKIND), intent(in) :: solcon  ! W/m2
    real(kind=RKIND), intent(in) :: gmt     ! decimal hours UTC
    real(kind=RKIND), intent(in) :: degrad  ! pi/180
    real(kind=RKIND), intent(out) :: rthratensw(col%nlev_valid)
    real(kind=RKIND), intent(out) :: rthratenlw(col%nlev_valid)
    ! Half-level (interface) flux profiles, W/m2, size nlev_valid+1,
    ! in the SAME top-to-bottom ERA5 order as everything else in this
    ! bridge (index 1 = TOA-most interface, index nlev_valid+1 =
    ! surface) -- index 1 is therefore the TOA flux and index
    ! nlev_valid+1 the surface flux; every index in between is the
    ! flux at that layer boundary ("nas camadas").
    real(kind=RKIND), intent(out) :: flux_sw_up(col%nlev_valid+1)
    real(kind=RKIND), intent(out) :: flux_sw_dn(col%nlev_valid+1)
    real(kind=RKIND), intent(out) :: flux_lw_up(col%nlev_valid+1)
    real(kind=RKIND), intent(out) :: flux_lw_dn(col%nlev_valid+1)

    type(single_level_type)   :: single_level
    type(thermodynamics_type) :: thermodynamics
    type(gas_type)            :: gas
    type(cloud_type)          :: cloud
    type(aerosol_type)        :: aerosol
    type(flux_type)           :: flux

    integer :: nlev, k
    real(jprb) :: tloctm, hrang, xxlat, cos_sza
    real(jprb) :: p_mid(col%nlev_valid), t_mid(col%nlev_valid)
    real(jprb) :: p_hl(col%nlev_valid+1), t_hl(col%nlev_valid+1)
    real(jprb) :: fzm, fzp
    real(jprb) :: scaling
    real(kind=RKIND), parameter :: gravity = 9.80616_RKIND, cp = 1004.5_RKIND

    if (.not. is_initialized) stop 'era5_ecrad_bridge: call era5_ecrad_init first'

    nlev = col%nlev_valid

    ! ------------------------------------------------------------
    ! Cosine of solar zenith angle -- ecRad, unlike rrtmg_swrad, does
    ! NOT compute this internally; the host must supply it. Same
    ! formula rrtmg_swrad uses (module_ra_rrtmg_sw.F), but with no
    ! xt24/radt forward-averaging term: this is an exact instantaneous
    ! angle for the ERA5 valid time (see the matching fix applied to
    ! rrtmg_offline_driver.F90's radt).
    ! ------------------------------------------------------------
    tloctm  = real(gmt,jprb) + real(col%lon,jprb)/15.0_jprb
    hrang   = 15.0_jprb*(tloctm-12.0_jprb)*real(degrad,jprb)
    xxlat   = real(col%lat,jprb)*real(degrad,jprb)
    cos_sza = sin(xxlat)*sin(real(declin,jprb)) + cos(xxlat)*cos(real(declin,jprb))*cos(hrang)

    ! ------------------------------------------------------------
    ! Allocate ecRad column types (ncol=1: this bridge is a
    ! column-at-a-time tool, matching rrtmg_offline_driver.F90)
    ! ------------------------------------------------------------
    call single_level%allocate(1, 1, 1, .false.)
    single_level%is_simple_surface = .true.
    single_level%cos_sza(1)          = max(cos_sza, 0.0_jprb)
    single_level%skin_temperature(1) = real(col%skt, jprb)
    single_level%sw_albedo(1,1)      = real(col%albedo, jprb)
    single_level%lw_emissivity(1,1)  = 0.98_jprb   ! ERA5 has no broadband
                                                     ! emissivity field,
                                                     ! same fallback as
                                                     ! rrtmg_offline_driver.F90

    ! DELIBERATE FIX (bug found empirically, see conversation log /
    ! project notes: only the ecRad path -- never RRTMG -- produced
    ! absurd values, concentrated near the surface, on real ERA5
    ! data): single_level%allocate() (radiation_single_level.F90)
    ! allocates single_level%iseed but never assigns it -- it is a
    ! genuinely uninitialized integer array until this point. It is
    ! read by the McICA cloud-generator solvers (solver_mcica_sw/lw)
    ! as the per-column random seed; ecRad's own source documents
    ! "it is expected that the host model generates this from the
    ! model time, longitude and latitude, in order that the model is
    ! deterministic" (radiation_single_level.F90). This bridge calls
    ! era5_ecrad_run_column once per column with ncol=1, so ecRad's
    ! own init_seed_simple() fallback (iseed(jcol)=jcol) would give
    ! EVERY column of the whole ERA5 grid/time range the same seed
    ! (1) -- defeating its purpose. Build the seed from lon/lat/time
    ! instead, matching ecRad's own recommendation.
    single_level%iseed(1) = int(abs( mod( &
         &   nint(col%lon*1000.0_jprb,jpim) &
         & + nint(col%lat*1000.0_jprb,jpim)*7919_jpim &
         & + nint(gmt*3600.0_jprb,jpim)*104729_jpim, &
         &   huge(1_jpim)/2_jpim ) ), kind(single_level%iseed))

    call thermodynamics%allocate(1, nlev, use_h2o_sat=.true.)

    ! ------------------------------------------------------------
    ! Half-level pressure/temperature. ERA5's own order (index 1 =
    ! TOA) is used directly. Top interface: log-pressure
    ! extrapolation beyond level 1. Interior interfaces:
    ! height-weighted (fzm/fzp) between adjacent full levels, same
    ! method as radiation_driver_monan_mpas.F90's
    ! fill_thermodynamics_from_mpas and rrtmg_offline_driver.F90.
    ! Surface interface (index nlev+1): the ACTUAL ERA5 surface
    ! pressure/skin temperature -- no extrapolation needed here,
    ! unlike the RRTMG path, since ERA5's own surface fields ARE the
    ! true lower boundary (RRTMG's kts=1 boundary had to be
    ! extrapolated from the lowest atmospheric level instead).
    ! ------------------------------------------------------------
    do k = 1, nlev
      p_mid(k) = real(col%p(k), jprb)
      t_mid(k) = real(col%t(k), jprb)
    end do

    ! top interface: extrapolate assuming it is as far above level 1
    ! (in log-pressure) as level 1 is above level 2, and linearly in
    ! temperature over that same step
    p_hl(1) = exp(2.0_jprb*log(p_mid(1)) - log(p_mid(2)))
    t_hl(1) = 2.0_jprb*t_mid(1) - t_mid(2)

    do k = 2, nlev
      fzm = real((col%z(k)-col%z(k-1)) / (col%z(min(k+1,nlev))-col%z(k-1)), jprb)
      if (k < nlev) then
        fzm = real((col%z(k)-col%z(k-1)) / (col%z(k+1)-col%z(k-1)), jprb)
        fzp = real((col%z(k+1)-col%z(k)) / (col%z(k+1)-col%z(k-1)), jprb)
        p_hl(k) = fzp*p_mid(k-1) + fzm*p_mid(k)
        t_hl(k) = fzp*t_mid(k-1) + fzm*t_mid(k)
      else
        ! second-to-last interface (k=nlev): interpolate between the
        ! two lowest full levels as usual
        p_hl(k) = 0.5_jprb*(p_mid(k-1)+p_mid(k))
        t_hl(k) = 0.5_jprb*(t_mid(k-1)+t_mid(k))
      end if
    end do

    p_hl(nlev+1) = real(col%sp, jprb)
    t_hl(nlev+1) = real(col%skt, jprb)

    thermodynamics%pressure_hl(1,:)    = p_hl
    thermodynamics%temperature_hl(1,:) = t_hl
    call thermodynamics%calc_saturation_wrt_liquid(1, 1)

    ! ------------------------------------------------------------
    ! Gases. Real ERA5 ozone used directly -- no CAM-climatology
    ! workaround needed here (see file header).
    ! ------------------------------------------------------------
    call gas%allocate(1, nlev)
    call gas%put(IH2O, IMassMixingRatio, reshape(real(col%qv(1:nlev),jprb), (/1,nlev/)))
    call gas%put(IO3,  IMassMixingRatio, reshape(real(col%o3(1:nlev),jprb), (/1,nlev/)))
    call gas%put_well_mixed(ICO2, IVolumeMixingRatio, 415.0e-6_jprb)
    call gas%put_well_mixed(IN2O, IVolumeMixingRatio, 332.0e-9_jprb)
    call gas%put_well_mixed(ICH4, IVolumeMixingRatio, 1921.0e-9_jprb)
    call gas%put_well_mixed(IO2,  IVolumeMixingRatio, 0.20944_jprb)

    ! ------------------------------------------------------------
    ! Clouds. ecRad's cloud_type has just two hydrometeor categories
    ! (liquid, ice); ERA5's four (clwc,ciwc,crwc,cswc) are combined
    ! liquid=clwc+crwc, ice=ciwc+cswc -- an approximation (rain/snow
    ! have different optical properties than cloud water/ice) but a
    ! reasonable one for a first offline comparison, and documented
    ! here rather than silently dropping rain/snow water content.
    ! ------------------------------------------------------------
    if (config%do_clouds) then
      call cloud%allocate(1, nlev, use_inhom_effective_size=.false.)
      do k = 1, nlev
        cloud%q_liq(1,k)    = real(col%qc(k)+col%qr(k), jprb)
        cloud%q_ice(1,k)    = real(col%qi(k)+col%qs(k), jprb)
        cloud%fraction(1,k) = real(col%cldfra(k), jprb)
        cloud%re_liq(1,k)   = 10.0e-6_jprb   ! ERA5 has no effective-radius
        cloud%re_ice(1,k)   = 40.0e-6_jprb   ! diagnostic; same fallback as
                                              ! rrtmg_offline_driver.F90
      end do
      call cloud%create_fractional_std(1, nlev, 0.0_jprb)
      call cloud%crop_cloud_fraction(1, 1, config%cloud_fraction_threshold, &
           &  config%cloud_mixing_ratio_threshold)

      ! DELIBERATE FIX (bug found empirically -- only the ecRad path,
      ! never RRTMG, produced non-physical values, concentrated near
      ! the surface, on real ERA5 data): cloud%allocate() eagerly
      ! allocates cloud%overlap_param (radiation_cloud.F90,
      ! allocate_cloud_arrays) but never assigns it -- it is genuinely
      ! uninitialized memory until this point. config.nam.sample's
      ! solvers (sw/lw_solver_name = "Tripleclouds") read
      ! cloud%overlap_param directly, with no internal fallback; it
      ! must be computed by the host from a decorrelation length via
      ! cloud%set_overlap_param(), which this bridge never called.
      ! 2000 m is a common default decorrelation length (same order of
      ! magnitude ecRad's own examples use); radiation_driver_monan_-
      ! mpas.F90 instead derives it per-column from latitude/day-of-
      ! year (NASA GMAO formula, report Sec. 4.2.2) -- a reasonable
      ! future refinement here, not required for this fix to be
      ! correct. thermodynamics%pressure_hl/temperature_hl (built
      ! above) already carry ERA5's own top-to-bottom order, which
      ! set_overlap_param detects and handles automatically.
      call cloud%set_overlap_param(thermodynamics, 2000.0_jprb)
    end if

    call flux%allocate(config, 1, 1, nlev)

    ! ------------------------------------------------------------
    ! Run ecRad and convert to K/s heating rates (same formula as
    ! radiation_driver_monan_mpas.F90 step 8, minus the exner
    ! division -- ERA5/this bridge works with actual temperature
    ! throughout, not potential temperature, so no theta conversion
    ! is needed here).
    ! ------------------------------------------------------------
    call driver_radiation_lw_sw(1, nlev, 1, 1, config, &
         &  single_level, thermodynamics, gas, cloud, aerosol, flux)

    do k = 1, nlev
      scaling = -real(gravity,jprb)/real(cp,jprb) &
           &  / (thermodynamics%pressure_hl(1,k+1) - thermodynamics%pressure_hl(1,k))
      rthratensw(k) = real( scaling * ( &
           &    (flux%sw_dn(1,k+1)-flux%sw_up(1,k+1)) &
           &  - (flux%sw_dn(1,k)  -flux%sw_up(1,k))   ), RKIND)
      rthratenlw(k) = real( scaling * ( &
           &    (flux%lw_dn(1,k+1)-flux%lw_up(1,k+1)) &
           &  - (flux%lw_dn(1,k)  -flux%lw_up(1,k))   ), RKIND)
    end do

    ! Interface flux profiles, same top-to-bottom order as everything
    ! else here -- index 1 (TOA) and index nlev+1 (surface) give the
    ! boundary ("topo"/"superfície") values directly; every other
    ! index is the flux at that internal layer boundary ("camadas").
    do k = 1, nlev+1
      flux_sw_up(k) = real(flux%sw_up(1,k), RKIND)
      flux_sw_dn(k) = real(flux%sw_dn(1,k), RKIND)
      flux_lw_up(k) = real(flux%lw_up(1,k), RKIND)
      flux_lw_dn(k) = real(flux%lw_dn(1,k), RKIND)
    end do

    call single_level%deallocate()
    call thermodynamics%deallocate()
    call gas%deallocate()
    if (config%do_clouds) call cloud%deallocate()
    call flux%deallocate()

  end subroutine era5_ecrad_run_column

end module era5_ecrad_bridge
