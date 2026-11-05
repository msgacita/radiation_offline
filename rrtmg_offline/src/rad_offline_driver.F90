! rad_offline_driver.F90
!
! Unified offline driver: reads &offline_control (namelist_offline.F90)
! to pick which radiation scheme to run -- 'rrtmg' (MPAS's own
! RRTMG-SW/RRTMG-LW, module_ra_rrtmg_sw.F/_lw.F, extracted verbatim
! from MPAS-Dev/MPAS-Model) or 'ecrad' (ECMWF's ecRad, via
! era5_ecrad_bridge.F90) -- and runs it column-by-column over the
! requested ERA5 time range, writing both schemes to the SAME CSV
! format so their output can be compared directly.
!
! This supersedes rrtmg_offline_driver.F90's command-line interface
! (kept in this project for anyone who only ever wants RRTMG and
! prefers not to write a namelist file) with a single entry point that
! mirrors how MONAN itself selects a radiation scheme -- a namelist
! switch (config_radt_sw_scheme/config_radt_lw_scheme in MONAN;
! radiation_scheme here) rather than a different driver subroutine per
! scheme.
!
! REQUIRED EXTERNAL DATA (see DADOS_ERA5.md and DADOS_ECRAD.md):
!   - if radiation_scheme='rrtmg': RRTMG_SW_DATA.DBL, RRTMG_LW_DATA.DBL
!     in the working directory
!   - if radiation_scheme='ecrad': the ecRad gas/cloud/aerosol optics
!     NetCDF files under the directory config.nam's directory_name
!     option points to (default "data/")
!   - one ERA5 pressure-level file and one ERA5 single-level file
!
! USAGE:
!   rad_offline_driver <namelist_file>
!
program rad_offline_driver

  use mpas_kind_types,    only : RKIND
  use mpas_derived_types, only : dm_info
  use mpas_atmphys_constants, only : cp, gravity
  use radconst_offline,  only : radconst
  use era5_bridge

  use namelist_offline, only : offline_control_type, read_offline_namelist

  use mpas_atmphys_rrtmg_swinit, only : rrtmg_initsw_forMPAS
  use mpas_atmphys_rrtmg_lwinit, only : rrtmg_initlw_forMPAS
  use module_ra_rrtmg_sw,        only : rrtmg_swrad
  use module_ra_rrtmg_lw,        only : rrtmg_lwrad

  use radiation_config,     only : config_type
  use era5_ecrad_bridge,    only : era5_ecrad_init, era5_ecrad_run_column
  use output_netcdf,        only : output_file_type, output_open, &
       &                           output_write_timestep, output_close

  implicit none

  type(offline_control_type) :: ctl
  character(len=512) :: namelist_file

  type(era5_file_type)   :: era5f
  type(era5_column_type) :: col
  type(config_type)      :: ecrad_config
  type(dm_info)           :: dminfo

  integer :: itime, nlev, k
  integer :: nlon, nlat, ilon_idx, ilat_idx, ilon, ilat, it_index
  type(output_file_type) :: outf
  real(kind=RKIND) :: julian, gmt, declin, solcon, degrad_local
  integer :: julday
  real(kind=RKIND), allocatable :: rthratensw(:), rthratenlw(:)
  real(kind=RKIND), allocatable :: flux_sw_up(:), flux_sw_dn(:)
  real(kind=RKIND), allocatable :: flux_sw_up_clear(:), flux_sw_dn_clear(:)
  real(kind=RKIND), allocatable :: flux_lw_up(:), flux_lw_dn(:)
  real(kind=RKIND), allocatable :: flux_lw_up_clear(:), flux_lw_dn_clear(:)

  ! ---- WRF/MPAS-style "1x1 tile" arrays for the RRTMG branch (same
  !      as rrtmg_offline_driver.F90 -- see that file for a more
  !      heavily commented version of this block) ----
  integer, parameter :: ids=1,ide=1, jds=1,jde=1
  integer, parameter :: ims=1,ime=1, jms=1,jme=1
  integer, parameter :: its=1,ite=1, jts=1,jte=1
  integer :: kds,kde, kms,kme, kts,kte,ilat_probe,ilon_probe,itime_probe
  logical :: warmed_up
  real(kind=RKIND), allocatable :: p3d(:,:,:), p8w(:,:,:), pi3d(:,:,:)
  real(kind=RKIND), allocatable :: t3d(:,:,:), t8w(:,:,:), dz8w(:,:,:)
  real(kind=RKIND), allocatable :: qv3d(:,:,:), qc3d(:,:,:), qr3d(:,:,:)
  real(kind=RKIND), allocatable :: qi3d(:,:,:), qs3d(:,:,:), qg3d(:,:,:)
  real(kind=RKIND), allocatable :: cldfra3d(:,:,:), o33d(:,:,:)
  real(kind=RKIND), allocatable :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)
  real(kind=RKIND), allocatable :: rre_cloud3d(:,:,:), rre_ice3d(:,:,:), rre_snow3d(:,:,:)
  real(kind=RKIND), allocatable :: reff_cloud(:), reff_ice(:)
  real(kind=RKIND), allocatable :: rthratensw3(:,:,:), rthratenlw3(:,:,:)
  real(kind=RKIND), allocatable :: swupflx(:,:,:), swupflxc(:,:,:)
  real(kind=RKIND), allocatable :: swdnflx(:,:,:), swdnflxc(:,:,:)
  real(kind=RKIND), allocatable :: lwupflx(:,:,:), lwupflxc(:,:,:)
  real(kind=RKIND), allocatable :: lwdnflx(:,:,:), lwdnflxc(:,:,:)
  real(kind=RKIND) :: tsk(ims:ime,jms:jme), albedo(ims:ime,jms:jme), emiss(ims:ime,jms:jme)
  real(kind=RKIND) :: xland(ims:ime,jms:jme), xice(ims:ime,jms:jme), snow(ims:ime,jms:jme)
  real(kind=RKIND) :: coszr(ims:ime,jms:jme), gsw(ims:ime,jms:jme), swcf(ims:ime,jms:jme)
  real(kind=RKIND) :: glw(ims:ime,jms:jme), olr(ims:ime,jms:jme), lwcf(ims:ime,jms:jme)
  real(kind=RKIND) :: xlat(ims:ime,jms:jme), xlong(ims:ime,jms:jme)
  real(kind=RKIND), parameter :: p0 = 100000.0_RKIND
  real(kind=RKIND), parameter :: Rd_over_cp = 287.0_RKIND/1004.5_RKIND
  real(kind=RKIND) :: z0,z1,z2,w1,w2,fzm,fzp

  ! =================================================================
  ! 0) Read namelist
  ! =================================================================
  if (command_argument_count() /= 1) then
    write(0,'(A)') 'usage: rad_offline_driver <namelist_file>'
    stop 1
  end if
  call get_command_argument(1, namelist_file)
  call read_offline_namelist(trim(namelist_file), ctl)

  write(*,'(A)') 'rad_offline_driver: radiation_scheme = '//trim(ctl%radiation_scheme)

  ! =================================================================
  ! 1) Initialize the selected scheme, ONCE
  ! =================================================================
  select case (trim(ctl%radiation_scheme))
  case ('rrtmg')
    call rrtmg_initsw_forMPAS(dminfo)
    call rrtmg_initlw_forMPAS(dminfo)
  case ('ecrad')
    call era5_ecrad_init(ecrad_config, trim(ctl%ecrad_namelist_path))
  end select

  call era5_open(trim(ctl%path_pressure_level), trim(ctl%path_single_level), era5f)

  ! -----------------------------------------------------------------
  ! WORKAROUND (not a root-cause fix -- see DADOS_ECRAD.md): the very
  ! first call to driver_radiation_lw_sw within a running process was
  ! observed, on at least one system/compiler combination, to return
  ! wildly wrong (unphysical, e.g. >0.1 K/s) longwave heating rates
  ! near the surface, while every SUBSEQUENT call in the same process
  ! -- even with nearly identical input data -- returned physically
  ! reasonable values. This was NOT reproduced with this project's own
  ! code compiled under gfortran with -fcheck=all -finit-real=snan
  ! (i.e. it is very unlikely to be an uninitialized variable in
  ! era5_ecrad_bridge.F90 itself), which points to some uninitialized
  ! or lazily-set internal state inside ecRad's own precompiled
  ! library (libradiation.a) that only becomes consistent after one
  ! full driver_radiation_lw_sw call has executed.
  !
  ! Rather than block on fully root-causing this inside ecRad's own
  ! (large, shared, validated-elsewhere) solver code, this "warm-up"
  ! call absorbs whatever that first-call effect is: it runs ecRad
  ! once, on the first requested column, and DISCARDS the result --
  ! every result actually written to the CSV output comes from a
  ! second-or-later call, which is the regime shown to be correct.
  ! -----------------------------------------------------------------
  if (trim(ctl%radiation_scheme) == 'ecrad') then
    warmed_up = .false.
    outer: do itime_probe = ctl%itime_start, ctl%itime_end
      do ilat_probe = ctl%ilat_start, ctl%ilat_end
        do ilon_probe = ctl%ilon_start, ctl%ilon_end
          call era5_read_column(era5f, ilon_probe, ilat_probe, itime_probe, col)
          if (col%nlev_valid >= 4) then
            allocate(rthratensw(col%nlev_valid), rthratenlw(col%nlev_valid))
            allocate(flux_sw_up(col%nlev_valid+1), flux_sw_dn(col%nlev_valid+1))
            allocate(flux_sw_up_clear(col%nlev_valid+1), flux_sw_dn_clear(col%nlev_valid+1))
            allocate(flux_lw_up(col%nlev_valid+1), flux_lw_dn(col%nlev_valid+1))
            allocate(flux_lw_up_clear(col%nlev_valid+1), flux_lw_dn_clear(col%nlev_valid+1))
            call era5_valid_time_to_julian_gmt(era5f, itime_probe, julian, gmt, julday)
            call radconst(declin, solcon, julian)
            call era5_ecrad_run_column(ecrad_config, col, declin, solcon, gmt, &
                 &  3.1415926_RKIND/180.0_RKIND, rthratensw, rthratenlw, &
                 &  flux_sw_up, flux_sw_dn, flux_lw_up, flux_lw_dn, &
                 &  flux_sw_up_clear, flux_sw_dn_clear, flux_lw_up_clear, flux_lw_dn_clear)
            deallocate(rthratensw, rthratenlw, flux_sw_up, flux_sw_dn, flux_sw_up_clear, flux_sw_dn_clear, flux_lw_up, flux_lw_dn, flux_lw_up_clear, flux_lw_dn_clear)
            warmed_up = .true.
          end if
          deallocate(col%p,col%t,col%z,col%qv,col%qc,col%qi,col%qr,col%qs,col%cldfra,col%o3)
          if (warmed_up) exit outer
        end do
      end do
    end do outer
    if (.not. warmed_up) then
      write(0,'(A)') 'rad_offline_driver: AVISO -- nenhuma coluna com nlev_valid>=4 '// &
           &          'encontrada para o aquecimento do ecRad; primeira chamada real '// &
           &          'pode retornar valores espurios (ver Secao 5.5 do relatorio)'
    else
      write(*,'(A)') 'rad_offline_driver: ecRad warm-up call done (result discarded)'
    end if
  end if


  ! Output NetCDF file is opened ONCE, up front, sized to the FULL
  ! requested grid (nlon x nlat) and the FULL input file's pressure
  ! levels (era5f%nlev, not any one column's nlev_valid -- see
  ! output_netcdf.F90's module header for why). This works regardless
  ! of how nlev_valid varies from column to column or hour to hour,
  ! since each column/hour is written using its OWN nlev_valid, with
  ! the untouched remainder of the fixed vertical dimension left at
  ! NetCDF's default fill value -- unlike an earlier single-column
  ! revision of this driver, which opened the file lazily and had to
  ! skip an hour outright if its nlev_valid differed from the first.
  nlon = ctl%ilon_end - ctl%ilon_start + 1
  nlat = ctl%ilat_end - ctl%ilat_start + 1
  call output_open(trim(ctl%path_output), nlon, nlat, era5f%nlev, &
       &  era5f%lon(ctl%ilon_start:ctl%ilon_end), &
       &  era5f%lat(ctl%ilat_start:ctl%ilat_end), &
       &  era5f%plev_hpa, era5f%time_units, trim(ctl%radiation_scheme), outf)

  ! =================================================================
  ! 2) Time loop (outer), area loop (inner: latitude, then longitude)
  ! =================================================================
  write(*,'(A,I0,A,I0)') 'rad_offline_driver: itime_start=', ctl%itime_start, &
       &  ' itime_end=', ctl%itime_end
  write(*,'(A,I0,A,I0,A,I0,A,I0)') 'rad_offline_driver: ilon=', ctl%ilon_start, &
       &  '..', ctl%ilon_end, ' ilat=', ctl%ilat_start, '..', ctl%ilat_end

  do itime = ctl%itime_start, ctl%itime_end
    write(*,'(A,I0)') 'rad_offline_driver: processing itime=', itime
    it_index = itime - ctl%itime_start + 1

    do ilat_idx = 1, nlat
      ilat = ctl%ilat_start + ilat_idx - 1
      do ilon_idx = 1, nlon
        ilon = ctl%ilon_start + ilon_idx - 1

        call era5_read_column(era5f, ilon, ilat, itime, col)
        nlev = col%nlev_valid
        ! Minimum is 4, not 2: RRTMG's own McICA subcolumn generator
        ! (mcica_subcol_gen, module_ra_rrtmg_sw.F/_lw.F) hard-requires
        ! at least 4 above-ground layers to seed its random number
        ! generator ("KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM
        ! FOUR LAYERS"), and stops the whole program outright if fewer
        ! are given -- not a recoverable per-column error from RRTMG's
        ! side, so this driver must never call it with nlev<4. Applied
        ! to both schemes for a consistent, physically sensible floor
        ! (2-3 valid levels is a degenerate column either way, e.g. a
        ! bad/fill-valued surface pressure at that point, not real
        ! atmosphere worth computing radiation for).
        if (nlev < 4) then
          write(0,'(A,I0,A,I0,A,I0,A,I0,A,F0.1,A,F0.2,A,F0.2)') 'skipping ilon=', ilon, ' ilat=', ilat, &
               &  ' itime=', itime, ' (only ', nlev, ' above-ground ERA5 levels, need >=4;', &
               &  ' col%sp(Pa)=', col%sp, ' lon=', col%lon, ' lat=', col%lat
          deallocate(col%p,col%t,col%z,col%qv,col%qc,col%qi,col%qr,col%qs,col%cldfra,col%o3)
          cycle
        end if

        call era5_valid_time_to_julian_gmt(era5f, itime, julian, gmt, julday)
        call radconst(declin, solcon, julian)
        degrad_local = 3.1415926_RKIND/180.0_RKIND

        allocate(rthratensw(nlev)  , rthratenlw(nlev))
	rthratensw=0.0;rthratenlw=0.0
        allocate(flux_sw_up(nlev+1), flux_sw_dn(nlev+1))
	flux_sw_up=0.0;flux_sw_dn=0.0
        allocate(flux_sw_up_clear(nlev+1), flux_sw_dn_clear(nlev+1))
	flux_sw_up_clear=0.0;flux_sw_dn_clear=0.0
        allocate(flux_lw_up(nlev+1), flux_lw_dn(nlev+1))
	flux_lw_up=0.0;flux_lw_dn=0.0
        allocate(flux_lw_up_clear(nlev+1), flux_lw_dn_clear(nlev+1))
	flux_lw_up_clear=0.0;flux_lw_dn_clear=0.0

    select case (trim(ctl%radiation_scheme))
    ! ---------------------------------------------------------------
    case ('ecrad')
    ! ---------------------------------------------------------------
      call era5_ecrad_run_column(ecrad_config, col, declin, solcon, gmt, &
           &  degrad_local, rthratensw, rthratenlw, &
           &  flux_sw_up, flux_sw_dn, flux_lw_up, flux_lw_dn, flux_sw_up_clear, flux_sw_dn_clear, flux_lw_up_clear, flux_lw_dn_clear)

      do k = 1, nlev
        if(ABS(rthratensw(k)) > 50.0) print*,'sw ',k,rthratensw(k), &
         col%p(k),col%t(k),col%z(k),col%qv(k),col%qc(k),col%qi(k),col%qr(k),col%qs(k),col%cldfra(k),col%o3(k) 
        if(ABS(rthratensw(k)) > 50.0) print*,'lw ',k,rthratensw(k), &
         col%p(k),col%t(k),col%z(k),col%qv(k),col%qc(k),col%qi(k),col%qr(k),col%qs(k),col%cldfra(k),col%o3(k) 
      end do
      do k = 1, nlev+1
        if(ABS(flux_lw_up(k)) > 5000.0) print*,'lwup ',k,flux_lw_up(k), &
         col%p(k),col%t(k),col%z(k),col%qv(k),col%qc(k),col%qi(k),col%qr(k),col%qs(k),col%cldfra(k),col%o3(k) 
 
        if(ABS(flux_lw_dn(k)) > 5000.0) print*,'lwdn ',k,flux_lw_dn(k), &
         col%p(k),col%t(k),col%z(k),col%qv(k),col%qc(k),col%qi(k),col%qr(k),col%qs(k),col%cldfra(k),col%o3(k) 

      end do

    ! ---------------------------------------------------------------
    case ('rrtmg')
    ! ---------------------------------------------------------------
      ! kme must be big enough for the OPTIONAL spectral flux
      ! diagnostics (swupflx/swdnflx/lwupflx/lwdnflx, dimension
      ! kms:kme+2 in the interface), which RRTMG fills out to
      ! "nlayers+1" -- nlayers = kte + max(nint(p_top_hPa/4),1),
      ! i.e. RRTMG synthesizes extra layers in fixed ~4 hPa steps from
      ! the model top down to ~0 hPa purely to populate these arrays
      ! (see the note at the rrtmg_swrad/rrtmg_lwrad calls below, and
      ! DADOS_ECRAD.md/README.md for the segfault this caused before
      ! it was accounted for). p_top_pa is estimated here with the same
      ! log-pressure extrapolation used for p8w's own top interface,
      ! computed from the raw ERA5 column BEFORE any array exists.
      block
        real(kind=RKIND) :: p_top_pa, nlayers_extra_r
        integer :: nlayers_extra
        p_top_pa = exp(2.0_RKIND*log(col%p(1)) - log(col%p(2)))
        nlayers_extra_r = anint((p_top_pa/100.0_RKIND)/4.0_RKIND)
        nlayers_extra = max(int(nlayers_extra_r), 1)
        kds=1; kde=nlev; kms=1; kme=nlev+1+nlayers_extra+2; kts=1; kte=nlev
      end block

      allocate(p3d(ims:ime,kms:kme,jms:jme), pi3d(ims:ime,kms:kme,jms:jme))
      allocate(t3d(ims:ime,kms:kme,jms:jme))
      allocate(qv3d(ims:ime,kms:kme,jms:jme), qc3d(ims:ime,kms:kme,jms:jme))
      allocate(qr3d(ims:ime,kms:kme,jms:jme), qi3d(ims:ime,kms:kme,jms:jme))
      allocate(qs3d(ims:ime,kms:kme,jms:jme), qg3d(ims:ime,kms:kme,jms:jme))
      allocate(cldfra3d(ims:ime,kms:kme,jms:jme), o33d(ims:ime,kms:kme,jms:jme))
      allocate(re_cloud(ims:ime,kms:kme,jms:jme), re_ice(ims:ime,kms:kme,jms:jme))
      allocate(re_snow(ims:ime,kms:kme,jms:jme))
      allocate(rre_cloud3d(ims:ime,kms:kme,jms:jme), rre_ice3d(ims:ime,kms:kme,jms:jme))
      allocate(rre_snow3d(ims:ime,kms:kme,jms:jme))
      rre_cloud3d = 0.0_RKIND
      rre_ice3d = 0.0_RKIND
      rre_snow3d = 0.0_RKIND
      allocate(dz8w(ims:ime,kms:kme,jms:jme))
      allocate(p8w(ims:ime,kms:kme,jms:jme), t8w(ims:ime,kms:kme,jms:jme))
      allocate(rthratensw3(ims:ime,kms:kme,jms:jme), rthratenlw3(ims:ime,kms:kme,jms:jme))
      allocate(swupflx(ims:ime,kms:kme+2,jms:jme), swupflxc(ims:ime,kms:kme+2,jms:jme))
      allocate(swdnflx(ims:ime,kms:kme+2,jms:jme), swdnflxc(ims:ime,kms:kme+2,jms:jme))
      allocate(lwupflx(ims:ime,kms:kme+2,jms:jme), lwupflxc(ims:ime,kms:kme+2,jms:jme))
      allocate(lwdnflx(ims:ime,kms:kme+2,jms:jme), lwdnflxc(ims:ime,kms:kme+2,jms:jme))
      swupflx=0.0_RKIND; swupflxc=0.0_RKIND; swdnflx=0.0_RKIND; swdnflxc=0.0_RKIND
      lwupflx=0.0_RKIND; lwupflxc=0.0_RKIND; lwdnflx=0.0_RKIND; lwdnflxc=0.0_RKIND

      do k = kts, kte
        p3d(1,k,1)      = col%p(nlev-k+1)
        t3d(1,k,1)      = col%t(nlev-k+1)
        qv3d(1,k,1)     = col%qv(nlev-k+1)
        qc3d(1,k,1)     = col%qc(nlev-k+1)
        qr3d(1,k,1)     = col%qr(nlev-k+1)
        qi3d(1,k,1)     = col%qi(nlev-k+1)
        qs3d(1,k,1)     = col%qs(nlev-k+1)
        qg3d(1,k,1)     = 0.0_RKIND
        cldfra3d(1,k,1) = col%cldfra(nlev-k+1)
        o33d(1,k,1)     = col%o3(nlev-k+1)
        pi3d(1,k,1)     = (p3d(1,k,1)/p0)**Rd_over_cp
        re_cloud(1,k,1) = 0.0_RKIND
        re_ice(1,k,1)   = 0.0_RKIND
        re_snow(1,k,1)  = 0.0_RKIND
      end do

      block
        real(kind=RKIND) :: zfull(kts:kte)
        do k = kts, kte
          zfull(k) = col%z(nlev-k+1)
        end do

        p8w(1,kts,1) = col%sp
        t8w(1,kts,1) = t3d(1,kts,1) + (t3d(1,kts,1)-t3d(1,kts+1,1)) &
             &       * (col%zsfc-zfull(kts)) / (zfull(kts)-zfull(kts+1))

        do k = kts+1, kte
          if (k < kte) then
            fzm = (zfull(k)-zfull(k-1)) / (zfull(k+1)-zfull(k-1))
            fzp = (zfull(k+1)-zfull(k)) / (zfull(k+1)-zfull(k-1))
            p8w(1,k,1) = fzm*p3d(1,k,1) + fzp*p3d(1,k-1,1)
            t8w(1,k,1) = fzm*t3d(1,k,1) + fzp*t3d(1,k-1,1)
          else
            z0 = zfull(kte) + (zfull(kte)-zfull(kte-1))
            z1 = zfull(kte); z2 = zfull(kte-1)
            w1 = (z0-z2)/(z1-z2); w2 = 1.0_RKIND-w1
            p8w(1,kte+1,1) = exp(w1*log(p3d(1,kte,1)) + w2*log(p3d(1,kte-1,1)))
            t8w(1,kte+1,1) = w1*t3d(1,kte,1) + w2*t3d(1,kte-1,1)
            p8w(1,k,1)   = 0.5_RKIND*(p3d(1,k,1)+p3d(1,k-1,1))
            t8w(1,k,1)   = 0.5_RKIND*(t3d(1,k,1)+t3d(1,k-1,1))
          end if
        end do

        do k = kts, kte
          dz8w(1,k,1) = (p8w(1,k,1)-p8w(1,k+1,1)) / (gravity*(p3d(1,k,1)/(287.0_RKIND*t3d(1,k,1))))
        end do
      end block

      tsk(1,1)    = col%skt
      albedo(1,1) = col%albedo
      emiss(1,1)  = 0.98_RKIND
      xland(1,1)  = 2.0_RKIND - min(max(col%lsm,0.0_RKIND),1.0_RKIND)
      xice(1,1)   = col%siconc
      snow(1,1)   = col%snow_kg_m2
      xlat(1,1)   = col%lat
      xlong(1,1)  = col%lon
      coszr(1,1)  = 0.0_RKIND
      gsw = 0.0_RKIND; swcf = 0.0_RKIND
      rthratensw3 = 0.0_RKIND; rthratenlw3 = 0.0_RKIND

      ! swupflx/swdnflx/lwupflx/lwdnflx (and their *c clear-sky
      ! counterparts) give the per-level spectral SW/LW up/down flux
      ! profile RRTMG computes internally. RRTMG extends the column
      ! with synthetic extra layers above the model top (fixed ~4 hPa
      ! increments down to ~0 hPa -- see "nlayers = kte +
      ! max(nint(pw1d(kte+1)/deltap),1)" in module_ra_rrtmg_sw.F/_lw.F)
      ! purely to fill these arrays out to index nlayers+1; kme was
      ! sized above (nlev+1+nlayers_extra+2) specifically to cover
      ! that, using the same formula computed from the raw ERA5 column
      ! before any array existed. Only indices kts:kte+1 (the physical
      ! column) are read back below -- the synthetic extra layers
      ! beyond kte+1 are discarded. A segmentation fault/heap
      ! corruption from writing past a too-small array was observed
      ! here before kme was sized this way (see DADOS_ECRAD.md).
      call rrtmg_swrad(p3d=p3d,p8w=p8w,pi3d=pi3d,t3d=t3d,t8w=t8w,dz8w=dz8w,        &
           &  qv3d=qv3d,qc3d=qc3d,qr3d=qr3d,qi3d=qi3d,qs3d=qs3d,qg3d=qg3d,          &
           &  cldfra3d=cldfra3d,o33d=o33d,tsk=tsk,albedo=albedo,                   &
           &  xland=xland,xice=xice,snow=snow,coszr=coszr,xtime=0.0_RKIND,gmt=gmt, &
           &  julday=julday,radt=60.0_RKIND,degrad=degrad_local,declin=declin,      &
           &  solcon=solcon,xlat=xlat,xlong=xlong,icloud=1,                        &
           &  cldovrlp=2,idcor=1,o3input=0,                                        &
           &  noznlevels=1,gsw=gsw,swcf=swcf,rthratensw=rthratensw3,               &
           &  has_reqc=0,has_reqi=0,has_reqs=0,re_cloud=re_cloud,                  &
           &  re_ice=re_ice,re_snow=re_snow,                                       &
           &  aer_opt=0,                                                           &
           &  swupflx=swupflx,swupflxc=swupflxc,swdnflx=swdnflx,swdnflxc=swdnflxc, &
           &  ids=ids,ide=ide,jds=jds,jde=jde,kds=kds,kde=kde,                     &
           &  ims=ims,ime=ime,jms=jms,jme=jme,kms=kms,kme=kme,                     &
           &  its=its,ite=ite,jts=jts,jte=jte,kts=kts,kte=kte)

      call rrtmg_lwrad(p3d=p3d,p8w=p8w,pi3d=pi3d,t3d=t3d,t8w=t8w,dz8w=dz8w,        &
           &  qv3d=qv3d,qc3d=qc3d,qr3d=qr3d,qi3d=qi3d,qs3d=qs3d,qg3d=qg3d,          &
           &  cldfra3d=cldfra3d,o33d=o33d,tsk=tsk,emiss=emiss,                     &
           &  xland=xland,xice=xice,snow=snow,xlat=xlat,julday=julday,icloud=1,    &
           &  cldovrlp=2,idcor=1,o3input=0,noznlevels=1,                           &
           &  glw=glw,olr=olr,lwcf=lwcf,rthratenlw=rthratenlw3,                    &
           &  has_reqc=0,has_reqi=0,has_reqs=0,re_cloud=re_cloud,                  &
           &  re_ice=re_ice,re_snow=re_snow,                                       &
           &  rre_cloud=rre_cloud3d,rre_ice=rre_ice3d,rre_snow=rre_snow3d,         &
           &  lwupflx=lwupflx,lwupflxc=lwupflxc,lwdnflx=lwdnflx,lwdnflxc=lwdnflxc,  &
           &  ids=ids,ide=ide,jds=jds,jde=jde,kds=kds,kde=kde,                     &
           &  ims=ims,ime=ime,jms=jms,jme=jme,kms=kms,kme=kme,                     &
           &  its=its,ite=ite,jts=jts,jte=jte,kts=kts,kte=kte)

      ! rrtmg's arrays are WRF-order (k=1 nearest surface); flip back
      ! to ERA5's own top-to-bottom order so both schemes' output uses
      ! the same level convention. Interfaces (nlev+1 values): ERA5
      ! index k maps to WRF index nlev+2-k (k=1/TOA <-> WRF kte+1;
      ! k=nlev+1/surface <-> WRF index 1) -- one more than the
      ! full-level flip (nlev-k+1) since there is one extra level.
      do k = 1, nlev
        if(ABS(rthratensw3(1,nlev-k+1,1)) > 50.0) print*,'sw ',k,rthratensw3(1,nlev-k+1,1)  
        if(ABS(rthratenlw3(1,nlev-k+1,1)) > 50.0) print*,'lw ',k,rthratenlw3(1,nlev-k+1,1)  
        rthratensw(k) = rthratensw3(1,nlev-k+1,1)
        rthratenlw(k) = rthratenlw3(1,nlev-k+1,1)
      end do
      do k = 1, nlev+1
        if(ABS(lwupflx(1,nlev-k+2,1)) > 5000.0) print*,'swup ',k,rthratensw3(1,nlev-k+2,1)  
        if(ABS(lwdnflx(1,nlev-k+2,1)) > 5000.0) print*,'swdn ',k,rthratenlw3(1,nlev-k+2,1)  
        flux_sw_up(k) = swupflx(1,nlev+2-k,1)
        flux_sw_dn(k) = swdnflx(1,nlev+2-k,1)
        flux_lw_up(k) = lwupflx(1,nlev+2-k,1)
        flux_lw_dn(k) = lwdnflx(1,nlev+2-k,1)
        flux_sw_up_clear(k) = swupflxc(1,nlev+2-k,1)
        flux_sw_dn_clear(k) = swdnflxc(1,nlev+2-k,1)
        flux_lw_up_clear(k) = lwupflxc(1,nlev+2-k,1)
        flux_lw_dn_clear(k) = lwdnflxc(1,nlev+2-k,1)
      end do
      allocate(reff_cloud(nlev), reff_ice(nlev))
      do k = 1, nlev
        reff_cloud(k) = rre_cloud3d(1,nlev-k+1,1)
        reff_ice(k) = rre_ice3d(1,nlev-k+1,1)
      end do

      deallocate(p3d,pi3d,t3d,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,cldfra3d,o33d)
      deallocate(re_cloud,re_ice,re_snow,dz8w,p8w,t8w,rthratensw3,rthratenlw3)
      deallocate(swupflx,swupflxc,swdnflx,swdnflxc,lwupflx,lwupflxc,lwdnflx,lwdnflxc)
      deallocate(rre_cloud3d,rre_ice3d,rre_snow3d)

    end select

    ! -------------------------------------------------------------
        ! Write output (same NetCDF file/variables regardless of scheme)
        ! -------------------------------------------------------------
        
        if (trim(ctl%radiation_scheme) == 'rrtmg') then
          call output_write_timestep(outf, ilon_idx, ilat_idx, it_index, nlev, &
               &  era5f%valid_time(itime), &
               &  col%t(1:nlev), col%qv(1:nlev), col%z(1:nlev), &
               &  col%qc(1:nlev), col%qi(1:nlev), col%qr(1:nlev), col%qs(1:nlev), &
               &  col%cldfra(1:nlev), col%o3(1:nlev), col%sp, col%skt, &
               &  rthratensw, rthratenlw, flux_sw_up, flux_sw_dn, flux_lw_up, flux_lw_dn, &
               &  flux_sw_up_clear, flux_sw_dn_clear, flux_lw_up_clear, flux_lw_dn_clear, &
               &  reliq=reff_cloud, reice=reff_ice)
        else
          call output_write_timestep(outf, ilon_idx, ilat_idx, it_index, nlev, &
               &  era5f%valid_time(itime), &
               &  col%t(1:nlev), col%qv(1:nlev), col%z(1:nlev), &
               &  col%qc(1:nlev), col%qi(1:nlev), col%qr(1:nlev), col%qs(1:nlev), &
               &  col%cldfra(1:nlev), col%o3(1:nlev), col%sp, col%skt, &
               &  rthratensw, rthratenlw, flux_sw_up, flux_sw_dn, flux_lw_up, flux_lw_dn, &
               &  flux_sw_up_clear, flux_sw_dn_clear, flux_lw_up_clear, flux_lw_dn_clear)
        end if

        deallocate(rthratensw, rthratenlw, flux_sw_up, flux_sw_dn, flux_sw_up_clear, flux_sw_dn_clear, flux_lw_up, flux_lw_dn, flux_lw_up_clear, flux_lw_dn_clear)
        if (allocated(reff_cloud)) deallocate(reff_cloud, reff_ice)
        deallocate(col%p,col%t,col%z,col%qv,col%qc,col%qi,col%qr,col%qs,col%cldfra,col%o3)

      end do
    end do
  end do

  call output_close(outf)
  call era5_close(era5f)

end program rad_offline_driver
