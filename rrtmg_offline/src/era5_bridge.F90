! era5_bridge.F90
!
! Reads ERA5 hourly pressure-level + single-level NetCDF files (as
! downloaded from the Copernicus Climate Data Store, CDS) and builds
! the plain Fortran arrays rrtmg_swrad()/rrtmg_lwrad() expect for one
! grid column at one valid time -- i.e. this is the offline,
! ERA5-driven replacement for MPAS's mpas_atmphys_interface.F +
! mpas_atmphys_driver_radiation_sw.F/_lw.F "*_from_MPAS" packing
! subroutines.
!
! Expected input files (CDS "reanalysis-era5-pressure-levels" and
! "reanalysis-era5-single-levels", NetCDF format):
!
!   Pressure-level file, variables (CDS short names):
!     t     - temperature                          [K]
!     q     - specific humidity                    [kg/kg]
!     z     - geopotential                         [m2/s2]
!     clwc  - specific cloud liquid water content   [kg/kg]
!     ciwc  - specific cloud ice water content      [kg/kg]
!     crwc  - specific rain water content           [kg/kg]
!     cswc  - specific snow water content           [kg/kg]
!     cc    - fraction of cloud cover               [0-1]
!     o3    - ozone mass mixing ratio                [kg/kg]
!   (all on dims (longitude,latitude,pressure_level,valid_time) or
!    (longitude,latitude,level,time) in netCDF's own storage order,
!    which is the REVERSE of Fortran's -- see read_column below)
!
!   Single-level file, variables:
!     sp      - surface pressure          [Pa]
!     skt     - skin temperature          [K]
!     fal     - forecast albedo           [0-1]
!     siconc  - sea ice area fraction     [0-1]  (older files may call
!               this "ci" instead -- both are tried)
!     sd      - snow depth (water equiv.) [m of water equivalent]
!     lsm     - land-sea mask             [0-1]  (invariant; some
!               single-level requests omit it -- see NOTE below)
!     z       - surface geopotential      [m2/s2] (invariant; same
!               name as the pressure-level geopotential, but this one
!               has no pressure_level/level dimension)
!
! NOTE on lsm/orography: these are INVARIANT fields; CDS treats them
! as ordinary time-varying variables that happen not to change, so
! request them for just one arbitrary hour/day if downloading them
! together with genuinely time-varying single-level fields, or
! download them separately (they are also available from the
! "reanalysis-era5-single-levels-monthly-means" invariant collection).
!
! NOTE on CDS NetCDF format drift: CDS changed its default netCDF
! conversion in November 2024. Both the dimension/coordinate names
! (valid_time/pressure_level vs. the older time/level) and the time
! units reference date changed. This module does NOT assume either
! convention: it probes for both dimension-name variants, and it
! parses the actual "<units> since <reference-date>" string from the
! time coordinate's attribute at run time rather than hard-coding an
! epoch. If you would rather have a stable, simpler format, you can
! also request data_format="netcdf_legacy" from the CDS API, which
! restores the pre-Nov-2024 time/level naming and the 1900-01-01
! epoch; this module handles that too.
!
module era5_bridge

  use mpas_kind_types, only : RKIND
  use netcdf

  implicit none
  private

  public :: era5_file_type
  public :: era5_open, era5_close
  public :: era5_column_type
  public :: era5_read_column
  public :: era5_valid_time_to_julian_gmt

  integer, parameter :: NOZNLEVELS_UNUSED = 0  ! this bridge does not use
                                                ! the CAM ozone-climatology
                                                ! path (o3input==2); real
                                                ! ERA5 ozone is supplied
                                                ! directly (see the
                                                ! deliberate o3input.ne.2
                                                ! modification documented
                                                ! in module_ra_rrtmg_sw.F/
                                                ! _lw.F in this same
                                                ! source tree)

  ! One open ERA5 file pair (pressure-level + single-level), with the
  ! bits of metadata (dimension ids/lengths, variable ids, grid
  ! coordinates, time units) needed to pull out individual columns.
  type :: era5_file_type
    integer :: ncid_pl = -1, ncid_sl = -1
    integer :: nlon = 0, nlat = 0, nlev = 0, ntime = 0
    real(kind=RKIND), allocatable :: lon(:), lat(:), plev_hpa(:)
    real(kind=RKIND), allocatable :: valid_time(:)   ! raw numeric time
                                                       ! coordinate values
    character(len=256) :: time_units = ''             ! e.g. "hours since
                                                       ! 1970-1-1 00:00:00"
    ! Some ERA5/CDS downloads store the pressure_level coordinate
    ! descending (surface-to-TOA, e.g. 1000,975,...,1 hPa) instead of
    ! the ascending (TOA-to-surface) order this whole bridge assumes
    ! throughout (era5_column_type below, era5_ecrad_bridge.F90,
    ! rad_offline_driver.F90's RRTMG branch). era5_open detects this
    ! once from plev_hpa's own on-disk order and reverses it in place;
    ! era5_read_column then reverses every level-dimension profile it
    ! reads (t/q/z/clwc/...) the same way, so nothing downstream of
    ! this module ever needs to know or care which order the source
    ! file actually used. Getting this wrong silently is NOT a
    ! cosmetic issue: it was the root cause of a real segfault (a
    ! bogus top-of-atmosphere pressure extrapolated from what was
    ! actually the near-surface level) -- see DADOS_ERA5.md.
    logical :: level_reversed = .false.
    ! Variable ids, pressure-level file (0 = not found/not requested)
    integer :: id_t=0, id_q=0, id_z=0, id_clwc=0, id_ciwc=0, id_crwc=0, &
         &     id_cswc=0, id_cc=0, id_o3=0
    ! Variable ids, single-level file
    integer :: id_sp=0, id_skt=0, id_fal=0, id_siconc=0, id_sd=0, &
         &     id_lsm=0, id_zsfc=0
  end type era5_file_type

  ! One column (one lon/lat/time) of ERA5 data, in TOP-TO-BOTTOM
  ! pressure order (plev(1) = highest altitude/lowest pressure,
  ! plev(nlev) = lowest altitude/highest pressure) -- i.e. exactly the
  ! order the CDS pressure_level coordinate is normally stored in.
  ! Levels at or below the surface (plev(k) >= sp) are still filled
  ! (with whatever CDS's extrapolated/undefined value is there) but
  ! flagged invalid via nlev_valid; rrtmg_offline_driver.F90 only
  ! passes the nlev_valid levels strictly above ground into
  ! rrtmg_swrad/rrtmg_lwrad.
  type :: era5_column_type
    integer :: nlev = 0, nlev_valid = 0
    real(kind=RKIND), allocatable :: p(:)      ! Pa
    real(kind=RKIND), allocatable :: t(:)      ! K
    real(kind=RKIND), allocatable :: z(:)      ! m (geopotential HEIGHT,
                                                ! already divided by g)
    real(kind=RKIND), allocatable :: qv(:), qc(:), qi(:), qr(:), qs(:) ! kg/kg
    real(kind=RKIND), allocatable :: cldfra(:) ! 0-1
    real(kind=RKIND), allocatable :: o3(:)     ! kg/kg (mass mixing ratio)
    real(kind=RKIND) :: sp = 0.0_RKIND         ! Pa
    real(kind=RKIND) :: skt = 0.0_RKIND        ! K
    real(kind=RKIND) :: albedo = 0.0_RKIND     ! 0-1
    real(kind=RKIND) :: siconc = 0.0_RKIND     ! 0-1
    real(kind=RKIND) :: snow_kg_m2 = 0.0_RKIND ! kg/m2 water equivalent
    real(kind=RKIND) :: lsm = 0.0_RKIND        ! 0-1 (land fraction)
    real(kind=RKIND) :: zsfc = 0.0_RKIND       ! m (surface height)
    real(kind=RKIND) :: lat = 0.0_RKIND        ! degrees north
    real(kind=RKIND) :: lon = 0.0_RKIND        ! degrees east
  end type era5_column_type

contains

  !=================================================================
  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(0,'(A)') 'era5_bridge NetCDF error at '//trim(where)//': '// &
           &          trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

  ! Look for a dimension under either of two candidate names (to
  ! handle the pre-/post-Nov-2024 CDS naming change); returns the
  ! dimid and length of whichever is found, or dimid=-1 if neither
  ! exists (caller decides whether that is fatal).
  subroutine find_dim(ncid, name_a, name_b, dimid, length)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name_a, name_b
    integer, intent(out) :: dimid, length
    integer :: status
    status = nf90_inq_dimid(ncid, name_a, dimid)
    if (status /= nf90_noerr) status = nf90_inq_dimid(ncid, name_b, dimid)
    if (status /= nf90_noerr) then
      dimid = -1; length = 0
      return
    end if
    call nc_check(nf90_inquire_dimension(ncid, dimid, len=length), &
         &        'find_dim inquire '//trim(name_a))
  end subroutine find_dim

  ! Variable id lookup that tolerates a missing/optional variable
  ! (returns 0, not fatal) -- used for fields this bridge can do
  ! without (e.g. siconc over an all-land domain).
  function find_var_optional(ncid, name) result(varid)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer :: varid, status
    status = nf90_inq_varid(ncid, name, varid)
    if (status /= nf90_noerr) varid = 0
  end function find_var_optional

  ! CDS/ERA5 pressure-level downloads use the NetCDF default
  ! missing_value sentinel (~9.969209968386869e36) for
  ! level/variable combinations that were not computed -- this is
  ! independent of whether the level is numerically above or below
  ! this column's surface pressure (nlev_valid, computed later from
  ! p(k) < sp, does NOT catch it: a level can be "above ground" and
  ! still carry this sentinel for clwc/ciwc/crwc/cswc). It is a huge
  ! POSITIVE number, so max(buf,0.0_RKIND) does not catch it either
  ! -- that clamp only guards against small negative noise. Left
  ! unmasked it is read as an ~1e33-1e36 kg/kg condensate load,
  ! which blows up cloud optical depth and diverges the two-stream
  ! flux solve level by level (huge, still-growing lw_up/lw_dn).
  ! Treat anything at or beyond this threshold as "no data", the
  ! same way an altogether-absent variable is already handled above
  ! (the id_* == 0 branches): 0.0 for a mixing ratio/fraction is the
  ! physically reasonable "nothing known here" default.
  elemental function mask_era5_missing(x) result(y)
    real(kind=RKIND), intent(in) :: x
    real(kind=RKIND) :: y
    real(kind=RKIND), parameter :: era5_missing_threshold = 1.0e30_RKIND
    if (abs(x) >= era5_missing_threshold) then
      y = 0.0_RKIND
    else
      y = x
    end if
  end function mask_era5_missing


  !=================================================================
  ! era5_open
  !
  ! Opens the pressure-level and single-level files and reads
  ! coordinates + variable ids. pressure_level values are converted
  ! to hPa here (double-checking the file's own "units" attribute
  ! rather than assuming); everything else stays in file units until
  ! era5_read_column.
  !=================================================================
  subroutine era5_open(path_pressure_level, path_single_level, f)

    character(len=*), intent(in) :: path_pressure_level, path_single_level
    type(era5_file_type), intent(out) :: f

    integer :: dimid_lon, dimid_lat, dimid_lev, dimid_time
    integer :: varid_lon, varid_lat, varid_lev, varid_time
    character(len=256) :: lev_units

    call nc_check(nf90_open(trim(path_pressure_level), nf90_nowrite, f%ncid_pl), &
         &        'open pressure-level file')
    call nc_check(nf90_open(trim(path_single_level), nf90_nowrite, f%ncid_sl), &
         &        'open single-level file')

    call find_dim(f%ncid_pl, 'longitude', 'lon', dimid_lon, f%nlon)
    call find_dim(f%ncid_pl, 'latitude',  'lat', dimid_lat, f%nlat)
    call find_dim(f%ncid_pl, 'pressure_level', 'level', dimid_lev, f%nlev)
    call find_dim(f%ncid_pl, 'valid_time', 'time', dimid_time, f%ntime)
    if (dimid_lev < 0)  stop 'era5_open: no pressure_level/level dimension found'
    if (dimid_time < 0) stop 'era5_open: no valid_time/time dimension found'

    allocate(f%lon(f%nlon), f%lat(f%nlat), f%plev_hpa(f%nlev), f%valid_time(f%ntime))

    call nc_check(nf90_inq_varid(f%ncid_pl, 'longitude', varid_lon), 'inq lon')
    call nc_check(nf90_get_var(f%ncid_pl, varid_lon, f%lon), 'get lon')
    call nc_check(nf90_inq_varid(f%ncid_pl, 'latitude', varid_lat), 'inq lat')
    call nc_check(nf90_get_var(f%ncid_pl, varid_lat, f%lat), 'get lat')

    if (nf90_inq_varid(f%ncid_pl, 'pressure_level', varid_lev) /= nf90_noerr) &
         &  call nc_check(nf90_inq_varid(f%ncid_pl, 'level', varid_lev), 'inq level')
    call nc_check(nf90_get_var(f%ncid_pl, varid_lev, f%plev_hpa), 'get level')
    lev_units = ''
    if (nf90_get_att(f%ncid_pl, varid_lev, 'units', lev_units) == nf90_noerr) then
      if (index(lev_units,'Pa') > 0 .and. index(lev_units,'hPa') == 0) then
        f%plev_hpa = f%plev_hpa / 100.0_RKIND   ! file was in Pa, not hPa
      end if
    end if

    ! Detect and normalize level order -- see the level_reversed
    ! comment on era5_file_type above for why this matters.
    if (f%nlev >= 2) then
      f%level_reversed = (f%plev_hpa(1) > f%plev_hpa(f%nlev))
      if (f%level_reversed) then
        f%plev_hpa = f%plev_hpa(f%nlev:1:-1)
        write(*,'(A)') 'era5_bridge: pressure_level was stored surface-to-TOA in this '// &
             &  'file -- reversing to TOA-to-surface order internally (see era5_bridge.F90).'
      end if
    end if

    if (nf90_inq_varid(f%ncid_pl, 'valid_time', varid_time) /= nf90_noerr) &
         &  call nc_check(nf90_inq_varid(f%ncid_pl, 'time', varid_time), 'inq time')
    call nc_check(nf90_get_var(f%ncid_pl, varid_time, f%valid_time), 'get time')
    call nc_check(nf90_get_att(f%ncid_pl, varid_time, 'units', f%time_units), &
         &        'get time units')

    f%id_t    = find_var_optional(f%ncid_pl, 't')
    f%id_q    = find_var_optional(f%ncid_pl, 'q')
    f%id_z    = find_var_optional(f%ncid_pl, 'z')
    f%id_clwc = find_var_optional(f%ncid_pl, 'clwc')
    f%id_ciwc = find_var_optional(f%ncid_pl, 'ciwc')
    f%id_crwc = find_var_optional(f%ncid_pl, 'crwc')
    f%id_cswc = find_var_optional(f%ncid_pl, 'cswc')
    f%id_cc   = find_var_optional(f%ncid_pl, 'cc')
    f%id_o3   = find_var_optional(f%ncid_pl, 'o3')
    if (f%id_t == 0) stop 'era5_open: required variable "t" not found in pressure-level file'
    if (f%id_q == 0) stop 'era5_open: required variable "q" not found in pressure-level file'
    if (f%id_z == 0) stop 'era5_open: required variable "z" not found in pressure-level file'

    f%id_sp     = find_var_optional(f%ncid_sl, 'sp')
    f%id_skt    = find_var_optional(f%ncid_sl, 'skt')
    f%id_fal    = find_var_optional(f%ncid_sl, 'fal')
    f%id_siconc = find_var_optional(f%ncid_sl, 'siconc')
    if (f%id_siconc == 0) f%id_siconc = find_var_optional(f%ncid_sl, 'ci')
    f%id_sd     = find_var_optional(f%ncid_sl, 'sd')
    f%id_lsm    = find_var_optional(f%ncid_sl, 'lsm')
    f%id_zsfc   = find_var_optional(f%ncid_sl, 'z')
    if (f%id_sp == 0)  stop 'era5_open: required variable "sp" not found in single-level file'
    if (f%id_skt == 0) stop 'era5_open: required variable "skt" not found in single-level file'

  end subroutine era5_open


  subroutine era5_close(f)
    type(era5_file_type), intent(inout) :: f
    if (f%ncid_pl >= 0) call nc_check(nf90_close(f%ncid_pl), 'close pressure-level file')
    if (f%ncid_sl >= 0) call nc_check(nf90_close(f%ncid_sl), 'close single-level file')
  end subroutine era5_close


  !=================================================================
  ! era5_read_column
  !
  ! Reads one (ilon,ilat,itime) column. NOTE on axis order: netCDF
  ! stores arrays in C/row-major order, so a file whose ncdump header
  ! lists dimensions as (valid_time, pressure_level, latitude,
  ! longitude) is read by nf90_get_var, when you ask for a single
  ! (lon,lat) point across all levels, with a start/count vector
  ! ordered (longitude,latitude,pressure_level,valid_time) -- i.e.
  ! REVERSED from the ncdump header order. This is standard
  ! netCDF-Fortran behaviour, not a bridge-specific quirk; it is
  ! implemented once here so nothing upstream needs to worry about it.
  !=================================================================
  subroutine era5_read_column(f, ilon, ilat, itime, col)

    type(era5_file_type), intent(in) :: f
    integer, intent(in) :: ilon, ilat, itime
    type(era5_column_type), intent(out) :: col

    integer :: start4(4), count4(4)
    integer :: start3(3), count3(3)
    real(kind=RKIND) :: buf(f%nlev)
    real(kind=RKIND) :: sbuf(1)

    col%nlev = f%nlev
    allocate(col%p(col%nlev), col%t(col%nlev), col%z(col%nlev), &
         &   col%qv(col%nlev), col%qc(col%nlev), col%qi(col%nlev), &
         &   col%qr(col%nlev), col%qs(col%nlev), col%cldfra(col%nlev), &
         &   col%o3(col%nlev))

    start4 = (/ ilon, ilat, 1, itime /)
    count4 = (/ 1, 1, f%nlev, 1 /)

    col%p = f%plev_hpa * 100.0_RKIND   ! hPa -> Pa, same profile every column
                                        ! (already reversed in era5_open if needed)

    call nc_check(nf90_get_var(f%ncid_pl, f%id_t, buf, start4, count4), 'read t')
    if (f%level_reversed) buf = buf(f%nlev:1:-1)
    col%t = buf

    call nc_check(nf90_get_var(f%ncid_pl, f%id_q, buf, start4, count4), 'read q')
    if (f%level_reversed) buf = buf(f%nlev:1:-1)
    col%qv = max(buf, 1.0e-12_RKIND)

    call nc_check(nf90_get_var(f%ncid_pl, f%id_z, buf, start4, count4), 'read z')
    if (f%level_reversed) buf = buf(f%nlev:1:-1)
    ! Height (m) from geopotential (m2/s2); reuse the RRTMG constants
    ! module's own gravity value for exact consistency with the rest
    ! of the driver.
    block
      use mpas_atmphys_constants, only : gravity
      col%z = buf / gravity
    end block

    if (f%id_clwc > 0) then
      call nc_check(nf90_get_var(f%ncid_pl, f%id_clwc, buf, start4, count4), 'read clwc')
      if (f%level_reversed) buf = buf(f%nlev:1:-1)
      col%qc = max(mask_era5_missing(buf),0.0_RKIND)
    else
      col%qc = 0.0_RKIND
    end if
    if (f%id_ciwc > 0) then
      call nc_check(nf90_get_var(f%ncid_pl, f%id_ciwc, buf, start4, count4), 'read ciwc')
      if (f%level_reversed) buf = buf(f%nlev:1:-1)
      col%qi = max(mask_era5_missing(buf),0.0_RKIND)
    else
      col%qi = 0.0_RKIND
    end if
    if (f%id_crwc > 0) then
      call nc_check(nf90_get_var(f%ncid_pl, f%id_crwc, buf, start4, count4), 'read crwc')
      if (f%level_reversed) buf = buf(f%nlev:1:-1)
      col%qr = max(mask_era5_missing(buf),0.0_RKIND)
    else
      col%qr = 0.0_RKIND
    end if
    if (f%id_cswc > 0) then
      call nc_check(nf90_get_var(f%ncid_pl, f%id_cswc, buf, start4, count4), 'read cswc')
      if (f%level_reversed) buf = buf(f%nlev:1:-1)
      col%qs = max(mask_era5_missing(buf),0.0_RKIND)
    else
      col%qs = 0.0_RKIND
    end if
    if (f%id_cc > 0) then
      call nc_check(nf90_get_var(f%ncid_pl, f%id_cc, buf, start4, count4), 'read cc')
      if (f%level_reversed) buf = buf(f%nlev:1:-1)
      col%cldfra = min(max(mask_era5_missing(buf),0.0_RKIND),1.0_RKIND)
    else
      ! Diagnose a 0/1 cloud flag from condensate if cc itself was not
      ! downloaded -- crude, but keeps the driver usable with a
      ! smaller ERA5 request.
      col%cldfra = merge(1.0_RKIND, 0.0_RKIND, (col%qc+col%qi+col%qr+col%qs) > 1.0e-8_RKIND)
    end if
    if (f%id_o3 > 0) then
      call nc_check(nf90_get_var(f%ncid_pl, f%id_o3, buf, start4, count4), 'read o3')
      if (f%level_reversed) buf = buf(f%nlev:1:-1)
      col%o3 = max(mask_era5_missing(buf),0.0_RKIND)
    else
      col%o3 = 0.0_RKIND
    end if

    start3 = (/ ilon, ilat, itime /)
    count3 = (/ 1, 1, 1 /)
    call nc_check(nf90_get_var(f%ncid_sl, f%id_sp, sbuf, start3, count3), 'read sp'); col%sp = sbuf(1)
    call nc_check(nf90_get_var(f%ncid_sl, f%id_skt, sbuf, start3, count3), 'read skt'); col%skt = sbuf(1)

    if (f%id_fal > 0) then
      call nc_check(nf90_get_var(f%ncid_sl, f%id_fal, sbuf, start3, count3), 'read fal'); col%albedo = sbuf(1)
    else
      col%albedo = 0.20_RKIND  ! generic fallback if forecast albedo was not requested
    end if
    if (f%id_siconc > 0) then
      call nc_check(nf90_get_var(f%ncid_sl, f%id_siconc, sbuf, start3, count3), 'read siconc'); col%siconc = sbuf(1)
    else
      col%siconc = 0.0_RKIND
    end if
    if (f%id_sd > 0) then
      call nc_check(nf90_get_var(f%ncid_sl, f%id_sd, sbuf, start3, count3), 'read sd'); col%snow_kg_m2 = sbuf(1)*1000.0_RKIND
    else
      col%snow_kg_m2 = 0.0_RKIND
    end if

    if (f%id_lsm > 0) then
      ! lsm is invariant: read at time index 1 regardless of itime, in
      ! case the caller only downloaded it for a single timestamp.
      call nc_check(nf90_get_var(f%ncid_sl, f%id_lsm, sbuf, &
           &  (/ilon,ilat,1/), (/1,1,1/)), 'read lsm'); col%lsm = sbuf(1)
    else
      col%lsm = 1.0_RKIND  ! default to "land" if not supplied
    end if
    if (f%id_zsfc > 0) then
      call nc_check(nf90_get_var(f%ncid_sl, f%id_zsfc, sbuf, &
           &  (/ilon,ilat,1/), (/1,1,1/)), 'read surface z')
      block
        use mpas_atmphys_constants, only : gravity
        col%zsfc = sbuf(1) / gravity
      end block
    else
      col%zsfc = 0.0_RKIND
    end if

    col%lon = f%lon(ilon)
    col%lat = f%lat(ilat)

    ! A pressure level is a valid ATMOSPHERIC level only if it is
    ! above the surface. CDS pressure-level fields still contain
    ! (extrapolated/unphysical) values below ground for
    ! high-topography columns, so this count must be redone per
    ! column, not assumed fixed.
    block
      integer :: k
      col%nlev_valid = 0
      do k = 1, col%nlev
        if (col%p(k) < col%sp) col%nlev_valid = col%nlev_valid + 1
      end do
    end block

  end subroutine era5_read_column


  !=================================================================
  ! era5_valid_time_to_julian_gmt
  !
  ! Parses the time coordinate's "units" attribute (e.g. "hours since
  ! 1970-1-1 00:00:00" or "hours since 1900-01-01 00:00:00") and the
  ! raw numeric value at index itime, and returns:
  !   julian - fractional day of year (1.0 = Jan 1st 00:00 UTC),
  !            exactly what radconst_offline's radconst() expects
  !   gmt    - Greenwich Mean Time, decimal hours (0.0-24.0)
  !   julday - integer day of year (for rrtmg_*rad's julday argument)
  ! Deliberately does not depend on any external date library --
  ! everything ERA5 needs is a proleptic Gregorian calendar, which a
  ! small hand-rolled day-count function handles exactly.
  !=================================================================
  subroutine era5_valid_time_to_julian_gmt(f, itime, julian, gmt, julday)

    type(era5_file_type), intent(in) :: f
    integer, intent(in) :: itime
    real(kind=RKIND), intent(out) :: julian, gmt
    integer, intent(out) :: julday

    integer :: ref_year, ref_month, ref_day, ref_hour, ref_min
    real(kind=RKIND) :: ref_sec, step_hours, total_hours
    integer :: cur_year, cur_month, cur_day
    real(kind=RKIND) :: hours_into_day
    character(len=32) :: unit_word

    call parse_time_units(f%time_units, unit_word, ref_year, ref_month, &
         &                 ref_day, ref_hour, ref_min, ref_sec)

    select case (trim(unit_word))
    case ('hours')
      step_hours = 1.0_RKIND
    case ('minutes')
      step_hours = 1.0_RKIND/60.0_RKIND
    case ('seconds')
      step_hours = 1.0_RKIND/3600.0_RKIND
    case ('days')
      step_hours = 24.0_RKIND
    case default
      stop 'era5_valid_time_to_julian_gmt: unrecognised time units (expected hours/minutes/seconds/days since ...)'
    end select

    total_hours = f%valid_time(itime)*step_hours &
         &      + real(ref_hour,RKIND) + real(ref_min,RKIND)/60.0_RKIND + ref_sec/3600.0_RKIND

    call add_days_from_hours(ref_year, ref_month, ref_day, total_hours, &
         &                    cur_year, cur_month, cur_day, hours_into_day)

    gmt    = hours_into_day
    julday = day_of_year(cur_year, cur_month, cur_day)
    julian = real(julday, RKIND) + hours_into_day/24.0_RKIND

  end subroutine era5_valid_time_to_julian_gmt


  ! Parses "<word> since YYYY-MM-DD[ HH:MM:SS]" (also accepts a "T"
  ! between date and time, and 1- or 2-digit month/day, both seen in
  ! real ERA5 files across CDS format versions).
  subroutine parse_time_units(units, unit_word, year, month, day, hour, minute, sec)
    character(len=*), intent(in) :: units
    character(len=*), intent(out) :: unit_word
    integer, intent(out) :: year, month, day, hour, minute
    real(kind=RKIND), intent(out) :: sec

    integer :: p_since, p_dash1, p_dash2, p_colon1, p_colon2, io
    character(len=256) :: rest

    hour = 0; minute = 0; sec = 0.0_RKIND

    p_since = index(units, 'since')
    if (p_since == 0) stop 'parse_time_units: "since" not found in time units string'
    read(units(1:p_since-1), *, iostat=io) unit_word
    if (io /= 0) stop 'parse_time_units: could not read the unit word before "since"'

    rest = adjustl(units(p_since+5:))
    p_dash1 = index(rest, '-')
    p_dash2 = index(rest(p_dash1+1:), '-') + p_dash1
    read(rest(1:p_dash1-1), *, iostat=io) year
    read(rest(p_dash1+1:p_dash2-1), *, iostat=io) month
    ! day runs until a space or 'T' or end of string
    block
      integer :: p_end
      p_end = scan(rest(p_dash2+1:), ' T')
      if (p_end == 0) then
        read(rest(p_dash2+1:), *, iostat=io) day
        return
      else
        read(rest(p_dash2+1:p_dash2+p_end-1), *, iostat=io) day
      end if
    end block

    p_colon1 = index(rest, ':')
    if (p_colon1 == 0) return   ! date only, no time-of-day part
    p_colon2 = index(rest(p_colon1+1:), ':')
    read(rest(p_colon1-2:p_colon1-1), *, iostat=io) hour
    if (p_colon2 == 0) then
      read(rest(p_colon1+1:), *, iostat=io) minute
    else
      p_colon2 = p_colon2 + p_colon1
      read(rest(p_colon1+1:p_colon2-1), *, iostat=io) minute
      read(rest(p_colon2+1:), *, iostat=io) sec
    end if
  end subroutine parse_time_units


  logical function is_leap(y)
    integer, intent(in) :: y
    is_leap = (mod(y,4)==0 .and. mod(y,100)/=0) .or. mod(y,400)==0
  end function is_leap

  integer function days_in_month(y, m)
    integer, intent(in) :: y, m
    integer, parameter :: dim(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    days_in_month = dim(m)
    if (m == 2 .and. is_leap(y)) days_in_month = 29
  end function days_in_month

  integer function day_of_year(y, m, d)
    integer, intent(in) :: y, m, d
    integer :: mm
    day_of_year = d
    do mm = 1, m-1
      day_of_year = day_of_year + days_in_month(y, mm)
    end do
  end function day_of_year

  ! Adds a (possibly large, possibly fractional) number of hours to a
  ! reference calendar date, returning the resulting date and the
  ! fractional hours remaining within that final day. Simple
  ! whole-day-at-a-time stepping is used rather than a closed-form
  ! Julian-day conversion, since ERA5 spans are always at most a few
  ! decades -- clarity over cleverness here.
  subroutine add_days_from_hours(y0, m0, d0, hours, y, m, d, hours_into_day)
    integer, intent(in) :: y0, m0, d0
    real(kind=RKIND), intent(in) :: hours
    integer, intent(out) :: y, m, d
    real(kind=RKIND), intent(out) :: hours_into_day

    real(kind=RKIND) :: remaining_hours
    integer :: whole_days

    y = y0; m = m0; d = d0
    whole_days = int(hours/24.0_RKIND)
    remaining_hours = hours - real(whole_days,RKIND)*24.0_RKIND

    do while (whole_days > 0)
      d = d + 1
      if (d > days_in_month(y,m)) then
        d = 1; m = m + 1
        if (m > 12) then
          m = 1; y = y + 1
        end if
      end if
      whole_days = whole_days - 1
    end do
    ! (ERA5 reference dates are always at/near the start of an epoch,
    ! so "hours" here is never negative in practice; a negative-hours
    ! branch is intentionally not implemented.)

    hours_into_day = remaining_hours

  end subroutine add_days_from_hours

end module era5_bridge
