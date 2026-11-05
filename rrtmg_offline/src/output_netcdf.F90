! output_netcdf.F90
!
! Writes the driver's output as a single NetCDF file per run, one
! (longitude, latitude, valid_time) slice per processed ERA5 column and
! timestep, mirroring the input file's own variable names/units/
! dimension conventions (era5_bridge.F90) as closely as makes sense,
! plus the radiative quantities computed by whichever scheme (RRTMG or
! ecRad) was run:
!
!   Dimensions: longitude(nlon), latitude(nlat),
!               pressure_level(nlev_total), half_level(nlev_total+1),
!               valid_time(unlimited)
!
!   nlon/nlat are the full requested grid (ilon_start:ilon_end,
!   ilat_start:ilat_end -- see namelist_offline.F90); a single column
!   is simply nlon=nlat=1, the driver's original mode of operation.
!
!   nlev_total is the FULL number of pressure levels in the ERA5 input
!   file (era5_file_type%nlev), not any one column's above-ground
!   nlev_valid (era5_column_type%nlev_valid) -- those can legitimately
!   differ from column to column (surface height) and from hour to
!   hour (surface pressure), so every column/timestep is written using
!   its OWN n_valid, filling indices 1:n_valid along the vertical
!   dimension(s) and leaving n_valid+1:nlev_total untouched. NetCDF's
!   default fill value is left in those untouched cells (this module
!   never disables fill mode), so a column with fewer valid levels than
!   another at the same timestep is simply padded, not skipped or
!   truncated -- unlike an earlier single-column revision of this
!   driver, which had to skip an entire hour outright if its
!   nlev_valid differed from the first-processed hour's.
!
!   Echoed input variables (full levels, same names as the ERA5
!   pressure-level file): t, q, z (here: geometric HEIGHT in m, not
!   geopotential -- see the units/long_name attribute -- since that is
!   what era5_bridge.F90 itself already converts to and is more
!   directly usable), clwc, ciwc, crwc, cswc, cc, o3. Plus two scalar
!   fields from the single-level input file: sp, skt.
!
!   Radiative output: rthratensw, rthratenlw (full levels, K/s) and
!   sw_up/sw_dn/lw_up/lw_dn (half levels, W/m2 -- index 1 is the
!   TOA-most interface, index n_valid+1 the surface for THAT column,
!   so the surface and TOA fluxes are simply the first/last WRITTEN
!   value of each column's own slice, and every index in between is
!   the flux at that internal layer boundary).
!
! A single global attribute "radiation_scheme" records which scheme
! (rrtmg/ecrad) produced the file.
module output_netcdf

  use mpas_kind_types, only : RKIND
  use netcdf

  implicit none
  private

  public :: output_file_type
  public :: output_open, output_write_timestep, output_close

  type :: output_file_type
    integer :: ncid = -1
    integer :: nlon = 0, nlat = 0
    integer :: nlev_total = 0   ! FIXED (whole-input-file) level count;
                                 ! see the module header for why this is
                                 ! not any one column's nlev_valid
    ! variable ids
    integer :: id_t, id_q, id_z, id_clwc, id_ciwc, id_crwc, id_cswc, &
         &     id_cc, id_o3, id_sp, id_skt
    integer :: id_rthratensw, id_rthratenlw
    integer :: id_sw_up, id_sw_dn, id_lw_up, id_lw_dn
    integer :: id_valid_time
  end type output_file_type

contains

  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(0,'(A)') 'output_netcdf error at '//trim(where)//': '// &
           &          trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check


  !=================================================================
  ! output_open
  !
  ! Creates the output file and defines all dimensions/variables.
  ! lons(nlon)/lats(nlat)/pressure_level_hpa(nlev_total) are written
  ! immediately (they never change); valid_time and all per-timestep
  ! variables are left empty (at NetCDF's default fill value) until
  ! output_write_timestep is called once per (column, processed ERA5
  ! hour) pair.
  !=================================================================
  subroutine output_open(path, nlon, nlat, nlev_total, lons, lats, &
       &  pressure_level_hpa, time_units, radiation_scheme, f)

    character(len=*), intent(in) :: path
    integer, intent(in) :: nlon, nlat, nlev_total
    real(kind=RKIND), intent(in) :: lons(nlon), lats(nlat)
    real(kind=RKIND), intent(in) :: pressure_level_hpa(nlev_total)
    character(len=*), intent(in) :: time_units
    character(len=*), intent(in) :: radiation_scheme
    type(output_file_type), intent(out) :: f

    integer :: dimid_lon, dimid_lat, dimid_lev, dimid_hlev, dimid_time
    integer :: varid_lon, varid_lat, varid_lev
    integer :: dims3(3)   ! (lon,lat,time) for scalar surface fields
    integer :: dims4(4)   ! (lon,lat,level,time) for profile fields

    f%nlon = nlon
    f%nlat = nlat
    f%nlev_total = nlev_total

    call nc_check(nf90_create(trim(path), NF90_CLOBBER, f%ncid), 'create')

    call nc_check(nf90_def_dim(f%ncid, 'longitude',     nlon, dimid_lon),  'def longitude')
    call nc_check(nf90_def_dim(f%ncid, 'latitude',      nlat, dimid_lat),  'def latitude')
    call nc_check(nf90_def_dim(f%ncid, 'pressure_level', nlev_total, dimid_lev), 'def pressure_level')
    call nc_check(nf90_def_dim(f%ncid, 'half_level',    nlev_total+1, dimid_hlev), 'def half_level')
    call nc_check(nf90_def_dim(f%ncid, 'valid_time', NF90_UNLIMITED, dimid_time), 'def valid_time')

    call nc_check(nf90_def_var(f%ncid, 'longitude', NF90_DOUBLE, dimid_lon, varid_lon), 'def var lon')
    call nc_check(nf90_put_att(f%ncid, varid_lon, 'units', 'degrees_east'), 'att lon')
    call nc_check(nf90_def_var(f%ncid, 'latitude', NF90_DOUBLE, dimid_lat, varid_lat), 'def var lat')
    call nc_check(nf90_put_att(f%ncid, varid_lat, 'units', 'degrees_north'), 'att lat')
    call nc_check(nf90_def_var(f%ncid, 'pressure_level', NF90_DOUBLE, dimid_lev, varid_lev), 'def var plev')
    call nc_check(nf90_put_att(f%ncid, varid_lev, 'units', 'hPa'), 'att plev')
    call nc_check(nf90_put_att(f%ncid, varid_lev, 'long_name', &
         &  'fixed input pressure-level grid; a given column/hour may only have '// &
         &  'real data on a leading subset of these -- see the fill value on t/q/etc.'), 'att plev long_name')

    call nc_check(nf90_def_var(f%ncid, 'valid_time', NF90_DOUBLE, dimid_time, f%id_valid_time), 'def var time')
    call nc_check(nf90_put_att(f%ncid, f%id_valid_time, 'units', trim(time_units)), 'att time')

    dims4 = (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /)
    call def_profile_var(f%ncid, 't',    dims4, 'K',      'temperature', f%id_t)
    call def_profile_var(f%ncid, 'q',    dims4, 'kg kg-1','specific_humidity', f%id_q)
    call def_profile_var(f%ncid, 'z',    dims4, 'm',      'height (geometric, not geopotential)', f%id_z)
    call def_profile_var(f%ncid, 'clwc', dims4, 'kg kg-1','specific_cloud_liquid_water_content', f%id_clwc)
    call def_profile_var(f%ncid, 'ciwc', dims4, 'kg kg-1','specific_cloud_ice_water_content', f%id_ciwc)
    call def_profile_var(f%ncid, 'crwc', dims4, 'kg kg-1','specific_rain_water_content', f%id_crwc)
    call def_profile_var(f%ncid, 'cswc', dims4, 'kg kg-1','specific_snow_water_content', f%id_cswc)
    call def_profile_var(f%ncid, 'cc',   dims4, '0-1',    'fraction_of_cloud_cover', f%id_cc)
    call def_profile_var(f%ncid, 'o3',   dims4, 'kg kg-1','ozone_mass_mixing_ratio', f%id_o3)
    call def_profile_var(f%ncid, 'rthratensw', dims4, 'K s-1', &
         &  'shortwave radiative heating rate (theta/T tendency, uncoupled)', f%id_rthratensw)
    call def_profile_var(f%ncid, 'rthratenlw', dims4, 'K s-1', &
         &  'longwave radiative heating rate (theta/T tendency, uncoupled)', f%id_rthratenlw)

    dims3 = (/ dimid_lon, dimid_lat, dimid_time /)
    call def_scalar_var(f%ncid, 'sp',  dims3, 'Pa', 'surface_pressure', f%id_sp)
    call def_scalar_var(f%ncid, 'skt', dims3, 'K',  'skin_temperature', f%id_skt)

    dims4 = (/ dimid_lon, dimid_lat, dimid_hlev, dimid_time /)
    call def_profile_var(f%ncid, 'sw_up', dims4, 'W m-2', &
         &  'upwelling shortwave flux at layer interfaces (index 1=TOA, last valid=surface)', f%id_sw_up)
    call def_profile_var(f%ncid, 'sw_dn', dims4, 'W m-2', &
         &  'downwelling shortwave flux at layer interfaces (index 1=TOA, last valid=surface)', f%id_sw_dn)
    call def_profile_var(f%ncid, 'lw_up', dims4, 'W m-2', &
         &  'upwelling longwave flux at layer interfaces (index 1=TOA, last valid=surface)', f%id_lw_up)
    call def_profile_var(f%ncid, 'lw_dn', dims4, 'W m-2', &
         &  'downwelling longwave flux at layer interfaces (index 1=TOA, last valid=surface)', f%id_lw_dn)

    call nc_check(nf90_put_att(f%ncid, NF90_GLOBAL, 'radiation_scheme', trim(radiation_scheme)), 'att scheme')
    call nc_check(nf90_put_att(f%ncid, NF90_GLOBAL, 'source', &
         &  'rad_offline_driver (RRTMG/ecRad forced by ERA5)'), 'att source')

    call nc_check(nf90_enddef(f%ncid), 'enddef')

    call nc_check(nf90_put_var(f%ncid, varid_lon, lons), 'put lon')
    call nc_check(nf90_put_var(f%ncid, varid_lat, lats), 'put lat')
    call nc_check(nf90_put_var(f%ncid, varid_lev, pressure_level_hpa), 'put plev')

  end subroutine output_open


  ! Small helpers so output_open's variable list above stays readable.
  ! Every profile/scalar variable declares an explicit _FillValue
  ! matching NF90_FILL_DOUBLE, the same value NetCDF's default fill
  ! mode (left on for this file -- see output_open's nf90_create
  ! call, no NF90_NOFILL) silently pads with whenever a column writes
  ! fewer than the full nlev_total/nlev_total+1 levels (i.e. any
  ! column where n_valid < nlev_total -- real topography, not a
  ! bug: see the n_valid comment on output_write_timestep above).
  ! Without this attribute, readers like ncview/cdo have no way to
  ! know that value means "not written" rather than real data, and
  ! either misreport min/max or fail outright trying to convert it.
  subroutine def_profile_var(ncid, name, dims4, units, long_name, varid)
    integer, intent(in) :: ncid, dims4(4)
    character(len=*), intent(in) :: name, units, long_name
    integer, intent(out) :: varid
    call nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dims4, varid), 'def var '//name)
    call nc_check(nf90_put_att(ncid, varid, 'units', units), 'att units '//name)
    call nc_check(nf90_put_att(ncid, varid, 'long_name', long_name), 'att long_name '//name)
    call nc_check(nf90_put_att(ncid, varid, '_FillValue', NF90_FILL_DOUBLE), 'att _FillValue '//name)
  end subroutine def_profile_var

  subroutine def_scalar_var(ncid, name, dims3, units, long_name, varid)
    integer, intent(in) :: ncid, dims3(3)
    character(len=*), intent(in) :: name, units, long_name
    integer, intent(out) :: varid
    call nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dims3, varid), 'def var '//name)
    call nc_check(nf90_put_att(ncid, varid, 'units', units), 'att units '//name)
    call nc_check(nf90_put_att(ncid, varid, 'long_name', long_name), 'att long_name '//name)
    call nc_check(nf90_put_att(ncid, varid, '_FillValue', NF90_FILL_DOUBLE), 'att _FillValue '//name)
  end subroutine def_scalar_var


  !=================================================================
  ! output_write_timestep
  !
  ! Writes ONE column's ONE timestep: position (ilon_idx, ilat_idx) in
  ! the output grid (1-based, i.e. 1..f%nlon / 1..f%nlat -- NOT the raw
  ! ERA5 file index, which may start elsewhere; the caller subtracts
  ! ilon_start/ilat_start -- see rad_offline_driver.F90) and
  ! it_index (1-based position along valid_time, i.e. itime-itime_start
  ! +1 -- computed by the caller, NOT auto-incremented here, since with
  ! an area the same it_index is written once per grid column, not
  ! once per call in a simple 1:1 sequence).
  !
  ! n_valid is THIS column's number of above-ground ERA5 levels
  ! (era5_column_type%nlev_valid); the data arrays are sized n_valid
  ! (mass levels) / n_valid+1 (interfaces), NOT f%nlev_total -- only
  ! indices 1:n_valid (resp. 1:n_valid+1) of the output file's fixed
  ! pressure_level/half_level dimension are written for this column;
  ! the rest keep NetCDF's default fill value (see the module header).
  !=================================================================
  subroutine output_write_timestep(f, ilon_idx, ilat_idx, it_index, &
       &  n_valid, valid_time_value, t, q, z, clwc, &
       &  ciwc, crwc, cswc, cc, o3, sp, skt, rthratensw, rthratenlw, &
       &  sw_up, sw_dn, lw_up, lw_dn)

    type(output_file_type), intent(inout) :: f
    integer, intent(in) :: ilon_idx, ilat_idx, it_index, n_valid
    real(kind=RKIND), intent(in) :: valid_time_value
    real(kind=RKIND), intent(in) :: t(n_valid), q(n_valid), z(n_valid)
    real(kind=RKIND), intent(in) :: clwc(n_valid), ciwc(n_valid)
    real(kind=RKIND), intent(in) :: crwc(n_valid), cswc(n_valid)
    real(kind=RKIND), intent(in) :: cc(n_valid), o3(n_valid)
    real(kind=RKIND), intent(in) :: sp, skt
    real(kind=RKIND), intent(in) :: rthratensw(n_valid), rthratenlw(n_valid)
    real(kind=RKIND), intent(in) :: sw_up(n_valid+1), sw_dn(n_valid+1)
    real(kind=RKIND), intent(in) :: lw_up(n_valid+1), lw_dn(n_valid+1)

    integer :: start4(4), count4(4), start3(3), count3(3), start1(1)

    ! Harmless to (re)write the same valid_time value once per grid
    ! column at the same it_index -- same value every time, and this
    ! keeps the bookkeeping here simple (no "have we already written
    ! this time index" state to track).
    start1 = (/ it_index /)
    call nc_check(nf90_put_var(f%ncid, f%id_valid_time, (/valid_time_value/), start1), 'put time')

    start4 = (/ ilon_idx, ilat_idx, 1, it_index /)
    count4 = (/ 1, 1, n_valid, 1 /)
    call nc_check(nf90_put_var(f%ncid, f%id_t,    reshape(t,   (/1,1,n_valid,1/)), start4, count4), 'put t')
    call nc_check(nf90_put_var(f%ncid, f%id_q,    reshape(q,   (/1,1,n_valid,1/)), start4, count4), 'put q')
    call nc_check(nf90_put_var(f%ncid, f%id_z,    reshape(z,   (/1,1,n_valid,1/)), start4, count4), 'put z')
    call nc_check(nf90_put_var(f%ncid, f%id_clwc, reshape(clwc,(/1,1,n_valid,1/)), start4, count4), 'put clwc')
    call nc_check(nf90_put_var(f%ncid, f%id_ciwc, reshape(ciwc,(/1,1,n_valid,1/)), start4, count4), 'put ciwc')
    call nc_check(nf90_put_var(f%ncid, f%id_crwc, reshape(crwc,(/1,1,n_valid,1/)), start4, count4), 'put crwc')
    call nc_check(nf90_put_var(f%ncid, f%id_cswc, reshape(cswc,(/1,1,n_valid,1/)), start4, count4), 'put cswc')
    call nc_check(nf90_put_var(f%ncid, f%id_cc,   reshape(cc,  (/1,1,n_valid,1/)), start4, count4), 'put cc')
    call nc_check(nf90_put_var(f%ncid, f%id_o3,   reshape(o3,  (/1,1,n_valid,1/)), start4, count4), 'put o3')
    call nc_check(nf90_put_var(f%ncid, f%id_rthratensw, reshape(rthratensw,(/1,1,n_valid,1/)), start4, count4), 'put rthratensw')
    call nc_check(nf90_put_var(f%ncid, f%id_rthratenlw, reshape(rthratenlw,(/1,1,n_valid,1/)), start4, count4), 'put rthratenlw')

    start3 = (/ ilon_idx, ilat_idx, it_index /)
    count3 = (/ 1, 1, 1 /)
    call nc_check(nf90_put_var(f%ncid, f%id_sp,  (/sp/),  start3, count3), 'put sp')
    call nc_check(nf90_put_var(f%ncid, f%id_skt, (/skt/), start3, count3), 'put skt')

    start4 = (/ ilon_idx, ilat_idx, 1, it_index /)
    count4 = (/ 1, 1, n_valid+1, 1 /)
    call nc_check(nf90_put_var(f%ncid, f%id_sw_up, reshape(sw_up,(/1,1,n_valid+1,1/)), start4, count4), 'put sw_up')
    call nc_check(nf90_put_var(f%ncid, f%id_sw_dn, reshape(sw_dn,(/1,1,n_valid+1,1/)), start4, count4), 'put sw_dn')
    call nc_check(nf90_put_var(f%ncid, f%id_lw_up, reshape(lw_up,(/1,1,n_valid+1,1/)), start4, count4), 'put lw_up')
    call nc_check(nf90_put_var(f%ncid, f%id_lw_dn, reshape(lw_dn,(/1,1,n_valid+1,1/)), start4, count4), 'put lw_dn')

  end subroutine output_write_timestep


  subroutine output_close(f)
    type(output_file_type), intent(inout) :: f
    call nc_check(nf90_close(f%ncid), 'close')
  end subroutine output_close

end module output_netcdf
