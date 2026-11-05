! namelist_offline.F90
!
! Reads the &offline_control namelist that selects which radiation
! scheme (RRTMG or ecRad) rad_offline_driver.F90 runs, and the file
! paths/grid range it needs. This is the offline-driver equivalent of
! MONAN's config_radt_sw_scheme/config_radt_lw_scheme namelist options
! (mpas_atmphys_driver_radiation_sw.F/_lw.F), reduced to a single
! "radiation_scheme" switch since this driver always runs SW and LW
! together for whichever scheme is chosen (see the header comments in
! era5_ecrad_bridge.F90 and rrtmg_offline_driver.F90 for why RRTMG and
! ecRad both compute SW+LW jointly).
!
! Example namelist file (see namelist.offline.sample in this project):
!
!   &offline_control
!     radiation_scheme     = 'ecrad'        ! 'rrtmg' or 'ecrad'
!     path_pressure_level  = 'era5_pressure_level.nc'
!     path_single_level    = 'era5_single_level.nc'
!     ilon_start           = 1              ! 1-based indices into the
!     ilon_end             = 1              ! ERA5 file's own
!     ilat_start           = 1              ! longitude/latitude
!     ilat_end             = 1              ! coordinate arrays
!     itime_start          = 1
!     itime_end            = 24
!     path_output          = 'output.nc'
!     ecrad_namelist_path  = 'config.nam'   ! only read if
!                                           ! radiation_scheme='ecrad'
!   /
!
! A single column (the only mode this driver had before area support
! was added) is just ilon_start=ilon_end and ilat_start=ilat_end. For
! backward compatibility with older namelist files, the singular
! "ilon"/"ilat" entries are still accepted and, if present, override
! BOTH the start and end of the corresponding range (i.e. "ilon=3" is
! shorthand for "ilon_start=3, ilon_end=3").
!
module namelist_offline

  implicit none
  private

  public :: offline_control_type
  public :: read_offline_namelist

  type :: offline_control_type
    character(len=16)  :: radiation_scheme    = 'rrtmg'
    character(len=512) :: path_pressure_level = ''
    character(len=512) :: path_single_level   = ''
    integer             :: ilon_start         = 1
    integer             :: ilon_end           = 1
    integer             :: ilat_start         = 1
    integer             :: ilat_end           = 1
    integer             :: itime_start        = 1
    integer             :: itime_end          = 1
    character(len=512) :: path_output         = 'output.nc'
    character(len=512) :: ecrad_namelist_path = 'config.nam'
  end type offline_control_type

contains

  subroutine read_offline_namelist(namelist_file, control)

    character(len=*), intent(in) :: namelist_file
    type(offline_control_type), intent(out) :: control

    character(len=16)  :: radiation_scheme
    character(len=512) :: path_pressure_level, path_single_level, path_output, ecrad_namelist_path
    integer :: ilon_start, ilon_end, ilat_start, ilat_end, itime_start, itime_end
    ! Deprecated singular shorthand -- 0 means "not set in this file"
    ! (a real 1-based index is always >=1, so 0 is a safe sentinel).
    integer :: ilon, ilat
    integer :: iunit, ios

    namelist /offline_control/ radiation_scheme, path_pressure_level, path_single_level, &
         &  ilon_start, ilon_end, ilat_start, ilat_end, ilon, ilat, &
         &  itime_start, itime_end, path_output, ecrad_namelist_path

    ! Defaults (in case the namelist file omits an optional entry)
    radiation_scheme     = control%radiation_scheme
    path_pressure_level  = control%path_pressure_level
    path_single_level    = control%path_single_level
    ilon_start            = control%ilon_start
    ilon_end              = control%ilon_end
    ilat_start            = control%ilat_start
    ilat_end              = control%ilat_end
    ilon = 0; ilat = 0
    itime_start          = control%itime_start
    itime_end            = control%itime_end
    path_output          = control%path_output
    ecrad_namelist_path  = control%ecrad_namelist_path

    open(newunit=iunit, file=trim(namelist_file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(0,'(A)') 'read_offline_namelist: could not open "'//trim(namelist_file)//'"'
      stop 1
    end if
    read(iunit, nml=offline_control, iostat=ios)
    if (ios /= 0) then
      write(0,'(A)') 'read_offline_namelist: error reading &offline_control from "'// &
           &          trim(namelist_file)//'" -- check for typos/missing commas'
      stop 1
    end if
    close(iunit)

    ! Deprecated singular shorthand: "ilon"/"ilat" override both ends
    ! of the corresponding range if present in the namelist file.
    if (ilon /= 0) then
      ilon_start = ilon; ilon_end = ilon
    end if
    if (ilat /= 0) then
      ilat_start = ilat; ilat_end = ilat
    end if

    select case (trim(radiation_scheme))
    case ('rrtmg', 'ecrad')
      ! ok
    case default
      write(0,'(A)') 'read_offline_namelist: radiation_scheme must be "rrtmg" or "ecrad", got "'// &
           &          trim(radiation_scheme)//'"'
      stop 1
    end select

    if (ilon_end < ilon_start .or. ilat_end < ilat_start) then
      write(0,'(A)') 'read_offline_namelist: ilon_end/ilat_end must be >= ilon_start/ilat_start'
      stop 1
    end if

    control%radiation_scheme    = radiation_scheme
    control%path_pressure_level = path_pressure_level
    control%path_single_level   = path_single_level
    control%ilon_start          = ilon_start
    control%ilon_end            = ilon_end
    control%ilat_start          = ilat_start
    control%ilat_end            = ilat_end
    control%itime_start         = itime_start
    control%itime_end           = itime_end
    control%path_output         = path_output
    control%ecrad_namelist_path = ecrad_namelist_path

  end subroutine read_offline_namelist

end module namelist_offline
