! mpas_stubs_offline.F90
!
! Minimal, serial (no-MPI) replacements for the handful of MPAS
! infrastructure modules that mpas_atmphys_rrtmg_swinit.F,
! mpas_atmphys_rrtmg_lwinit.F, module_ra_rrtmg_sw.F, and
! module_ra_rrtmg_lw.F "use". This is the COMPLETE list of MPAS
! dependencies those four files have (confirmed by grepping the
! original MPAS-Model source; everything else RRTMG needs -- all the
! rrsw_kgXX/rrlw_kgXX absorption-coefficient modules, mcica_subcol_gen_*,
! MersenneTwister, rrtmg_sw_rad, rrtmg_lw_rad, etc. -- is bundled
! inside those same two module_ra_rrtmg_*.F files and needs no stub at
! all).
!
! None of this changes any radiative-transfer physics: it only
! replaces (a) MPI broadcast calls with no-ops, since an offline
! single-process driver has nothing to broadcast to, and (b) MPAS's
! own logging/unit-management/error-handling helpers with trivial
! plain-Fortran equivalents.
!
! DELIBERATELY NOT reproducing MPAS's real mpas_dmpar/mpas_log/
! mpas_io_units (those are large, general-purpose modules); this file
! only implements the exact narrow subset actually called by the four
! files above.

module mpas_kind_types
  implicit none
  integer, parameter :: RKIND   = selected_real_kind(12)
  integer, parameter :: StrKIND = 512
end module mpas_kind_types


module mpas_derived_types
  implicit none
  ! Real MPAS's dm_info carries full MPI communicator/rank/nprocs
  ! bookkeeping; here it is reduced to just the one field the RRTMG
  ! init files actually read (dminfo%my_proc_id), always 0 (see
  ! mpas_dmpar%IO_NODE below).
  type :: dm_info
    integer :: my_proc_id = 0
    integer :: nprocs     = 1
  end type dm_info
end module mpas_derived_types


module mpas_dmpar
  use mpas_derived_types, only : dm_info
  implicit none
  integer, parameter :: IO_NODE = 0
contains
  ! No-op broadcasts: a serial offline driver has only one "rank", so
  ! whatever was just read by that rank already IS the broadcast
  ! result. These exist only so the DM_BCAST_INTEGER/MACRO/REAL
  ! preprocessor macros used throughout mpas_atmphys_rrtmg_swinit.F /
  ! _lwinit.F still expand to something callable.
  subroutine mpas_dmpar_bcast_int(dminfo, i)
    type(dm_info), intent(in) :: dminfo
    integer, intent(inout) :: i
  end subroutine mpas_dmpar_bcast_int

  subroutine mpas_dmpar_bcast_real(dminfo, r)
    use mpas_kind_types, only : RKIND
    type(dm_info), intent(in) :: dminfo
    real(kind=RKIND), intent(inout) :: r
  end subroutine mpas_dmpar_bcast_real

  subroutine mpas_dmpar_bcast_reals(dminfo, n, r)
    use mpas_kind_types, only : RKIND
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: n
    real(kind=RKIND), dimension(*), intent(inout) :: r
  end subroutine mpas_dmpar_bcast_reals
end module mpas_dmpar


module mpas_atmphys_constants
  ! Physical constants module_ra_rrtmg_sw.F/_lw.F actually import
  ! ("use mpas_atmphys_constants, only : cp, g=>gravity") -- values
  ! taken from MPAS-Model's own mpas_atmphys_constants.F, so results
  ! match exactly what MPAS itself would compute.
  use mpas_kind_types, only : RKIND
  implicit none
  real(kind=RKIND), parameter :: cp      = 1004.5_RKIND
  real(kind=RKIND), parameter :: gravity = 9.80616_RKIND
end module mpas_atmphys_constants


module mpas_atmphys_utilities
  implicit none
contains
  subroutine physics_error_fatal(mess)
    character(len=*), intent(in) :: mess
    write(0,'(A)') 'FATAL: '//trim(mess)
    stop 1
  end subroutine physics_error_fatal
end module mpas_atmphys_utilities


module mpas_log
  implicit none
contains
  ! Real mpas_log_write supports a small templating language ($i, $r,
  ! $c... substituted from intArgs/realArgs/charArgs). The RRTMG init
  ! files only ever pass zero or one integer, so that one case is
  ! handled directly; anything else just prints the raw message.
  subroutine mpas_log_write(mess, intArgs, realArgs, charArgs)
    character(len=*), intent(in) :: mess
    integer, dimension(:), intent(in), optional :: intArgs
    real,    dimension(:), intent(in), optional :: realArgs
    character(len=*), dimension(:), intent(in), optional :: charArgs

    if (present(intArgs) .and. len_trim(mess) > 0) then
      write(*,'(A,1x,I0)') trim(mess), intArgs(1)
    else if (len_trim(mess) > 0) then
      write(*,'(A)') trim(mess)
    end if
  end subroutine mpas_log_write
end module mpas_log


module mpas_io_units
  implicit none
contains
  subroutine mpas_new_unit(unit_num, unformatted)
    integer, intent(out) :: unit_num
    logical, intent(in), optional :: unformatted
    integer :: i
    logical :: is_open
    do i = 20, 999
      inquire(unit=i, opened=is_open)
      if (.not. is_open) then
        unit_num = i
        return
      end if
    end do
    unit_num = -1
  end subroutine mpas_new_unit

  subroutine mpas_release_unit(unit_num)
    integer, intent(in) :: unit_num
    ! Nothing to release explicitly; the unit is simply closed by the
    ! caller (mpas_atmphys_rrtmg_swinit.F/_lwinit.F already call
    ! close(rrtmg_unit) themselves before this).
  end subroutine mpas_release_unit
end module mpas_io_units
