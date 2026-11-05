! radconst_offline.F90
!
! Faithful, unmodified-math port of subroutine radconst() from
! MPAS-Model's src/core_atmosphere/physics/
! mpas_atmphys_driver_radiation_sw.F (itself borrowed from WRF).
! Computes the solar declination angle and the eccentricity-corrected
! solar constant for a given day of year. Pure math, no MPAS
! dependency at all in the original -- copied here verbatim (only
! wrapped in a module and given its two module-level constants,
! solcon_0/degrad/dpd, explicitly as arguments/parameters instead of
! "use mpas_atmphys_constants").
!
! rrtmg_swrad() itself then computes the actual cosine of solar
! zenith angle (coszrs) internally from declin + latitude + hour angle
! (see module_ra_rrtmg_sw.F around "coszrs = sin(xxlat)*sin(declin) +
! cos(xxlat)*cos(declin)*cos(hrang)") -- the offline driver does not
! need to (and should not) duplicate that calculation; it only needs
! to supply declin/solcon (from here) plus gmt/xtime/julday/xlat/xlong
! so rrtmg_swrad can compute coszrs itself, exactly as MPAS does.
module radconst_offline

  use mpas_kind_types, only : RKIND

  implicit none
  private
  public :: radconst

  real(kind=RKIND), parameter :: degrad   = 3.1415926_RKIND/180.0_RKIND
  real(kind=RKIND), parameter :: dpd      = 360.0_RKIND/365.0_RKIND
  real(kind=RKIND), parameter :: solcon_0 = 1370.0_RKIND  ! W/m2, MPAS's own value
                                                           ! (mpas_atmphys_constants.F
                                                           ! has an alternative
                                                           ! 1365.0 commented out)

contains

  subroutine radconst(declin, solcon, julian)

    real(kind=RKIND), intent(in)  :: julian   ! day of year, 1.0-366.0 (fractional
                                               ! within the day is fine, e.g.
                                               ! 45.5 for noon on day 45)
    real(kind=RKIND), intent(out) :: declin   ! solar declination, radians
    real(kind=RKIND), intent(out) :: solcon   ! solar constant, W/m2

    real(kind=RKIND) :: obecl, sinob, sxlong, arg, djul, rjul, eccfac

    declin = 0.0_RKIND
    solcon = 0.0_RKIND

    ! obliquity = 23.5 degrees
    obecl = 23.5_RKIND*degrad
    sinob = sin(obecl)

    ! longitude of the sun from vernal equinox
    if (julian >= 80.0_RKIND) then
      sxlong = dpd*(julian-80.0_RKIND)
    else
      sxlong = dpd*(julian+285.0_RKIND)
    end if
    sxlong = sxlong*degrad
    arg    = sinob*sin(sxlong)
    declin = asin(arg)

    ! solar constant eccentricity factor (Paltridge and Platt 1976)
    djul   = julian*360.0_RKIND/365.0_RKIND
    rjul   = djul*degrad
    eccfac = 1.000110_RKIND + 0.034221_RKIND*cos(rjul) + 0.001280_RKIND*sin(rjul) &
         &   + 0.000719_RKIND*cos(2*rjul) + 0.000077_RKIND*sin(2*rjul)
    solcon = solcon_0*eccfac

  end subroutine radconst

end module radconst_offline
