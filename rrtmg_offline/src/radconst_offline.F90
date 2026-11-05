! radconst_offline.F90
!
! Computes solar declination, the equation-of-time-corrected cosine of
! the solar zenith angle, and the eccentricity-corrected solar constant.
module radconst_offline

  use mpas_kind_types, only : RKIND

  implicit none
  private
  public :: radconst

  real(kind=RKIND), parameter :: degrad   = 3.1415926_RKIND/180.0_RKIND
  real(kind=RKIND), parameter :: solcon_0 = 1366.1_RKIND  ! W/m2, MPAS's own value
                                                           ! (mpas_atmphys_constants.F
                                                           ! has an alternative
                                                           ! 1365.0 commented out)

!   real(kind=RKIND), parameter :: solcon_0 = 1370.0_RKIND  ! W/m2, MPAS's own value
!                                                            ! (mpas_atmphys_constants.F
!                                                            ! has an alternative
!                                                            ! 1365.0 commented out)
contains

  subroutine radconst(declin, solcon, coszrs, julian, xlat, xlong, gmt, xtime, radt)

    real(kind=RKIND), intent(in)  :: julian   ! day of year, 1.0-366.0 (fractional
                                               ! within the day is fine, e.g.
                                               ! 45.5 for noon on day 45)
    real(kind=RKIND), intent(in)  :: xlat     ! latitude, degrees
    real(kind=RKIND), intent(in)  :: xlong    ! longitude, degrees
    real(kind=RKIND), intent(in)  :: gmt      ! GMT offset, hours
    real(kind=RKIND), intent(in)  :: xtime    ! local time, minutes
    real(kind=RKIND), intent(in)  :: radt     ! radiation time step, minutes
    real(kind=RKIND), intent(out) :: declin   ! solar declination, radians
    real(kind=RKIND), intent(out) :: solcon   ! solar constant, W/m2
    real(kind=RKIND), intent(out) :: coszrs   ! cosine of solar zenith angle

    real(kind=RKIND) :: djul, rjul, eccfac
    real(kind=RKIND) :: hour, gamma, eot, angle

    declin = 0.0_RKIND
    solcon = 0.0_RKIND
    coszrs = 0.0_RKIND

    ! julian already includes the UTC fraction of the day. xtime and
    ! half the radiation step advance that instant to the RRTMG midpoint.
    hour = gmt + (xtime+radt*0.5)/60.0
    gamma = 2.0_RKIND*acos(-1.0_RKIND)/365.0_RKIND * &
          &  (julian-1.0_RKIND+(xtime+radt*0.5)/1440.0_RKIND)
    ! NOAA fractional-year solar declination and equation of time (minutes):
    ! https://gml.noaa.gov/grad/solcalc/solareqns.PDF
    declin = 0.006918_RKIND - 0.399912_RKIND*cos(gamma) + 0.070257_RKIND*sin(gamma) &
          & - 0.006758_RKIND*cos(2.0_RKIND*gamma) + 0.000907_RKIND*sin(2.0_RKIND*gamma) &
          & - 0.002697_RKIND*cos(3.0_RKIND*gamma) + 0.001480_RKIND*sin(3.0_RKIND*gamma)
    eot = 229.18_RKIND * (0.000075_RKIND + 0.001868_RKIND*cos(gamma) - 0.032077_RKIND*sin(gamma) &
          &  - 0.014615_RKIND*cos(2.0_RKIND*gamma) - 0.040849_RKIND*sin(2.0_RKIND*gamma))
    angle = 15.0_RKIND*(hour+xlong/15.0_RKIND+eot/60.0_RKIND-12.0_RKIND)*degrad
    coszrs = sin(xlat*degrad)*sin(declin) &
          & + cos(xlat*degrad)*cos(declin)*cos(angle)
    ! solar constant eccentricity factor (Paltridge and Platt 1976)
    djul   = julian*360.0_RKIND/365.0_RKIND
    rjul   = djul*degrad
    eccfac = 1.000110_RKIND + 0.034221_RKIND*cos(rjul) + 0.001280_RKIND*sin(rjul) &
         &   + 0.000719_RKIND*cos(2*rjul) + 0.000077_RKIND*sin(2*rjul)
    solcon = solcon_0*eccfac

  end subroutine radconst

end module radconst_offline
