! rrtmg_offline_driver.F90
!
! Standalone, offline driver for MPAS's RRTMG-SW/RRTMG-LW physics
! (module_ra_rrtmg_sw.F / module_ra_rrtmg_lw.F, extracted verbatim
! from MPAS-Dev/MPAS-Model -- see the header comments in those two
! files and in era5_bridge.F90 for exactly what was and was not
! changed), driven entirely by ERA5 hourly pressure-level data instead
! of an MPAS forecast.
!
! REQUIRED EXTERNAL DATA (not included, must be supplied by the user
! -- see the project README):
!   RRTMG_SW_DATA.DBL, RRTMG_LW_DATA.DBL
!     RRTMG absorption look-up tables, the same ones MPAS/WRF need.
!     Must be present in the working directory this program is run
!     from (rrtmg_initsw_forMPAS/rrtmg_initlw_forMPAS open them by
!     that exact relative name).
!   One ERA5 pressure-level file and one ERA5 single-level file
!   (NetCDF), covering the lon/lat/time range you want -- see the
!   variable list documented at the top of era5_bridge.F90.
!
! USAGE:
!   rrtmg_offline_driver <pressure_level.nc> <single_level.nc> \
!       <ilon> <ilat> <itime_start> <itime_end> <output.csv>
!   (1-based NetCDF indices; run once per grid column of interest --
!   this is a column-by-column research/diagnostic tool, not a
!   reimplementation of MPAS's 2-D tiled physics driver)
!
program rrtmg_offline_driver

  use mpas_kind_types,    only : RKIND
  use mpas_derived_types, only : dm_info
  use mpas_atmphys_constants, only : cp, gravity
  use radconst_offline,  only : radconst
  use era5_bridge

  use mpas_atmphys_rrtmg_swinit, only : rrtmg_initsw_forMPAS
  use mpas_atmphys_rrtmg_lwinit, only : rrtmg_initlw_forMPAS
  use module_ra_rrtmg_sw,        only : rrtmg_swrad
  use module_ra_rrtmg_lw,        only : rrtmg_lwrad

  implicit none

  ! ---- command-line arguments ----
  character(len=512) :: path_pl, path_sl, path_out
  character(len=32)  :: arg
  integer :: ilon, ilat, itime_start, itime_end, itime

  ! ---- ERA5 ----
  type(era5_file_type)   :: era5f
  type(era5_column_type) :: col

  ! ---- WRF/MPAS-style "1x1 tile" dimensions RRTMG expects ----
  integer, parameter :: ids=1,ide=1, jds=1,jde=1
  integer, parameter :: ims=1,ime=1, jms=1,jme=1
  integer, parameter :: its=1,ite=1, jts=1,jte=1
  integer :: kds,kde, kms,kme, kts,kte
  integer :: nlev   ! = kte-kts+1, the number of ATMOSPHERIC (above
                     ! ground) ERA5 levels for the current column;
                     ! recomputed every column/time since it can
                     ! differ (see era5_bridge.F90's nlev_valid note)

  ! ---- arrays passed to rrtmg_swrad/rrtmg_lwrad, dimensioned
  !      (ims:ime,kms:kme,jms:jme) = (1,nlev,1) or (1,nlev+1,1) or
  !      (1,nlev+2,1) as appropriate, allocated fresh each column
  !      since nlev varies. Kept as allocatables rather than automatic
  !      arrays sized by a fixed max level count so a below-ground-
  !      heavy column (e.g. over the Andes/Himalaya/Tibetan Plateau)
  !      never silently overflows a fixed bound. ----
  real(kind=RKIND), allocatable :: p3d(:,:,:), p8w(:,:,:), pi3d(:,:,:)
  real(kind=RKIND), allocatable :: t3d(:,:,:), t8w(:,:,:), dz8w(:,:,:)
  real(kind=RKIND), allocatable :: qv3d(:,:,:), qc3d(:,:,:), qr3d(:,:,:)
  real(kind=RKIND), allocatable :: qi3d(:,:,:), qs3d(:,:,:), qg3d(:,:,:)
  real(kind=RKIND), allocatable :: cldfra3d(:,:,:), o33d(:,:,:)
  real(kind=RKIND), allocatable :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)
  real(kind=RKIND), allocatable :: rthratensw(:,:,:), rthratenlw(:,:,:)
  ! (spectral flux diagnostics swupflx/swdnflx/lwupflx/lwdnflx removed --
  ! see note at the rrtmg_swrad/rrtmg_lwrad calls below)

  real(kind=RKIND) :: tsk(ims:ime,jms:jme), albedo(ims:ime,jms:jme)
  real(kind=RKIND) :: emiss(ims:ime,jms:jme)
  real(kind=RKIND) :: xland(ims:ime,jms:jme), xice(ims:ime,jms:jme)
  real(kind=RKIND) :: snow(ims:ime,jms:jme)
  real(kind=RKIND) :: coszr(ims:ime,jms:jme), gsw(ims:ime,jms:jme), swcf(ims:ime,jms:jme)
  real(kind=RKIND) :: glw(ims:ime,jms:jme), olr(ims:ime,jms:jme), lwcf(ims:ime,jms:jme)
  real(kind=RKIND) :: xlat(ims:ime,jms:jme), xlong(ims:ime,jms:jme)

  ! ---- astronomical / time ----
  real(kind=RKIND) :: declin, solcon, julian, gmt, radt, degrad_local
  integer :: julday

  ! ---- misc ----
  type(dm_info) :: dminfo
  integer :: ounit, k
  real(kind=RKIND), parameter :: p0 = 100000.0_RKIND  ! Pa, reference pressure for exner
  real(kind=RKIND), parameter :: Rd_over_cp = 287.0_RKIND/1004.5_RKIND
  real(kind=RKIND) :: z0, z1, z2, w1, w2, fzm, fzp

  ! =================================================================
  ! 0) Parse command line
  ! =================================================================
  if (command_argument_count() /= 7) then
    write(0,'(A)') 'usage: rrtmg_offline_driver <pressure_level.nc> <single_level.nc> '// &
         &          '<ilon> <ilat> <itime_start> <itime_end> <output.csv>'
    stop 1
  end if
  call get_command_argument(1, path_pl)
  call get_command_argument(2, path_sl)
  call get_command_argument(3, arg); read(arg,*) ilon
  call get_command_argument(4, arg); read(arg,*) ilat
  call get_command_argument(5, arg); read(arg,*) itime_start
  call get_command_argument(6, arg); read(arg,*) itime_end
  call get_command_argument(7, path_out)

  ! =================================================================
  ! 1) Initialize RRTMG (loads RRTMG_SW_DATA.DBL / RRTMG_LW_DATA.DBL
  !    from the working directory) -- ONCE, before the time loop.
  ! =================================================================
  call rrtmg_initsw_forMPAS(dminfo)
  call rrtmg_initlw_forMPAS(dminfo)

  call era5_open(trim(path_pl), trim(path_sl), era5f)

  open(newunit=ounit, file=trim(path_out), status='replace', action='write')
  write(ounit,'(A)') 'itime,julian,gmt,level_k,z_m,p_pa,rthratensw_K_s,rthratenlw_K_s'

  ! =================================================================
  ! 2) Time loop -- one column, many hours
  ! =================================================================
  do itime = itime_start, itime_end

    call era5_read_column(era5f, ilon, ilat, itime, col)
    nlev = col%nlev_valid
    if (nlev < 2) then
      write(0,'(A,I0)') 'skipping itime=', itime, ' (fewer than 2 above-ground ERA5 levels -- check ilon/ilat)'
      cycle
    end if

    kds=1; kde=nlev; kms=1; kme=nlev+1; kts=1; kte=nlev

    allocate(p3d(ims:ime,kms:kme,jms:jme), pi3d(ims:ime,kms:kme,jms:jme))
    allocate(t3d(ims:ime,kms:kme,jms:jme))
    allocate(qv3d(ims:ime,kms:kme,jms:jme), qc3d(ims:ime,kms:kme,jms:jme))
    allocate(qr3d(ims:ime,kms:kme,jms:jme), qi3d(ims:ime,kms:kme,jms:jme))
    allocate(qs3d(ims:ime,kms:kme,jms:jme), qg3d(ims:ime,kms:kme,jms:jme))
    allocate(cldfra3d(ims:ime,kms:kme,jms:jme), o33d(ims:ime,kms:kme,jms:jme))
    allocate(re_cloud(ims:ime,kms:kme,jms:jme), re_ice(ims:ime,kms:kme,jms:jme))
    allocate(re_snow(ims:ime,kms:kme,jms:jme))
    allocate(dz8w(ims:ime,kms:kme,jms:jme))
    allocate(p8w(ims:ime,kms:kme,jms:jme), t8w(ims:ime,kms:kme,jms:jme))
    allocate(rthratensw(ims:ime,kms:kme,jms:jme), rthratenlw(ims:ime,kms:kme,jms:jme))
    ! swupflx/swdnflx/lwupflx/lwdnflx intentionally NOT allocated/passed
    ! (optional diagnostics, unused here) -- see note at the calls below.

    ! -------------------------------------------------------------
    ! 3) Full-level (mass point) fields: ERA5 columns are stored
    !    top-to-bottom (index 1 = highest altitude); RRTMG/WRF-style
    !    arrays here are bottom-to-top (kts=1 = nearest the surface),
    !    matching MPAS convention -- see the note in era5_bridge.F90.
    !    ERA5 index for WRF level k (k=kts..kte) is (nlev-k+1).
    ! -------------------------------------------------------------
    do k = kts, kte
      p3d(1,k,1)      = col%p(nlev-k+1)
      t3d(1,k,1)      = col%t(nlev-k+1)
      qv3d(1,k,1)     = col%qv(nlev-k+1)
      qc3d(1,k,1)     = col%qc(nlev-k+1)
      qr3d(1,k,1)     = col%qr(nlev-k+1)
      qi3d(1,k,1)     = col%qi(nlev-k+1)
      qs3d(1,k,1)     = col%qs(nlev-k+1)
      qg3d(1,k,1)     = 0.0_RKIND   ! ERA5/IFS has no graupel category
      cldfra3d(1,k,1) = col%cldfra(nlev-k+1)
      o33d(1,k,1)     = col%o3(nlev-k+1)
      pi3d(1,k,1)     = (p3d(1,k,1)/p0)**Rd_over_cp
      ! re_cloud/re_ice/re_snow: ERA5 does not provide these; passing
      ! has_reqc=has_reqi=has_reqs=0 below tells RRTMG to diagnose its
      ! own default effective radii from water content instead of
      ! using these arrays, so they are allocated but never read --
      ! zero-filled only so they are not uninitialized memory.
      re_cloud(1,k,1) = 0.0_RKIND
      re_ice(1,k,1)   = 0.0_RKIND
      re_snow(1,k,1)  = 0.0_RKIND
    end do

    ! -------------------------------------------------------------
    ! 4) Half-level (interface) pressure/temperature/thickness.
    !    Interior interfaces: height-weighted interpolation between
    !    adjacent full levels (same fzm/fzp method used in MPAS's own
    !    mpas_atmphys_interface.F -- see the ecRad/MONAN coupling work
    !    in this project for the derivation). Surface and top
    !    interfaces: linear/log-pressure extrapolation, also the same
    !    method MPAS itself uses.
    ! -------------------------------------------------------------
    block
      real(kind=RKIND) :: zfull(kts:kte)  ! height of WRF-order full levels
      do k = kts, kte
        zfull(k) = col%z(nlev-k+1)
      end do

      ! surface interface (k=kts): linear extrapolation to zsfc using
      ! the two lowest full levels
      p8w(1,kts,1) = col%sp
      t8w(1,kts,1) = t3d(1,kts,1) + (t3d(1,kts,1)-t3d(1,kts+1,1)) &
           &       * (col%zsfc-zfull(kts)) / (zfull(kts)-zfull(kts+1))

      ! interior interfaces
      do k = kts+1, kte
        fzm = (zfull(k)-zfull(k-1)) / max(zfull(min(k+1,kte))-zfull(k-1), 1.0_RKIND)
        if (k < kte) then
          fzm = (zfull(k)-zfull(k-1)) / (zfull(k+1)-zfull(k-1))
          fzp = (zfull(k+1)-zfull(k)) / (zfull(k+1)-zfull(k-1))
          p8w(1,k,1) = fzm*p3d(1,k,1) + fzp*p3d(1,k-1,1)
          t8w(1,k,1) = fzm*t3d(1,k,1) + fzp*t3d(1,k-1,1)
        else
          ! top interface (k=kte, i.e. kme=kte+1 slot): log-pressure
          ! extrapolation beyond the topmost level, same style as
          ! MPAS's own top-of-model treatment
          z0 = zfull(kte) + (zfull(kte)-zfull(kte-1))
          z1 = zfull(kte); z2 = zfull(kte-1)
          w1 = (z0-z2)/(z1-z2); w2 = 1.0_RKIND-w1
          p8w(1,kme,1) = exp(w1*log(p3d(1,kte,1)) + w2*log(p3d(1,kte-1,1)))
          t8w(1,kme,1) = w1*t3d(1,kte,1) + w2*t3d(1,kte-1,1)
          p8w(1,k,1)   = 0.5_RKIND*(p3d(1,k,1)+p3d(1,k-1,1))
          t8w(1,k,1)   = 0.5_RKIND*(t3d(1,k,1)+t3d(1,k-1,1))
        end if
      end do

      do k = kts, kte
        dz8w(1,k,1) = (p8w(1,k,1)-p8w(1,k+1,1)) / (gravity*(p3d(1,k,1)/(287.0_RKIND*t3d(1,k,1))))
      end do
    end block

    ! -------------------------------------------------------------
    ! 5) Surface / astronomical inputs
    ! -------------------------------------------------------------
    tsk(1,1)    = col%skt
    albedo(1,1) = col%albedo
    emiss(1,1)  = 0.98_RKIND   ! typical broadband surface emissivity;
                               ! ERA5 does not provide this directly
    xland(1,1)  = 2.0_RKIND - min(max(col%lsm,0.0_RKIND),1.0_RKIND)  ! WRF: 1=land, 2=water
    xice(1,1)   = col%siconc
    snow(1,1)   = col%snow_kg_m2
    xlat(1,1)   = col%lat
    xlong(1,1)  = col%lon
    coszr(1,1)  = 0.0_RKIND   ! recomputed internally by rrtmg_swrad
    gsw = 0.0_RKIND; swcf = 0.0_RKIND

    call era5_valid_time_to_julian_gmt(era5f, itime, julian, gmt, julday)
    call radconst(declin, solcon, julian)
    degrad_local = 3.1415926_RKIND/180.0_RKIND
    radt = 0.0_RKIND     ! minutes. IMPORTANT: rrtmg_swrad computes
                         ! coszrs using xt24=mod(xtime+radt*0.5,1440),
                         ! i.e. it centers the zenith angle on
                         ! xtime+radt/2 (WRF/MPAS's forward-averaging
                         ! convention between infrequent radiation
                         ! calls). Since this offline driver already
                         ! supplies gmt for the EXACT hour of each ERA5
                         ! timestep and wants an exact instantaneous
                         ! angle (not a forward average), radt must be
                         ! 0.0 here -- a nonzero value (as in an
                         ! earlier revision of this file) silently
                         ! shifts every coszrs by radt/2 minutes.

    rthratensw = 0.0_RKIND
    rthratenlw = 0.0_RKIND

    ! -------------------------------------------------------------
    ! 6) Run RRTMG-SW and RRTMG-LW
    !
    ! swupflx/swdnflx/lwupflx/lwdnflx (and their *c clear-sky
    ! counterparts) are OPTIONAL diagnostic outputs, deliberately not
    ! requested below: RRTMG internally extends the column with
    ! synthetic extra layers above the model top (fixed ~4 hPa
    ! increments down to ~0 hPa -- see "nlayers = kte +
    ! max(nint(pw1d(kte+1)/deltap),1)" in module_ra_rrtmg_sw.F/_lw.F)
    ! purely to fill these per-level spectral flux arrays, whose
    ! required size therefore depends on each column's model-top
    ! pressure, not on kme+2 as their declared interface suggests.
    ! Since this driver never reads these arrays (only rthratensw/
    ! rthratenlw and gsw/swcf/glw/olr/lwcf), omitting them avoids that
    ! sizing problem entirely -- a segmentation fault/heap corruption
    ! from writing past a too-small array was observed here before
    ! this fix.
    ! -------------------------------------------------------------
    call rrtmg_swrad(p3d=p3d,p8w=p8w,pi3d=pi3d,t3d=t3d,t8w=t8w,dz8w=dz8w,        &
         &  qv3d=qv3d,qc3d=qc3d,qr3d=qr3d,qi3d=qi3d,qs3d=qs3d,qg3d=qg3d,          &
         &  cldfra3d=cldfra3d,o33d=o33d,tsk=tsk,albedo=albedo,                   &
         &  xland=xland,xice=xice,snow=snow,coszr=coszr,xtime=0.0_RKIND,gmt=gmt, &
         &  julday=julday,radt=radt,degrad=degrad_local,declin=declin,          &
         &  solcon=solcon,xlat=xlat,xlong=xlong,icloud=1,                        &
         &  cldovrlp=2,idcor=1,o3input=0,                                        &
         &  noznlevels=1,gsw=gsw,swcf=swcf,rthratensw=rthratensw,                &
         &  has_reqc=0,has_reqi=0,has_reqs=0,re_cloud=re_cloud,                  &
         &  re_ice=re_ice,re_snow=re_snow,                                       &
         &  aer_opt=0,                                                           &
         &  ids=ids,ide=ide,jds=jds,jde=jde,kds=kds,kde=kde,                     &
         &  ims=ims,ime=ime,jms=jms,jme=jme,kms=kms,kme=kme,                     &
         &  its=its,ite=ite,jts=jts,jte=jte,kts=kts,kte=kte)

    call rrtmg_lwrad(p3d=p3d,p8w=p8w,pi3d=pi3d,t3d=t3d,t8w=t8w,dz8w=dz8w,        &
         &  qv3d=qv3d,qc3d=qc3d,qr3d=qr3d,qi3d=qi3d,qs3d=qs3d,qg3d=qg3d,          &
         &  cldfra3d=cldfra3d,o33d=o33d,tsk=tsk,emiss=emiss,                     &
         &  xland=xland,xice=xice,snow=snow,xlat=xlat,julday=julday,icloud=1,    &
         &  cldovrlp=2,idcor=1,o3input=0,noznlevels=1,                           &
         &  glw=glw,olr=olr,lwcf=lwcf,rthratenlw=rthratenlw,                     &
         &  has_reqc=0,has_reqi=0,has_reqs=0,re_cloud=re_cloud,                  &
         &  re_ice=re_ice,re_snow=re_snow,                                       &
         &  ids=ids,ide=ide,jds=jds,jde=jde,kds=kds,kde=kde,                     &
         &  ims=ims,ime=ime,jms=jms,jme=jme,kms=kms,kme=kme,                     &
         &  its=its,ite=ite,jts=jts,jte=jte,kts=kts,kte=kte)

    ! -------------------------------------------------------------
    ! 7) Write output
    ! -------------------------------------------------------------
    do k = kts, kte
      write(ounit,'(I0,2(",",F0.4),",",I0,4(",",ES14.6))') &
           &  itime, julian, gmt, k, nint(col%z(nlev-k+1))*1.0_RKIND, &
           &  p3d(1,k,1), rthratensw(1,k,1), rthratenlw(1,k,1)
    end do

    deallocate(p3d,pi3d,t3d,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,cldfra3d,o33d)
    deallocate(re_cloud,re_ice,re_snow,dz8w,p8w,t8w,rthratensw,rthratenlw)
    deallocate(col%p,col%t,col%z,col%qv,col%qc,col%qi,col%qr,col%qs,col%cldfra,col%o3)

  end do

  close(ounit)
  call era5_close(era5f)

end program rrtmg_offline_driver
