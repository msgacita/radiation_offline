program test_era5_bridge
  use mpas_kind_types, only : RKIND
  use era5_bridge
  implicit none

  type(era5_file_type)   :: f
  type(era5_column_type) :: col
  real(kind=RKIND) :: julian, gmt
  integer :: julday, itime, k

  call era5_open('../data/test_pressure_level.nc', '../data/test_single_level.nc', f)

  print '(A,I0,A,I0,A,I0,A,I0)', 'nlon=', f%nlon, ' nlat=', f%nlat, &
       &  ' nlev=', f%nlev, ' ntime=', f%ntime
  print '(A)', 'time_units="'//trim(f%time_units)//'"'

  do itime = 1, f%ntime
    call era5_read_column(f, 1, 1, itime, col)
    call era5_valid_time_to_julian_gmt(f, itime, julian, gmt, julday)
    print '(A,I0,A,F10.3,A,F6.3,A,I0,A,I0)', &
         &  'itime=', itime, ' julian=', julian, ' gmt=', gmt, &
         &  ' julday=', julday, ' nlev_valid=', col%nlev_valid
    print '(A,F10.1,A,F6.2)', '  sp=', col%sp, ' skt=', col%skt
    do k = 1, col%nlev
      print '(A,I0,A,F8.1,A,F7.2,A,F9.1,A,L1)', &
           &  '  k=', k, ' p=', col%p(k), ' t=', col%t(k), &
           &  ' z=', col%z(k), ' valid=', col%p(k) < col%sp
    end do
    deallocate(col%p,col%t,col%z,col%qv,col%qc,col%qi,col%qr,col%qs,col%cldfra,col%o3)
  end do

  call era5_close(f)

end program test_era5_bridge
