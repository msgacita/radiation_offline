 ! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!PK #define DM_BCAST_INTEGER(A) call mpas_dmpar_bcast_int(dminfo,A)
!PK #define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)
!PK #define DM_BCAST_REAL(A) call mpas_dmpar_bcast_real(dminfo,A)

!=================================================================================================================
!PK  module mpas_atmphys_ecrad_swinit
!PK  use mpas_dmpar
!PK  use mpas_kind_types
!PK  use mpas_derived_types
!PK  use mpas_atmphys_constants
!PK  use mpas_atmphys_utilities
!PK  use mpas_log, only : mpas_log_write
!PK  use mpas_io_units

!wrf physics
!PK use module_ra_ecrad_sw

 implicit none
 private
 public:: ecrad_initsw_forMPAS


!Main initialization module for the ecrad short wave radiation code, based on the MPI decomposition
!used in MPAS. It replaces the initialization in module_ra_rrtgm_sw.F.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_ecrad_swinit:
! -----------------------------------------
! ecrad_swinit_forMPAS: call ecrad_swlookuptable from subroutine init_radiation_sw.
! ecrad_swlookuptable : read and broadcast all input data on MPAS nodes using MPI decomposition.
! ecrad_swinit        : added initialization specific to ecrad_sw.


 contains



!=================================================================================================================
 subroutine ecrad_initsw_forMPAS(dminfo)
!=================================================================================================================
!input arguments:
 type(dm_info):: dminfo

!-----------------------------------------------------------------------------------------------------------------

!read in absorption coefficients and other data:
!PK   call ecrad_swlookuptable(dminfo)

!Perform g-point reduction and other initializations: specific heat of dry air (cp) used in
!flux to heating rate conversion factor.
!PK  call ecrad_sw_ini(cp)

 end subroutine ecrad_initsw_forMPAS


!=================================================================================================================
 end module mpas_atmphys_ecrad_swinit
!=================================================================================================================
