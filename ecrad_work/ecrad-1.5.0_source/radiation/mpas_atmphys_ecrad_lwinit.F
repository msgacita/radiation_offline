! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!PK #define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)

!=================================================================================================================
 module mpas_atmphys_ecrad_lwinit
 
!PK  use mpas_dmpar
!PK  use mpas_kind_types
!PK  use mpas_derived_types
!PK  use mpas_atmphys_constants
!PK  use mpas_atmphys_utilities
!PK  use mpas_log, only : mpas_log_write
!PK  use mpas_io_units

!wrf physics
!pk use module_ra_ecrad_lw

 implicit none
 private
 public:: ecrad_initlw_forMPAS


!Main initialization module for the ecrad long wave radiation code, based on the MPI decomposition
!used in MPAS. It replaces the initialization in module_ra_rrtgm_lw.F.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_ecrad_lwinit:
! -----------------------------------------
! ecrad_lwinit_forMPAS: call ecrad_lwlookuptable from subroutine init_radiation_lw.
! ecrad_lwlookuptable : read and broadcast all input data on MPAS nodes using MPI decomposition.
 

 contains

!=================================================================================================================
 subroutine ecrad_initlw_forMPAS(dminfo)
!=================================================================================================================

!input arguments:
 type(dm_info):: dminfo

!-----------------------------------------------------------------------------------------------------------------

!read in absorption coefficients and other data:
!pk call ecrad_lwlookuptable(dminfo)

!Perform g-point reduction and other initializations: specific heat of dry air (cp) used in
!flux to heating rate conversion factor.
!pk call ecrad_lw_ini(cp)

 end subroutine ecrad_initlw_forMPAS

!=================================================================================================================
 end module mpas_atmphys_ecrad_lwinit
!=================================================================================================================
