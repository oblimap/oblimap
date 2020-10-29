! File name: oblimap_par_im_to_gcm_program.f90
!
! Copyright (C) 2020 Thomas Reerink.
!
! This file is distributed under the terms of the
! GNU General Public License.
!
! This file is part of OBLIMAP 2.2
!
! See Reerink et al. (2010,2016) for OBLIMAP's scientific documentation:
!  https://www.geosci-model-dev.net/3/13/2010/
!  https://www.geosci-model-dev.net/9/4111/2016/
!
! The OBLIMAP User Guide (Reerink, 2016) can be found at:
!  https://github.com/oblimap/oblimap/blob/master/documentation/
!
! The OBLIMAP code can be downloaded by:
!  svn checkout https://svn.science.uu.nl/repos/project.oblimap
! or from OBLIMAP's Github by:
!  git clone https://github.com/oblimap/oblimap
!
! OBLIMAP is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! OBLIMAP is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with OBLIMAP. If not, see https://www.gnu.org/licenses/
!
!
! OBLIMAP is maintained by:
!
! Thomas Reerink
! Institute for Marine and Atmospheric Research Utrecht (IMAU)
! Utrecht University
! Princetonplein 5
! 3584 CC Utrecht
! The Netherlands
!
! email: tjreerink@gmail.com
!

PROGRAM oblimap_par_im_to_gcm_program
  USE oblimap_configuration_module, ONLY: C, PAR, initialize_config_variables, oblimap_licence
  USE oblimap_im_to_gcm_mapping_module, ONLY: oblimap_im_to_gcm_mapping
  use mpi_f08
  use mpi_helpers_mod
  IMPLICIT NONE

  INTEGER :: ierror

  ! Output: ierror
  CALL MPI_Init(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, PAR%rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, PAR%nprocs, ierror)

  call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, PAR%rank_shared, MPI_INFO_NULL, PAR%shared_comm)
  call MPI_COMM_RANK(PAR%shared_comm, PAR%rank_shared, ierror)
  call MPI_COMM_SIZE(PAR%shared_comm, PAR%nshared_procs, ierror)

  call MPI_Comm_split(MPI_COMM_WORLD, PAR%rank_shared, 0, PAR%inter_comm)
  call MPI_Comm_rank(PAR%inter_comm, PAR%rank_inter)
  call MPI_Comm_size(PAR%inter_comm, PAR%n_nodes)

  ! Read the configuration file and initialization of the struckt C%:
  CALL initialize_config_variables()

  ! FIXME with 2D decomposition
  IF(PAR%nprocs > C%NX) THEN
   CALL MPI_Finalize(ierror)
   IF(PAR%rank == 0) WRITE(UNIT=*, FMT='(/2A, I4, A, I4/)') C%OBLIMAP_ERROR, ' You are using ', PAR%nprocs, ' processors, which is too much for the current implementated domain decomposition. Lower that number to at least: ', C%NLON
   STOP
  END IF

  ! 1D decomposition
#ifndef CONFIG_X_DECOMP
!! by default, decompose in X direction
#define CONFIG_X_DECOMP
#endif

#ifdef CONFIG_X_DECOMP
  call decompose(C%NX, PAR%n_nodes, PAR%rank_inter, PAR%node_nx0, PAR%node_nx1)
  call decompose(PAR%node_nx0, PAR%node_nx1, PAR%nshared_procs, PAR%rank_shared, PAR%nx0, PAR%nx1)

  call decompose(C%NLON, PAR%n_nodes, PAR%rank_inter, PAR%node_nlon0, PAR%node_nlon1)
  call decompose(PAR%node_nlon0, PAR%node_nlon1, PAR%nshared_procs, PAR%rank_shared, PAR%nlon0, PAR%nlon1)
#else
  PAR%nx0 = 1
  PAR%nx1 = C%NX

  PAR%nlon0 = 1
  PAR%nlon1 = C%NLON
#endif

#ifdef CONFIG_Y_DECOMP
  call decompose(C%NY, PAR%n_nodes, PAR%rank_inter, PAR%node_ny0, PAR%node_ny1)
  call decompose(PAR%node_ny0, PAR%node_ny1, PAR%nshared_procs, PAR%rank_shared, PAR%ny0, PAR%ny1)

  call decompose(C%NLAT, PAR%n_nodes, PAR%rank_inter, PAR%node_nlat0, PAR%node_nlat1)
  call decompose(PAR%node_nlat0, PAR%node_nlat1, PAR%nshared_procs, PAR%rank_shared, PAR%nlat0, PAR%nlat1)
#else
  PAR%ny0 = 1
  PAR%ny1 = C%NY

  PAR%nlat0 = 1
  PAR%nlat1 = C%NLAT
#endif

#if (defined CONFIG_X_DECOMP && defined CONFIG_Y_DECOMP)
#error "CONFIG_X_DECOMP and CONFIG_Y_DECOMP are not compatible (yet)"
#endif

  IF(PAR%rank == 0) THEN
   ! Output: -
   CALL oblimap_licence('oblimap_gcm_to_im_program')

   WRITE(UNIT=*,FMT='(1(A, I4)/)') '  OBLIMAP-PAR(shared) runs with: number_of_processors  = ', PAR%nprocs
   WRITE(UNIT=*,FMT='(2(A, I4)/)') ', NLON = ', C%NLON, ', max_nr_of_lines_per_partition_block = ', PAR%nlon1-PAR%nlon0, ', load unbalance = ', PAR%nprocs * (PAR%nlon1-PAR%nlon0) - C%nlon
   WRITE(UNIT=*,FMT='(2(A, I4)/)') ', NLAT = ', C%NLAT, ', max_nr_of_columns_per_partition_block = ', PAR%nlat1-PAR%nlat0, ', load unbalance = ', PAR%nprocs * (PAR%nlat1-PAR%nlat0) - C%NLAT
  END IF

  WRITE(UNIT=*,FMT='(3I4,3(A, I4))') PAR%rank, PAR%rank_inter, PAR%rank_shared, ', NLON = ', C%NLON, ', max_nr_of_lines_per_partition_block = ', PAR%nlon1-PAR%nlon0, ', load unbalance = ', PAR%nprocs * (PAR%nlon1-PAR%nlon0) - C%NLON
  WRITE(UNIT=*,FMT='(3I4,3(A, I4))') PAR%rank, PAR%rank_inter, PAR%rank_shared, ', NLAT = ', C%NLAT, ', max_nr_of_columns_per_partition_block = ', PAR%nlat1-PAR%nlat0, ', load unbalance = ', PAR%nprocs * (PAR%nlat1-PAR%nlat0) - C%NLAT

  ! Calling the oblimap_im_to_gcm_mapping :
  CALL oblimap_im_to_gcm_mapping()

  ! Output: ierror
  CALL MPI_Finalize(ierror)

END PROGRAM oblimap_par_im_to_gcm_program
