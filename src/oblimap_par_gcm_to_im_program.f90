! File name: oblimap_gcm_to_im_program.f90
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

PROGRAM oblimap_gcm_to_im_program
  USE oblimap_configuration_module, ONLY: C, PAR, initialize_config_variables, oblimap_licence
  USE oblimap_gcm_to_im_mapping_module, ONLY: oblimap_gcm_to_im_mapping
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

  ! Read the configuration file and initialization of the struckt C%:
  CALL initialize_config_variables()

  ! FIXME with 2D decomposition
  IF(PAR%nprocs > C%NX) THEN
   CALL MPI_Finalize(ierror)
   IF(PAR%rank_shared == 0) WRITE(UNIT=*, FMT='(/2A, I4, A, I4/)') C%OBLIMAP_ERROR, ' You are using ', PAR%nprocs, ' processors, which is too much for the current implementated domain decomposition. Lower that number to at least: ', C%NX
   STOP
  END IF

  ! 1D decomposition
  PAR%ny0 = 1
  PAR%ny1 = C%NY

  call decompose(C%NX, PAR%nshared_procs, PAR%rank_shared, PAR%nx0, PAR%nx1)

  PAR%nlat0 = 1
  PAR%nlat1 = C%NLAT

  call decompose(C%NLON, PAR%nshared_procs, PAR%rank_shared, PAR%nlon0, PAR%nlon1)

  write(*,*) PAR%nshared_procs, PAR%rank_shared, PAR%nlon0, PAR%nlon1

  IF(PAR%rank == 0) THEN
   ! Output: -
   CALL oblimap_licence('oblimap_gcm_to_im_program')

   WRITE(UNIT=*,FMT='(4(A, I4)/)') '  OBLIMAP-PAR(shared) runs with: number_of_processors  = ', PAR%nprocs , ', NX = ', C%NX, ', max_nr_of_lines_per_partition_block = ', PAR%nx1-PAR%Nx0, ', load unbalance = ', PAR%nprocs * (PAR%nx1-PAR%Nx0) - C%NX
  END IF

  ! Calling the oblimap_gcm_to_im_mapping :
  CALL oblimap_gcm_to_im_mapping()

  ! Output: ierror
  CALL MPI_Finalize(ierror)

END PROGRAM oblimap_gcm_to_im_program
