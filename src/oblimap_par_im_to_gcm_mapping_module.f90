! File name: oblimap_im_to_gcm_mapping_module.f90
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

MODULE oblimap_im_to_gcm_mapping_module

CONTAINS
  SUBROUTINE oblimap_im_to_gcm_mapping()
    ! This program reads several field variables from a rectangular Ice Model (IM) netcdf file which
    ! are projected by an inverse oblique stereographic projection on a Global Circulation Model (GCM)
    ! grid. The IM grid point coordinates are projected on the GCM grid, where in general they fall
    ! irregular and between the GCM grid points. Therefore the field values defined on the IM grid
    ! points have to be interpolated at each GCM grid point, with help of the nearby projected IM
    ! grid points. Two interpolation methods are available in this program both based on a distance
    ! weigthing Shepard technique. One method, the 'quadrant method', searches within each quadrant
    ! around each GCM grid point the nearest projected IM point and interpolates with help of this
    ! four points to obtain the field value at such a GCM grid point. Another method, the 'radius
    ! method', searches all projected IM points within a certain radius and interpolates with help
    ! of these points to obtain the field value at such a GCM grid point. Finally the projected and
    ! interpolated GCM fields are written to a GCM netcdf file. The GCM points which were not affected
    ! by the mapping (= projection + interpolation) remain their values as before.
    !
    ! The GCM geographical coordinates are in longitude (lon) and latitude (lat) in degrees, while the
    ! IM rectangular coordinates are in x and y in meters. The (lon,lat) coordinates are defined in the
    ! curved spherical surface S, and the (x,y) coordinates in the flat surface S'. For a more extended
    ! description of the projection and the interpolation method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    !
    ! In the NAMELIST or config file, several options can be changed without compiling the program.
    ! For example the center of the area of interest (the center of the projected area) can be specified
    ! by setting the:
    !   lamda_M_config = 319  (value for the case of Greenland)
    !   phi_M_config   =  72  (value for the case of Greenland)
    ! in the config file. The exact (inverse) oblique stereographic projection can be chosen with
    !   alpha_stereographic_config = 7.1 (value for the case of Greenland)
    ! The center of the IM grid will coincide with (lamda_M_config,phi_M_config), and the extensions of
    ! the IM grid are determined by the IM grid spacings C%dx and C%dy and the IM grid sizes C%NX and C%NY.
    !
    !
    ! For the mapping from IM to GCM we make use of initial_gcm_field_*, mapped_gcm_field_*, and gcm_field_*.
    !
    ! initial_gcm_field_*
    !   Contains the initial GCM fields. This fields are used because a large part of each GCM field is
    !   unaffected by the to and fro mapping because the IM grid is only local.
    !   The unaffected (mapping_participation_mask(i,j) == 0) field values of the initial fields are
    !   afterwards merged with the to and fro mapped values. Besides in the 'scan' phase the longitude
    !   and latitude coordinates of the GCM grid points have to be read from the GCM file which contains
    !   these initial_gcm_field_* fields.
    !
    ! mapped_gcm_field_*
    !   Contains the GCM fields which are mapped (= projected + interpolated) from the IM grid on to the GCM
    !   grid, these fields only have non-zero values on the grid points (i,j) with:
    !    mapping_participation_mask(i,j) == 1.
    !
    ! gcm_field_*
    !   Contains the to and fro mapped GCM fields merged with the initial_gcm_field_*.
    !
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type, check_directory_existence, PAR
    USE oblimap_read_and_write_module, ONLY: oblimap_netcdf_file_type, oblimap_open_netcdf_file, initialize_im_coordinates, create_netcdf_for_gcm_grid, &
          oblimap_read_netcdf_fields, oblimap_write_netcdf_fields, oblimap_close_netcdf_file, reduce_dummy_dimensions
    USE oblimap_scan_contributions_module, ONLY: make_mapping_participation_mask, scan_with_radius_method_im_to_gcm, scan_with_quadrant_method_im_to_gcm, &
          shifting_center_im_grid, determining_scan_parameters
    USE oblimap_mapping_module, ONLY: oblimap_ddo_type, oblimap_read_sid_file, oblimap_mapping, oblimap_deallocate_ddo
    use mpi_helpers_mod
    use mpi_f08
    IMPLICIT NONE

    ! Local variables:
    !REAL(dp), DIMENSION(                          C%NLON,C%NLAT                            ) :: lon_gcm
    !REAL(dp), DIMENSION(                          C%NLON,C%NLAT                            ) :: lat_gcm
    real(DP), dimension(:,:), pointer :: lon_gcm, lon_gcm_ ! global shared pointer and local one
    real(DP), dimension(:,:), pointer :: lat_gcm, lat_gcm_

    !REAL(dp), DIMENSION(                          C%NX  ,C%NY                              ) :: x_coordinates_of_im_grid_points  ! The x-coordinates of the IM points in S'
    !REAL(dp), DIMENSION(                          C%NX  ,C%NY                              ) :: y_coordinates_of_im_grid_points  ! The y-coordinates of the IM points in S'
    real(DP), dimension(:,:), pointer :: x_coordinates_of_im_grid_points, x_coordinates_of_im_grid_points_ ! global shared pointer and local one
    real(DP), dimension(:,:), pointer :: y_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points_
    integer, dimension(:,:), pointer :: mapping_participation_mask, mapping_participation_mask_

    LOGICAL,  DIMENSION(:,:,:,:), pointer :: mask_of_invalid_contributions, mask_of_invalid_contributions_    ! For each field and for each layer a mask represents the invalid contributions (like e.g. missing values) of the GCM grid points
    REAL(dp), DIMENSION(:,:,:,:), pointer :: gcm_field, gcm_field_
    REAL(dp), DIMENSION(:,:,:,:), pointer :: im_field, im_field_
    !REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers) :: initial_gcm_field
    !REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers) :: mapped_gcm_field
    REAL(dp), DIMENSION(:,:,:,:), pointer :: initial_gcm_field, initial_gcm_field_
    REAL(dp), DIMENSION(:,:,:,:), pointer :: mapped_gcm_field, mapped_gcm_field_
    REAL(dp)                                                                                 :: time                             ! The time value for the considered record
    INTEGER                                                                                  :: field_counter                    ! The counter in the loop over the field numbers
    INTEGER                                                                                  :: record_counter                   ! The counter over the time records
    INTEGER                                                                                  :: layer_counter                    ! The counter over the vertical layers
    TYPE(oblimap_netcdf_file_type)                                                           :: im_netcdf_file
    TYPE(oblimap_netcdf_file_type)                                                           :: initial_gcm_netcdf_file
    TYPE(oblimap_netcdf_file_type)                                                           :: created_gcm_netcdf_file
    TYPE(oblimap_ddo_type)                                                                   :: oblimap_ddo                      ! The DDO containing all the scanned contributions
    TYPE(oblimap_scan_parameter_type)                                                        :: advised_scan_parameter
    type(MPI_Win) :: lon_gcm_win, lat_gcm_win
    type(MPI_Win) :: x_coordinates_of_im_grid_points_win
    type(MPI_Win) :: y_coordinates_of_im_grid_points_win
    type(MPI_Win) :: mapping_participation_mask_win
    type(MPI_Win) :: mask_of_invalid_contributions_win
    type(MPI_Win) :: gcm_field_win
    type(MPI_Win) :: im_field_win
    type(MPI_Win) :: mapped_gcm_field_win
    type(MPI_Win) :: initial_gcm_field_win
    ! FIXME for testing
    integer :: p, row, col

    ! TODO allocate lon_gcm
    ! TODO allocate lat_gcm
    call alloc_shared( C%NLON, PAR%io_in_nlon0, PAR%io_in_nlon1 &
                     , C%NLAT, PAR%io_in_nlat0, PAR%io_in_nlat1 &
                     , lon_gcm, lon_gcm_ &
                     , lon_gcm_win &
                     , PAR%shared_comm)

    call alloc_shared( C%NLON, PAR%io_in_nlon0, PAR%io_in_nlon1 &
                     , C%NLAT, PAR%io_in_nlat0, PAR%io_in_nlat1 &
                     , lat_gcm, lat_gcm_ &
                     , lat_gcm_win &
                     , PAR%shared_comm)

    call alloc_shared( C%NX, PAR%io_in_nx0, PAR%io_in_nx1 &
                     , C%NY, PAR%io_in_ny0, PAR%io_in_ny1 &
                     , x_coordinates_of_im_grid_points &
                     , x_coordinates_of_im_grid_points_ &
                     , x_coordinates_of_im_grid_points_win, PAR%shared_comm)

    call alloc_shared( C%NX, PAR%io_in_nx0, PAR%io_in_nx1 &
                     , C%NY, PAR%io_in_ny0, PAR%io_in_ny1 &
                     , y_coordinates_of_im_grid_points &
                     , y_coordinates_of_im_grid_points_ &
                     , y_coordinates_of_im_grid_points_win, PAR%shared_comm)

    call alloc_shared( C%NLON, PAR%io_in_nlon0, PAR%io_in_nlon1 &
                     , C%NLAT, PAR%io_in_nlat0, PAR%io_in_nlat1 &
                     , mapping_participation_mask, mapping_participation_mask_ &
                     , mapping_participation_mask_win &
                     , PAR%shared_comm)

    ! Check whether the directory in the path of C%gcm_created_filename exists:
    CALL check_directory_existence(C%gcm_created_filename)

    ! Opening the initial GCM netcdf file, and reading the longitude and latitude coordinates of the initial GCM grid:
    ! Output: initial_gcm_netcdf_file, lon_gcm, lat_gcm
    CALL oblimap_open_netcdf_file(C%gcm_input_filename, C%number_of_mapped_fields, C%gcm_field_name, C%NLON, C%NLAT, lon_gcm, lat_gcm, .TRUE., nc = initial_gcm_netcdf_file)

    IF(C%choice_projection_method == 'rotation_projection') THEN
     lon_gcm_ = lon_gcm_ - C%shift_x_coordinate_rotation_projection
     lat_gcm_ = lat_gcm_ - C%shift_y_coordinate_rotation_projection
    ELSE
     ! It is required that all angles are in the 0 - 360 degree range:
     WHERE(lon_gcm_ <    0._dp) lon_gcm_ = lon_gcm_ + 360._dp
     WHERE(lon_gcm_ >= 360._dp) lon_gcm_ = lon_gcm_ - 360._dp
    END IF

    ! Output: im_netcdf_file
    CALL oblimap_open_netcdf_file(C%im_input_filename, C%number_of_mapped_fields, C%im_field_name, C%NX, C%NY, nc = im_netcdf_file)

    ! Determine the x, y coordinates of the IM grid points in plane S'
    ! Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    CALL initialize_im_coordinates(C%NX, C%NY, C%dx, C%dy, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)

    ! In/Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    IF(C%enable_shift_im_grid) CALL shifting_center_im_grid(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)

    ! This routine has to be called prior to make_mapping_participation_mask, for the case that alpha is determined automatically.
    ! Output: advised_scan_parameter
    IF(C%scanning_mode) CALL determining_scan_parameters('im-to-gcm', lon_gcm, lat_gcm, advised_scan_parameter)

    ! The IM grid is defined on just a small part of the world (e.g. Antarctica, Greenland etc.).
    ! While mapping the fields at the IM plane to the GCM grid, only those GCM points participate in
    ! the mapping which lay within the area covered by the IM grid. Therefore at those GCM grid points
    ! which did not participate in the mapping the field values of the original pre-mapped GCM file
    ! are taken (if not ignored by the C%ignore_reading_pre_mapped_fields option):
    !  gcm_field(mapping_participation_mask == 0) = gcm_field_pre_mapped(mapping_participation_mask == 0)
    ! So a mask is created for the GCM points which are involved (mask == 1) in the mapping:
    ! Output: mapping_participation_mask
    CALL make_mapping_participation_mask(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, mapping_participation_mask)


    ! Ahead of the mapping, the scanning has to be done. In case that the file C%sid_filename exists,
    ! the scanning phase can be omitted because this file contains the scanned projection data. This file contains for each
    ! target grid point the coordinates and the distances between the contributing points and the target point. These
    ! coordinates and distances are used for the field interpolation during the mapping phase.
    ! The points which are projected most nearby the target point contribute relative to the sheppard weighing depending on
    ! their distance. There are two methods: the 'radius method' and the 'quadrant method' which can be used to select the
    ! contributing points. This 'scan' part in fact contains the most technical and CPU consuming part of OBLIMAP, covering
    ! the projection and the selection of points for the interpolation (and keeping the relative distances of each projected
    ! point relative to the target point).
    IF(C%scanning_mode) THEN
     IF(C%choice_quadrant_method) THEN
      ! Output: the C%sid_filename is file created
      CALL scan_with_quadrant_method_im_to_gcm(mapping_participation_mask, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
     ELSE
      ! Output: the C%sid_filename file is created
      CALL scan_with_radius_method_im_to_gcm(mapping_participation_mask, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
     END IF
    END IF


    ! Reading the contributions of the scanned projection data into the Dynamic Data Object (DDO):
    ! Output: oblimap_ddo
    CALL oblimap_read_sid_file(C%sid_filename, oblimap_ddo)

    ! TODO find a way to load the mapped points according to the parallel decomposition
    ! this assume that the mapped points are sorted by paralllel proc (true by concatenate the sid file)
    oblimap_ddo%total_mapped_points0=oblimap_ddo%total_mapped_points
    oblimap_ddo%total_mapped_points1=1
    do p = 1, oblimap_ddo%total_mapped_points
      row = oblimap_ddo%row_index_destination(p)
      col = oblimap_ddo%column_index_destination(p)

      if ( row >= PAR%nlon0 .and. row <= PAR%nlon1 .and. &
           col >= PAR%nlat0 .and. col <= PAR%nlat1) then
           oblimap_ddo%total_mapped_points0 = min(oblimap_ddo%total_mapped_points0, p)
           oblimap_ddo%total_mapped_points1 = max(oblimap_ddo%total_mapped_points1, p)
      endif
    end do
    !!! FIXME read only indices within nx0:nx1 and ny0:ny1
    !!do p = oblimap_ddo%total_mapped_points0, oblimap_ddo%total_mapped_points1
    !!  row = oblimap_ddo%row_index_destination(p)
    !!  col = oblimap_ddo%column_index_destination(p)

    !!  if ( row < PAR%node_nlon0 .or. row > PAR%node_nlon1 .or. &
    !!       col < PAR%node_nlat0 .or. col > PAR%node_nlat1) then
    !!       write(*,'(a,i4,a,i8,a,i8,a)') 'WARNING: on ', PAR%rank, ' point (', row,' , ', col, ') will not be written well if multi-node run'
    !!  endif
    !!end do

    ! Output: -
    CALL create_netcdf_for_gcm_grid(lon_gcm, lat_gcm, im_netcdf_file, created_gcm_netcdf_file)

    call alloc_shared( C%NX, PAR%io_in_nx0, PAR%io_in_nx1 &
                     , C%NY, PAR%io_in_ny0, PAR%io_in_ny1 &
                     , C%number_of_vertical_layers, 1, C%number_of_vertical_layers &
                     , C%number_of_mapped_fields, 1, C%number_of_mapped_fields &
                     , mask_of_invalid_contributions, mask_of_invalid_contributions_ &
                     , mask_of_invalid_contributions_win &
                     , PAR%shared_comm)
    call alloc_shared( C%NLON, PAR%nlon0, PAR%nlon1 &
                     , C%NLAT, PAR%nlat0, PAR%nlat1 &
                     , C%number_of_vertical_layers, 1, C%number_of_vertical_layers &
                     , C%number_of_mapped_fields, 1, C%number_of_mapped_fields &
                     , gcm_field, gcm_field_ &
                     , gcm_field_win &
                     , PAR%shared_comm)
    call alloc_shared( C%NLON, PAR%io_in_nlon0, PAR%io_in_nlon1 &
                     , C%NLAT, PAR%io_in_nlat0, PAR%io_in_nlat1 &
                     , C%number_of_vertical_layers, 1, C%number_of_vertical_layers &
                     , C%number_of_mapped_fields, 1, C%number_of_mapped_fields &
                     , initial_gcm_field, initial_gcm_field_ &
                     , initial_gcm_field_win &
                     , PAR%shared_comm)
    call alloc_shared( C%NLON, PAR%nlon0, PAR%nlon1 &
                     , C%NLAT, PAR%nlat0, PAR%nlat1 &
                     , C%number_of_vertical_layers, 1, C%number_of_vertical_layers &
                     , C%number_of_mapped_fields, 1, C%number_of_mapped_fields &
                     , mapped_gcm_field, mapped_gcm_field_ &
                     , mapped_gcm_field_win &
                     , PAR%shared_comm)
    call alloc_shared( C%NX, PAR%io_in_nx0, PAR%io_in_nx1 &
                     , C%NY, PAR%io_in_ny0, PAR%io_in_ny1 &
                     , C%number_of_vertical_layers, 1, C%number_of_vertical_layers &
                     , C%number_of_mapped_fields, 1, C%number_of_mapped_fields &
                     , im_field, im_field_ &
                     , im_field_win &
                     , PAR%shared_comm)

    ! From the IM netcdf file we read the IM variables which will be mapped:
    DO record_counter = 0, C%im_record_range(2) - C%im_record_range(1)

     ! Output: initial_gcm_field
     CALL oblimap_read_netcdf_fields(initial_gcm_netcdf_file, C%gcm_record_range(1) + record_counter, initial_gcm_field)

     ! Output: im_field, time
     CALL oblimap_read_netcdf_fields(im_netcdf_file, C%im_record_range(1) + record_counter, im_field, time)
     CALL MPI_Barrier(PAR%shared_comm)

     ! Determine the mask_of_invalid_contributions based on the invalid values for the specified field:
     mask_of_invalid_contributions_ = .FALSE.
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
       IF(C%masked_fields(field_counter)) WHERE(im_field_(:,:,layer_counter,C%field_which_determines_invalid_value_mask(field_counter)) == C%invalid_input_value(field_counter)) &
        mask_of_invalid_contributions_(:,:,layer_counter,field_counter) = .TRUE.
     END DO
     END DO

     ! Rescaling each field by multiplication with a im_to_gcm_factor and by adding a im_to_gcm_shift (in case the units differ):
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
       WHERE(im_field_(:,:,layer_counter,field_counter) /= C%invalid_input_value(field_counter)) &
               im_field_(:,:,layer_counter,field_counter) = C%field_factor(field_counter) * im_field_(:,:,layer_counter,field_counter) + C%field_shift(field_counter)
     END DO
     END DO

     ! The IM extended fields are mapped (= projected + interpolated) on the GCM grid. For each target grid point the coordinates and
     ! the relative distances of the nearest projected points are stored in the C%sid_filename file, these are
     ! used here to map the fields:
     ! Output: mapped_gcm_field
     CALL MPI_Barrier(PAR%shared_comm)
     CALL oblimap_mapping(oblimap_ddo, C%NX, C%NY, C%NLON, C%NLAT, mask_of_invalid_contributions, im_field, mapped_gcm_field)

     ! Compose the fields which are given back to GCM. Part of these fields (the coordinates with mapping_participation_mask == 0) are
     ! the same as the original GCM fields: initial_gcm_field_*, and the rest of them equal the mapped values as in mapped_gcm_field_*.
     CALL MPI_Barrier(PAR%shared_comm)
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
       gcm_field_(:,:,layer_counter,field_counter) =       mapping_participation_mask_(:,:)  * mapped_gcm_field_ (:,:,layer_counter,field_counter) &
                                                    + (1 - mapping_participation_mask_(:,:)) * initial_gcm_field_(:,:,layer_counter,field_counter)
       WHERE(gcm_field_(:,:,layer_counter,field_counter) == C%invalid_input_value(field_counter)) gcm_field_(:,:,layer_counter,field_counter) = C%invalid_output_value(field_counter)
     END DO
     END DO

     IF(C%oblimap_message_level > 1) THEN
      DO field_counter = 1, C%number_of_mapped_fields
      DO layer_counter = 1, C%number_of_vertical_layers
       IF(C%masked_fields(field_counter)) THEN
        ! This estimate makes only sense for the IM -> GCM mapping direction, because a rectangular grid is assumed in this estimate.
        WRITE(UNIT=*, FMT='(/A, F8.3, 3A, I5, A, F8.3, A/, A/)') ' An optimal alpha_stereographic_config =', &
         C%radians_to_degrees * ASIN(SQRT((COUNT(mask_of_invalid_contributions(C%field_which_determines_invalid_value_mask(field_counter),:,:,layer_counter)) * C%dx * C%dy) / (2._dp * C%pi)) / C%earth_radius), &
         ' degrees, based on the masked IM grid area for the variable "', TRIM(C%im_field_name(field_counter)), '" at layer = ', layer_counter, '. This run uses alpha_stereographic_config = ', C%radians_to_degrees * C%alpha_stereographic, ' degrees.', &
         ' Note that this estimate of an optimal alpha_stereographic_config can be provided with the reverse IM -> GCM mapping only, but concerns both mapping directions.'
       END IF
      END DO
      END DO
     END IF

     IF(C%im_record_range(2) - C%im_record_range(1) > 0) &
      WRITE(UNIT=*, FMT='(A, I4, A)') ' Time record ', 1 + record_counter, ' is written by OBLIMAP.'

     ! Finally the mapped fields gcm_field are written to an output file:
     ! Output: -
     CALL MPI_Barrier(PAR%shared_comm)
     CALL oblimap_write_netcdf_fields(created_gcm_netcdf_file, 1 + record_counter, gcm_field, time)
    END DO

    ! Output: -
    CALL oblimap_close_netcdf_file(im_netcdf_file)
    CALL oblimap_close_netcdf_file(initial_gcm_netcdf_file)
    CALL oblimap_close_netcdf_file(created_gcm_netcdf_file)

    IF(C%reduce_dummy_dimensions) CALL reduce_dummy_dimensions(C%gcm_created_filename, C%number_of_mapped_fields, C%gcm_field_name, C%NLON, C%NLAT)

    IF(C%oblimap_message_level > 0 .and. PAR%rank_shared == 0) THEN
     WRITE(UNIT=*, FMT='(A)')
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
      IF(C%masked_fields(field_counter)) WRITE(UNIT=*, FMT='(A, E24.16, A, A20, A, I3)') &
       ' Using masked mapping for the invalid value = ', C%invalid_input_value(field_counter), ' in the field ', TRIM(C%im_field_name(field_counter)), '   for layer = ', layer_counter
     END DO
     END DO
    END IF

    ! Finishing message:
    WRITE(UNIT=*, FMT='(/3A/2A/)') ' Finished! The file  ', TRIM(C%gcm_created_filename), '  is created. Which can be viewed by:', '  ncview ', TRIM(C%gcm_created_filename)

    ! Output: -
    CALL oblimap_deallocate_ddo(oblimap_ddo)
    !deallocate(mask_of_invalid_contributions)
    !deallocate(gcm_field)
    !deallocate( im_field)
    call MPI_Win_free(mask_of_invalid_contributions_win)
    call MPI_Win_free(gcm_field_win)
    call MPI_Win_free(im_field_win)
    call MPI_Win_free(mapped_gcm_field_win)
    call MPI_Win_free(initial_gcm_field_win)

    call MPI_Win_free(lon_gcm_win)
    call MPI_Win_free(lat_gcm_win)
    call MPI_Win_free(x_coordinates_of_im_grid_points_win)
    call MPI_Win_free(y_coordinates_of_im_grid_points_win)
    call MPI_Win_free(mapping_participation_mask_win)

  END SUBROUTINE oblimap_im_to_gcm_mapping

END MODULE oblimap_im_to_gcm_mapping_module
