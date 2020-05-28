! Compile with:
!  module load c/intel/64/15.0.2; module load openmpi/intel/2.0.0
!  make all
! Run with:
!  ./long-local-run.csh fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example

PROGRAM fortran_mpi_win_allocate_shared_example_program

  USE mpi
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
  IMPLICIT NONE
  
  INTEGER                                :: ierr, status(MPI_STATUS_SIZE)
  INTEGER                                :: process_rank, number_of_processes
  
  INTEGER,  DIMENSION(100, 100)          :: data_array
  INTEGER,  DIMENSION(:, :), POINTER     :: p1
  INTEGER(KIND=MPI_ADDRESS_KIND)         :: windowsize
  INTEGER                                :: disp_unit
  INTEGER                                :: win
  TYPE(C_PTR)                            :: baseptr
      
  ! Split off slave programs
  CALL MPI_INIT(ierr)
  
  ! Get rank of current process and total number of processes
  CALL MPI_COMM_RANK(       MPI_COMM_WORLD, process_rank, ierr)
  CALL MPI_COMM_SIZE(       MPI_COMM_WORLD, number_of_processes, ierr)

  IF (process_rank==0) WRITE(*,'(A)') ' Program fortran_mpi_example starts initializing .....'
  IF (process_rank==0) WRITE(*,'(A,I2,A)') ' MPI master: created ', number_of_processes-1, ' slave programs.'
  
  ! ==========
  ! Code for master:
  IF (process_rank==0) THEN
    ! Fill in some values in data_array:
    data_array        = 3
    data_array(13,37) = 42
    
    WRITE(*,'(A,I5,A)')    ' MPI master: SUM(data_array) = ', SUM(data_array)
    WRITE(*,'(A,I2,A)')    ' MPI master: data_array(13,37) = ', data_array(13,37)
    
    ! Send the array to the slave program
    WRITE(*,'(A)')         ' MPI master: sending data_array to slave...'
    CALL MPI_SEND(data_array, 10000, MPI_INTEGER, 1, 1, MPI_COMM_WORLD, ierr)
    
    ! Allocate some shared memory space with an associated window object:
    windowsize  = 100*100*4_MPI_ADDRESS_KIND
    disp_unit   = 100*4
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p1, [100,100])
    
    ! copy the values from "data_array" to the newly allocated memory space.
    p1 = data_array
    
    ! Send the window to this memory space to the slave.
    CALL MPI_SEND( win, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
  END IF
  ! End of code for master
  ! ==========
  
  ! ==========
  ! Code for slaves:
  IF (process_rank/=0) THEN
  
    ! Receive data_array from master
    CALL MPI_RECV(data_array, 10000, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    WRITE(*,'(A)')         ' MPI  slave: received data_array from master.'
    WRITE(*,'(A,I5,A)')    ' MPI  slave: SUM(data_array) = ', SUM(data_array)
    WRITE(*,'(A,I2,A)')    ' MPI  slave: data_array(13,37) = ', data_array(13,37)
    
    ! Allocate some shared memory space with an associated window object:
    windowsize  = 0_MPI_ADDRESS_KIND
    disp_unit   = 1
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    ! Use MPI_RECV to receive the window to the master's allocated memory space
    CALL MPI_RECV( win, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space
    CALL C_F_POINTER(baseptr, p1, [100,100])
    
    ! See what's in there
    WRITE(*,'(A)')         ' MPI  slave: reading data from master memory space.'
    WRITE(*,'(A,I5,A)')    ' MPI  slave: SUM(p1) = ', SUM(p1)
    WRITE(*,'(A,I2,A)')    ' MPI  slave: p1(13,37) = ', p1(13,37)
    
  END IF
  ! End of code for slaves
  ! ==========
  
  ! Finalize all MPI processes
  CALL MPI_FINALIZE(ierr)

  IF (process_rank==0) WRITE(*,'(A)') ' Program fortran_mpi_example has finished!'

END PROGRAM fortran_mpi_win_allocate_shared_example_program

