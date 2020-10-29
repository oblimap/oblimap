module mpi_helpers_mod
  use mpi_f08
  use oblimap_configuration_module
  implicit none

  private

  public :: decompose
  interface decompose
    module procedure :: decompose_range
    module procedure :: decompose_size
  end interface decompose

  public :: alloc_shared
  interface alloc_shared
    module procedure alloc_shared_1d_dp
    module procedure alloc_shared_1d_i
    module procedure alloc_shared_2d_dp
    module procedure alloc_shared_2d_i

    module procedure alloc_shared_4d_l
    module procedure alloc_shared_4d_dp
  end interface alloc_shared

contains
  !! https://arxiv.org/pdf/1804.09536.pdf
  !! balanced block decomposition
  !! decompose N amongst M procs for proc p:
  !! n: number of elements
  !! s: starting index
  subroutine decompose_range(S, N, M, p, n0, n1)
    implicit none
    integer, intent(in) :: S, N, M, p
    integer, intent(inout) :: n0, n1
    integer :: n_, q_, r_

     q_ = (N-S+1) / M ;
     r_ = mod((N-S+1), M) ! N - (INT(N/M) * M) ;
     if (r_ > p) then
             n_ = q_ + 1
             n0 = q_ * p + p
     else
             n_ = q_
             n0 = q_ * p + r_
     endif
     n1 = n0 + n_ + S -1
     n0 = n0 + S !from n0:n1
  end subroutine decompose_range

  subroutine decompose_size(N, M, p, n0, n1)
    implicit none
    integer, intent(in) :: N, M, p
    integer, intent(inout) :: n0, n1

    call decompose_range(1, N, M, p, n0, n1)
  end subroutine decompose_size

  subroutine alloc_shared_1d_dp( nx, x0, x1 &
                               , a_sm, a_, a_win, shared_comm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: nx, x0, x1
    real(dp), dimension(:), pointer, intent(out) :: a_sm, a_
    type(MPI_Win), intent(out) :: a_win
    type(MPI_Comm), intent(in) :: shared_comm

    type(MPI_Info) :: info
    type(C_PTR) :: a_sm_c

    integer :: rank_shared, rank
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, real_extent, size_sm
    integer :: sizeofreal
    integer(KIND=MPI_ADDRESS_KIND) :: size_

    call MPI_Comm_rank(shared_comm, rank_shared)

    size_sm = 0
    if(rank_shared == 0) size_sm = (nx)

    call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lowerbound, real_extent)
    sizeofreal = real_extent

    size_ = size_sm*sizeofreal
    info = MPI_INFO_NULL
    call MPI_Win_allocate_shared(size_, sizeofreal, info, shared_comm, a_sm_c, a_win)
    rank = 0
    call MPI_Win_shared_query(a_win, rank, size_sm, sizeofreal, a_sm_c)
    size_sm = size_sm/sizeofreal
    call c_f_pointer(a_sm_c, a_sm, (/ nx /))

    a_ => a_sm(x0:x1)

  end subroutine alloc_shared_1d_dp

  subroutine alloc_shared_2d_dp( nx, x0, x1 &
                               , ny, y0, y1 &
                               , a_sm, a_, a_win, shared_comm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: nx, x0, x1
    integer, intent(in) :: ny, y0, y1
    real(dp), dimension(:,:), pointer, intent(out) :: a_sm, a_
    type(MPI_Win), intent(out) :: a_win
    type(MPI_Comm), intent(in) :: shared_comm

    type(MPI_Info) :: info
    type(C_PTR) :: a_sm_c

    integer :: rank_shared, rank
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, real_extent, size_sm
    integer :: sizeofreal
    integer :: ii, ij
    integer(KIND=MPI_ADDRESS_KIND) :: size_

    call MPI_Comm_rank(shared_comm, rank_shared)

    size_sm = 0
    if(rank_shared == 0) size_sm = (nx)*(ny)

    call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lowerbound, real_extent)
    sizeofreal = real_extent
    !! USE mpi_f08
    !! MPI_Win_allocate_shared(size, disp_unit, info, comm, baseptr, win, ierror)
    !!     USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
    !!     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
    !!     INTEGER, INTENT(IN) :: disp_unit
    !!     TYPE(MPI_Info), INTENT(IN) :: info
    !!     TYPE(MPI_Comm), INTENT(IN) :: comm
    !!     TYPE(C_PTR), INTENT(OUT) :: baseptr
    !!     TYPE(MPI_Win), INTENT(OUT) :: win
    !!     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
    size_ = size_sm*sizeofreal
    info = MPI_INFO_NULL
    call MPI_Win_allocate_shared(size_, sizeofreal, info, shared_comm, a_sm_c, a_win)
    rank = 0
    call MPI_Win_shared_query(a_win, rank, size_sm, sizeofreal, a_sm_c)
    size_sm = size_sm/sizeofreal
    call c_f_pointer(a_sm_c, a_sm, (/ nx, ny /))

    a_ => a_sm(x0:x1 &
              ,y0:y1)

   !if(rank_shared == 0) a_sm = -1
   !call MPI_Barrier(shared_comm, ier)

   !!a_ = 0 !rank_shared+10
   !!call MPI_Barrier(shared_comm, ier)

   !! if(rank_shared == 0) then
   !!   do ii=1,nx
   !!     do ij=1,ny
   !!       write(*,'(e14.6$)') a_sm(ii,ij)
   !!     end do
   !!     print*, ''
   !!   end do
   !!   print*, 'done'
   !! endif


  end subroutine alloc_shared_2d_dp

  subroutine alloc_shared_1d_i( nx, x0, x1 &
                              , a_sm, a_, a_win, shared_comm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: nx, x0, x1
    integer, dimension(:), pointer, intent(out) :: a_sm, a_
    type(MPI_Win), intent(out) :: a_win
    type(MPI_Comm), intent(in) :: shared_comm

    type(MPI_Info) :: info
    type(C_PTR) :: a_sm_c

    integer :: rank_shared, rank
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, real_extent, size_sm
    integer :: sizeofreal
    integer :: ii, ij
    integer(KIND=MPI_ADDRESS_KIND) :: size_

    call MPI_Comm_rank(shared_comm, rank_shared)

    size_sm = 0
    if(rank_shared == 0) size_sm = (nx)

    call MPI_Type_get_extent(MPI_INTEGER, lowerbound, real_extent)
    sizeofreal = real_extent

    size_ = size_sm*sizeofreal
    info = MPI_INFO_NULL
    call MPI_Win_allocate_shared(size_, sizeofreal, info, shared_comm, a_sm_c, a_win)
    rank = 0
    call MPI_Win_shared_query(a_win, rank, size_sm, sizeofreal, a_sm_c)
    size_sm = size_sm/sizeofreal
    call c_f_pointer(a_sm_c, a_sm, (/ nx /))

    a_ => a_sm(x0:x1)

  end subroutine alloc_shared_1d_i

  subroutine alloc_shared_2d_i( nx, x0, x1 &
                              , ny, y0, y1 &
                              , a_sm, a_, a_win, shared_comm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: nx, x0, x1
    integer, intent(in) :: ny, y0, y1
    integer, dimension(:,:), pointer, intent(out) :: a_sm, a_
    type(MPI_Win), intent(out) :: a_win
    type(MPI_Comm), intent(in) :: shared_comm

    type(MPI_Info) :: info
    type(C_PTR) :: a_sm_c

    integer :: rank_shared, rank
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, real_extent, size_sm
    integer :: sizeofreal
    integer :: ii, ij
    integer(KIND=MPI_ADDRESS_KIND) :: size_

    call MPI_Comm_rank(shared_comm, rank_shared)

    size_sm = 0
    if(rank_shared == 0) size_sm = (nx)*(ny)

    call MPI_Type_get_extent(MPI_INTEGER, lowerbound, real_extent)
    sizeofreal = real_extent

    size_ = size_sm*sizeofreal
    info = MPI_INFO_NULL
    call MPI_Win_allocate_shared(size_, sizeofreal, info, shared_comm, a_sm_c, a_win)
    rank = 0
    call MPI_Win_shared_query(a_win, rank, size_sm, sizeofreal, a_sm_c)
    size_sm = size_sm/sizeofreal
    call c_f_pointer(a_sm_c, a_sm, (/ nx, ny /))

    a_ => a_sm(x0:x1 &
              ,y0:y1)

  end subroutine alloc_shared_2d_i

  subroutine alloc_shared_4d_l ( nx, x0, x1 &
                               , ny, y0, y1 &
                               , nz, z0, z1 &
                               , nk, k0, k1 &
                               , a_sm, a_, a_win, shared_comm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: nx, x0, x1
    integer, intent(in) :: ny, y0, y1
    integer, intent(in) :: nz, z0, z1
    integer, intent(in) :: nk, k0, k1
    logical, dimension(:,:,:,:), pointer, intent(out) :: a_sm, a_
    type(MPI_Win), intent(out) :: a_win
    type(MPI_Comm), intent(in) :: shared_comm

    type(MPI_Info) :: info
    type(C_PTR) :: a_sm_c

    integer :: rank_shared, rank
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, real_extent, size_sm
    integer :: sizeofreal
    integer :: ii, ij
    integer(KIND=MPI_ADDRESS_KIND) :: size_

    call MPI_Comm_rank(shared_comm, rank_shared)

    size_sm = 0
    if(rank_shared == 0) size_sm = (nx)*(ny)*(nz)*(nk)

    call MPI_Type_get_extent(MPI_LOGICAL, lowerbound, real_extent)
    sizeofreal = real_extent

    size_ = size_sm*sizeofreal
    info = MPI_INFO_NULL
    call MPI_Win_allocate_shared(size_, sizeofreal, info, shared_comm, a_sm_c, a_win)
    rank = 0
    call MPI_Win_shared_query(a_win, rank, size_sm, sizeofreal, a_sm_c)
    size_sm = size_sm/sizeofreal
    call c_f_pointer(a_sm_c, a_sm, (/ nx, ny, nz, nk /))

    a_ => a_sm(x0:x1 &
              ,y0:y1 &
              ,z0:z1 &
              ,k0:k1 &
              )

  end subroutine alloc_shared_4d_l

  subroutine alloc_shared_4d_dp( nx, x0, x1 &
                               , ny, y0, y1 &
                               , nz, z0, z1 &
                               , nk, k0, k1 &
                               , a_sm, a_, a_win, shared_comm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: nx, x0, x1
    integer, intent(in) :: ny, y0, y1
    integer, intent(in) :: nz, z0, z1
    integer, intent(in) :: nk, k0, k1
    real(dp), dimension(:,:,:,:), pointer, intent(out) :: a_sm, a_
    type(MPI_Win), intent(out) :: a_win
    type(MPI_Comm), intent(in) :: shared_comm

    type(MPI_Info) :: info
    type(C_PTR) :: a_sm_c

    integer :: rank_shared, rank
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, real_extent, size_sm
    integer :: sizeofreal
    integer :: ii, ij
    integer(KIND=MPI_ADDRESS_KIND) :: size_

    call MPI_Comm_rank(shared_comm, rank_shared)

    size_sm = 0
    if(rank_shared == 0) size_sm = (nx)*(ny)*(nz)*(nk)

    call MPI_Type_get_extent(MPI_DOUBLE_PRECISION, lowerbound, real_extent)
    sizeofreal = real_extent

    size_ = size_sm*sizeofreal
    info = MPI_INFO_NULL
    call MPI_Win_allocate_shared(size_, sizeofreal, info, shared_comm, a_sm_c, a_win)
    rank = 0
    call MPI_Win_shared_query(a_win, rank, size_sm, sizeofreal, a_sm_c)
    size_sm = size_sm/sizeofreal
    call c_f_pointer(a_sm_c, a_sm, (/ nx, ny, nz, nk /))

    a_ => a_sm(x0:x1 &
              ,y0:y1 &
              ,z0:z1 &
              ,k0:k1 &
              )

  end subroutine alloc_shared_4d_dp

end module mpi_helpers_mod
