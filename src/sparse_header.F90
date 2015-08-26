module sparse_header

  use constants, only: NULL_COLUMN, ZERO

  implicit none

!===============================================================================
! SRARSECSRCOMPLEX stores indices and values of a sparse matrix stored in
! compressed sparse row (CSR) format. This format has the non-zeros of the
! matrix stored in contiguous memory locations. The 'indices' array stores the
! column indices of the elements of the non-zero array and the 'row_ptr' array
! stores stores the indices in the indices and data arrays corresponding to
! beginning of each row. The total number of elements required is 2*n_nonzero +
! n + 1.
!===============================================================================

  type SparseCsrComplex
    ! Information about size of matrix
    integer :: m         ! number of rows
    integer :: n         ! number of columns
    integer :: nnz ! number of non-zero elements

    ! Storage of locations and values
    integer,    allocatable :: indptr(:) ! indices within values array
    integer,    allocatable :: indices(:) ! column indices with non-zeros
    complex(8), allocatable :: data(:)  ! non-zero values in matrix
  contains
    procedure :: init => sparse_csr_init_complex
    procedure :: expand => sparse_csr_expand_complex
  end type SparseCsrComplex

contains

!===============================================================================
! SPARSE_CSR_INIT allocates space for a compressed space row based on the size
! of the matrix and the expected number of non-zeros.
!===============================================================================

  subroutine sparse_csr_init_complex(A, m, n, nnz)

    class(SparseCsrComplex) :: A
    integer, intent(in)     :: m
    integer, intent(in)     :: n
    integer, intent(in)     :: nnz

    ! Set integer constants
    A % m = m
    A % n = n
    A % nnz = nnz

    ! If any components are already allocated, remove that allocation
    if (allocated(A % indptr)) deallocate(A % indptr)
    if (allocated(A % indices)) deallocate(A % indices)
    if (allocated(A % data)) deallocate(A % data)

    ! Allocate components
    allocate(A % indptr(m + 1))
    allocate(A % indices(nnz))
    allocate(A % data(nnz))

    ! Initialize component arrays
    A % indptr = 0
    A % indices = NULL_COLUMN
    A % data = ZERO

  end subroutine sparse_csr_init_complex

!===============================================================================
! SPARSE_CSR_EXPAND adds extra space to a row in a sparse matrix. This is used
! in sparse factorization when fill-in elements are added.
!===============================================================================

  subroutine sparse_csr_expand_complex(A, i, space)

    class(SparseCsrComplex) :: A     ! sparse matrix
    integer, intent(in)     :: i     ! row to expand
    integer, intent(in)     :: space ! number of items to add

    integer :: n      ! original size of A%indices
    integer :: i_val  ! index in indices/data
    integer :: i_val2 ! another index in indices/data
    integer,    allocatable :: indices(:) ! expanded copy of indices
    complex(8), allocatable :: data(:)  ! expanded copy of data

    ! Get original size of indices
    n = size(A % indices)

    ! Allocate temporary arrays
    allocate(indices(n + space))
    allocate(data(n + space))

    ! Get pointers to start of rows i and i+1
    i_val = A % indptr(i)
    i_val2 = A % indptr(i+1)

    ! Copy rows 1, i-1
    indices(1 : i_val-1) = A % indices(1 : i_val-1)
    data(1 : i_val-1)  = A % data(1 : i_val-1)

    ! Copy row i
    indices(i_val : i_val2-1) = A % indices(i_val : i_val2-1)
    data(i_val : i_val2-1)  = A % data(i_val : i_val2-1)

    ! Initialize extra space in row i
    indices(i_val2 : i_val2+space-1) = -1
    data(i_val2 : i_val2+space-1)   = ZERO

    ! Copy rows i+1, n
    indices(i_val2+space:) = A % indices(i_val2:)
    data(i_val2+space:)  = A % data(i_val2:)

    ! Move allocation back
    call move_alloc(FROM=indices, TO=A%indices)
    call move_alloc(FROM=data, TO=A%data)

    ! Adjust row pointers
    A % indptr(i+1:) = A % indptr(i+1:) + space

    ! Increase number of non-zeros
    A % nnz = A % nnz + space

  end subroutine sparse_csr_expand_complex

end module sparse_header
