module depletion_header

  use constants, only: NULL_COLUMN, ZERO, NONE

  implicit none

!===============================================================================
! DEPLETIONSTEP
!===============================================================================

  type :: DepletionStep
    integer :: units        ! MWd/kgU or days
    real(8) :: value        ! user-provided value
    real(8) :: time  = NONE ! length of depletion step in seconds
    real(8) :: power = 1.0  ! power fraction for this step
  end type DepletionStep

!===============================================================================
! SPARSEMATRIX
!===============================================================================

  type SparseMatrix
    ! Information about size of matrix
    integer :: m         ! number of rows
    integer :: n         ! number of columns
    integer :: n_nonzero ! number of non-zero elements

    ! Storage of locations and values
    integer,    allocatable :: row_ptr(:) ! indices within values array
    integer,    allocatable :: columns(:) ! column indices with non-zeros
    complex(8), allocatable :: values(:)  ! non-zero values in matrix
  end type SparseMatrix

contains

  subroutine sparse_matrix_init(A, m, n, n_nonzero)

    type(SparseMatrix), intent(inout) :: A
    integer,            intent(in)    :: m
    integer,            intent(in)    :: n
    integer,            intent(in)    :: n_nonzero

    ! Set integer constants
    A % m = m
    A % n = n
    A % n_nonzero = n_nonzero

    ! If any components are already allocated, remove that allocation
    if (allocated(A % row_ptr)) deallocate(A % row_ptr)
    if (allocated(A % columns)) deallocate(A % columns)
    if (allocated(A % values)) deallocate(A % values)

    ! Allocate components
    allocate(A % row_ptr(m + 1))
    allocate(A % columns(n_nonzero))
    allocate(A % values(n_nonzero))

    ! Initialize component arrays
    A % row_ptr = 0
    A % columns = NULL_COLUMN
    A % values = ZERO

  end subroutine sparse_matrix_init

end module depletion_header
