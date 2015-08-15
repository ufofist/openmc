module matrix_coo_header

  use list_header, only: ListInt, ListReal

  implicit none
  private

  type, public :: CooMatrix
    integer :: m            ! number of rows
    integer :: n            ! number of columns
    integer :: nnz = 0      ! number of non-zeros
    type(ListInt) :: row_   ! temp list to construct matrix 1 value at a time
    type(ListInt) :: col_   ! temp list to construct matrix 1 value at a time
    type(ListReal) :: data_ ! temp list to construct matrix 1 value at a time
    integer, allocatable :: row(:)  ! row indices
    integer, allocatable :: col(:)  ! column indices
    real(8), allocatable :: data(:) ! value of matrix at (row, col)
    logical :: frozen = .false.     ! whether the matrix has been frozen
  contains
    procedure :: initialize => coo_matrix_initialize
    procedure :: set => coo_matrix_set
    procedure :: freeze => coo_matrix_freeze
    procedure :: to_dense => coo_matrix_to_dense
    procedure :: to_csr => coo_matrix_to_csr
  end type CooMatrix

contains

!===============================================================================
! COO_MATRIX_SET
!===============================================================================

  subroutine coo_matrix_set(this, row, col, data)
    class(CooMatrix), intent(inout) :: this
    integer, intent(in) :: row
    integer, intent(in) :: col
    real(8), intent(in) :: data

    call this%row_%append(row)
    call this%col_%append(col)
    call this%data_%append(data)
  end subroutine coo_matrix_set

  subroutine coo_matrix_initialize(this, row, col, data, dims)
    class(CooMatrix), intent(inout) :: this
    real(8), allocatable, intent(in) :: row(:)
    real(8), allocatable, intent(in) :: col(:)
    real(8), allocatable, intent(in) :: data(:)
    integer, intent(in) :: dims(2)

    ! Set dimensions of matrix
    this%m = dims(1)
    this%n = dims(2)

    ! Allocate arrays for storing data
    this%nnz = size(row)
    allocate(this%row(this%nnz))
    allocate(this%col(this%nnz))
    allocate(this%data(this%nnz))

    ! Copy data
    this%row(:) = row
    this%col(:) = col
    this%data(:) = data

    ! Freeze the matrix
    this%frozen = .true.
  end subroutine coo_matrix_initialize

  subroutine coo_matrix_freeze(this, dims)
    class(CooMatrix), intent(inout) :: this
    integer, intent(in) :: dims(2)

    integer :: i

    ! Set dimensions of matrix
    this%m = dims(1)
    this%n = dims(2)

    ! Allocate arrays for storing data
    this%nnz = this%row_%size()
    allocate(this%row(this%nnz))
    allocate(this%col(this%nnz))
    allocate(this%data(this%nnz))

    ! Copy from lists to arrays
    do i = 1, this%nnz
      this%row(i) = this%row_%get_item(i)
      this%col(i) = this%row_%get_item(i)
      this%data(i) = this%row_%get_item(i)
    end do

    ! Clear lists
    call this%row_%clear()
    call this%col_%clear()
    call this%data_%clear()

    ! Freeze the matrix
    this%frozen = .true.
  end subroutine coo_matrix_freeze

  subroutine coo_matrix_to_dense(this, matrix)
    class(CooMatrix), intent(inout) :: this
    real(8), allocatable, intent(inout) :: matrix(:,:)

    integer :: i

    do i = 1, this%nnz
      matrix(this%row(i), this%col(i)) = this%data(i)
    end do
  end subroutine coo_matrix_to_dense

  subroutine coo_matrix_to_csr(this)
    class(CooMatrix), intent(inout) :: this
  end subroutine coo_matrix_to_csr

end module matrix_coo_header
