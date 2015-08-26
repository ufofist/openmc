module matrix_csr_header

  implicit none
  private

  type, public :: MatrixCSRReal
    integer :: m            ! number of rows
    integer :: n            ! number of columns
    integer :: nnz = 0      ! number of non-zeros
    integer, allocatable :: indptr(:)  ! cumulative number of non-zeros by row
    integer, allocatable :: indices(:) ! column indices
    real(8), allocatable :: data(:)    ! matrix value vector
  contains
    generic :: initialize => csrreal_initialize_from_dense, &
         csrreal_initialize_from_arrays
    procedure, private :: csrreal_initialize_from_dense
    procedure, private :: csrreal_initialize_from_arrays
  end type MatrixCSRReal

contains

  subroutine csrreal_initialize_from_dense(this, matrix)
    class(MatrixCSRReal), intent(inout) :: this
    real(8), allocatable :: matrix(:,:)

    integer :: i, j ! row and column indices
    integer :: m, n ! number of rows and columns
    integer :: nnz  ! number of non-zeros
    integer :: innz ! indices for non-zeros

    ! Set size of matrix
    m = size(matrix, 1)
    n = size(matrix, 2)
    this%m = m
    this%n = n

    ! Determine number of non-zeros
    nnz = count(abs(matrix) > 0.0_8)
    this%nnz = nnz

    ! Allocate arrays
    allocate(this%indptr(m + 1))
    allocate(this%indices(nnz))
    allocate(this%data(nnz))

    ! Construct sparse matrix data
    innz = 1
    do i = 1, m
      this%indptr(i) = innz
      do j = 1, n
        if (abs(matrix(i,j)) > 0.0_8) then
          this%indices(innz) = j
          this%data(innz) = matrix(i,j)
          innz = innz + 1
        end if
      end do
    end do
    this%indptr(i) = innz

  end subroutine csrreal_initialize_from_dense

  subroutine csrreal_initialize_from_arrays(this, indptr, indices, data, dims)
    class(MatrixCSRReal), intent(inout) :: this
    integer, allocatable :: indptr(:)
    integer, allocatable :: indices(:)
    real(8), allocatable :: data(:)
    integer :: dims(2)

    this%m = dims(1)
    this%n = dims(2)
    this%nnz = size(indices)

    ! Allocate and copy data
    allocate(this%indptr(this%m + 1), SOURCE=indptr)
    allocate(this%indices(this%nnz), SOURCE=indices)
    allocate(this%data(this%nnz), SOURCE=data)

  end subroutine csrreal_initialize_from_arrays

end module matrix_csr_header
