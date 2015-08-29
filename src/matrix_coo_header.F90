module matrix_coo_header

  use matrix_csr_header, only: MatrixCSRReal
  use stl_vector, only: VectorInt, VectorReal

  implicit none
  private

  type, public :: MatrixCOOReal
    integer :: m             ! number of rows
    integer :: n             ! number of columns
    integer :: nnz = 0       ! number of non-zeros
    type(VectorInt) :: row   ! row indices
    type(VectorInt) :: col   ! column indices
    type(VectorReal) :: data ! value of matrix at (row, col)
  contains
    generic :: initialize => &
         initialize_arrays_real, &
         initialize_csr_real, &
         initialize_dense_real, &
         initialize_empty_real
    procedure, private :: initialize_arrays_real
    procedure, private :: initialize_csr_real
    procedure, private :: initialize_dense_real
    procedure, private :: initialize_empty_real
    procedure :: set => set_real
    procedure :: sort => sort_real
    procedure :: to_dense => to_dense_real
    procedure :: to_csr => to_csr_real
  end type MatrixCOOReal

contains

!===============================================================================
! Implementation of MatrixCOOReal
!===============================================================================

  subroutine initialize_arrays_real(this, row, col, data, dims)
    class(MatrixCOOReal), intent(inout) :: this
    real(8), allocatable, intent(in) :: row(:)
    real(8), allocatable, intent(in) :: col(:)
    real(8), allocatable, intent(in) :: data(:)
    integer, intent(in) :: dims(2)

    ! Set dimensions of matrix
    this%m = dims(1)
    this%n = dims(2)

    ! Clear existing data
    call this%row%clear()
    call this%row%shrink_to_fit()
    call this%col%clear()
    call this%col%shrink_to_fit()
    call this%data%clear()
    call this%data%shrink_to_fit()

    ! Allocate arrays for storing data
    this%nnz = size(row)
    call this%row%reserve(this%nnz)
    call this%col%reserve(this%nnz)
    call this%data%reserve(this%nnz)

    ! Copy data
    this%row%data(:) = row
    this%col%data(:) = col
    this%data%data(:) = data
  end subroutine initialize_arrays_real

  subroutine initialize_csr_real(this, csr)
    class(MatrixCOOReal), intent(inout) :: this
    type(MatrixCSRReal), intent(in) :: csr

    integer :: i, j, k

    ! Set size of matrix
    this%m = csr%m
    this%n = csr%n

    ! Preallocate arrays
    call this%row%reserve(csr%nnz)
    call this%col%reserve(csr%nnz)
    call this%data%reserve(csr%nnz)

    ! Copy data from CSR matrix
    do i = 1, csr%m
      do k = csr%indptr(i), csr%indptr(i+1) - 1
        j = csr%indices(k)
        call this%set(i, j, csr%data(k))
      end do
    end do
  end subroutine initialize_csr_real

  subroutine initialize_dense_real(this, matrix)
    class(MatrixCOOReal), intent(inout) :: this
    real(8), intent(in), allocatable :: matrix(:,:)

    integer :: i
    integer :: j

    this%m = size(matrix, 1)
    this%n = size(matrix, 2)

    do i = 1, this%m
      do j = 1, this%n
        if (abs(matrix(i,j)) > 0.0_8) then
          call this%set(i, j, matrix(i,j))
        end if
      end do
    end do
  end subroutine initialize_dense_real

  subroutine initialize_empty_real(this, m, n)
    class(MatrixCOOReal), intent(inout) :: this
    integer, intent(in) :: m
    integer, intent(in) :: n

    ! Set dimensions of matrix
    this%m = m
    this%n = n
  end subroutine initialize_empty_real

  subroutine set_real(this, row, col, data)
    class(MatrixCOOReal), intent(inout) :: this
    integer, intent(in) :: row
    integer, intent(in) :: col
    real(8), intent(in) :: data

    call this%row%push_back(row)
    call this%col%push_back(col)
    call this%data%push_back(data)
    this%nnz = this%nnz + 1
  end subroutine set_real

  subroutine sort_real(this)
    class(MatrixCOOReal), intent(inout) :: this

    if (this%nnz > 0) call quicksort_real(this, 1, this%nnz)
  end subroutine sort_real

  recursive subroutine quicksort_real(this, lo, hi)
    class(MatrixCOOReal), intent(inout) :: this
    integer, intent(in) :: lo
    integer, intent(in) :: hi

    integer :: p

    if (lo < hi) then
      p = partition_real(this, lo, hi)
      call quicksort_real(this, lo, p)
      call quicksort_real(this, p + 1, hi)
    end if
  end subroutine quicksort_real

  function partition_real(this, lo, hi) result(j)
    class(MatrixCOOReal), intent(inout) :: this
    integer, intent(in) :: lo
    integer, intent(in) :: hi
    integer :: j

    integer :: i
    integer :: pivot
    integer :: test
    integer :: tmp_row, tmp_col
    real(8) :: tmp_data

    ! Since we need to sort both rows and columns, we can use the value (row*N +
    ! col) where N is the number of columns as the basis for sorting

    pivot = this%row%data(lo)*this%n + this%col%data(lo)
    i = lo - 1
    j = hi + 1
    do
      do
        j = j - 1
        test = this%row%data(j)*this%n + this%col%data(j)
        if (test <= pivot) exit
      end do
      do
        i = i + 1
        test = this%row%data(i)*this%n + this%col%data(i)
        if (test >= pivot) exit
      end do
      if (i < j) then
        ! Swap i with j
        tmp_row = this%row%data(j)
        tmp_col = this%col%data(j)
        tmp_data = this%data%data(j)
        this%row%data(j) = this%row%data(i)
        this%col%data(j) = this%col%data(i)
        this%data%data(j) = this%data%data(i)
        this%row%data(i) = tmp_row
        this%col%data(i) = tmp_col
        this%data%data(i) = tmp_data
      else
        exit
      end if
    end do
  end function partition_real

  subroutine to_dense_real(this, matrix)
    class(MatrixCOOReal), intent(inout) :: this
    real(8), allocatable, intent(out) :: matrix(:,:)

    integer :: i

    do i = 1, this%nnz
      matrix(this%row%data(i), this%col%data(i)) = this%data%data(i)
    end do
  end subroutine to_dense_real

  subroutine to_csr_real(this, csr)
    class(MatrixCOOReal), intent(inout) :: this
    type(MatrixCSRReal), intent(out) :: csr

    integer :: i
    integer :: row

    ! Make sure COO data is sorted
    call this%sort()

    ! Set dimensions and number of non-zeros
    csr%m = this%m
    csr%n = this%n
    csr%nnz = this%nnz

    ! Allocate arrays for CSR
    allocate(csr%indptr(this%m + 1))
    allocate(csr%indices(this%nnz))
    allocate(csr%data(this%nnz))

    ! Create indptr array
    i = 1
    do row = 1, this%m
      do while (this%row%data(i) < row)
        i = i + 1
      end do
      csr%indptr(row) = i
    end do
    csr%indptr(this%m + 1) = this%nnz + 1

    ! Since columns are already sorted, we can directly copy values for indices
    ! and data
    csr%indices(:) = this%col%data(1:this%nnz)
    csr%data(:) = this%data%data(1:this%nnz)

  end subroutine to_csr_real

end module matrix_coo_header
