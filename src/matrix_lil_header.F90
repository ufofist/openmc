module matrix_lil_header

  use matrix_csr_header, only: MatrixCSRReal
  use stl_vector, only: VectorInt, VectorReal

  implicit none
  private

  type, public :: MatrixLILReal
    integer :: m
    integer :: n
    integer :: nnz = 0
    type(VectorInt), allocatable :: rows(:)
    type(VectorReal), allocatable :: data(:)
  contains
    generic :: initialize => &
         initialize_csr_real, &
         initialize_dense_real, &
         initialize_empty_real
    procedure, private :: initialize_csr_real
    procedure, private :: initialize_dense_real
    procedure, private :: initialize_empty_real
    procedure :: set => set_real
    procedure :: sort => sort_real
    procedure :: to_dense => to_dense_real
    procedure :: to_csr => to_csr_real
  end type MatrixLILReal

contains

  subroutine initialize_csr_real(this, csr)
    class(MatrixLILReal), intent(inout) :: this
    type(MatrixCSRReal), intent(in) :: csr

    integer :: i, j, k
    integer :: nnz

    ! Set size of matrix
    this%m = csr%m
    this%n = csr%n

    ! Allocate vectors
    allocate(this%rows(this%m))
    allocate(this%data(this%m))

    ! Copy data from CSR matrix
    do i = 1, csr%m
      ! Preallocate indices/data
      nnz = csr%indptr(i+1) - csr%indptr(i)
      call this%rows(i)%reserve(nnz)
      call this%data(i)%reserve(nnz)

      ! Add each row of data
      do k = csr%indptr(i), csr%indptr(i+1) - 1
        j = csr%indices(k)
        call this%set(i, j, csr%data(k))
      end do
    end do
  end subroutine initialize_csr_real

  subroutine initialize_dense_real(this, matrix)
    class(MatrixLILReal), intent(inout) :: this
    real(8), intent(in), allocatable :: matrix(:,:)

    integer :: i
    integer :: j

    ! Set size of matrix
    this%m = size(matrix, 1)
    this%n = size(matrix, 2)

    ! Allocate vectors
    allocate(this%rows(this%m))
    allocate(this%data(this%m))

    do i = 1, this%m
      do j = 1, this%n
        if (abs(matrix(i,j)) > 0.0_8) then
          call this%set(i, j, matrix(i,j))
        end if
      end do
    end do
  end subroutine initialize_dense_real

  subroutine initialize_empty_real(this, m, n)
    class(MatrixLILReal), intent(inout) :: this
    integer, intent(in) :: m
    integer, intent(in) :: n

    ! Set dimensions of matrix
    this%m = m
    this%n = n

    ! Allocate vectors
    allocate(this%rows(this%m))
    allocate(this%data(this%m))
  end subroutine initialize_empty_real

  subroutine set_real(this, row, col, data)
    class(MatrixLILReal), intent(inout) :: this
    integer, intent(in) :: row
    integer, intent(in) :: col
    real(8), intent(in) :: data

    call this%rows(row)%push_back(col)
    call this%data(row)%push_back(data)
    this%nnz = this%nnz + 1
  end subroutine set_real

  subroutine sort_real(this)
    class(MatrixLILReal), intent(inout) :: this

    integer :: i, j, n
    integer :: row
    integer :: col
    real(8) :: data

    ! Perform an insertion sort on each row
    do row = 1, this%m
      n = this%rows(row)%size()
      if (n > 1) then
        do j = 2, n
          col = this%rows(row)%data(j)
          data = this%data(row)%data(j)
          i = j - 1
          do while (i > 0 .and. this%rows(row)%data(i) > col)
            this%rows(row)%data(i+1) = this%rows(row)%data(i)
            this%data(row)%data(i+1) = this%data(row)%data(i)
            i = i - 1
          end do
          this%rows(row)%data(i+1) = col
          this%data(row)%data(i+1) = data
        end do
      end if
    end do
  end subroutine sort_real

  subroutine to_dense_real(this, matrix)
    class(MatrixLILReal), intent(inout) :: this
    real(8), allocatable, intent(out) :: matrix(:,:)

    integer :: i, j, k

    do i = 1, this%m
      do k = 1, this%rows(i)%size()
        j = this%rows(i)%data(k)
        matrix(i,j) = this%data(i)%data(k)
      end do
    end do
  end subroutine to_dense_real

  subroutine to_csr_real(this, csr)
    class(MatrixLILReal), intent(inout) :: this
    type(MatrixCSRReal), intent(out) :: csr

    integer :: innz, i, k

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
    innz = 1
    do i = 1, this%m
      csr%indptr(i) = innz
      do k = 1, this%rows(i)%size()
        csr%indices(innz) = this%rows(i)%data(k)
        csr%data(innz) = this%data(i)%data(k)
        innz = innz + 1
      end do
    end do
    csr%indptr(this%m + 1) = innz
  end subroutine to_csr_real


end module matrix_lil_header
