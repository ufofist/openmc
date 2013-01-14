module depletion_header

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

end module depletion_header
