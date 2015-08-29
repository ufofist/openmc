module depletion

  use constants
  use depletion_header
  use eigenvalue,        only: run_eigenvalue
  use global
  use matrix_csr_header, only: MatrixCSRReal, MatrixCSRComplex
  use matrix_lil_header, only: MatrixLILReal
  use output,            only: header
  use sparse_header,     only: SparseCsrComplex
  use stl_vector,        only: VectorInt

contains

!===============================================================================
! RUN_DEPLETION
!===============================================================================

  subroutine run_depletion()

    integer :: i

    if (master) call header("DEPLETION SIMULATION", level=1)

    do i = 1, n_depletion_steps
      ! call run_eigenvalue()
    end do

  end subroutine run_depletion

!===============================================================================
! SOLVE_CRAM solves the matrix exponential using the Chebyshev rational
! approximation method. This method is described in M. Pusa and J. Leppanen,
! "Computing the Matrix Exponential in Burnup Calculations," Nucl. Sci. Eng.,
! 164, p. 140--150 (2010) as well as M. Pusa, "Rational Approximations to the
! Matrix Exponential in Burnup Calculations," Nucl. Sci. Eng., 169, p. 155--167
! (2011).
!===============================================================================

  subroutine solve_cram(A, x0, x, t)

    type(MatrixCSRReal), intent(in) :: A    ! burnup matrix
    complex(8), intent(in)  :: x0(:)        ! starting concentration vector
    complex(8), intent(out) :: x(:)         ! end-of-step concentration vector
    real(8),    intent(in)  :: t            ! time interval

    integer :: i     ! row index
    integer :: j     ! column index
    integer :: k     ! loop index for CRAM order
    integer :: i_val ! index in indices/data
    integer :: n     ! size of solution vector
    complex(8), allocatable :: b(:) ! right-hand side
    complex(8), allocatable :: y(:) ! solution to linear system
    type(MatrixCSRReal)    :: fill  ! fill-in matrix
    type(MatrixCSRComplex) :: lu    ! matrix used during LU decomposition

    ! Allocate arrays
    n = size(x0)
    allocate(b(n))
    allocate(y(n))

    ! Perform symbolic factorization on the burnup matrix for LU decomposition
    call symbolic_factorization(A, fill)

    ! Multiply elements in LU by time interval
    fill%data(:) = fill%data(:) * t

    ! Compute first term of solution
    x = CRAM_ALPHA0 * x0

    do k = 1, CRAM_ORDER/2
      ! Copy fill matrix
      lu = fill

      ! Calculate right hand side
      b = CRAM_ALPHA(k) * x0

      ! Subtract theta on diagonal terms
      ROWS: do i = 1, n
        COLUMNS: do i_val = lu%indptr(i), lu%indptr(i+1) - 1
          ! Get index of column
          j = lu%indices(i_val)

          ! Subtract theta from diagonal term
          if (i == j) lu%data(i_val) = lu%data(i_val) - CRAM_THETA(k)
        end do COLUMNS
      end do ROWS

      ! Perform gaussian elimination
      call numerical_elimination(lu, b, y)

      ! Add to solution
      x = x + TWO * y
    end do

  end subroutine solve_cram

!===============================================================================
! SYMBOLIC_FACTORIZATION computes the non-zero structure of the fill-in matrix
! resulting from Gaussian elimination on a sparse matrix 'A'. The algorithm used
! here is the FILL2 algorithm from D. J. Rose and R. E. Tarjan, "Algorithmic
! aspects of vertex elimination on directed graphs," SIAM J. Appl. Math, 40,
! pp. 176--197 (1978).
!===============================================================================

  subroutine symbolic_factorization(matrix, fill)

    type(MatrixCSRReal), intent(in)    :: matrix
    type(MatrixCSRReal), intent(inout) :: fill

    integer :: i         ! row index
    integer :: i_val     ! index in indices/data
    integer :: j         ! column index
    integer :: k         ! another column index
    integer :: n         ! number of columns
    integer, allocatable :: a(:) ! temporary fill list (one row)
    type(VectorInt) :: Omega
    type(MatrixLILReal) :: fill_

    n = matrix%n ! Copy size of matrix

    ! Allocate Omega and fill_row
    allocate(a(n))
    call Omega%reserve(n)

    ! Create LIL sparse matrix using data from the burnup matrix
    call fill_%initialize(matrix)

    ! ==========================================================================
    ! SYMBOLIC DECOMPOSITION

    ROWS: do i = 2, n
      ! Indicate that the Omega set is empty
      call Omega%clear()

      ! Set the vector 'a' to zero
      a(:) = 0

      ! ========================================================================
      ! Find non-zeros in row i in lower triangular part of original matrix

      COLUMNS_IN_ROW_I: do i_val = matrix%indptr(i), matrix%indptr(i+1) - 1
        ! Get index of column
        j = matrix%indices(i_val)

        ! Save positions of non-zero index
        a(j) = 1

        ! If non-zero is in lower triangular part, add it to Omega
        if (j < i) call Omega%push_back(j)
      end do COLUMNS_IN_ROW_I

      ! ========================================================================
      ! Add fill-in where necessary

      ! Loop while the list of neighbor rows is empty
      NEIGHBOR_LIST: do while (Omega%size() > 0)
        ! Get last item on Omega list
        j = Omega%data(Omega%size())
        call Omega%pop_back()

        ! Now loop over the columns in row j
        COLUMNS_IN_ROW_J: do i_val = 1, fill_%rows(j)%size() ! fill%indptr(j), fill%indptr(j+1) - 1
          ! Get index of column
          k = fill_%rows(j)%data(i_val)

          ! If column k is in upper triangular part and the k-th column in row i
          ! is a zero entry, then (i,k) should be added to the fill-in
          if (j < k .and. a(k) == 0) then
            ! Add (i,k) to fill-in matrix
            call fill_%set(i, k, ZERO)

            ! Save k to temporary fill list
            a(k) = 1

            ! If j < k < i, there is a prospect for a longer path to node i
            ! through the nodes j and k, and therefore the node k should be
            ! added to Omega
            if (k < i) call Omega%push_back(k)
          end if
        end do COLUMNS_IN_ROW_J
      end do NEIGHBOR_LIST
    end do ROWS

    ! Convert to CSR format
    call fill_%to_csr(fill)

  end subroutine symbolic_factorization

!===============================================================================
! NUMERICAL_ELIMINATION zeros out entries in the matrix below the diagonal using
! Gaussian elimination. A recent paper -- M. Pusa and J. Leppanen, "An efficient
! implementation of the Chebyshev Rational Approximation Method for solving the
! burnup equations," Proc. PHYSOR 2012, Knoxville, TN, Apr. 15-20, 2012 --
! argues that partial pivoting is not required (based on properties of the
! burnup matrix). The actual algorithm for elimination on the sparse matrix is
! the CELIMINATE algorithm from R. E. Tarjan, "Graph Theory and Gaussian
! Elimination," Stanford STAN-CS-75-526, November 1975.
!===============================================================================

  subroutine numerical_elimination(fill, b, x)

    type(MatrixCSRComplex), intent(inout) :: fill ! fill matrix
    complex(8),         intent(inout) :: b(:) ! right-hand side
    complex(8),         intent(out)   :: x(:) ! solution

    integer :: i          ! row index
    integer :: i_val      ! index in indices/data
    integer :: i_val2     ! another index in indices/data
    integer :: j          ! column index
    integer :: k          ! another column index
    integer :: n          ! number of columns
    complex(8) :: fill_ij ! (i,j) element in fill matrix
    complex(8) :: fill_jj ! (j,j) element in fill matrix
    complex(8) :: fill_jk ! (j,k) element in fill matrix
    complex(8) :: frac    ! F_ij / F_jj
    complex(8) :: sum     ! temporary sum for back substitution
    complex(8), allocatable :: v(:)
    complex(8), allocatable :: diag(:)

    ! Size of matrix
    n = fill%n

    ! Allocate arrays for storing a single row and the diagonal
    allocate(v(n))
    allocate(diag(n))

    ! Initialize v and diag
    v(:) = ZERO
    diag(:) = ZERO

    ! Since diag(1) is not set below (loop only covers rows 2 and above),
    ! manually set it here
    if (fill%indices(1) == 1) diag(1) = fill%data(1)

    ROWS: do i = 2, n

      ! Copy row i to vector v and save diagonal
      do i_val = fill%indptr(i), fill%indptr(i+1) - 1
        j = fill%indices(i_val)

        ! Copy value into v vector
        v(j) = fill%data(i_val)

        ! Save diagonal
        if (i == j) diag(j) = v(j)
      end do

      COL_IN_ROW_I: do i_val = fill%indptr(i), fill%indptr(i+1) - 1
        j = fill%indices(i_val)

        if (j < i) then
          fill_ij = v(j)
          fill_jj = diag(j)

          ! TODO: Check norm of diagonal term

          ! Update right-hand side term in row i
          frac = fill_ij/fill_jj
          b(i) = b(i) - frac*b(j)

          ! Update matrix elements
          COL_IN_ROW_J: do i_val2 = fill%indptr(j), fill%indptr(j+1) - 1
            k = fill%indices(i_val2)

            ! We only need to update the terms
            if (j < k) then
              fill_jk = fill%indices(i_val2)
              v(k) = v(k) - frac*fill_jk
            end if
          end do COL_IN_ROW_J
        end if

      end do COL_IN_ROW_I

      ! ========================================================================
      ! UPDATE ROW I

      do i_val = fill%indptr(i), fill%indptr(i+1) - 1
        j = fill%indices(i_val)

        ! Update diagonal
        if (i == j) diag(j) = v(j)

        ! Update element (i,j)
        fill%data(i_val) = v(j)
      end do

    end do ROWS

    !===========================================================================
    ! BACK_SUBSTITUTION

    do i = n, 1, -1
      ! Initialize sum
      sum = ZERO

      do i_val = fill%indptr(i), fill%indptr(i+1) - 1
        ! Get index of column
        j = fill%indices(i_val)

        ! Add A_ij * x_j to the sum -- only use upper triangular portion
        if (j > i) sum = sum + fill%data(i_val) * x(j)
      end do

      ! TODO: check norm of diagonal

      ! Calculate i-th element of solution
      x(i) = (b(i) - sum)/diag(i)
    end do

  end subroutine numerical_elimination

end module depletion
