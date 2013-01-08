module endf_header

  implicit none

!===============================================================================
! ENDFCONT represents a control record
!===============================================================================

  type EndfCont
    real(8) :: rvalues(2) ! C1, C2
    integer :: ivalues(4) ! L1, L2, N1, N2
    integer :: MAT
    integer :: MF
    integer :: MT
    integer :: NS
  end type EndfCont

!===============================================================================
! ENDFLIST represents a list record containing a series of numbers
!===============================================================================

  type EndfList
    real(8) :: rvalues(2) ! C1, C2
    integer :: ivalues(4) ! L1, L2, NPL, N2
    integer :: MAT
    integer :: MF
    integer :: MT
    integer :: NS
    real(8), allocatable :: B(:)
  end type EndfList

!===============================================================================
! TAB1 represents a one-dimensional interpolable function 
!===============================================================================

  type Tab1
    real(8) :: rvalues(2) ! C1, C2
    integer :: ivalues(2) ! L1, L2

    integer :: n_regions = 0       ! # of interpolation regions
    integer, allocatable :: nbt(:) ! values separating interpolation regions
    integer, allocatable :: int(:) ! interpolation scheme
    integer :: n_pairs             ! # of pairs of (x,y) values
    real(8), allocatable :: x(:)   ! values of abscissa
    real(8), allocatable :: y(:)   ! values of ordinate
  end type Tab1

end module endf_header
