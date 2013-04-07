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
    
    ! Type-Bound procedures
    contains
      procedure :: clear => Tab1_clear ! deallocates a Tab1 Object.
  end type Tab1
  
  contains
  
!===============================================================================
! TAB1_CLEAR deallocates the items in Tab1
!===============================================================================

    subroutine tab1_clear(this)
      
      class(Tab1), intent(inout) :: this ! The Tab1 to clear
      
      if (allocated(this % nbt)) &
           deallocate(this % nbt, this % int)
        
      if (allocated(this % x)) &
           deallocate(this % x, this % y)
        
    end subroutine tab1_clear

!===============================================================================
! ENDFDECAYMODE
!===============================================================================

  type EndfDecayMode
    integer :: type
    integer :: state
    real(8) :: Q_value
    real(8) :: branching_ratio
  end type EndfDecayMode

!===============================================================================
! ENDFDECAY
!===============================================================================

  type EndfDecay
    integer :: zzaaam
    real(8) :: awr
    real(8) :: lambda
    real(8) :: energy
    integer :: n_modes
    type(EndfDecayMode), allocatable :: modes(:)
  end type EndfDecay

!===============================================================================
! FISSIONPRODUCT
!===============================================================================

  type FissionProduct
    integer :: zzaaam
    real(8) :: yield
  end type FissionProduct

!===============================================================================
! ENDFYIELD
!===============================================================================

  type EndfYield
    real(8) :: energy
    integer :: n_products
    type(FissionProduct), allocatable :: products(:)
  end type EndfYield

!===============================================================================
! ENDFFISSIONYIELD
!===============================================================================

  type EndfFissionYield
    integer :: zzaaam
    integer :: n_energy
    type(EndfYield), allocatable :: yield(:)
  end type EndfFissionYield

end module endf_header
