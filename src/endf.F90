module endf

  use constants
  use endf_header
  use string, only: to_str

contains

!===============================================================================
! REACTION_NAME gives the name of the reaction for a given MT value
!===============================================================================

  function reaction_name(MT) result(string)

    integer, intent(in) :: MT
    character(20)       :: string

    select case (MT)
    case (TOTAL_XS)
      string = '(n,total)'
    case (ELASTIC)
      string = '(n,elastic)'
    case (N_LEVEL)
      string = '(n,level)'
    case (N_2ND)
      string = '(n,2nd)'
    case (N_2N)
      string = '(n,2n)'
    case (N_3N)
      string = '(n,3n)'
    case (N_FISSION)
      string = '(n,fission)'
    case (N_F)
      string = '(n,f)'
    case (N_NF)
      string = '(n,nf)'
    case (N_2NF)
      string = '(n,2nf)'
    case (N_NA)
      string = '(n,na)'
    case (N_N3A)
      string = '(n,n3a)'
    case (N_2NA)
      string = '(n,2na)'
    case (N_3NA)
      string = '(n,3na)'
    case (N_NP)
      string = '(n,np)'
    case (N_N2A)
      string = '(n,n2a)'
    case (N_2N2A)
      string = '(n,2n2a)'
    case (N_ND)
      string = '(n,nd)'
    case (N_NT)
      string = '(n,nt)'
    case (N_N3HE)
      string = '(n,nHe-3)'
    case (N_ND2A)
      string = '(n,nd2a)'
    case (N_NT2A)
      string = '(n,nt2a)'
    case (N_4N)
      string = '(n,4n)'
    case (N_3NF)
      string = '(n,3nf)'
    case (N_2NP)
      string = '(n,2np)'
    case (N_3NP)
      string = '(n,3np)'
    case (N_N2P)
      string = '(n,n2p)'
    case (N_NPA)
      string = '(n,npa)'
    case (N_N1 : N_N40)
      string = '(n,n' // trim(to_str(MT-50)) // ')'
    case (N_NC)
      string = '(n,nc)'
    case (N_DISAPPEAR)
      string = '(n,disappear)'
    case (N_GAMMA)
      string = '(n,gamma)'
    case (N_P)
      string = '(n,p)'
    case (N_D)
      string = '(n,d)'
    case (N_T)
      string = '(n,t)'
    case (N_3HE)
      string = '(n,3He)'
    case (N_A)
      string = '(n,a)'
    case (N_2A)
      string = '(n,2a)'
    case (N_3A)
      string = '(n,3a)'
    case (N_2P)
      string = '(n,2p)'
    case (N_PA)
      string = '(n,pa)'
    case (N_T2A)
      string = '(n,t2a)'
    case (N_D2A)
      string = '(n,d2a)'
    case (N_PD)
      string = '(n,pd)'
    case (N_PT)
      string = '(n,pt)'
    case (N_DA)
      string = '(n,da)'
    case (201)
      string = '(n,Xn)'
    case (202)
      string = '(n,Xgamma)'
    case (203)
      string = '(n,Xp)'
    case (204)
      string = '(n,Xd)'
    case (205)
      string = '(n,Xt)'
    case (206)
      string = '(n,X3He)'
    case (207)
      string = '(n,Xa)'
    case (444)
      string = '(damage)'
    case (600 : 648)
      string = '(n,p' // trim(to_str(MT-600)) // ')'
    case (649)
      string = '(n,pc)'
    case (650 : 698)
      string = '(n,d' // trim(to_str(MT-650)) // ')'
    case (699)
      string = '(n,dc)'
    case (700 : 748)
      string = '(n,t' // trim(to_str(MT-700)) // ')'
    case (749)
      string = '(n,tc)'
    case (750 : 798)
      string = '(n,3He' // trim(to_str(MT-750)) // ')'
    case (799)
      string = '(n,3Hec)'
    case (800 : 848)
      string = '(n,a' // trim(to_str(MT-800)) // ')'
    case (849)
      string = '(n,tc)'
    case default
      string = 'MT=' // trim(to_str(MT))
    end select

  end function reaction_name

!===============================================================================
! READ_CONT_RECORD
!===============================================================================

  subroutine read_cont_record(fh, cont)

    integer,        intent(in)    :: fh
    type(EndfCont), intent(inout) :: cont

    read(fh,'(2G11.0,4I11,I4,I2,I3,I5)') cont % rvalues, cont % ivalues, &
         cont % MAT, cont % MF, cont % MT, cont % NS

  end subroutine read_cont_record

!===============================================================================
! READ_LIST_RECORD
!===============================================================================

  subroutine read_list_record(fh, list)

    integer, intent(in) :: fh
    type(EndfList), intent(inout) :: list

    integer :: i ! implied do index
    integer :: n ! length of list

    ! Read first line and determine list length
    read(fh,'(2G11.0,4I11,I4,I2,I3,I5)') list % rvalues, list % ivalues, &
         list % MAT, list % MF, list % MT, list % NS
    n = list % ivalues(3)

    ! Allocate array for list
    if (allocated(list%B)) deallocate(list%B)
    allocate(list%B(n))

    ! Read list
    read(fh,'(6G11.0)') (list%B(i), i=1,n)

  end subroutine read_list_record

!===============================================================================
! READ_TAB1_RECORD
!===============================================================================

  subroutine read_tab1_record(fh, t)

    integer,    intent(in)    :: fh
    type(Tab1), intent(inout) :: t

    integer :: i ! implied do loop index

    ! read first list with information on number of regions and pairs
    read(fh,'(2G11.0,4I11,I4,I2,I3,I5)') t % rvalues, t % ivalues, &
         t % n_regions, t % n_pairs

    ! deallocate any allocatable components
    if (allocated(t%nbt)) deallocate(t%nbt)
    if (allocated(t%int)) deallocate(t%int)
    if (allocated(t%x)) deallocate(t%x)
    if (allocated(t%y)) deallocate(t%y)

    ! allocate components
    allocate(t%nbt(t%n_regions))
    allocate(t%int(t%n_regions))
    allocate(t%x(t%n_pairs))
    allocate(t%y(t%n_pairs))

    ! read interpolation schemes, then (x,y) pairs
    read(fh,'(6I11)') (t%nbt(i), t%int(i), i=1,t%n_regions)
    read(fh,'(6G11.0)') (t%x(i), t%y(i), i=1,t%n_pairs)

  end subroutine read_tab1_record

!===============================================================================
! IS_FISSION determines if a given MT number is that of a fission event. This
! accounts for aggregate fission (MT=18) as well as partial fission reactions.
!===============================================================================

  function is_fission(MT) result(fission_event)

    integer, intent(in) :: MT
    logical             :: fission_event

    if (MT == N_FISSION .or. MT == N_F .or. MT == N_NF .or. MT == N_2NF & 
         .or. MT == N_3NF) then
      fission_event = .true.
    else
      fission_event = .false.
    end if

  end function is_fission

!===============================================================================
! IS_SCATTER determines if a given MT number is that of a scattering event
!===============================================================================

  function is_scatter(MT) result(scatter_event)

    integer, intent(in) :: MT
    logical             :: scatter_event

    if (MT < 100) then
      if (MT == N_FISSION .or. MT == N_F .or. MT == N_NF .or. MT == N_2NF & 
           .or. MT == N_3NF) then
        scatter_event = .false.
      else
        scatter_event = .true.
      end if
    else
      scatter_event = .false.
    end if

  end function is_scatter

end module endf
