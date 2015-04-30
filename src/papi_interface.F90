module papi_interface

  use, intrinsic :: ISO_C_BINDING

  use global, only: master, n_procs, n_particles, n_batches, gen_per_batch
#ifdef MPI
  use mpi
#endif

#ifdef PAPI
  include "f90papi.h"

  character(PAPI_MAX_STR_LEN), allocatable :: event_names(:)
  integer(C_INT), save        :: eventset = PAPI_NULL
  integer(C_LONG_LONG), save  :: values(PAPI_MAX_HWCTRS)

!$omp threadprivate(eventset, values)

contains

!===============================================================================
! PAPI_INITIALIZE
!===============================================================================

  subroutine papi_initialize()

#ifdef _OPENMP
    use omp_lib
#endif

    integer        :: i
    integer(C_INT) :: check

    ! Initialize PAPI library --- this should be outside of a parallel region
    check = PAPI_VER_CURRENT
    call PAPIF_library_init(check)
    if (check /= PAPI_VER_CURRENT) then
      call PAPIF_perror("PAPIF_library_init")
      stop
    end if

!$omp parallel
!$omp critical
#ifdef _OPENMP
    ! Initialize thread support
    call PAPIF_thread_init(omp_get_thread_num, check)
    if (check /= PAPI_OK) then
      call PAPIF_perror("PAPIF_thread_init")
      stop
    end if
#endif

    ! Create event set
    call PAPIF_create_eventset(eventset, check)
    if (check /= PAPI_OK) then
      call PAPIF_perror("PAPIF_create_eventset")
      stop
    end if

    ! Add events to event set
    do i = 1, size(event_names)
      ! Add event
      call PAPIF_add_named_event(eventset, event_names(i), check)
      if (check /= PAPI_OK) then
        call PAPIF_perror("PAPIF_add_named_event (" // trim(event_names(i)) // ")")
        stop
      end if
    end do
!$omp end critical
!$omp end parallel

  end subroutine papi_initialize

!===============================================================================
! PAPI_START_COUNTING
!===============================================================================

  subroutine papi_start_counting()

    integer(C_INT) :: check

!$omp parallel
!$omp critical
    ! Start counting events
    call PAPIF_start(eventset, check)
    if (check /= PAPI_OK) then
      call PAPIF_perror("PAPIF_start")
      stop
    end if
!$omp end critical
!$omp end parallel

  end subroutine papi_start_counting

!===============================================================================
! PAPI_STOP_COUNTING
!===============================================================================

  subroutine papi_stop_counting()

#ifdef _OPENMP
    use omp_lib
#endif
    use global

    integer :: i, j, n
    real(8) :: l3_tcm = 0., l3_tca = 0., flop = 0.
    real(8) :: l2_dca = 0., l2_dcm = 0., l2_ica = 0., l2_icm = 0.
    real(8) :: l1_dcm = 0., l1_icm = 0.
    real(8) :: tlb_im = 0., tlb_dm = 0.
    real(8) :: stl_icy = 0., br_msp = 0.
    real(8) :: l2_pending = 0., l2_stalls = 0.
    real(8) :: vec_sp = 0., vec_dp = 0.
#ifdef MPI
    integer(C_LONG_LONG), allocatable :: values_process(:)
#endif
    real(8) :: seconds
    integer(C_INT) :: n_events
    integer(C_INT) :: check
    integer(C_INT) :: event_code
    integer(C_INT) :: count, flags
    integer(C_LONG_LONG), allocatable :: values_sum(:)
    character(PAPI_MAX_STR_LEN), allocatable :: symbol(:)
    character(PAPI_MAX_STR_LEN), allocatable :: long_descr(:)
    character(PAPI_MAX_STR_LEN), allocatable :: short_descr(:)
    character(PAPI_MAX_STR_LEN), allocatable :: event_note(:)

    n_events = size(event_names)
    allocate(values_sum(n_events), symbol(n_events), long_descr(n_events), &
         short_descr(n_events), event_note(n_events))

    values_sum(:) = 0.

!$omp parallel
    ! Stop counting events
    call PAPIF_stop(eventset, values, check)
    if (check /= PAPI_OK) then
      call PAPIF_perror("PAPIF_stop")
      stop
    end if

!$omp critical
    do i = 1, n_events
      call PAPIF_event_name_to_code(event_names(i), event_code, check)
      if (check /= PAPI_OK) then
        call PAPIF_perror("PAPIF_event_name_to_code")
        stop
      end if

      call PAPIF_get_event_info(event_code, symbol(i), long_descr(i), short_descr(i), &
           count, event_note(i), flags, check)
      if (check /= PAPI_OK) then
        call PAPIF_perror("PAPIF_get_event_info")
        stop
      end if

      ! Accumulate values
      values_sum(i) = values_sum(i) + values(i)
    end do
!$omp end critical
!$omp end parallel

    ! Determine wall-clock time. This is used rather than cycles to determine
    ! bandwidth (latter may be inaccurate due to turboboost)
    seconds = time_inactive%elapsed + time_active%elapsed

#ifdef MPI
    allocate(values_process(n_procs*n_events))
    call MPI_GATHER(values_sum, n_events, MPI_INTEGER8, values_process, n_events, &
         MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_err)
    if (master) then
      values_sum(:) = 0.
      do j = 0, n_procs - 1
        do i = 1, n_events
          values_sum(i) = values_sum(i) + values_process(j*n_events + i)
        end do
      end do
    end if
#endif

    if (master) then
      do i = 1, n_events
        select case (symbol(i))
        case ("PAPI_FP_INS")
          flop = values_sum(i)
        case ("PAPI_L3_TCA")
          l3_tca = values_sum(i)
        case ("PAPI_L3_TCM")
          l3_tcm = values_sum(i)
        case ("PAPI_L2_ICA")
          l2_ica = values_sum(i)
        case ("PAPI_L2_ICM")
          l2_icm = values_sum(i)
        case ("PAPI_L2_DCA")
          l2_dca = values_sum(i)
        case ("PAPI_L2_DCM")
          l2_dcm = values_sum(i)
        case ("PAPI_L1_ICM")
          l1_icm = values_sum(i)
        case ("PAPI_L1_DCM")
          l1_dcm = values_sum(i)
        case ("PAPI_TLB_IM")
          tlb_im = values_sum(i)
        case ("PAPI_TLB_DM")
          tlb_dm = values_sum(i)
        case ("PAPI_STL_ICY")
          stl_icy = values_sum(i)
        case ("CYCLE_ACTIVITY:CYCLES_L2_PENDING")
          l2_pending = values_sum(i)
        case ("CYCLE_ACTIVITY:STALLS_L2_PENDING")
          l2_stalls = values_sum(i)
        case ("PAPI_BR_MSP")
          br_msp = values_sum(i)
        case ("PAPI_VEC_SP")
          vec_sp = values_sum(i)
        case ("PAPI_VEC_DP")
          vec_dp = values_sum(i)
        end select
        print *, values_sum(i), trim(symbol(i))
      end do

      n = n_particles*n_batches*gen_per_batch
      if (flop > 0.) then
        print *, "  MFLOPs            = ", flop/(1e6*seconds)
      end if
      if (l3_tcm > 0.) then
        print *, "  Bandwidth (MiB/s) = ", l3_tcm*64/(1024**2*seconds)
      end if
      if (l3_tca > 0.) then
        print *, "  L3 accesses*      = ", real(l3_tca)/n
        if (l3_tcm > 0.) print *, "  L3 miss rate      = ", real(l3_tcm)/real(l3_tca)
      end if
      if (l2_ica > 0.) then
        print *, "  L2I accesses*     = ", real(l2_ica)/n
        if (l2_icm > 0.) print *, "  L2I miss rate     = ", real(l2_icm)/real(l2_ica)
      end if
      if (l2_dca > 0.) then
        print *, "  L2D accesses*     = ", real(l2_dca)/n
        if (l2_dcm > 0.) print *, "  L2D miss rate     = ", real(l2_dcm)/real(l2_dca)
      end if
      if (l1_icm > 0.) then
        print *, "  L1I misses*     = ", real(l1_icm)/n
      end if
      if (l1_dcm > 0.) then
        print *, "  L1D misses*     = ", real(l1_dcm)/n
      end if
      if (tlb_im > 0.) then
        print *, "  TLBI misses*      = ", real(tlb_im)/n
      end if
      if (tlb_dm > 0.) then
        print *, "  TLBD misses*      = ", real(tlb_dm)/n
      end if
      if (stl_icy > 0.) then
        print *, "  No-issue cycles*  = ", real(stl_icy)/n
      end if
      if (l2_pending > 0.) then
        print *, "  Cycle w/ L2 pend* = ", real(l2_pending)/n
      end if
      if (l2_stalls > 0.) then
        print *, "  Stalls L2 pend*   = ", real(l2_stalls)/n
      end if
      if (br_msp > 0.) then
        print *, "  Mispredict brnch* = ", real(br_msp)/n
      end if
      if (vec_sp > 0.) then
        print *, "  Single vec/SIMD*  = ", real(vec_sp)/n
      end if
      if (vec_dp > 0.) then
        print *, "  Double vec/SIMD*  = ", real(vec_dp)/n
      end if
    end if

#ifdef MPI
    deallocate(values_process)
#endif
    deallocate(values_sum, symbol, long_descr, short_descr, event_note)

  end subroutine papi_stop_counting

#endif

end module papi_interface
