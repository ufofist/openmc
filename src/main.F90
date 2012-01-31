program main

  use constants
  use global
  use finalize,        only: finalize_run
  use initialize,      only: initialize_run
  use intercycle,      only: shannon_entropy, calculate_keff, synchronize_bank
  use output,          only: write_message, header
  use particle_header, only: Particle
  use plot,            only: run_plot
  use physics,         only: transport
  use random_lcg,      only: set_particle_seed
  use source,          only: get_source_particle
  use string,          only: to_str
  use tally,           only: synchronize_tallies 
  use timing,          only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

  implicit none

  ! set up problem
  call initialize_run()

  ! start problem
  if (plotting) then
     call run_plot()
  else
     call run_problem()
  end if

  call calculate_leakage()

  ! finalize run
  call finalize_run()
  
contains

!===============================================================================
! RUN_PROBLEM encompasses all the main logic where iterations are performed over
! the cycles and histories.
!===============================================================================

  subroutine run_problem()

    integer                 :: i_cycle     ! cycle index
    integer                 :: final_stage
    integer(8)              :: i_particle  ! history index
    type(Particle), pointer :: p => null()

    if (master) call header("BEGIN SIMULATION", level=1)

    tallies_on = .false.
    call timer_start(time_inactive)

    ! Display column titles
    if (entropy_on) then
       message = " Cycle   k(cycle)   Entropy         Average k"
       call write_message(1)
       message = " =====   ========   =======    ==================="
       call write_message(1)
    else
       message = " Cycle   k(cycle)          Average k"
       call write_message(1)
       message = " =====   ========     ==================="
       call write_message(1)
    end if

    ! ==========================================================================
    ! LOOP OVER CYCLES
    CYCLE_LOOP: do i_cycle = 1, n_cycles

       ! Start timer for computation
       call timer_start(time_compute)

       message = "Simulating cycle " // trim(to_str(i_cycle)) // "..."
       call write_message(8)

       ! Set all tallies to zero
       n_bank = 0

       ! =======================================================================
       ! LOOP OVER HISTORIES
       HISTORY_LOOP: do

          ! grab source particle from bank
          p => get_source_particle()
          if ( .not. associated(p) ) then
             ! no particles left in source bank
             exit HISTORY_LOOP
          end if

          ! set random number seed
          i_particle = (i_cycle-1)*n_particles + p % id
          call set_particle_seed(i_particle)
          
          ! set particle trace
          trace = .false.
          if (i_cycle == trace_cycle .and. &
               p % id == trace_particle) trace = .true.

          ! transport particle
          call transport(p)

       end do HISTORY_LOOP

       ! Accumulate time for computation
       call timer_stop(time_compute)

       ! =======================================================================
       ! WRAP UP FISSION BANK AND COMPUTE TALLIES, KEFF, ETC

       ! Start timer for inter-cycle synchronization
       call timer_start(time_intercycle)

       ! Collect tallies
       if (tallies_on) then
          call timer_start(time_ic_tallies)
          call synchronize_tallies()
          call timer_stop(time_ic_tallies)
       end if

       ! Calculate shannon entropy
       if (entropy_on) call shannon_entropy()

       ! Distribute fission bank across processors evenly
       call synchronize_bank(i_cycle)

       ! Collect results and statistics
       call calculate_keff(i_cycle)

       ! print cycle information

       if (tallies_on) then
#ifdef MPI
          call MPI_REDUCE(last_stage, final_stage, 1, MPI_INTEGER, MPI_MAX, &
               0, MPI_COMM_WORLD, mpi_err)
#else
          final_stage = last_stage
#endif
          final_stage_count(final_stage) = final_stage_count(final_stage) + 1
       end if

       ! Turn tallies on once inactive cycles are complete
       if (i_cycle == n_inactive) then
          tallies_on = .true.
          call timer_stop(time_inactive)
          call timer_start(time_active)
       end if

       ! Stop timer for inter-cycle synchronization
       call timer_stop(time_intercycle)

    end do CYCLE_LOOP

    call timer_stop(time_active)

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_problem

!===============================================================================
! CALCULATE_LEAKAGE
!===============================================================================

  subroutine calculate_leakage()

    integer :: i, j, k, m
#ifdef MPI
    integer :: n
#endif
    real(8) :: src
    real(8) :: leak
    real(8) :: leak_fraction

#ifdef MPI
    ! If running in parallel, we first need to combine the values from all
    ! processors
    n = lmesh_nx * lmesh_ny * lmesh_nz * MAX_STAGES
    if (master) then
       call MPI_REDUCE(MPI_IN_PLACE, starting_source, n, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
       call MPI_REDUCE(MPI_IN_PLACE, leakage, n, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
       call MPI_REDUCE(MPI_IN_PLACE, last_stage, 1, MPI_INTEGER, MPI_MAX, &
            0, MPI_COMM_WORLD, mpi_err)
    else
       call MPI_REDUCE(starting_source, starting_source, n, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
       call MPI_REDUCE(leakage, leakage, n, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
       call MPI_REDUCE(last_stage, last_stage, 1, MPI_INTEGER, MPI_MAX, &
            0, MPI_COMM_WORLD, mpi_err)
    end if
#endif

    if (master) then
       ! open files for writing
       open(UNIT=50, FILE='source.out', STATUS='replace', ACTION='write')
       open(UNIT=51, FILE='leakage.out', STATUS='replace', ACTION='write')
       open(UNIT=52, FILE='leakage_fraction.out', STATUS='replace', ACTION='write')
       open(UNIT=53, FILE='stage_counts.out', STATUS='replace', ACTION='write')

       ! Write dimension of mesh on each file
       write(50,*) lmesh_nx, lmesh_ny, lmesh_nz, last_stage
       write(51,*) lmesh_nx, lmesh_ny, lmesh_nz, last_stage
       write(52,*) lmesh_nx, lmesh_ny, lmesh_nz, last_stage

       ! Loop over all mesh cells
       do m = 1, last_stage
          do k = 1, lmesh_nz
             do j = 1, lmesh_ny
                do i = 1, lmesh_nx
                   ! calculate fraction of source sites in each mesh
                   src = starting_source(i,j,k,m)/(n_particles*(n_cycles - n_inactive))

                   ! calculate absolute leakage
                   leak = leakage(i,j,k,m)/(n_particles*(n_cycles - n_inactive))

                   ! calculate leakage fraction
                   if (starting_source(i,j,k,m) > ZERO) then
                      leak_fraction = leakage(i,j,k,m)/starting_source(i,j,k,m)
                   else
                      leak_fraction = ZERO
                   end if

                   ! write values to each file
                   write(50,*) i,j,k,m,src
                   write(51,*) i,j,k,m,leak
                   write(52,*) i,j,k,m,leak_fraction
                end do
             end do
          end do
       end do

       ! Write stage counts
       do m = 1, MAX_STAGES
          write(53,*) m, final_stage_count(m)
       end do

       ! close files
       close(UNIT=50)
       close(UNIT=51)
       close(UNIT=52)
       close(UNIT=53)
    end if

  end subroutine calculate_leakage

end program main
