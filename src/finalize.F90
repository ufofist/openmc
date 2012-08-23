module finalize

  use global
  use output,         only: print_runtime
  use source,         only: write_source_binary
  use tally,          only: write_tallies, tally_statistics, statistics_score
  use timing,         only: timer_start, timer_stop

#ifdef MPI
  use state_point,    only: server_create_state_point
#endif

#ifdef HDF5
  use hdf5_interface, only: hdf5_write_results, hdf5_close_output
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics, writing out tallies, and writing hdf5 output.
!===============================================================================

  subroutine finalize_run()

    ! Start finalization timer
    call timer_start(time_finalize)

    if (run_mode /= MODE_PLOTTING) then
       if (use_servers) then
          if (server) then
             if (compute_rank == 0) then
                ! Get k_batch and entropy from master processor
                call MPI_RECV(k_batch, n_batches, MPI_REAL8, 0, &
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                call MPI_RECV(entropy, n_batches, MPI_REAL8, 0, &
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                call MPI_RECV(global_tallies, n_global_tallies, &
                     MPI_TALLYSCORE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                     mpi_err)
             end if
             
             call server_create_state_point()
             call statistics_score(server_scores)
          elseif (master) then
             ! Send k_batch, entropy, and global tallies to first server
             call MPI_SEND(k_batch, n_batches, MPI_REAL8, support_ratio - 1, &
                  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
             call MPI_SEND(entropy, n_batches, MPI_REAL8, support_ratio - 1, &
                  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
             call MPI_SEND(global_tallies, n_global_tallies, MPI_TALLYSCORE, &
                  support_ratio - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                  mpi_err)

             call statistics_score(global_tallies)
          end if
       else
          ! Calculate statistics for tallies and write to tallies.out
          call tally_statistics()
          if (master) call write_tallies()

          ! Write out binary source
          if (write_source) call write_source_binary()
       end if
    end if

    ! stop timers and show timing statistics
    call timer_stop(time_finalize)
    call timer_stop(time_total)
    if (master .and. (run_mode /= MODE_PLOTTING)) call print_runtime()

#ifdef HDF5
    ! Write time statistics to HDF5 output 
    if (master) then
       call hdf5_write_results()
       call hdf5_close_output()
    end if
#endif

    ! deallocate arrays
    call free_memory()
    
#ifdef MPI
    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(mpi_err)
#endif

  end subroutine finalize_run

end module finalize
