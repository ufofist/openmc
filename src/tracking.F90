module tracking

  use constants
  use cross_section,      only: calculate_xs
  use error,              only: fatal_error, warning
  use geometry_header,    only: cells
  use geometry,           only: find_cell, distance_to_boundary, cross_surface, &
                                cross_lattice, check_cell_overlap
  use output,             only: write_message
  use material_header,    only: materials
  use message_passing
  use mgxs_header
  use nuclide_header
  use particle_header,    only: LocalCoord, Particle
  use physics,            only: collision
  use physics_mg,         only: collision_mg
  use random_lcg,         only: prn
  use settings
  use simulation_header
  use string,             only: to_str
  use tally_header
  use tally,              only: score_analog_tally, score_tracklength_tally, &
                                score_collision_tally, score_surface_current, &
                                score_track_derivative, score_surface_tally, &
                                score_collision_derivative, zero_flux_derivs
  use track_output,       only: initialize_particle_track, write_particle_track, &
                                add_particle_track, finalize_particle_track

  implicit none

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport(p)

    type(Particle), intent(inout) :: p

    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! sampled distance to collision
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: tau_hat                ! Optical depth used in newton's method
    real(8) :: optical_depth          ! Optical depth integral
    real(8), allocatable :: xs_t(:)   ! The total cross section along flight path
    real(8) :: PNC                    ! Probilitiy of no collision
    real(8) :: xyz_orig(3)            ! Original xyz coordinate of the
                                      ! of the particle

    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      call write_message("Simulating Particle " // trim(to_str(p % id)))
    end if

    ! Initialize number of events to zero
    n_event = 0

    ! Add paricle's starting weight to count for normalizing tallies later
!$omp atomic
    total_weight = total_weight + p % wgt

    ! Force calculation of cross-sections by setting last energy to zero
    if (run_CE) then
      micro_xs % last_E = ZERO
    end if

    ! Prepare to write out particle track.
    if (p % write_track) then
      call initialize_particle_track()
    endif

    ! Every particle starts with no accumulated flux derivative.
    if (active_tallies % size() > 0) call zero_flux_derivs()

    EVENT_LOOP: do
      ! Store pre-collision particle properties
      p % last_wgt = p % wgt
      p % last_E   = p % E
      p % last_uvw = p % coord(1) % uvw
      p % last_xyz = p % coord(1) % xyz

      ! If the cell hasn't been determined based on the particle's location,
      ! initiate a search for the current cell. This generally happens at the
      ! beginning of the history and again for any secondary particles
      if (p % coord(p % n_coord) % cell == NONE) then
        call find_cell(p, found_cell)
        if (.not. found_cell) then
          call fatal_error("Could not locate particle " // trim(to_str(p % id)))
        end if

        ! set birth cell attribute
        if (p % cell_born == NONE) p % cell_born = p % coord(p % n_coord) % cell
      end if

      ! Write particle track.
      if (p % write_track) call write_particle_track(p)

      if (check_overlaps) call check_cell_overlap(p)

      ! Calculate microscopic and macroscopic cross sections
      if (run_CE) then
        ! If the material is the same as the last material and the temperature
        ! hasn't changed, we don't need to lookup cross sections again.
        if (p % material /= p % last_material .or. &
             p % sqrtkT /= p % last_sqrtkT) call calculate_xs(p)
      else
        ! Since the MGXS can be angle dependent, this needs to be done
        ! After every collision for the MGXS mode
        if (p % material /= MATERIAL_VOID) then
          ! Update the temperature index
          call macro_xs(p % material) % obj % find_temperature(p % sqrtkT)
          ! Get the data
          call macro_xs(p % material) % obj % calculate_xs(p % g, &
               p % coord(p % n_coord) % uvw, material_xs)
        else
          material_xs % total      = ZERO
          material_xs % absorption = ZERO
          material_xs % nu_fission = ZERO
        end if

        ! Finally, update the particle group while we have already checked for
        ! if multi-group
        p % last_g = p % g
      end if

      ! Find the distance to the nearest boundary
      call distance_to_boundary(p, d_boundary, surface_crossed, &
           lattice_translation, next_level)

      ! Sample a distance to collision
      if (material_xs % total == ZERO) then
        d_collision = INFINITY
      else if (materials(p % material) % continuous_num_density) then
        ! Continuous material tracking
        if (.not. allocated(xs_t)) allocate(xs_t(num_intervals+1))

        ! Integrate along complete neutron flight path
        call simpsons_path_integration(p, optical_depth, d_boundary, xs_t, .false., 0)

        ! Sample the collision for analytic density
        PNC = exp(-optical_depth)

        if (prn() <= PNC) then
          ! No Collision
          d_collision = INFINITY
        else
          ! Collision

          ! Sample optical depth
          tau_hat = -log(ONE - (ONE - PNC) * prn())

          ! Get the flight distance for sampled optical depth
          call estimate_flight_distance(xs_t, d_boundary, tau_hat, d_collision)

          ! Move particle to the point of collision so we can make sure that
          ! cross sections are updated for later tallying of keff and user tallies
          xyz_orig = p % coord(p % n_coord) % xyz
          call move_particle_coord(p % coord(p % n_coord), d_collision)
          call calculate_xs(p)

          ! Move particle back
          p % coord(p % n_coord) % xyz = xyz_orig
        end if
      else
        d_collision = -log(prn()) / material_xs % total
      end if

      ! Select smaller of the two distances
      distance = min(d_boundary, d_collision)

      ! Advance particle
      do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
      end do

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) then
        call score_tracklength_tally(p, distance)
      end if

      ! Score track-length estimate of k-eff
      if (run_mode == MODE_EIGENVALUE) then
        global_tally_tracklength = global_tally_tracklength + p % wgt * &
             distance * material_xs % nu_fission
      end if

      ! Score flux derivative accumulators for differential tallies.
      if (active_tallies % size() > 0) call score_track_derivative(p, distance)

      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        if (next_level > 0) p % n_coord = next_level

        ! Saving previous cell data
        do j = 1, p % n_coord
          p % last_cell(j) = p % coord(j) % cell
        end do
        p % last_n_coord = p % n_coord

        p % coord(p % n_coord) % cell = NONE
        if (any(lattice_translation /= 0)) then
          ! Particle crosses lattice boundary
          p % surface = NONE
          call cross_lattice(p, lattice_translation)
          p % event = EVENT_LATTICE
        else
          ! Particle crosses surface
          p % surface = surface_crossed

          call cross_surface(p)
          p % event = EVENT_SURFACE
        end if
        ! Score cell to cell partial currents
        if(active_surface_tallies % size() > 0) call score_surface_tally(p)
      else
        ! ====================================================================
        ! PARTICLE HAS COLLISION

        ! Score collision estimate of keff
        if (run_mode == MODE_EIGENVALUE) then
          global_tally_collision = global_tally_collision + p % wgt * &
               material_xs % nu_fission / material_xs % total
        end if

        ! score surface current tallies -- this has to be done before the collision
        ! since the direction of the particle will change and we need to use the
        ! pre-collision direction to figure out what mesh surfaces were crossed

        if (active_current_tallies % size() > 0) call score_surface_current(p)

        ! Clear surface component
        p % surface = NONE

        if (run_CE) then
          call collision(p)
        else
          call collision_mg(p)
        end if

        ! Score collision estimator tallies -- this is done after a collision
        ! has occurred rather than before because we need information on the
        ! outgoing energy for any tallies with an outgoing energy filter
        if (active_collision_tallies % size() > 0) call score_collision_tally(p)
        if (active_analog_tallies % size() > 0) call score_analog_tally(p)

        ! Reset banked weight during collision
        p % n_bank   = 0
        p % wgt_bank = ZERO
        p % n_delayed_bank(:) = 0

        ! Reset fission logical
        p % fission = .false.

        ! Save coordinates for tallying purposes
        p % last_xyz_current = p % coord(1) % xyz

        ! Set last material to none since cross sections will need to be
        ! re-evaluated
        p % last_material = NONE

        ! Set all uvws to base level -- right now, after a collision, only the
        ! base level uvws are changed
        do j = 1, p % n_coord - 1
          if (p % coord(j + 1) % rotated) then
            ! If next level is rotated, apply rotation matrix
            p % coord(j + 1) % uvw = matmul(cells(p % coord(j) % cell) % &
                 rotation_matrix, p % coord(j) % uvw)
          else
            ! Otherwise, copy this level's direction
            p % coord(j + 1) % uvw = p % coord(j) % uvw
          end if
        end do

        ! Score flux derivative accumulators for differential tallies.
        if (active_tallies % size() > 0) call score_collision_derivative(p)
      end if

      ! If particle has too many events, display warning and kill it
      n_event = n_event + 1
      if (n_event == MAX_EVENTS) then
        if (master) call warning("Particle " // trim(to_str(p%id)) &
             &// " underwent maximum number of events.")
        p % alive = .false.
      end if

      ! Check for secondary particles if this particle is dead
      if (.not. p % alive) then
        if (p % n_secondary > 0) then
          call p % initialize_from_source(p % secondary_bank(p % n_secondary), &
                                          run_CE, energy_bin_avg)
          p % n_secondary = p % n_secondary - 1
          n_event = 0

          ! Enter new particle in particle track file
          if (p % write_track) call add_particle_track()
        else
          exit EVENT_LOOP
        end if
      end if
    end do EVENT_LOOP

    ! Finish particle track output.
    if (p % write_track) then
      call write_particle_track(p)
      call finalize_particle_track(p)
    endif

  end subroutine transport

!===============================================================================
! MOVE_PARTICLE_COORD moves particle along integration path
!===============================================================================

  subroutine move_particle_coord(coord,ds)

    type(LocalCoord), intent(inout) :: coord
    real(8), intent(in) :: ds          ! Path length to move particle

    ! Move particle along path
    coord % xyz(1) = coord % xyz(1) + ds*coord % uvw(1);
    coord % xyz(2) = coord % xyz(2) + ds*coord % uvw(2);
    coord % xyz(3) = coord % xyz(3) + ds*coord % uvw(3);

  end subroutine move_particle_coord

!===============================================================================
! SIMPSONS_PATH_INTEGRATION uses Simpson's rule to integrate the total cross
! section along the particle path
!===============================================================================

  subroutine simpsons_path_integration(p, optical_depth, distance,  &
                                      xs_t, dbg_file, it_num)

    type(Particle), intent(inout) :: p      ! Particle of interest
    real(8), intent(inout) :: optical_depth ! Total optical depth integration
    real(8), intent(in) :: distance         ! Total distance to integrate
    logical, intent(in) :: dbg_file         ! Debug flag to write to file
    real(8) :: ds                           ! Differential path length
    real(8) :: new_temp                     ! Temperature at current position
    character(len=90) :: format             ! Format of the dbg output
    character(len=90) :: fname              ! Filename to print total cross
    integer :: i                            ! Iterator for intervals
    integer :: it_num                       ! This is needed only for output
    type(LocalCoord) :: coord               ! Coordinate level of
                                            ! continuous transport cell
                                            ! path
    real(8), intent(inout) :: xs_t(:)

    ! We need to get a the lowest coordinate level
    ! This is an assumption of the method.
    coord = p % coord(p % n_coord)

    ! Calculate differential path length
    ds = distance / real(num_intervals,8)

    if (dbg_file) then
       write(fname, '( "particle_",I0.5,"_fc_",I0.5,"_it_",I0.5,".out" )' ) &
            p%id, dbg_file_counter, it_num
       write(*,*) 'Opening file', fname
       open(unit = 1001, file=fname)
       format = '( E21.10, ",",  E21.10, ",", E21.10, ",", E21.10, ",", E21.10, ",", E21.10, ",", E21.10)'
       !            x             y          z             temp          xs_t
       ! or first row u           v          w           energy           ds
       write(1001, format) p % coord(p % n_coord) % uvw(1), p % coord(p % n_coord) % uvw(2), &
            p % coord(p % n_coord) % uvw(3), p % E, ds
    endif

    do i=1,(num_intervals+1)
       ! Recalculate the cross section
       call calculate_xs(p)
       ! Save the total cross section
       xs_t(i) = material_xs % total

       if(dbg_file) then
          write(*,*) 'String being written to file ', fname, ':'
          write(*,*) p % coord(p % n_coord) % xyz(1),  p % coord(p % n_coord) % xyz(2), &
               p % coord(p % n_coord) % xyz(3), new_temp, xs_t(i)
          write(1001, format) p % coord(p % n_coord) % xyz(1), p % coord(p % n_coord) % xyz(2), &
               p % coord(p % n_coord) % xyz(3), new_temp, xs_t(i)
       endif
       ! Move particle along path
       call move_particle_coord(p % coord(p % n_coord),ds)
    enddo

    p % coord(p % n_coord) = coord

    optical_depth = 0.0

    do i=1,(num_intervals-1),2
       ! Add to the integral
       optical_depth = optical_depth + &
            2.0_8 * ds / 6.0_8 * ( xs_t(i) + 4.0_8 * xs_t(i+1) + xs_t(i+2))
    enddo

    if(dbg_file) close(1001)

  end subroutine simpsons_path_integration

  subroutine estimate_flight_distance(xs_t, distance, tau_hat, s)

    real(8), intent(inout) :: s              ! The estimated flight distance
    real(8), intent(in) :: distance          ! Total distance to integrate
    real(8), intent(in)  :: xs_t(:)          ! Total cross sections at each point on path
    real(8), intent(in)  :: tau_hat          ! The sampled optical depth
    real(8) :: ds                            ! Spacing between points.  Assumed to be constant.
    real(8) :: delta_tau_hat                 ! The delta tau hat that we are solving for in the bin
    real(8) :: optical_depth
    integer :: index
    logical :: tau_overrun                   ! Indicates if the search for the correct bin failed
    real(8) :: a,b,c,d                       ! Polynomial expansion values
    real(8) :: m                             ! Slope of a linear guess for quadratic polynomial
    real(8) :: dds                           ! The differential lenght in the particular pin

    ! Calculate differential path length
    ds = distance / real(num_intervals,8)

    ! Set variables for loop
    optical_depth = 0.0
    index = 1
    tau_overrun = .FALSE.

    ! Find the index of the cell that goes over the sampled path integration
    do while(optical_depth <= tau_hat)
       if (index > (num_intervals - 1)) then
          tau_overrun = .TRUE.
          optical_depth = tau_hat
       else
          optical_depth = optical_depth + &
               2.0_8 * ds / 6.0_8 * ( xs_t(index) + &
               4.0_8 * xs_t(index+1) + xs_t(index+2))
          index = index + 2
       endif
    enddo

    ! Subtract off the last delta optical depth that we added
    index = index - 2
    optical_depth = optical_depth - &
         2.0_8 * ds / 6.0_8 * ( xs_t(index) + &
         4.0_8 * xs_t(index+1) + xs_t(index+2))
    ! Calculate the delta in the optical depth
    delta_tau_hat = tau_hat - optical_depth

    if (tau_overrun) then
       write(*,*) 'WARNING: The search for the optical depth bin failed'
       write(*,*) 'tau_hat = ', tau_hat
       write(*,*) 'optical_depth = ', optical_depth
       s = distance
    else

       ! Determine the coefficients for a second order expansion in that bin
       ! sigma_t = ax^2 + bx + c
       ! Note that we are shifting the portion of the curve so that the first
       ! point lies at zero.  This makes integration much easier.
       c = xs_t(index)
       b = ( 2.0_8*xs_t(index+1) - 0.5_8*xs_t(index+2) - 1.5_8*c ) / (ds)
       a = ( xs_t(index+2) - c - 2.0_8*ds*b ) / (4.0_8*ds*ds)
       ! Define quantities for solution of the cubic equation
       b = b / 2.0_8
       a = a / 3.0_8
       d = -delta_tau_hat
       ! If the polynomial is not truly cubic, the soluiton will blow up
       ! Choose between integrated polynomial of cubic,quadratic, and linear
       if ( abs(a) > 1E-10) then
          call get_cubic_root(a,b,c,d,0.0_8,2.0_8*ds,dds)
          s = ds * (index-1) + dds

          !write(*,*) 'a = ', a
          !write(*,*) 'b = ', b
          !write(*,*) 'c = ', c
          !write(*,*) 'd = ', d
          !write(*,*) 'ds = ', ds
          !write(*,*) 'index = ', index
          !write(*,*) 'xs_t(1) = ', xs_t(index)
          !write(*,*) 'xs_t(2) = ', xs_t(index+1)
          !write(*,*) 'xs_t(3) = ', xs_t(index+2)
          !write(*,*) 's = ', s
          !write(*,*) 'distance = ', distance
          !write(*,*) 'optical_depth = ', optical_depth
          !write(*,*) 'tau_hat = ', tau_hat
          !write(*,*) 'delta_tau_hat = ', delta_tau_hat
          !write(*,*) 'dds = ', dds
          !call fatal_error('WARNING: non-real roots detected for cubic equation')
       else if ( abs(b) > 1E-10 ) then
          ! Get the best guest for root based on a linear fit
          m = ( xs_t(index+2) - xs_t(index) ) / (2.0_8 * ds)
          call get_quadratic_root(0.5*m, c, d, 0.0_8, 2.0*ds, dds)
          s = ds * (index-1) + dds
          !write(*,*) 'a = ', 0.5*m
          !write(*,*) 'b = ', c
          !write(*,*) 'c = ', d
          !write(*,*) 'ds = ', ds
          !write(*,*) 'index = ', index
          !write(*,*) 'xs_t(1) = ', xs_t(index)
          !write(*,*) 'xs_t(2) = ', xs_t(index+1)
          !write(*,*) 'xs_t(3) = ', xs_t(index+2)
          !write(*,*) 's = ', s
          !write(*,*) 'distance = ', distance
          !write(*,*) 'optical_depth = ', optical_depth
          !write(*,*) 'tau_hat = ', tau_hat
          !write(*,*) 'delta_tau_hat = ', delta_tau_hat
          !write(*,*) 'dds = ', dds
       else
          ! This is the case where the cross section is constant.
          s = ds * (index-1) - d/c
          ! Put some checks in for now
          if (s < (ds * (index-1)) .OR. s > (ds * (index+1)) ) then
             write(*,*) 'a = ', a
             write(*,*) 'b = ', b
             write(*,*) 'c = ', c
             write(*,*) 'ds = ', ds
             write(*,*) 'xs_t(1) = ', xs_t(index)
             write(*,*) 'xs_t(2) = ', xs_t(index+1)
             write(*,*) 'xs_t(3) = ', xs_t(index+2)
             write(*,*) 's = ', s
             write(*,*) 'distance = ', distance
             write(*,*) 'optical_depth = ', optical_depth
             write(*,*) 'tau_hat = ', tau_hat
             write(*,*) 'delta_tau_hat = ', delta_tau_hat
             write(*,*) 'index = ' , index
             write(*,*) 'ds*(index-1) = ', ds * (index-1)
             write(*,*) 'ds*(index) = ', ds * (index)
             call fatal_error('s for the constant cross section case lies out of bounds')
          endif
       endif
    endif
    if (s > distance .OR. s < 0.0) then
       write(*,*) 'a = ', a
       write(*,*) 'b = ', b
       write(*,*) 'c = ', c
       write(*,*) 'd = ', d
       write(*,*) 'ds = ', ds
       write(*,*) 'index = ', index
       write(*,*) 'ds*(index-1)', ds*(index-1)
       write(*,*) 'xs_t(1) = ', xs_t(index)
       write(*,*) 'xs_t(2) = ', xs_t(index+1)
       write(*,*) 'xs_t(3) = ', xs_t(index+2)
       write(*,*) 's = ', s
       write(*,*) 'distance = ', distance
       write(*,*) 'optical_depth = ', optical_depth
       write(*,*) 'tau_hat = ', tau_hat
       write(*,*) 'delta_tau_hat = ', delta_tau_hat
       write(*,*) 'dds = ', dds
       call fatal_error('The predicted distance is greater than distance to boundary or less than zero')
    endif

  end subroutine estimate_flight_distance

  subroutine get_quadratic_root(a, b, c, lower_b, upper_b, root)

    ! This function returns the first quadratic root in the interval
    ! [lower_b, upper_b]

    real(8), intent(in) :: a,b,c    ! Coefficients in ax^2 + bx + c
    real(8), intent(in) :: lower_b  ! Lower bound of search
    real(8), intent(in) :: upper_b  ! Upper bound of search
    real(8), intent(inout) :: root  ! The root of the equation
    real(8) :: q                    ! Variable for numeric solution
    real(8) :: C_1, C_2             ! The two roots
    real(8) :: inner_sqrt           ! Temporary variabel to check if real
    real(8) :: e1, e2, e3, e4       ! Errors when root goes over bound
    real(8) :: min_error            ! Min error used to determie which root
                                    ! should be assigned

    inner_sqrt = b*b - 4.0_8*a*c

    if( inner_sqrt < 0.0) then
       write(*,*) 'a = ', a
       write(*,*) 'b = ', b
       write(*,*) 'c = ', c
       write(*,*) 'inner_sqrt = ', inner_sqrt
       call fatal_error('Non real roots for quadratic equation')
    else if (inner_sqrt >= 0) then

       q = -1.0_8 / 2.0_8 * ( b + sign(1.0_8,b) * sqrt(inner_sqrt) )
       C_1 = q / a
       C_2 = c / q

       if (C_1 >= lower_b .AND. C_1 <= upper_b) then
          root = C_1
       else if (C_2 >= lower_b .AND. C_2 <= upper_b) then
          root = C_2
       else
          e1 = abs(C_1 - lower_b)
          e2 = abs(C_1 - upper_b)
          e3 = abs(C_2 - lower_b)
          e4 = abs(C_2 - upper_b)
          min_error = min(e1, e2, e3, e4)
          if ( abs(e1) == min_error .OR. abs(e3) == min_error) then
             root = lower_b
          else if ( abs(e2) == min_error .OR. abs(e4) == min_error) then
             root = upper_b
          else
             call fatal_error('Invalid case for handling quadratic bounds overload')
          endif

          write(*,*) 'C_1 = ', C_1
          write(*,*) 'C_2 = ', C_2
          write(*,*) 'Lower bound = ', lower_b
          write(*,*) 'Upper bound = ', upper_b
          write(*,*) 'C_1 - lower_bound = ', C_1 - lower_b
          write(*,*) 'C_2 - lower_bound = ', C_2 - lower_b
          write(*,*) 'C_1 - upper_bound = ', C_1 - upper_b
          write(*,*) 'C_2 - upper_bound = ', C_2 - upper_b
          write(*,*) 'a = ', a
          write(*,*) 'b = ', b
          write(*,*) 'c = ', c
          write(*,*) 'ERROR:  Neither quadratic roots satisifed the criteria'
          call warning('Quadratic roots not within 1E-5 of bounds')
       endif
    endif

  end subroutine get_quadratic_root

  subroutine get_cubic_root(a, b, c, d, lower_b, upper_b, root)

    real(8), intent(in) :: a,b,c,d  ! Coefficients in ax^3 + bx^2 + cx + d = 0
    real(8), intent(in) :: lower_b  ! Lower bound of search
    real(8), intent(in) :: upper_b  ! Upper bound of search
    real(8), intent(inout) :: root  ! The root of the equation
    real(8) :: aa, bb, cc           ! Coefficients in x^3 + aa x^2 + bb x + cc = 0
    real(8) :: Q, R, dd, ee, theta  ! Temporary variables needed in root finding

    ! Calculate coefficients as they appear in the numerical recipes book
    aa = b / a
    bb = c / a
    cc = d / a

    Q = (aa*aa - 3.0_8 * bb) / (9.0_8)
    R = (2.0_8*aa*aa*aa - 9.0_8*aa*bb + 27.0_8*cc) / (54.0_8)

    if (R*R < Q*Q*Q) then
       theta = acos(R/sqrt(Q*Q*Q))
       root = -2.0_8 * sqrt(Q) * cos(theta/3.0_8) - aa / 3.0_8
       if (root < lower_b .OR. root > upper_b) then
          root = -2.0_8 * sqrt(Q) * cos( (theta + 2.0_8*PI)/3.0_8) - aa / 3.0_8
          if (root < lower_b .OR. root > upper_b) then
             root = -2.0_8 * sqrt(Q) * cos( (theta - 2.0_8*PI)/3.0_8) - aa / 3.0_8
             if (root < lower_b .OR. root > upper_b) then
                if (root < lower_b) then
                   !write(*,*) 'Setting to lower bound'
                   root = lower_b
                else
                   root = upper_b
                endif

                write(*,*) 'root = ', root
                write(*,*) 'lower_b = ', lower_b
                write(*,*) 'upper_b = ', upper_b
                write(*,*) 'a = ', a
                write(*,*) 'b = ', b
                write(*,*) 'c = ', c
                write(*,*) 'd = ', d
                call warning('Acceptable cubic root not found')
             endif
          endif
       endif
    else
       if (R >= 0.0) then
          dd = -(abs(R) + sqrt( abs(R*R-Q*Q*Q) ))**(1.0_8/3.0_8)
       else
          dd = (abs(R) + sqrt( abs(R*R-Q*Q*Q) ))**(1.0_8/3.0_8)
       endif

       if (dd == 0.0) then
          ee = 0.0
       else
          ee = Q/dd
       endif

       root = (dd + ee) - aa / 3.0_8
       if (root < lower_b .OR. root > upper_b) then

          if (root < lower_b) then
             !write(*,*) 'Setting to lower bound'
             root = lower_b
          else
             root = upper_b
          endif

          write(*,*) 'root = ', root
          write(*,*) 'lower_b = ', lower_b
          write(*,*) 'upper_b = ', upper_b
          write(*,*) 'a = ', a
          write(*,*) 'b = ', b
          write(*,*) 'c = ', c
          write(*,*) 'd = ', d
          call warning('Acceptable cubic root not found')
       endif

    endif
  end subroutine get_cubic_root

end module tracking
