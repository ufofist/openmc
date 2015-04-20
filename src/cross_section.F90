module cross_section

  use ace_header,      only: Nuclide, SAlphaBeta, Reaction, UrrData
  use constants
  use energy_grid,     only: grid_method, log_spacing
  use error,           only: fatal_error
  use fission,         only: nu_total
  use global
  use list_header,     only: ListElemInt
  use material_header, only: Material
  use particle_header, only: Particle
  use random_lcg,      only: prn
  use search,          only: binary_search

  implicit none
  save

  integer :: union_grid_index
!$omp threadprivate(union_grid_index)

contains

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle), intent(inout) :: p

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in mat % i_sab_nuclides
    real(8) :: atom_density  ! atom density of a nuclide
    logical :: check_sab     ! should we check for S(a,b) table?
    type(Material), pointer :: mat ! current material

    ! Set all material macroscopic cross sections to zero
    p % material_xs % total          = ZERO
    p % material_xs % elastic        = ZERO
    p % material_xs % absorption     = ZERO
    p % material_xs % fission        = ZERO
    p % material_xs % nu_fission     = ZERO
    p % material_xs % kappa_fission  = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    mat => materials(p % material)

    ! Find energy index on global or material unionized grid
    if (grid_method == GRID_MAT_UNION) &
         call find_energy_index(p % E, p % material)

    ! Determine if this material has S(a,b) tables
    check_sab = (mat % n_sab > 0)

    ! Initialize position in i_sab_nuclides
    j = 1

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
      ! ========================================================================
      ! CHECK FOR S(A,B) TABLE

      i_sab = 0

      ! Check if this nuclide matches one of the S(a,b) tables specified -- this
      ! relies on i_sab_nuclides being in sorted order
      if (check_sab) then
        if (i == mat % i_sab_nuclides(j)) then
          ! Get index in sab_tables
          i_sab = mat % i_sab_tables(j)

          ! If particle energy is greater than the highest energy for the S(a,b)
          ! table, don't use the S(a,b) table
          if (p % E > sab_tables(i_sab) % threshold_inelastic) i_sab = 0

          ! Increment position in i_sab_nuclides
          j = j + 1

          ! Don't check for S(a,b) tables if there are no more left
          if (j > mat % n_sab) check_sab = .false.
        end if
      end if

      ! ========================================================================
      ! CALCULATE MICROSCOPIC CROSS SECTION

      ! Determine microscopic cross sections for this nuclide
      i_nuclide = mat % nuclide(i)

      ! Calculate microscopic cross section for this nuclide
      if (p % E /= p % micro_xs(i_nuclide) % last_E) then
        call calculate_nuclide_xs(p, i_nuclide, i_sab, i)
      else if (i_sab /= p % micro_xs(i_nuclide) % last_index_sab) then
        call calculate_nuclide_xs(p, i_nuclide, i_sab, i)
      end if

      ! ========================================================================
      ! ADD TO MACROSCOPIC CROSS SECTION

      ! Copy atom density of nuclide in material
      atom_density = mat % atom_density(i)

      ! Add contributions to material macroscopic total cross section
      p % material_xs % total = p % material_xs % total + &
           atom_density * p % micro_xs(i_nuclide) % total

      ! Add contributions to material macroscopic scattering cross section
      p % material_xs % elastic = p % material_xs % elastic + &
           atom_density * p % micro_xs(i_nuclide) % elastic

      ! Add contributions to material macroscopic absorption cross section
      p % material_xs % absorption = p % material_xs % absorption + &
           atom_density * p % micro_xs(i_nuclide) % absorption

      ! Add contributions to material macroscopic fission cross section
      p % material_xs % fission = p % material_xs % fission + &
           atom_density * p % micro_xs(i_nuclide) % fission

      ! Add contributions to material macroscopic nu-fission cross section
      p % material_xs % nu_fission = p % material_xs % nu_fission + &
           atom_density * p % micro_xs(i_nuclide) % nu_fission

      ! Add contributions to material macroscopic energy release from fission
      p % material_xs % kappa_fission = p % material_xs % kappa_fission + &
           atom_density * p % micro_xs(i_nuclide) % kappa_fission
    end do

  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(p, i_nuclide, i_sab, i_nuc_mat)

    type(Particle), intent(inout) :: p
    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    integer, intent(in) :: i_nuc_mat ! index into nuclides array for a material
    integer :: i_grid ! index on nuclide energy grid
    integer :: i_low  ! lower logarithmic mapping index
    integer :: i_high ! upper logarithmic mapping index
    integer :: u      ! index into logarithmic mapping array
    real(8) :: f             ! interp factor on nuclide energy grid
    type(Nuclide),  pointer :: nuc
    type(Material), pointer :: mat

    ! Set pointer to nuclide and material
    nuc => nuclides(i_nuclide)
    mat => materials(p % material)

    ! Determine index on nuclide energy grid
    select case (grid_method)
    case (GRID_MAT_UNION)

      i_grid = mat % nuclide_grid_index(i_nuc_mat, union_grid_index)

    case (GRID_LOGARITHM)
      ! Determine the energy grid index using a logarithmic mapping to reduce
      ! the energy range over which a binary search needs to be performed

      if (p % E < nuc % energy(1)) then
        i_grid = 1
      elseif (p % E > nuc % energy(nuc % n_grid)) then
        i_grid = nuc % n_grid - 1
      else
        ! Determine bounding indices based on which equal log-spaced interval
        ! the energy is in
        u = int(log(p % E/1.0e-11_8)/log_spacing)
        i_low  = nuc % grid_index(u)
        i_high = nuc % grid_index(u + 1) + 1

        ! Perform binary search over reduced range
        i_grid = binary_search(nuc % energy(i_low:i_high), &
             i_high - i_low + 1, p % E) + i_low - 1
      end if

    case (GRID_NUCLIDE)
      ! Perform binary search on the nuclide energy grid in order to determine
      ! which points to interpolate between

      if (p % E <= nuc % energy(1)) then
        i_grid = 1
      elseif (p % E > nuc % energy(nuc % n_grid)) then
        i_grid = nuc % n_grid - 1
      else
        i_grid = binary_search(nuc % energy, nuc % n_grid, p % E)
      end if

    end select

    ! check for rare case where two energy points are the same
    if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1

    ! calculate interpolation factor
    f = (p % E - nuc%energy(i_grid))/(nuc%energy(i_grid+1) - nuc%energy(i_grid))

    p % micro_xs(i_nuclide) % index_grid    = i_grid
    p % micro_xs(i_nuclide) % interp_factor = f

    ! Initialize sab treatment to false
    p % micro_xs(i_nuclide) % index_sab   = NONE
    p % micro_xs(i_nuclide) % elastic_sab = ZERO

    ! Initialize URR probability table treatment to false
    p % micro_xs(i_nuclide) % use_ptable  = .false.

    ! Initialize nuclide cross-sections to zero
    p % micro_xs(i_nuclide) % fission    = ZERO
    p % micro_xs(i_nuclide) % nu_fission = ZERO
    p % micro_xs(i_nuclide) % kappa_fission  = ZERO

    ! Calculate microscopic nuclide total cross section
    p % micro_xs(i_nuclide) % total = (ONE - f) * nuc % total(i_grid) &
         + f * nuc % total(i_grid+1)

    ! Calculate microscopic nuclide elastic cross section
    p % micro_xs(i_nuclide) % elastic = (ONE - f) * nuc % elastic(i_grid) &
         + f * nuc % elastic(i_grid+1)

    ! Calculate microscopic nuclide absorption cross section
    p % micro_xs(i_nuclide) % absorption = (ONE - f) * nuc % absorption( &
         i_grid) + f * nuc % absorption(i_grid+1)

    if (nuc % fissionable) then
      ! Calculate microscopic nuclide total cross section
      p % micro_xs(i_nuclide) % fission = (ONE - f) * nuc % fission(i_grid) &
           + f * nuc % fission(i_grid+1)

      ! Calculate microscopic nuclide nu-fission cross section
      p % micro_xs(i_nuclide) % nu_fission = (ONE - f) * nuc % nu_fission( &
           i_grid) + f * nuc % nu_fission(i_grid+1)

      ! Calculate microscopic nuclide kappa-fission cross section
      ! The ENDF standard (ENDF-102) states that MT 18 stores
      ! the fission energy as the Q_value (fission(1))
      p % micro_xs(i_nuclide) % kappa_fission = &
           nuc % reactions(nuc % index_fission(1)) % Q_value * &
           p % micro_xs(i_nuclide) % fission
    end if

    ! If there is S(a,b) data for this nuclide, we need to do a few
    ! things. Since the total cross section was based on non-S(a,b) data, we
    ! need to correct it by subtracting the non-S(a,b) elastic cross section and
    ! then add back in the calculated S(a,b) elastic+inelastic cross section.

    if (i_sab > 0) call calculate_sab_xs(p, i_nuclide, i_sab)

    ! if the particle is in the unresolved resonance range and there are
    ! probability tables, we need to determine cross sections from the table

    if (urr_ptables_on .and. nuc % urr_present) then
      if (p % E > nuc % urr_data % energy(1) .and. &
           p % E < nuc % urr_data % energy(nuc % urr_data % n_energy)) then
        call calculate_urr_xs(p, i_nuclide)
      end if
    end if

    p % micro_xs(i_nuclide) % last_E = p % E
    p % micro_xs(i_nuclide) % last_index_sab = i_sab

  end subroutine calculate_nuclide_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace
! whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(p, i_nuclide, i_sab)

    type(Particle), intent(inout) :: p
    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array

    integer :: i_grid    ! index on S(a,b) energy grid
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section
    type(SAlphaBeta), pointer :: sab

    ! Set flag that S(a,b) treatment should be used for scattering
    p % micro_xs(i_nuclide) % index_sab = i_sab

    ! Get pointer to S(a,b) table
    sab => sab_tables(i_sab)

    ! Get index and interpolation factor for inelastic grid
    if (p % E < sab % inelastic_e_in(1)) then
      i_grid = 1
      f = ZERO
    else
      i_grid = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, p % E)
      f = (p % E - sab%inelastic_e_in(i_grid)) / &
           (sab%inelastic_e_in(i_grid+1) - sab%inelastic_e_in(i_grid))
    end if

    ! Calculate S(a,b) inelastic scattering cross section
    inelastic = (ONE - f) * sab % inelastic_sigma(i_grid) + &
         f * sab % inelastic_sigma(i_grid + 1)

    ! Check for elastic data
    if (p % E < sab % threshold_elastic) then
      ! Determine whether elastic scattering is given in the coherent or
      ! incoherent approximation. For coherent, the cross section is
      ! represented as P/E whereas for incoherent, it is simply P

      if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
        if (p % E < sab % elastic_e_in(1)) then
          ! If energy is below that of the lowest Bragg peak, the elastic
          ! cross section will be zero
          elastic = ZERO
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, p % E)
          elastic = sab % elastic_P(i_grid) / p % E
        end if
      else
        ! Determine index on elastic energy grid
        if (p % E < sab % elastic_e_in(1)) then
          i_grid = 1
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, p % E)
        end if

        ! Get interpolation factor for elastic grid
        f = (p % E - sab%elastic_e_in(i_grid))/(sab%elastic_e_in(i_grid+1) - &
             sab%elastic_e_in(i_grid))

        ! Calculate S(a,b) elastic scattering cross section
        elastic = (ONE - f) * sab % elastic_P(i_grid) + &
             f * sab % elastic_P(i_grid + 1)
      end if
    else
      ! No elastic data
      elastic = ZERO
    end if

    ! Correct total and elastic cross sections
    p % micro_xs(i_nuclide) % total = p % micro_xs(i_nuclide) % total - &
         p % micro_xs(i_nuclide) % elastic + inelastic + elastic
    p % micro_xs(i_nuclide) % elastic = inelastic + elastic

    ! Store S(a,b) elastic cross section for sampling later
    p % micro_xs(i_nuclide) % elastic_sab = elastic

  end subroutine calculate_sab_xs

!===============================================================================
! CALCULATE_URR_XS determines cross sections in the unresolved resonance range
! from probability tables
!===============================================================================

  subroutine calculate_urr_xs(p, i_nuclide)

    type(Particle), intent(inout) :: p
    integer, intent(in) :: i_nuclide ! index into nuclides array

    integer :: i            ! loop index
    integer :: i_energy     ! index for energy
    integer :: i_low        ! band index at lower bounding energy
    integer :: i_up         ! band index at upper bounding energy
    integer :: same_nuc_idx ! index of same nuclide
    real(8) :: f            ! interpolation factor
    real(8) :: r            ! pseudo-random number
    real(8) :: elastic      ! elastic cross section
    real(8) :: capture      ! (n,gamma) cross section
    real(8) :: fission      ! fission cross section
    real(8) :: inelastic    ! inelastic cross section
    logical :: same_nuc     ! do we know the xs for this nuclide at this energy?
    type(UrrData),  pointer :: urr
    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn

    p % micro_xs(i_nuclide) % use_ptable = .true.

    ! get pointer to probability table
    nuc => nuclides(i_nuclide)
    urr => nuc % urr_data

    ! determine energy table
    i_energy = 1
    do
      if (p % E < urr % energy(i_energy + 1)) exit
      i_energy = i_energy + 1
    end do

    ! determine interpolation factor on table
    f = (p % E - urr % energy(i_energy)) / &
         (urr % energy(i_energy + 1) - urr % energy(i_energy))

    ! sample probability table using the cumulative distribution

    ! if we're dealing with a nuclide that we've previously encountered at
    ! this energy but a different temperature, use the original random number to
    ! preserve correlation of temperature in probability tables
    same_nuc = .false.
    do i = 1, nuc % nuc_list % size()
      if (p % E /= ZERO .and. p % E == p % micro_xs(nuc % nuc_list % get_item(i)) % last_E) then
        same_nuc = .true.
        same_nuc_idx = i
        exit
      end if
    end do

    if (same_nuc) then
      r = p % micro_xs(nuc % nuc_list % get_item(same_nuc_idx)) % last_prn
    else
      r = prn()
      p % micro_xs(i_nuclide) % last_prn = r
    end if

    i_low = 1
    do
      if (urr % prob(i_energy, URR_CUM_PROB, i_low) > r) exit
      i_low = i_low + 1
    end do
    i_up = 1
    do
      if (urr % prob(i_energy + 1, URR_CUM_PROB, i_up) > r) exit
      i_up = i_up + 1
    end do

    ! determine elastic, fission, and capture cross sections from probability
    ! table
    if (urr % interp == LINEAR_LINEAR) then
      elastic = (ONE - f) * urr % prob(i_energy, URR_ELASTIC, i_low) + &
           f * urr % prob(i_energy + 1, URR_ELASTIC, i_up)
      fission = (ONE - f) * urr % prob(i_energy, URR_FISSION, i_low) + &
           f * urr % prob(i_energy + 1, URR_FISSION, i_up)
      capture = (ONE - f) * urr % prob(i_energy, URR_N_GAMMA, i_low) + &
           f * urr % prob(i_energy + 1, URR_N_GAMMA, i_up)
    elseif (urr % interp == LOG_LOG) then
      ! Get logarithmic interpolation factor
      f = log(p % E / urr % energy(i_energy)) / &
           log(urr % energy(i_energy + 1) / urr % energy(i_energy))

      ! Calculate elastic cross section/factor
      elastic = ZERO
      if (urr % prob(i_energy, URR_ELASTIC, i_low) > ZERO .and. &
          urr % prob(i_energy + 1, URR_ELASTIC, i_up) > ZERO) then
        elastic = exp((ONE - f) * log(urr % prob(i_energy, URR_ELASTIC, &
             i_low)) + f * log(urr % prob(i_energy + 1, URR_ELASTIC, &
             i_up)))
      end if

      ! Calculate fission cross section/factor
      fission = ZERO
      if (urr % prob(i_energy, URR_FISSION, i_low) > ZERO .and. &
          urr % prob(i_energy + 1, URR_FISSION, i_up) > ZERO) then
        fission = exp((ONE - f) * log(urr % prob(i_energy, URR_FISSION, &
             i_low)) + f * log(urr % prob(i_energy + 1, URR_FISSION, &
             i_up)))
      end if

      ! Calculate capture cross section/factor
      capture = ZERO
      if (urr % prob(i_energy, URR_N_GAMMA, i_low) > ZERO .and. &
          urr % prob(i_energy + 1, URR_N_GAMMA, i_up) > ZERO) then
        capture = exp((ONE - f) * log(urr % prob(i_energy, URR_N_GAMMA, &
             i_low)) + f * log(urr % prob(i_energy + 1, URR_N_GAMMA, &
             i_up)))
      end if
    end if

    ! Determine treatment of inelastic scattering
    inelastic = ZERO
    if (urr % inelastic_flag > 0) then
      ! Get pointer to inelastic scattering reaction
      rxn => nuc % reactions(nuc % urr_inelastic)

      ! Get index on energy grid and interpolation factor
      i_energy = p % micro_xs(i_nuclide) % index_grid
      f = p % micro_xs(i_nuclide) % interp_factor

      ! Determine inelastic scattering cross section
      if (i_energy >= rxn % threshold) then
        inelastic = (ONE - f) * rxn % sigma(i_energy - rxn%threshold + 1) + &
             f * rxn % sigma(i_energy - rxn%threshold + 2)
      end if
    end if

    ! Multiply by smooth cross-section if needed
    if (urr % multiply_smooth) then
      elastic = elastic * p % micro_xs(i_nuclide) % elastic
      capture = capture * (p % micro_xs(i_nuclide) % absorption - &
           p % micro_xs(i_nuclide) % fission)
      fission = fission * p % micro_xs(i_nuclide) % fission
    end if

    ! Check for negative values
    if (elastic < ZERO) elastic = ZERO
    if (fission < ZERO) fission = ZERO
    if (capture < ZERO) capture = ZERO

    ! Set elastic, absorption, fission, and total cross sections. Note that the
    ! total cross section is calculated as sum of partials rather than using the
    ! table-provided value
    p % micro_xs(i_nuclide) % elastic = elastic
    p % micro_xs(i_nuclide) % absorption = capture + fission
    p % micro_xs(i_nuclide) % fission = fission
    p % micro_xs(i_nuclide) % total = elastic + inelastic + capture + fission

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      p % micro_xs(i_nuclide) % nu_fission = nu_total(nuc, p % E) * &
           p % micro_xs(i_nuclide) % fission
    end if

  end subroutine calculate_urr_xs

!===============================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid at a certain
! energy
!===============================================================================

  subroutine find_energy_index(E, i_mat)

    real(8), intent(in) :: E       ! energy of particle
    integer, intent(in) :: i_mat   ! material index
    type(Material), pointer :: mat ! pointer to current material

    mat => materials(i_mat)

    ! if the energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (E <= mat % e_grid(1)) then
      union_grid_index = 1
    elseif (E > mat % e_grid(mat % n_grid)) then
      union_grid_index = mat % n_grid - 1
    else
      union_grid_index = binary_search(mat % e_grid, mat % n_grid, E)
    end if

  end subroutine find_energy_index

!===============================================================================
! 0K_ELASTIC_XS determines the microscopic 0K elastic cross section
! for a given nuclide at the trial relative energy used in resonance scattering
!===============================================================================

  function elastic_xs_0K(E, nuc) result(xs_out)

    type(Nuclide), pointer :: nuc    ! target nuclide at temperature
    integer                :: i_grid ! index on nuclide energy grid
    real(8)                :: f      ! interp factor on nuclide energy grid
    real(8), intent(inout) :: E      ! trial energy
    real(8)                :: xs_out ! 0K xs at trial energy

    ! Determine index on nuclide energy grid
    if (E < nuc % energy_0K(1)) then
      i_grid = 1
    elseif (E > nuc % energy_0K(nuc % n_grid_0K)) then
      i_grid = nuc % n_grid_0K - 1
    else
      i_grid = binary_search(nuc % energy_0K, nuc % n_grid_0K, E)
    end if

    ! check for rare case where two energy points are the same
    if (nuc % energy_0K(i_grid) == nuc % energy_0K(i_grid+1)) then
      i_grid = i_grid + 1
    end if

    ! calculate interpolation factor
    f = (E - nuc % energy_0K(i_grid)) &
      & / (nuc % energy_0K(i_grid + 1) - nuc % energy_0K(i_grid))

    ! Calculate microscopic nuclide elastic cross section
    xs_out = (ONE - f) * nuc % elastic_0K(i_grid) &
      & + f * nuc % elastic_0K(i_grid + 1)

  end function elastic_xs_0K

end module cross_section
