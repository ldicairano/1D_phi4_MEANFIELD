program main_phi4
  use mod_kinds,              only: dp
  use mod_types_phi4,         only: Phi4Params, IOParams, RNGState, Phi4State, Phi4Obs
  use mod_cli_phi4,           only: parse_cli_phi4, cli_msg_restart_already_complete
  use mod_input_phi4,         only: read_input_file
  use mod_random,             only: rng_init
  use mod_phi4_init,          only: allocate_state, init_configuration
  use mod_metropolis_phi4,    only: sweep_iface
  use mod_ensemble_dispatch,  only: select_sweep
  use mod_io_observables_phi4,only: reset_obs, update_iface, reinflate_iface, select_obs_ops
  use mod_io_append_dispatch, only: select_append, append_iface
  use mod_io_restart_phi4,    only: read_restart, write_restart
  implicit none

!   ██╗██████╗      ███╗   ███╗███████╗    ██████╗██╗  ██╗██╗██╗  ██╗
! ████║██╔══███╗    ████╗ ████║██╔════╝    ██╔═██║██║  ██║██║██║  ██║
! ╚═██║██║  ███║    ██╔████╔██║█████╗      ██████║███████║██║███████║
!   ██║██║  ███║    ██║╚██╔╝██║██╔══╝      ██╔═══╝██╔══██║██║╚════██║
!   ██║██████╔╝     ██║ ╚═╝ ██║██║         ██║    ██║  ██║██║     ██║
!   ╚═╝╚═════╝      ╚═╝     ╚═╝╚═╝         ╚═╝    ╚═╝  ╚═╝╚═╝     ╚═╝
!
!                                  ϕ^4   (MEAN-FIELD, 1D)
!                                   Luxembourg, 25 Feb 2026
!  //============================================================================\\
!  |                        1D  MEAN-FIELD  ϕ^4  MODEL                              |
!  |                     (CHAIN  +  ALL-TO-ALL  TERM)                               |
!  |                    MICROCANONICAL  MONTE  CARLO  (MMC)                         |
!  |                    (MOMENTA INTEGRATED OUT, ϕ-ONLY)                            |
!  \\===========================================================================//

!                            FROM SUMMER 2024
! Rome,       15 Aug 2024 (ZN_THEORY)
! Luxembourg, 19 Sep 2024 (SINE-GORDON)
! Luxembourg, 26 Nov 2024 (phi^4)
! Luxembourg, 18 Jun 2025 (real phi^4)

! ------------------------------------------------
!  Hamiltonian (typical form):
!
! H = sum_i [ p_i^2/2 ] +
!     sum_i [ (lambda/4!) phi_i^4 - (mu^2/2) phi_i^2 ] +
!     (2J/(N-1)) * ( N*sum_i phi_i^2 - (sum_i phi_i)^2 )    (all-to-all mean-field)
!      NOTE: in code par%mu is mu^2 (paper notation).
!
!  Microcanonical configuration weight after integrating the momenta:
!
!    w(ϕ) ∝ [E_tot - V(ϕ)]^(N/2 - 1)  Θ(E_tot - V(ϕ))

  type(Phi4Params) :: par
  type(IOParams)   :: io
  type(RNGState)   :: rng
  type(Phi4State)  :: st
  type(Phi4Obs)    :: obs

  procedure(sweep_iface),     pointer :: sweep     => null()
  procedure(append_iface),    pointer :: append    => null()
  procedure(update_iface),    pointer :: update    => null()
  procedure(reinflate_iface), pointer :: reinflate => null()

  character(len=256) :: input_file
  integer :: n_samp_cli, L_cli, n_steps_cli, n_jump_cli, n_realiz_cli, n_rest_cli
  integer :: step, step0, n_print, n_sweep
  integer :: n_MC
  real(dp) :: ekin
  real(dp) :: t_in, t_out, time_tot, time_config_in, time_config_out

  time_tot = 0.0_dp
  CALL CPU_TIME(time_config_in)

  ! ---------------------------
  ! CLI parsing
  !   ./phi4_mmc [-i input.inp] n_samp L n_steps n_jump n_realiz n_rest
  ! ---------------------------
  call parse_cli_phi4(input_file, n_samp_cli, L_cli, n_steps_cli, n_jump_cli, n_realiz_cli, n_rest_cli)

  ! ---------------------------
  ! Read input file (namelist)
  ! ---------------------------
  call read_input_file(trim(input_file), par, io, rng)

  ! ---------------------------
  ! Override from CLI (legacy-compatible)
  ! ---------------------------
  par%L       = L_cli
  par%N       = par%L

  ! ---------------------------
  ! Auto-normalize mean-field coupling using runtime N
  ! Target: -(1/(4N)) (sum_i q_i)^2
  ! In this code the MF piece is:  -(2*coup/(N-1)) (sum_i q_i)^2
  ! => choose coup = (N-1)/(8N)
  ! ---------------------------
  if (par%auto_coup) then
    if (par%N <= 1) stop "ERROR: auto_coup requires N>1"
    par%coup = real(par%N - 1, dp) / (8.0_dp * real(par%N, dp))
  end if

  par%n_steps = n_steps_cli
  par%n_jump  = n_jump_cli
  n_print     = par%n_print
  ! External run identifier (trajectory id)
  par%n_realiz = n_realiz_cli
  
  ! Segment index for chained restarts
  par%n_rest   = n_rest_cli

  ! --- legacy sampling id + energy-per-site ---
  par%n_samp = n_samp_cli


  if (trim(par%ensemble) == "micro") then
    par%en_in = 0.01_dp * real(par%n_samp, dp)
    par%Etot  = real(par%N, dp) * par%en_in
  else if (trim(par%ensemble) == "canon") then
    par%beta_mc = 0.00001_dp * real(par%n_samp, dp)
  else
    stop "ERROR: ensemble must be 'micro' or 'canon'"
  end if

  ! ---------------------------
  ! RNG + allocate
  ! ---------------------------
  call rng_init(rng, par%seed)
  call allocate_state(par, st)

  ! ---------------------------
  ! Select sweep  
  ! ---------------------------
  call select_sweep(par, sweep)

  ! ---------------------------
  ! Select subroutine for printing  
  ! --------------------------- 
  call select_append(par, append)

  ! ---------------------------
  ! Select subroutine for updating observable 
  ! --------------------------- 
  call select_obs_ops(par, update, reinflate)
  ! ---------------------------
  ! Realizations loop
  ! ---------------------------

  ! initialize the observables
    call reset_obs(obs)
    n_MC = 0
    step0 = 0
    
    ! Set initial conditions or restart the trajectory
    if (par%n_rest == 1) then
      call init_configuration(par, st, rng)
      n_MC = 0
    else
      call read_restart(io, par, st, obs, rng, step0, par%n_realiz)   ! auto: reads only if n_rest>1
      call reinflate(obs)  !
      n_MC = obs%n_MC
    end if


    ! Check if the total length has been already reached
    if (par%n_rest > 1) then
      if (step0 >= par%n_steps) then
        call cli_msg_restart_already_complete(par, step0)
        stop 0
      end if
    end if


    call cpu_time(time_config_out)
    time_tot = time_config_out - time_config_in



    ! ---- INTEGRATION ----
    do step = step0 + 1, par%n_steps
      call cpu_time(t_in)

      ! ---- MONTE CARLO STEP ----
      DO n_sweep = 1, par%n_jump
         call sweep(par, st, rng, obs)
      END DO
      ! ---------------------
      call cpu_time(t_out)
      time_tot = time_tot + abs(t_out - t_in)

      ! --- early stop close to walltime ---
      if (time_tot >= 0.95_dp * par%cluster_time) then
        obs%n_MC = n_MC
        call write_restart(io, par, st, obs, rng, step, par%n_realiz)
        exit
      end if
      n_MC = n_MC + 1
      ekin = st%K

      !    --- measurement schedule ---
      call update(par, obs, n_MC, ekin, st%V, st%M)

      !    --- printing schedule ---
      if (mod(n_MC, n_print) == 0) then
        call append(io, par, obs, time_tot, n_MC)
      end if

    end do

    obs%n_MC = n_MC
    call write_restart(io, par, st, obs, rng, par%n_steps, par%n_realiz)

end program main_phi4
