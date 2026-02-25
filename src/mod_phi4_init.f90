module mod_phi4_init
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State, RNGState, Phi4Obs
  use mod_random,     only: rng_uniform, rng_gauss
  use mod_phi4_model, only: compute_potential, compute_kinetic
  use mod_metropolis_phi4, only: mmc_sweep_micro, mmc_sweep_canon
  implicit none
contains

  subroutine allocate_state(par, st)
    type(Phi4Params), intent(in) :: par
    type(Phi4State),  intent(inout) :: st
    integer :: k

    allocate(st%phi(par%N))
    if (trim(par%ensemble) == "canon") then
      allocate(st%pi(par%N))
      st%pi = 0.0_dp
    end if

    allocate(st%index(par%N))
    do k = 1, par%N
       st%index(k) = k
    end do
  end subroutine allocate_state


  subroutine init_configuration(par, st, rng)
    ! Dispatcher: choose micro vs canon once here (not in the hot loop)
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng

    if (trim(par%ensemble) == "micro") then
      call init_configuration_micro(par, st, rng)
    else
      call init_configuration_canon(par, st, rng)
    end if
  end subroutine init_configuration


  subroutine init_configuration_canon(par, st, rng)
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng
    integer :: i, KK
    real(dp) :: s, inv_sqrt_beta
    type(Phi4Params) :: p
    type(Phi4Obs)    :: obs

    if (.not. allocated(st%pi)) then
      write(*,*) "ERROR: canonical init requires st%pi allocated."
      stop
    end if
    if (par%beta_mc <= 0.0_dp) then
      write(*,*) "ERROR: canonical init requires beta_mc > 0."
      stop
    end if

    ! Initialize phi (Gaussian)
    s = 0.5_dp
    do i = 1, par%N
      st%phi(i) = s * rng_gauss(rng)
    end do
    call compute_potential(par, st)

    ! Initialize momenta pi ~ N(0, 1/beta_mc)
    inv_sqrt_beta = 1.0_dp / sqrt(par%beta_mc)
    do i = 1, par%N
      st%pi(i) = rng_gauss(rng) * inv_sqrt_beta
    end do
    call compute_kinetic(par, st)

    ! Tune delta to target acceptance window
    call tune_delta_canon(par, st, rng)

    ! Short warm-up with tuned delta
    p = par
    obs = Phi4Obs()
    do KK = 1, par%n_therm
      call mmc_sweep_canon(p, st, rng, obs)
    end do

    call compute_potential(par, st)
    call compute_kinetic(par, st)
  end subroutine init_configuration_canon


  subroutine tune_delta_canon(par, st, rng)
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng

    integer :: LL, KK
    real(dp) :: delta_try, acc
    type(Phi4Params) :: p
    type(Phi4Obs)    :: obs
    real(dp), parameter :: acc_lo = 0.51_dp
    real(dp), parameter :: acc_hi = 0.56_dp

    p = par
    acc = 0.0_dp
    delta_try = par%delta

    do LL = 1, 1000
      delta_try = 0.1_dp * real(LL, dp)
      p%delta = delta_try

      obs = Phi4Obs()
      do KK = 1, 10
        call mmc_sweep_canon(p, st, rng, obs)
      end do

      if (obs%n_att > 0) acc = real(obs%n_acc, dp) / real(obs%n_att, dp)
      if (acc <= acc_hi .and. acc > acc_lo) exit
    end do

    par%delta = delta_try
  end subroutine tune_delta_canon


  subroutine init_configuration_micro(par, st, rng)
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng

    integer :: i, LL, KK
    real(dp) :: fraz, r
    real(dp) :: delta_try, acc
    type(Phi4Params) :: p
    type(Phi4Obs)    :: obs
    real(dp), parameter :: acc_lo = 0.51_dp
    real(dp), parameter :: acc_hi = 0.56_dp

    ! Uniform initialization (legacy-style)
    fraz = 2.0_dp
    do i = 1, par%N
      r = rng_uniform(rng)
      st%phi(i) = fraz * (r - 0.5_dp)
    end do
    call compute_potential(par, st)

    ! Ensure E - V > 0
    if (par%Etot - st%V <= 0.0_dp) then
      fraz = 0.5_dp
      do i = 1, par%N
        r = rng_uniform(rng)
        st%phi(i) = fraz * (r - 0.5_dp)
      end do
      call compute_potential(par, st)
    end if
    if (par%Etot - st%V <= 0.0_dp) then
      st%phi = 0.0_dp
      call compute_potential(par, st)
    end if

    ! Tune delta (micro)
    p = par
    acc = 0.0_dp
    delta_try = par%delta
    do LL = 1, 1000
      delta_try = 0.1_dp * real(LL, dp)
      p%delta = delta_try
      obs = Phi4Obs()
      do KK = 1, 10
        call mmc_sweep_micro(p, st, rng, obs)
      end do
      if (obs%n_att > 0) acc = real(obs%n_acc, dp) / real(obs%n_att, dp)
      if (acc <= acc_hi .and. acc > acc_lo) exit
    end do

    par%delta = delta_try
    !write(*,'(A,F10.4,A,F10.4)') " [MICRO INIT] tuned delta = ", par%delta, "  acc ~ ", acc

    ! Warm-up
    p = par
    obs = Phi4Obs()
    do KK = 1, par%n_therm
      call mmc_sweep_micro(p, st, rng, obs)
    end do

    call compute_potential(par, st)

    ! Keep kinetic energy available in the state (micro: K = E - V)
    st%K = par%Etot - st%V
    if (allocated(st%pi)) st%pi = 0.0_dp
  end subroutine init_configuration_micro

end module mod_phi4_init