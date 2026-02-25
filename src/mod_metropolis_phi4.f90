module mod_metropolis_phi4
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State, Phi4Obs, RNGState
  use mod_random,     only: rng_uniform
  use mod_phi4_model, only: deltaV_local_fast
  implicit none

  public :: sweep_iface
  public :: mmc_sweep_canon, mmc_sweep_micro

  abstract interface
    subroutine sweep_iface(par, st, rng, obs)
      import :: Phi4Params, Phi4State, RNGState, Phi4Obs
      type(Phi4Params), intent(in)    :: par
      type(Phi4State),  intent(inout) :: st
      type(RNGState),   intent(inout) :: rng
      type(Phi4Obs),    intent(inout) :: obs
    end subroutine sweep_iface
  end interface

contains

  pure real(dp) function micro_ratio(par, Vold, Vnew) result(r)
    type(Phi4Params), intent(in) :: par
    real(dp), intent(in) :: Vold, Vnew
    real(dp) :: aold, anew, expo
    expo = 0.5_dp*real(par%N,dp) - 1.0_dp
    aold = par%Etot - Vold
    anew = par%Etot - Vnew
    if (aold <= 0.0_dp .or. anew <= 0.0_dp) then
      r = 0.0_dp
    else
      r = (anew/aold)**expo
    end if
  end function micro_ratio

  pure subroutine update_macros_accept(par, st, phi_old, phi_new)
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
    real(dp), intent(in) :: phi_old, phi_new
    real(dp) :: d1, d2, d4, o2, n2

    d1 = phi_new - phi_old
    o2 = phi_old*phi_old
    n2 = phi_new*phi_new
    d2 = n2 - o2
    d4 = n2*n2 - o2*o2

    st%S1 = st%S1 + d1
    st%S2 = st%S2 + d2
    st%S4 = st%S4 + d4

    st%M    = st%S1 / real(par%N, dp)
    st%phi2 = st%S2 / real(par%N, dp)
    st%phi4 = st%S4 / real(par%N, dp)
  end subroutine update_macros_accept


  subroutine mmc_sweep_canon(par, st, rng, obs)
    ! Canonical Metropolis on full Hamiltonian H = K(pi) + V(phi).
    ! One sweep:
    !   1) Fisher–Yates shuffle of site indices (uniform random permutation)
    !   2) Visit sites in shuffled order (exactly one attempt per site)
    ! Optimization:
    !   - avoid calling exp() for downhill moves (dH <= 0)
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng
    type(Phi4Obs),    intent(inout) :: obs
  
    integer :: k, i, j, tmp
    real(dp) :: phi_old, phi_new, dV, Vnew
    real(dp) :: p_old, p_new, dK
    real(dp) :: dH, u
  
    ! ---------------------------
    ! Fisher–Yates shuffle (in-place)
    ! ---------------------------
    do k = par%N, 2, -1
      u = rng_uniform(rng)                     ! u in [0,1)
      j = 1 + int(u * real(k, dp))             ! j in [1,k]
      tmp = st%index(j)
      st%index(j) = st%index(k)
      st%index(k) = tmp
    end do
  
    ! ---------------------------
    ! Metropolis updates in permuted order
    ! ---------------------------
    do k = 1, par%N
      i = st%index(k)
  
      phi_old = st%phi(i)
      phi_new = phi_old + par%delta*(2.0_dp*rng_uniform(rng) - 1.0_dp)
  
      call deltaV_local_fast(par, st%S1, st%phi(i), phi_new, dV)
      Vnew = st%V + dV
  
      p_old = st%pi(i)
      p_new = p_old + par%delta*(2.0_dp*rng_uniform(rng) - 1.0_dp)
      dK = 0.5_dp*(p_new*p_new - p_old*p_old)
  
      dH = dV + dK
  
      obs%n_att = obs%n_att + 1
  
      ! Accept if downhill; otherwise accept with exp(-beta_mc * dH)
      if (dH <= 0.0_dp) then
        st%phi(i) = phi_new
        st%V      = Vnew
        call update_macros_accept(par, st, phi_old, phi_new)
  
        st%pi(i) = p_new
        st%K     = st%K + dK
  
        obs%n_acc = obs%n_acc + 1
      else
        u = rng_uniform(rng)
        if (u < exp(-par%beta_mc*dH)) then
          st%phi(i) = phi_new
          st%V      = Vnew
          call update_macros_accept(par, st, phi_old, phi_new)
  
          st%pi(i) = p_new
          st%K     = st%K + dK
  
          obs%n_acc = obs%n_acc + 1
        end if
      end if
    end do
  
    obs%acc_rate = real(obs%n_acc,dp)/real(obs%n_att,dp)
  end subroutine mmc_sweep_canon


  subroutine mmc_sweep_micro(par, st, rng, obs)
    use mod_kinds,      only: dp
    use mod_random,     only: rng_uniform
    use mod_phi4_model, only: deltaV_local_fast
    implicit none
  
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng
    type(Phi4Obs),    intent(inout) :: obs
  
    integer  :: k, i, j, tmp
    real(dp) :: phi_old, phi_new, dV, Vnew
    real(dp) :: rr, u, dW, expo
  
    ! ---- exponent (microcanonical configurational weight) ----
    ! This corresponds to alpha = N_kin/2 - 1. If you use a different N_dof,
    ! replace the next line accordingly.
    expo = 0.5_dp*real(par%N,dp) - 1.0_dp
  
    ! ---------------------------
    ! Fisher–Yates shuffle (in-place)
    ! ---------------------------
    do k = par%N, 2, -1
      rr = rng_uniform(rng)                 ! rr in [0,1)
      j  = 1 + int(rr * real(k, dp))        ! j in [1,k]
      tmp = st%index(j)
      st%index(j) = st%index(k)
      st%index(k) = tmp
    end do
  
    ! ---------------------------
    ! Microcanonical MMC sweep (permuted order)
    ! ---------------------------
    do k = 1, par%N
      i = st%index(k)
  
      phi_old = st%phi(i)
      phi_new = phi_old + par%delta*(2.0_dp*rng_uniform(rng) - 1.0_dp)
  
      call deltaV_local_fast(par, st%S1, phi_old, phi_new, dV)
      Vnew = st%V + dV
  
      obs%n_att = obs%n_att + 1
  
      dW = log(par%Etot - Vnew) - log(par%Etot - st%V)
      u  = rng_uniform(rng)

      ! Hard constraint: Etot must exceed the new potential
      if (par%Etot > Vnew .AND. u < exp(expo*dW)) then
          st%phi(i) = phi_new
          st%V      = Vnew
          call update_macros_accept(par, st, phi_old, phi_new)
          obs%n_acc = obs%n_acc + 1
      end if
  
    end do
  
    obs%acc_rate = real(obs%n_acc,dp)/real(obs%n_att,dp)
  
    ! Keep kinetic energy available in the state (micro: K = E - V)
    st%K = par%Etot - st%V
  end subroutine mmc_sweep_micro

end module mod_metropolis_phi4