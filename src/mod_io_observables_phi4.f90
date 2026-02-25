module mod_io_observables_phi4
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4Obs, IOParams
  implicit none
  private
  public :: reset_obs
  public :: append_iface
  public :: update_iface, reinflate_iface, select_obs_ops
  public :: update_obs_micro, update_obs_canon
  public :: reinflate_micro, reinflate_canon


  abstract interface
  subroutine append_iface(io, par, obs, time_tot, n_MC)
    import :: IOParams, Phi4Params, Phi4Obs, dp
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4Obs),    intent(in) :: obs
    real(dp),         intent(in) :: time_tot
    integer,          intent(in) :: n_MC
  end subroutine append_iface
  end interface

  abstract interface
  subroutine update_iface(par, obs, n_MC, ekin, pot, mag)
    import :: Phi4Params, Phi4Obs, dp
    type(Phi4Params), intent(in)    :: par
    type(Phi4Obs),    intent(inout) :: obs
    integer,          intent(in)    :: n_MC
    real(dp),         intent(in)    :: ekin, pot, mag
  end subroutine
  end interface

  abstract interface
    subroutine reinflate_iface(obs)
    import :: Phi4Obs
    type(Phi4Obs), intent(inout) :: obs
    end subroutine
  end interface

contains


subroutine select_obs_ops(par, update, reinflate)
  type(Phi4Params), intent(in) :: par
  procedure(update_iface),    pointer, intent(out) :: update
  procedure(reinflate_iface), pointer, intent(out) :: reinflate

  select case (trim(par%ensemble))
  case ("micro")
    update    => update_obs_micro
    reinflate => reinflate_micro
  case ("canon")
    update    => update_obs_canon
    reinflate => reinflate_canon
  case default
    stop "ERROR: ensemble must be 'micro' or 'canon' (select_obs_ops)"
  end select
end subroutine



  subroutine reset_obs(obs)
    type(Phi4Obs), intent(inout) :: obs
    obs = Phi4Obs()
  end subroutine reset_obs


  subroutine update_obs_micro(par, obs, n_MC, ekin, pot, mag)
    type(Phi4Params), intent(in)    :: par
    type(Phi4Obs),    intent(inout) :: obs
    integer,          intent(in)    :: n_MC
    real(dp),         intent(in)    :: ekin, pot, mag
  
    real(dp) :: Neff, kin_per_dof
    real(dp) :: a1, a2, a3
  
    Neff = real(par%N - 1, dp)
    kin_per_dof = ekin / Neff
  
    ! --- accumulators (micro-style, kinetic moments) ---
    obs%av_mag  = obs%av_mag  + abs(mag) / real(par%N, dp)
    obs%av_mag2 = obs%av_mag2 + mag*mag
    obs%av_mag4 = obs%av_mag4 + mag*mag*mag*mag
  
    obs%av_pot  = obs%av_pot + pot / real(par%N, dp)
    obs%av_kin  = obs%av_kin + ekin / real(par%N, dp)
  
    obs%av_kin_1 = obs%av_kin_1 + 1.0_dp/kin_per_dof
    obs%av_kin_2 = obs%av_kin_2 + 1.0_dp/(kin_per_dof*kin_per_dof)
    obs%av_kin_3 = obs%av_kin_3 + 1.0_dp/(kin_per_dof*kin_per_dof*kin_per_dof)
  
    obs%av_en   = obs%av_en + (ekin + pot) / real(par%N, dp)
  
    ! --- running means ---
    obs%c_mag  = obs%av_mag  / real(n_MC, dp)
    obs%c_mag2 = obs%av_mag2 / real(n_MC, dp)
    obs%c_mag4 = obs%av_mag4 / real(n_MC, dp)
  
    obs%binder = 1.0_dp - obs%c_mag4/(3.0_dp*obs%c_mag2*obs%c_mag2)
  
    obs%c_pot = obs%av_pot / real(n_MC, dp)
    obs%c_kin = obs%av_kin / real(n_MC, dp)
  
    obs%c_kin_1 = obs%av_kin_1 / real(n_MC, dp)
    obs%c_kin_2 = obs%av_kin_2 / real(n_MC, dp)
    obs%c_kin_3 = obs%av_kin_3 / real(n_MC, dp)
  
    obs%c_en = obs%av_en / real(n_MC, dp)
  
    ! --- microcanonical beta and entropy derivatives (legacy formulas) ---
    a1 = (0.5_dp - 1.0_dp/Neff)
    a2 = (0.5_dp - 2.0_dp/Neff)
    a3 = (0.5_dp - 3.0_dp/Neff)
  
    obs%beta = a1 * obs%c_kin_1
  
    obs%ders_2 = Neff * ( a1*a2*obs%c_kin_2 - (a1*a1)*(obs%c_kin_1*obs%c_kin_1) )
  
    obs%ders_3 = Neff*Neff * ( &
         a1*a2*a3*obs%c_kin_3 &
       - 3.0_dp*(a1*a1)*a2*obs%c_kin_2*obs%c_kin_1 &
       + 2.0_dp*(a1*a1*a1)*(obs%c_kin_1*obs%c_kin_1*obs%c_kin_1) )
  
    ! canonical-only fields: keep harmlessly zero
    obs%c_en2 = 0.0_dp
    obs%cv    = 0.0_dp
  
    obs%n_MC = n_MC
  end subroutine update_obs_micro


  subroutine update_obs_canon(par, obs, n_MC, ekin, pot, mag)
    type(Phi4Params), intent(in)    :: par
    type(Phi4Obs),    intent(inout) :: obs
    integer,          intent(in)    :: n_MC
    real(dp),         intent(in)    :: ekin, pot, mag
  
    real(dp) :: e
  
    e = (ekin + pot) / real(par%N, dp)
  
    ! --- accumulators (canonical: energy moments) ---
    obs%av_mag  = obs%av_mag  + abs(mag) / real(par%N, dp)
    obs%av_mag2 = obs%av_mag2 + mag*mag
    obs%av_mag4 = obs%av_mag4 + mag*mag*mag*mag
  
    obs%av_pot  = obs%av_pot + pot / real(par%N, dp)
    obs%av_kin  = obs%av_kin + ekin / real(par%N, dp)
  
    obs%av_en   = obs%av_en  + e
    obs%av_en2  = obs%av_en2 + e*e
  
    ! --- running means ---
    obs%c_mag  = obs%av_mag  / real(n_MC, dp)
    obs%c_mag2 = obs%av_mag2 / real(n_MC, dp)
    obs%c_mag4 = obs%av_mag4 / real(n_MC, dp)
    obs%binder = 1.0_dp - obs%c_mag4/(3.0_dp*obs%c_mag2*obs%c_mag2)
  
    obs%c_pot = obs%av_pot / real(n_MC, dp)
    obs%c_kin = obs%av_kin / real(n_MC, dp)
  
    obs%c_en  = obs%av_en  / real(n_MC, dp)
    obs%c_en2 = obs%av_en2 / real(n_MC, dp)
  
    ! canonical Cv per site
    obs%cv = (par%beta_mc*par%beta_mc) * real(par%N, dp) * (obs%c_en2 - obs%c_en*obs%c_en)
  
    ! micro-only fields: keep harmlessly zero
    obs%c_kin_1 = 0.0_dp
    obs%c_kin_2 = 0.0_dp
    obs%c_kin_3 = 0.0_dp
    obs%beta    = 0.0_dp
    obs%ders_2  = 0.0_dp
    obs%ders_3  = 0.0_dp
  
    obs%n_MC = n_MC
  end subroutine update_obs_canon


  subroutine reinflate_micro(obs)
    type(Phi4Obs), intent(inout) :: obs
    real(dp) :: n
  
    n = real(obs%n_MC, dp)
    if (obs%n_MC <= 0) return
  
    obs%av_mag   = obs%c_mag  * n
    obs%av_mag2  = obs%c_mag2 * n
    obs%av_mag4  = obs%c_mag4 * n
  
    obs%av_pot   = obs%c_pot  * n
    obs%av_kin   = obs%c_kin  * n
  
    obs%av_kin_1 = obs%c_kin_1 * n
    obs%av_kin_2 = obs%c_kin_2 * n
    obs%av_kin_3 = obs%c_kin_3 * n
  
    obs%av_en    = obs%c_en   * n
  
    ! binder/beta/derivatives will be recomputed by update_obs_micro
  end subroutine reinflate_micro


  subroutine reinflate_canon(obs)
    type(Phi4Obs), intent(inout) :: obs
    real(dp) :: n
  
    n = real(obs%n_MC, dp)
    if (obs%n_MC <= 0) return
  
    obs%av_mag   = obs%c_mag  * n
    obs%av_mag2  = obs%c_mag2 * n
    obs%av_mag4  = obs%c_mag4 * n
  
    obs%av_pot   = obs%c_pot  * n
    obs%av_kin   = obs%c_kin  * n
  
    obs%av_en    = obs%c_en   * n
    obs%av_en2   = obs%c_en2  * n
  
    ! cv will be recomputed by update_obs_canon
  end subroutine reinflate_canon




end module mod_io_observables_phi4