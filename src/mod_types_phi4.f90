module mod_types_phi4
  use mod_kinds, only: dp
  implicit none

  type :: Phi4Params
    integer :: L = 0
    integer :: N = 0

    ! --- legacy sampling id (for obs_*.dat naming) ---
    integer  :: n_samp = 0
    real(dp) :: en_in  = 0.0_dp   ! energy per site (legacy)
    real(dp) :: Etot   = 0.0_dp   ! total energy (micro)

    ! --- run control ---
    integer :: n_steps  = 0
    integer :: n_jump   = 1
    integer :: n_therm  = 0
    integer :: n_realiz = 1   ! trajectory id (CLI)
    integer :: n_rest = 1     ! segment id (CLI): 1=fresh, >1 restart from n_rest-1
    integer :: n_print = 1000   ! print/append frequency in units of measurements (n_MC)
    ! Time / walltime control (from &system)
    real(dp) :: cluster_time = 172800.0_dp   ! walltime in seconds (default 48h)

    ! --- model params ---
    real(dp) :: coup   = 1.0_dp
    real(dp) :: mu     = 1.0_dp   ! 2.0_dp
    real(dp) :: lambda = 6.0_dp   ! 3.0_dp/5.0_dp
    ! If .true., compute coup(J) from N after CLI overrides (see main).
    logical  :: auto_coup = .false.
    
    ! --- ensemble ---
    character(len=16) :: ensemble = "micro"  ! "micro" or "canon"
    real(dp) :: beta_mc = 0.0_dp             ! canonical input inverse temperature

    real(dp) :: delta = 0.1_dp   ! proposal step; will be tuned in init

    ! --- RNG ---
    integer :: seed = 12345
  end type Phi4Params


  type :: Phi4State
    real(dp), allocatable :: phi(:)
    real(dp), allocatable :: pi(:)   ! canonical momenta (allocated when needed)

    integer, allocatable :: index(:)  ! array for random selection in monte carlo loop

    real(dp) :: V = 0.0_dp
    real(dp) :: K = 0.0_dp   ! kinetic energy (canonical) or K=Etot-V (micro)

    ! Raw sums for O(1) macro updates (phi only)
    real(dp) :: S1 = 0.0_dp
    real(dp) :: S2 = 0.0_dp
    real(dp) :: S4 = 0.0_dp

    ! Normalized macros (phi only)
    real(dp) :: M    = 0.0_dp
    real(dp) :: phi2 = 0.0_dp
    real(dp) :: phi4 = 0.0_dp
  end type Phi4State


  type :: Phi4Obs
    integer :: n_acc = 0
    integer :: n_att = 0
    real(dp) :: acc_rate = 0.0_dp

    integer :: n_MC = 0

    real(dp) :: av_mag   = 0.0_dp
    real(dp) :: av_mag2  = 0.0_dp
    real(dp) :: av_mag4  = 0.0_dp

    real(dp) :: av_pot   = 0.0_dp
    real(dp) :: av_kin   = 0.0_dp
    real(dp) :: av_kin_1 = 0.0_dp
    real(dp) :: av_kin_2 = 0.0_dp
    real(dp) :: av_kin_3 = 0.0_dp

    real(dp) :: av_en    = 0.0_dp
    real(dp) :: av_en2   = 0.0_dp

    real(dp) :: c_mag   = 0.0_dp
    real(dp) :: c_mag2  = 0.0_dp
    real(dp) :: c_mag4  = 0.0_dp
    real(dp) :: binder  = 0.0_dp

    real(dp) :: c_pot   = 0.0_dp
    real(dp) :: c_kin   = 0.0_dp
    real(dp) :: c_kin_1 = 0.0_dp
    real(dp) :: c_kin_2 = 0.0_dp
    real(dp) :: c_kin_3 = 0.0_dp

    real(dp) :: c_en    = 0.0_dp
    real(dp) :: c_en2   = 0.0_dp
    real(dp) :: cv      = 0.0_dp

    ! Microcanonical observable beta and derivatives (from kinetic moments)
    real(dp) :: beta   = 0.0_dp
    real(dp) :: ders_2 = 0.0_dp
    real(dp) :: ders_3 = 0.0_dp

    real(dp) :: av_m   = 0.0_dp
    real(dp) :: av_m2  = 0.0_dp
    real(dp) :: c_m    = 0.0_dp
    real(dp) :: c_m2   = 0.0_dp
    real(dp) :: chi    = 0.0_dp

  end type Phi4Obs


  type :: IOParams
    character(len=256) :: out_dir = "out"
    character(len=256) :: restart_dir = "."
  end type IOParams


  type :: RNGState
    integer :: seed = 12345
  end type RNGState

end module mod_types_phi4