module mod_io_restart_phi4
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State, Phi4Obs, IOParams, RNGState
  implicit none
  private
  public :: read_restart, write_restart
  public :: read_last_nmc_from_dat
  integer, parameter :: RESTART_VERSION = 20260225

contains

subroutine ensure_dir(path)
  character(len=*), intent(in) :: path
  integer :: istat
  if (len_trim(path) == 0) return
  call execute_command_line("mkdir -p " // trim(path), exitstat=istat)
  if (istat /= 0) then
    write(*,*) "ERROR: cannot create directory: ", trim(path)
    stop
  end if
end subroutine ensure_dir

  subroutine restart_fname(io, par, ir, nrest, fname)
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    integer,          intent(in) :: ir, nrest
    character(len=*), intent(out):: fname
    character(len=256) :: base

    ! restart_<ensemble>_<n_samp>_<ir>_<N>_<n_rest>.bin
    write(base,'("restart_",A,"_",I0,"_",I0,"_",I0,"_",I0,".bin")') trim(par%ensemble), par%n_samp, ir, par%N, nrest
    fname = trim(io%restart_dir)//"/"//trim(base)
  end subroutine restart_fname


  subroutine write_restart(io, par, st, obs, rng, step, ir, nrest_out)
    integer, intent(in), optional :: nrest_out
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4State),  intent(in) :: st
    type(Phi4Obs),    intent(in) :: obs
    type(RNGState),   intent(in) :: rng
    integer,          intent(in) :: step
    integer,          intent(in) :: ir

    integer :: nrest_write
    integer :: u
    character(len=512) :: fname

    nrest_write = par%n_rest
    if (present(nrest_out)) nrest_write = nrest_out

    call ensure_dir(io%restart_dir)
    call restart_fname(io, par, ir, nrest_write, fname)

    open(newunit=u, file=trim(fname), access="stream", form="unformatted", status="replace", action="write")

    write(u) RESTART_VERSION

    ! Minimal params needed to guarantee compatibility
    write(u) par%L, par%N
    write(u) par%coup, par%mu, par%lambda
    write(u) par%ensemble
    write(u) par%Etot, par%beta_mc
    write(u) par%delta
    write(u) par%n_samp, par%en_in
    write(u) nrest_write

    write(u) step
    write(u) rng%seed

    ! You said you prefer small restart: store only n_MC and c_* (means)
    write(u) obs%n_acc, obs%n_att, obs%acc_rate
    write(u) obs%n_MC

    write(u) obs%c_mag, obs%c_mag2, obs%c_mag4, obs%binder
    write(u) obs%c_pot, obs%c_kin, obs%c_kin_1, obs%c_kin_2, obs%c_kin_3
    write(u) obs%c_en,  obs%c_en2
    write(u) obs%beta, obs%ders_2, obs%ders_3

    ! State
    write(u) st%V, st%K
    write(u) st%S1, st%S2, st%S4
    write(u) st%M, st%phi2, st%phi4
    write(u) st%phi
    if (allocated(st%pi)) write(u) st%pi

    close(u)
  end subroutine write_restart


  subroutine read_restart(io, par, st, obs, rng, step0, ir)
    type(IOParams),   intent(in)    :: io
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(Phi4Obs),    intent(inout) :: obs
    type(RNGState),   intent(inout) :: rng
    integer,          intent(out)   :: step0
    integer,          intent(in)    :: ir

    integer :: u, ios, ver
    integer :: Lr, Nr, nrestr
    real(dp) :: coupr, mur, lambdar
    character(len=16) :: ensr
    real(dp) :: Etotr, betamcr, deltar
    integer  :: nsampr
    real(dp) :: eninr
    character(len=512) :: fname

    step0 = 0

    ! Automatic rule:
    ! - n_rest == 1 : no restart file is read
    ! - n_rest >  1 : read from n_rest-1
    if (par%n_rest <= 1) return

    call restart_fname(io, par, ir, par%n_rest-1, fname)

    open(newunit=u, file=trim(fname), access="stream", form="unformatted", status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: restart file not found for n_rest-1."
      write(*,*) "  expected: ", trim(fname)
      stop
    end if

    read(u) ver
    if (ver /= RESTART_VERSION) then
      write(*,*) "ERROR: restart version mismatch. file=", ver, " expected=", RESTART_VERSION
      stop
    end if

    read(u) Lr, Nr
    read(u) coupr, mur, lambdar
    read(u) ensr
    read(u) Etotr, betamcr
    read(u) deltar
    read(u) nsampr, eninr
    read(u) nrestr

    ! Basic consistency checks
    if (nrestr /= par%n_rest-1) then
      write(*,*) "ERROR: restart mismatch (n_rest). file:", nrestr, " expected:", par%n_rest-1
      stop
    end if
    if (Lr /= par%L .or. Nr /= par%N) then
      write(*,*) "ERROR: restart mismatch (L,N). file:", Lr, Nr, " input:", par%L, par%N
      stop
    end if
    if (abs(coupr-par%coup) > 1.0e-12_dp .or. abs(mur-par%mu) > 1.0e-12_dp .or. abs(lambdar-par%lambda) > 1.0e-12_dp) then
      write(*,*) "ERROR: restart mismatch (model params)."
      stop
    end if
    if (trim(ensr) /= trim(par%ensemble)) then
      write(*,*) "ERROR: restart mismatch (ensemble)."
      stop
    end if

    ! Restore runtime parameters that must continue seamlessly
    par%Etot    = Etotr
    par%beta_mc = betamcr
    par%delta   = deltar
    par%n_samp  = nsampr
    par%en_in   = eninr
    ! par%n_rest stays the current segment (we do not overwrite it)

    read(u) step0
    read(u) rng%seed

    read(u) obs%n_acc, obs%n_att, obs%acc_rate
    read(u) obs%n_MC

    read(u) obs%c_mag, obs%c_mag2, obs%c_mag4, obs%binder
    read(u) obs%c_pot, obs%c_kin, obs%c_kin_1, obs%c_kin_2, obs%c_kin_3
    read(u) obs%c_en,  obs%c_en2
    read(u) obs%beta, obs%ders_2, obs%ders_3

    read(u) st%V, st%K
    read(u) st%S1, st%S2, st%S4
    read(u) st%M, st%phi2, st%phi4
    read(u) st%phi
    if (allocated(st%pi)) read(u) st%pi

    close(u)

    ! Ensure normalized macros match sums
    st%M    = st%S1 / real(par%N, dp)
    st%phi2 = st%S2 / real(par%N, dp)
    st%phi4 = st%S4 / real(par%N, dp)

  end subroutine read_restart

subroutine read_last_nmc_from_dat(io, par, ir, nrest, nmc_last, found)
  use mod_kinds,      only: dp
  use mod_types_phi4, only: IOParams, Phi4Params
  implicit none
  type(IOParams),   intent(in)  :: io
  type(Phi4Params), intent(in)  :: par
  integer,          intent(in)  :: ir, nrest
  integer,          intent(out) :: nmc_last
  logical,          intent(out) :: found

  character(len=512) :: fname, line
  integer :: u, ios, ios2
  real(dp) :: t_dummy
  integer  :: nmc_tmp

  found    = .false.
  nmc_last = 0

  ! Nome IDENTICO a quello usato in append (vedi mod_io_append_dispatch.f90)
  if (trim(par%ensemble) == "micro") then
    write(fname,'(A,"/obs_micro_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, nrest, par%N
  else if (trim(par%ensemble) == "canon") then
    write(fname,'(A,"/obs_canon_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, nrest, par%N
  else
    return
  end if

  open(newunit=u, file=trim(fname), status="old", action="read", iostat=ios)
  if (ios /= 0) return

  do
    read(u,'(A)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle

    ! Nel .dat la 2a colonna è SEMPRE n_MC (sia micro che canon)
    read(line, *, iostat=ios2) t_dummy, nmc_tmp
    if (ios2 == 0) then
      nmc_last = nmc_tmp
      found = .true.
    end if
  end do

  close(u)
end subroutine read_last_nmc_from_dat

end module mod_io_restart_phi4
