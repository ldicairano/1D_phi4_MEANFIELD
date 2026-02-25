module mod_io_append_dispatch
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4Obs, IOParams
  implicit none
  private
  public :: append_iface, select_append
  public :: append_observables_micro, append_observables_canon

  abstract interface
  subroutine append_iface(io, par, obs, time_tot, n_MC)
    use mod_kinds,      only: dp
    use mod_types_phi4, only: IOParams, Phi4Params, Phi4Obs
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4Obs),    intent(in) :: obs
    real(dp),         intent(in) :: time_tot
    integer,          intent(in) :: n_MC
  end subroutine append_iface
end interface

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


  subroutine select_append(par, append)
    type(Phi4Params), intent(in) :: par
    procedure(append_iface), pointer, intent(out) :: append

    select case (trim(par%ensemble))
    case ("micro")
      append => append_observables_micro
    case ("canon")
      append => append_observables_canon
    case default
      write(*,*) "ERROR: ensemble must be 'micro' or 'canon' (select_append)."
      stop
    end select
  end subroutine select_append


  subroutine append_observables_micro(io, par, obs, time_tot, n_MC)
    use mod_kinds,      only: dp
    use mod_types_phi4, only: IOParams, Phi4Params, Phi4Obs
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4Obs),    intent(in) :: obs
    real(dp),         intent(in) :: time_tot
    integer,          intent(in) :: n_MC
    character(len=256) :: fname
    integer :: u

    ! Columns:
    !  1 time_tot
    !  2 n_MC
    !  3 en_in
    !  4 delta
    !  5 c_en
    !  6 c_en2
    !  7 cv
    !  8 c_kin
    !  9 c_pot
    ! 10 c_kin_1
    ! 11 c_kin_2
    ! 12 c_kin_3
    ! 13 beta
    ! 14 ders_2
    ! 15 ders_3
    ! 16 c_mag
    ! 17 c_mag2
    ! 18 c_mag4
    ! 19 binder

    call ensure_dir(io%out_dir)

    write(fname,'(A,"/obs_micro_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
    trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="unknown", action="write")
  
    write(u,'(ES14.6,1X,I10,1X,ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
      time_tot, n_MC, par%en_in, par%delta, &
      obs%c_en, &
      obs%c_kin, obs%c_pot, &
      obs%c_kin_1, obs%c_kin_2, obs%c_kin_3, &
      obs%beta, obs%ders_2, obs%ders_3, &
      obs%c_mag, obs%c_mag2, obs%c_mag4, obs%binder
  
    close(u)
  end subroutine append_observables_micro

  subroutine append_observables_canon(io, par, obs, time_tot, n_MC)
    use mod_kinds,      only: dp
    use mod_types_phi4, only: IOParams, Phi4Params, Phi4Obs
    implicit none
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4Obs),    intent(in) :: obs
    real(dp),         intent(in) :: time_tot
    integer,          intent(in) :: n_MC
  
    character(len=256) :: fname
    integer :: u

    call ensure_dir(io%out_dir)

    ! Same naming convention as micro: obs_<n_samp>_<n_realiz>_<n_rest>_<N>.dat
    write(fname,'(A,"/obs_canon_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
    trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="unknown", action="write")
  
    ! Columns (canonical):
    !  1 time_tot
    !  2 n_MC
    !  3 beta_mc
    !  4 delta
    !  5 c_en
    !  6 c_en2
    !  7 cv
    !  8 c_kin
    !  9 c_pot
    ! 10 c_mag
    ! 11 c_mag2
    ! 12 c_mag4
    ! 13 binder
    write(u,'(ES14.6,1X,I10,1X,ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,'// &
            'ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6)') &
      time_tot, n_MC, par%beta_mc, par%delta, &
      obs%c_en, obs%c_en2, obs%cv, &
      obs%c_kin, obs%c_pot, &
      obs%c_mag, obs%c_mag2, obs%c_mag4, obs%binder
  
    close(u)
  end subroutine append_observables_canon

end module mod_io_append_dispatch
