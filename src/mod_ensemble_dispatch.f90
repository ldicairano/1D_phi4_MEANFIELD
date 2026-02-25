module mod_ensemble_dispatch
  use mod_types_phi4,      only: Phi4Params
  use mod_metropolis_phi4, only: sweep_iface, mmc_sweep_canon, mmc_sweep_micro
  implicit none
  private
  public :: select_sweep

contains

  subroutine select_sweep(par, sweep)
    type(Phi4Params), intent(in) :: par
    procedure(sweep_iface), pointer, intent(out) :: sweep

    select case (adjustl(trim(par%ensemble)))
    case ("micro")
      sweep => mmc_sweep_micro
    case ("canon")
      sweep => mmc_sweep_canon
    case default
      write(*,*) "ERROR: unknown ensemble = ", trim(par%ensemble)
      stop
    end select
  end subroutine select_sweep

end module mod_ensemble_dispatch