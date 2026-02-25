module mod_cli_phi4
  implicit none
contains

  subroutine parse_cli_phi4(input_file, n_samp, L, n_steps, n_jump, n_realiz, n_rest)
    character(len=*), intent(inout) :: input_file
    integer, intent(out) :: n_samp, L, n_steps, n_jump, n_realiz, n_rest

    integer :: nargs, ioff
    character(len=256) :: arg

    input_file = "input.inp"
    nargs = command_argument_count()

    if (nargs < 6) then
      call print_usage()
      stop 1
    end if

    ioff = 0
    call get_command_argument(1, arg)
    if (trim(arg) == "-i") then
      if (nargs < 8) then
        write(*,*) "ERROR: with -i you must provide:"
        write(*,*) "  ./phi4_mmc -i input.inp n_samp L n_steps n_jump n_realiz n_rest"
        stop 1
      end if
      call get_command_argument(2, input_file)
      ioff = 2
    end if

    ! now the 6 mandatory numeric args follow (either from 1..6 or 3..8)
    call get_command_argument(1+ioff, arg); read(arg,*) n_samp
    call get_command_argument(2+ioff, arg); read(arg,*) L
    call get_command_argument(3+ioff, arg); read(arg,*) n_steps
    call get_command_argument(4+ioff, arg); read(arg,*) n_jump
    call get_command_argument(5+ioff, arg); read(arg,*) n_realiz
    call get_command_argument(6+ioff, arg); read(arg,*) n_rest

    if (L <= 0) stop "ERROR: L must be > 0"
    if (n_steps <= 0) stop "ERROR: n_steps must be > 0"
    if (n_jump <= 0) stop "ERROR: n_jump must be > 0"
    if (n_realiz <= 0) stop "ERROR: n_realiz must be > 0"
    if (n_rest < 0) stop "ERROR: n_rest must be >= 0"
  end subroutine parse_cli_phi4

  subroutine print_usage()
    write(*,*) "Usage:"
    write(*,*) "  ./phi4_mmc -i input.inp n_samp L n_steps n_jump n_realiz n_rest"
    write(*,*) "or"
    write(*,*) "  ./phi4_mmc n_samp L n_steps n_jump n_realiz n_rest   (uses input.inp)"
    write(*,*) ""
    write(*,*) "Meaning (compatible with legacy code):"
    write(*,*) "  n_samp  : energy index -> en_in = 0.01*n_samp ; Etot = N*en_in"
    write(*,*) "  L       : lattice size (N = L here)"
    write(*,*) "  n_steps : total MMC steps"
    write(*,*) "  n_jump  : measurement stride"
    write(*,*) "  n_realiz: number of realizations"
    write(*,*) "  n_rest  : restart_every (0 disables)"
  end subroutine print_usage


  subroutine cli_msg_restart_already_complete(par, step0)
    use, intrinsic :: iso_fortran_env, only: output_unit
    use mod_types_phi4, only: Phi4Params
    implicit none
    type(Phi4Params), intent(in) :: par
    integer,          intent(in) :: step0
  
    write(output_unit,'(A)') "------------------------------------------------------------"
    write(output_unit,'(A)') "INFO: restart requested but the run is already complete."
    write(output_unit,'(A,A)') "  ensemble              = ", trim(par%ensemble)
    write(output_unit,'(A,I0)') "  n_rest (requested)    = ", par%n_rest
    write(output_unit,'(A,I0)') "  step0 (from restart)  = ", step0
    write(output_unit,'(A,I0)') "  n_steps (input max)   = ", par%n_steps
    write(output_unit,'(A)') "  Action: exiting now (nothing to simulate)."
    write(output_unit,'(A)') "------------------------------------------------------------"
    flush(output_unit)
  end subroutine cli_msg_restart_already_complete

end module mod_cli_phi4
