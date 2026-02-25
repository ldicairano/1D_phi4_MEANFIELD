module mod_random
  use mod_kinds,      only: dp
  use mod_types_phi4, only: RNGState
  implicit none
contains

  subroutine rng_init(rng, seed_in)
    type(RNGState), intent(inout) :: rng
    integer, intent(in) :: seed_in
    rng%seed = max(1, seed_in)
  end subroutine rng_init

  function rng_uniform(rng) result(u)
    type(RNGState), intent(inout) :: rng
    real(dp) :: u
    integer, parameter :: a = 1103515245, c = 12345
    integer :: x
    rng%seed = a*rng%seed + c
    x = iand(rng%seed, int(z'7fffffff'))
    u = real(x, dp) / real(int(z'7fffffff'), dp)
    if (u <= 0.0_dp) u = 1.0e-16_dp
  end function rng_uniform

  function rng_gauss(rng) result(g)
    type(RNGState), intent(inout) :: rng
    real(dp) :: g
    real(dp) :: u1, u2
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    g = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*acos(-1.0_dp)*u2)
  end function rng_gauss

end module mod_random
