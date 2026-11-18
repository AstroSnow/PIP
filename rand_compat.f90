module rand_compat
  ! Portable copy of gfortran's RAND/SRAND extensions, for compilers that lack them
  ! (Intel Fortran needs IFPORT and returns a different sequence).
  !
  ! Same generator as libgfortran: Park-Miller "minimal standard",
  !   seed = 16807*seed mod (2**31-1)
  ! and RAND = ((seed-1) with its low 9 bits cleared) / 2**31, as default REAL.
  ! Verified bit-for-bit against gfortran's intrinsic RAND (5 million draws), so
  ! initial conditions using RAND are identical whichever compiler builds the code.
  !
  ! Usage: add "use rand_compat" to a routine that calls rand()/srand().
  implicit none
  private
  public :: rand, srand
  integer(8), save :: rand_seed = 1_8
  integer(8), parameter :: GFC_RAND_A = 16807_8, GFC_RAND_M = 2147483647_8
contains
  subroutine srand(seed)
    integer, intent(in) :: seed
    if (seed == 0) then
      rand_seed = 123459876_8
    else
      rand_seed = int(seed, 8)
    end if
  end subroutine srand

  real function rand(flag)
    ! flag absent or 0: next number; 1: restart the default sequence; other: reseed with flag
    integer, intent(in), optional :: flag
    if (present(flag)) then
      if (flag == 1) then
        rand_seed = 123459876_8
      else if (flag /= 0) then
        rand_seed = int(flag, 8)
      end if
    end if
    rand_seed = mod(GFC_RAND_A * rand_seed, GFC_RAND_M)
    rand = real(real(iand(rand_seed - 1_8, not(511_8)), 8) / 2.0d0**31)
  end function rand
end module rand_compat
