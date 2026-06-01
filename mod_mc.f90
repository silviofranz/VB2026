!======================================================================
! mod_mc: Monte Carlo core — THE PARALLELIZATION TARGET.
!
! Contains the energy-change functions and the Metropolis sweep.
! All data is owned by other modules; this module is pure computation.
!
! For OpenMP parallelisation over disorder samples, the arrays owned
! by mod_spins, mod_replica, mod_graph, and mod_random will need to
! be made thread-private (or passed as arguments).
!
! Depends on: mod_params, mod_random, mod_graph, mod_spins, mod_replica.
!======================================================================
module mod_mc
  use mod_params
  use mod_random,  only: idummy, ran3
  use mod_graph,   only: nn, coup, conn, nn1, coup1
  use mod_spins,   only: s, s1, ener, ener1, rf, oh, oh1, oh2
  use mod_replica, only: beta, perm, exp_table
  implicit none

contains

  !--------------------------------------------------------------------
  ! Energy change for flipping spin x in replica r (copy s).
  ! dE = 2 * s(r,x) * local_field(r,x)
  !--------------------------------------------------------------------
  pure function dh(r, x) result(de)
    integer, intent(in) :: r, x
    integer :: de
    integer :: y

    de = 0
    do y = 1, conn(x)
      de = de + coup(x,y) * s(r, nn(x,y))
    end do
    de = 2 * s(r,x) * de
  end function dh

  !--------------------------------------------------------------------
  ! Energy change for flipping spin x in replica r (copy s1).
  !--------------------------------------------------------------------
  pure function dh1(r, x) result(de)
    integer, intent(in) :: r, x
    integer :: de
    integer :: y

    de = 0
    do y = 1, conn(x)
      de = de + coup1(x,y) * s1(r, nn1(x,y))
    end do
    de = 2 * s1(r,x) * de
  end function dh1

  !--------------------------------------------------------------------
  ! One Metropolis sweep: for each replica, attempt l random single-spin
  ! flips on both copies s and s1 independently.
  !--------------------------------------------------------------------
  subroutine monte()
    integer :: r, y, x
    real    :: rr
    integer :: ddh

    do r = 1, rep
      do y = 1, l
        rr = ran3(idummy)
        x  = min(int(rr*l) + 1, l)

        ! --- copy s: flip and update energy + overlap ---
        ddh = dh(r, x)
        if (ddh <= 0 .or. ran3(idummy) < exp_table(r, ddh)) then
          ! update overlap before flipping sign
          if (x <= l/2) then
            oh1(r) = oh1(r) - 2 * s(r,x) * s1(r,x)
          else
            oh2(r) = oh2(r) - 2 * s(r,x) * s1(r,x)
          end if
          oh(r)   = oh1(r) + oh2(r)
          s(r,x)  = -s(r,x)
          ener(r) = ener(r) + real(ddh)
        end if

        ! --- copy s1: flip and update energy + overlap ---
        ddh = dh1(r, x)
        if (ddh <= 0 .or. ran3(idummy) < exp_table(r, ddh)) then
          ! update overlap before flipping sign
          if (x <= l/2) then
            oh1(r) = oh1(r) - 2 * s(r,x) * s1(r,x)
          else
            oh2(r) = oh2(r) - 2 * s(r,x) * s1(r,x)
          end if
          oh(r)    = oh1(r) + oh2(r)
          s1(r,x)  = -s1(r,x)
          ener1(r) = ener1(r) + real(ddh)
        end if
      end do
    end do
  end subroutine monte

end module mod_mc
