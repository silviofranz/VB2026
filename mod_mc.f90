!======================================================================
! mod_mc: 16-replica 16-bit Bit-Sliced Monte Carlo core.
!======================================================================
module mod_mc
  use mod_params
  use mod_random,  only: idummy, fast_ran_64
  use mod_graph,   only: nn, coup, conn, nn1, coup1
  use mod_spins,   only: ms_s, ms_s1, update_all_metrics
  use mod_replica, only: T_slice
  implicit none

contains

  subroutine monte()
    integer :: x, y, s, h, i, de_int, z_i
    integer(8) :: h0, h1, h2, h3, h4, carry, contrib
    integer(8) :: R(0:15), s_mask, h_mask, C, eq, lt, flipped, flip_mask
    integer(8) :: m_bit(0:4, 0:1)

    ! UPDATING ms_s
    do x = 1, l
       h0=0; h1=0; h2=0; h3=0; h4=0
       z_i = conn(x)
       do y = 1, z_i
          if (coup(x, y) == 1) then
             contrib = ms_s(nn(x, y))
          else
             contrib = NOT(ms_s(nn(x, y)))
          end if
          carry = h0; h0 = ieor(h0, contrib); contrib = iand(carry, contrib)
          carry = h1; h1 = ieor(h1, contrib); contrib = iand(carry, contrib)
          carry = h2; h2 = ieor(h2, contrib); contrib = iand(carry, contrib)
          carry = h3; h3 = ieor(h3, contrib); contrib = iand(carry, contrib)
          h4 = ieor(h4, contrib)
       end do

       do i = 0, 15
          R(i) = fast_ran_64(idummy)
       end do

       m_bit(0, 1) = h0; m_bit(0, 0) = NOT(h0)
       m_bit(1, 1) = h1; m_bit(1, 0) = NOT(h1)
       m_bit(2, 1) = h2; m_bit(2, 0) = NOT(h2)
       m_bit(3, 1) = h3; m_bit(3, 0) = NOT(h3)
       m_bit(4, 1) = h4; m_bit(4, 0) = NOT(h4)

       flip_mask = 0_8
       do s = 0, 1
          if (s == 1) then
             s_mask = ms_s(x)
          else
             s_mask = NOT(ms_s(x))
          end if
          s_mask = iand(s_mask, 65535_8) ! Mask for 16 replicas
          if (s_mask == 0_8) cycle

          do h = 0, z_i
             h_mask = iand(s_mask, m_bit(0, iand(h, 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(1, iand(ishft(h, -1), 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(2, iand(ishft(h, -2), 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(3, iand(ishft(h, -3), 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(4, iand(ishft(h, -4), 1)))
             if (h_mask == 0_8) cycle

             de_int = 4 * h - 2 * z_i
             if (s == 0) de_int = -de_int

             if (de_int <= 0) then
                flipped = h_mask
             else
                C = 0_8
                do i = 0, 15
                   eq = NOT(ieor(R(i), T_slice(i, de_int)))
                   lt = iand(NOT(R(i)), T_slice(i, de_int))
                   C = ior(lt, iand(eq, C))
                end do
                flipped = iand(h_mask, C)
             end if
             flip_mask = ior(flip_mask, flipped)
          end do
       end do
       ms_s(x) = ieor(ms_s(x), flip_mask)
    end do

    ! UPDATING ms_s1
    do x = 1, l
       h0=0; h1=0; h2=0; h3=0; h4=0
       z_i = conn(x)
       do y = 1, z_i
          if (coup1(x, y) == 1) then
             contrib = ms_s1(nn1(x, y))
          else
             contrib = NOT(ms_s1(nn1(x, y)))
          end if
          carry = h0; h0 = ieor(h0, contrib); contrib = iand(carry, contrib)
          carry = h1; h1 = ieor(h1, contrib); contrib = iand(carry, contrib)
          carry = h2; h2 = ieor(h2, contrib); contrib = iand(carry, contrib)
          carry = h3; h3 = ieor(h3, contrib); contrib = iand(carry, contrib)
          h4 = ieor(h4, contrib)
       end do

       do i = 0, 15
          R(i) = fast_ran_64(idummy)
       end do

       m_bit(0, 1) = h0; m_bit(0, 0) = NOT(h0)
       m_bit(1, 1) = h1; m_bit(1, 0) = NOT(h1)
       m_bit(2, 1) = h2; m_bit(2, 0) = NOT(h2)
       m_bit(3, 1) = h3; m_bit(3, 0) = NOT(h3)
       m_bit(4, 1) = h4; m_bit(4, 0) = NOT(h4)

       flip_mask = 0_8
       do s = 0, 1
          if (s == 1) then
             s_mask = ms_s1(x)
          else
             s_mask = NOT(ms_s1(x))
          end if
          s_mask = iand(s_mask, 65535_8)
          if (s_mask == 0_8) cycle

          do h = 0, z_i
             h_mask = iand(s_mask, m_bit(0, iand(h, 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(1, iand(ishft(h, -1), 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(2, iand(ishft(h, -2), 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(3, iand(ishft(h, -3), 1)))
             if (h_mask == 0_8) cycle
             h_mask = iand(h_mask, m_bit(4, iand(ishft(h, -4), 1)))
             if (h_mask == 0_8) cycle

             de_int = 4 * h - 2 * z_i
             if (s == 0) de_int = -de_int

             if (de_int <= 0) then
                flipped = h_mask
             else
                C = 0_8
                do i = 0, 15
                   eq = NOT(ieor(R(i), T_slice(i, de_int)))
                   lt = iand(NOT(R(i)), T_slice(i, de_int))
                   C = ior(lt, iand(eq, C))
                end do
                flipped = iand(h_mask, C)
             end if
             flip_mask = ior(flip_mask, flipped)
          end do
       end do
       ms_s1(x) = ieor(ms_s1(x), flip_mask)
    end do

    call update_all_metrics()

  end subroutine monte

end module mod_mc
