!======================================================================
! mod_spins: spin configurations, energies, random fields — OpenMP-safe.
! All physics arrays are threadprivate: each thread works on its own
! independent spin configuration.
!======================================================================
module mod_spins
  use mod_params
  use mod_random, only: idummy1, ran3
  use mod_graph,  only: nn, coup, conn, nn1, coup1
  implicit none

  integer, allocatable :: s(:,:)        ! primary spin copy
  integer, allocatable :: s1(:,:)       ! secondary spin copy (for overlap)
  real, allocatable    :: ener(:)       ! energy of copy s
  real, allocatable    :: ener1(:)      ! energy of copy s1
  real, allocatable    :: rf(:)         ! quenched random fields
  real    :: sigma       = 0.0   ! random-field amplitude
  real    :: sigma_campo = 0.0   ! legacy (output only)
  real    :: aux                 ! auxiliary accumulator

  ! Running overlap sums (maintained incrementally — O(1) update per spin flip)
  integer, allocatable :: oh(:)         ! total overlap = oh1 + oh2
  integer, allocatable :: oh1(:)        ! left-half overlap  sum_i s(r,i)*s1(r,i), i<=l/2
  integer, allocatable :: oh2(:)        ! right-half overlap sum_i s(r,i)*s1(r,i), i>l/2

  !$omp threadprivate(s, s1, ener, ener1, rf, aux, oh, oh1, oh2)
  ! sigma, sigma_campo: set once in main, shared read-only

contains

  subroutine alloc_spins()
    if (.not. allocated(s)) then
      allocate(s(rep,l), s1(rep,l))
      allocate(ener(rep), ener1(rep), rf(l))
      allocate(oh(rep), oh1(rep), oh2(rep))
    end if
  end subroutine alloc_spins

  subroutine init_s()
    integer :: i, r
    do i = 1, l
      rf(i) = 1.0
      if (ran3(idummy1) >= 0.5) rf(i) = -1.0
      rf(i) = sigma * rf(i)
    end do
    do r = 1, rep
      do i = 1, l
        s(r,i)  = 1;  if (ran3(idummy1) >= 0.5) s(r,i)  = -1
        s1(r,i) = 1;  if (ran3(idummy1) >= 0.5) s1(r,i) = -1
      end do
      ener(r)  = 0.0
      ener1(r) = 0.0
    end do
    call energia()
  end subroutine init_s

  ! Compute oh1, oh2, oh from scratch — call after init_s or after a rewiring step
  subroutine init_overlaps()
    integer :: i, r
    do r = 1, rep
      oh1(r) = 0;  oh2(r) = 0
      do i = 1, l/2
        oh1(r) = oh1(r) + s(r,i) * s1(r,i)
      end do
      do i = l/2+1, l
        oh2(r) = oh2(r) + s(r,i) * s1(r,i)
      end do
      oh(r) = oh1(r) + oh2(r)
    end do
  end subroutine init_overlaps

  subroutine energia()
    integer :: i, n, k, r
    real    :: hh, hh1
    do r = 1, rep
      ener(r) = 0.0;  ener1(r) = 0.0
    end do
    do r = 1, rep
      do i = 1, l
        hh = 0.0;  hh1 = 0.0
        do n = 1, conn(i)
          hh  = hh  + coup(i,n) * s(r,  nn(i,n))
          hh1 = hh1 + coup1(i,n) * s1(r, nn1(i,n))
        end do
        ener(r)  = ener(r)  + s(r,i)  * hh  * 0.5 + s(r,i)  * rf(i)
        ener1(r) = ener1(r) + s1(r,i) * hh1 * 0.5 + s1(r,i) * rf(i)
      end do
      ener(r)  = -ener(r)
      ener1(r) = -ener1(r)
    end do
  end subroutine energia

end module mod_spins
