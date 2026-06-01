!======================================================================
! main.f90 (TER version)
!
! Progressive rewiring: for EACH disorder sample, the two graphs start
! identical and graph2 (nn1/coup1) is rewired ONE step at a time.
! The MC is run between each rewiring step so we track how the overlap
! between s (living on graph1) and s1 (living on graph2) evolves as
! the two graphs gradually diverge.
!
! Parallelization: outer loop over samples is parallelized with OpenMP.
! Each thread owns its own graph, spins, and RNG (all threadprivate).
!
! Output line format:
!   energia is s_diff=  <is>  <s_diff>  <t>  <E/N> <q^2> <qL^2> <qR^2> ...
!======================================================================
program vb_mc
  use omp_lib
  use mod_params
  use mod_sim,          only: is, delta
  use mod_random,       only: idummy, idummy1, gasdev_iset, xoro_init
  use mod_spins,        only: sigma, sigma_campo, ener, ener1
  use mod_replica,      only: beta, beta_c, alloc_replica
  implicit none

  integer :: r, tid, isamp

  ! Timestamp vars
  character(len=14) :: strong
  integer :: values(8)
  character(len=3), parameter :: months(12) = &
       ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

  ! ---- Simulation parameters ----
  idummy  = idummy_seed
  idummy1 = idummy1_seed
  sigma   = 0.0

  beta_c = log( (1.0 + 1.0/sqrt(real(connett)-1.0)) / &
                (1.0 - 1.0/sqrt(real(connett)-1.0)) ) / 2.0
  print *, 'beta_c = ', beta_c
  print *, 'Running with', omp_get_max_threads(), 'OpenMP threads'

  ! Generate timestamp for output file
  call date_and_time(values=values)
  write(strong, '(i2.2,a,a,i2.2,a,i2.2,a,i2.2)') &
    values(3), '-', months(values(2)), &
    values(5), ':', values(6), ':', values(7)

  open(unit=40, status='unknown', file=trim(strong)//'-VB-dat.dat')

  write(40,*) '# tw=', tw, ' sample=', sample, ' sigma=', sigma
  write(40,*) '# sigma_campo=', sigma_campo, ' l=', l, ' connett=', connett
  write(40,*) '# idummy=', idummy, ' idummy1=', idummy1, ' nt=', nt

  ! Temperature ladder (shared read-only)
  call alloc_replica()
  do r = 1, rep
    beta(r) = beta_c * 2.0 * (1.0 - (real(r)-1.0)/15.0)
  end do
  do r = 1, rep
    write(40,*) '# temp, rep:', 1.0/beta(r), r, 1.0/beta_c
  end do
  call flush(40)

  ! Initialise per-thread RNGs with unique seeds
  !$omp parallel private(tid)
    tid     = omp_get_thread_num()
    idummy  = idummy_seed  - tid * 98765
    idummy1 = idummy1_seed - tid * 54321
    gasdev_iset = 0
    xoro_init = 0
  !$omp end parallel

  !--------------------------------------------------------------------
  ! Outer loop: disorder samples — PARALLELIZED
  !--------------------------------------------------------------------
  !$omp parallel do schedule(dynamic) default(shared)
  do isamp = 1, sample
    is = isamp
    call prog()
  end do
  !$omp end parallel do

  call flush(40)

end program vb_mc


!======================================================================
! prog(): one disorder sample — progressive rewiring.
! 1. Generate base graph; copy to graph2 (identical at start).
! 2. Initialize BOTH spin copies on their respective graphs.
! 3. For s_diff = 0 .. sdiff_max:
!       if s_diff > 0: apply one additional swap to graph2
!       run tw MC steps
!       measure and output overlaps
!======================================================================
subroutine prog()
  use mod_params
  use mod_sim,          only: is
  use mod_random,       only: idummy
  use mod_graph,        only: genera_grafo, add_one_swap, s_diff, reset_graph, alloc_graph
  use mod_spins,        only: s, s1, ener, ener1, aux, init_s, init_overlaps, &
                               oh, oh1, oh2, energia, alloc_spins
  use mod_replica,      only: beta, beta_c, perm, acc, count, &
                               permutazione, scambio, init_exp_table, alloc_replica
  use mod_mc,           only: monte
  implicit none

  integer :: t, r, ja, conta, meas_time, totale, step, ii, n_step, s_diff_old
  real(8) :: enemed(rep), ove(rep), ove1(rep), ove2(rep)
  real(8) :: eee, fff

  call alloc_graph()
  call alloc_spins()
  call alloc_replica()

  ! Temperature ladder reset for this sample
  do r = 1, rep
    beta(r) = beta_c * 2.0 * (1.0 - (real(r)-1.0)/15.0)
  end do
  call init_exp_table()

  call permutazione()
  call genera_grafo()    ! sets nn/coup AND nn1/coup1 = nn/coup, s_diff=0
  call init_s()          ! initializes both s and s1 on the respective graphs
  call init_overlaps()   ! compute initial identity-graph overlaps
  call energia()         ! redundant but safe, sets ener/ener1 base values

  !----------------------------------------------------------------
  ! Progressive rewiring loop
  ! This loop continues until the graph is totally rewired.
  !----------------------------------------------------------------
  do step = 0, 1000
    ! Incremental rewiring: we do NOT reset graph2.
    ! We add a certain number of NEW swaps to the current state of nn1.
    if (step == 0) then
       n_step = 0
    else
       n_step = int(1.14**step)
    end if

    s_diff_old = s_diff

    if (n_step > 0) then
       do ii = 1, n_step
          call add_one_swap()
       end do
       ! Re-sync energies with the new topology nn1/coup1
       call energia()
    end if

    ! Stop if the graph is fully rewired (no new successful swaps could be applied)
    if (n_step > 0 .and. s_diff == s_diff_old) exit

    ! --- Relaxation / equilibration ---
    aux = 0.0;  meas_time = 1;  conta = 0
    enemed = 0.0;  acc = 0.0;  count = 0.0
    ove = 0.0;  ove1 = 0.0;  ove2 = 0.0

    do t = 1, tw
      conta = conta + 1
      call monte()

      eee = 0.0;  fff = 0.0
      do r = 1, rep
        enemed(r) = enemed(r) + ener(perm(r))
        ! oh, oh1, oh2 are maintained incrementally in monte() — just read them
        eee    = eee + ener(r) / l
        fff    = fff + ener(r)**2
      end do
      do r = 1, rep
        ove(r)  = ove(r)  + (real(oh(perm(r)))      /l)**2
        ove1(r) = ove1(r) + (real(oh1(perm(r)))*2.0 /l)**2
        ove2(r) = ove2(r) + (real(oh2(perm(r)))*2.0 /l)**2
      end do
      call scambio()

      ! Log at doubly-spaced measurement times within this relaxation
      if (t == meas_time) then
        meas_time = meas_time * 2
        !$omp critical (write_unit40)
        write(40,'(a,i8,3i8,100f20.6)') &
          'energia is s_diff=', is, step,s_diff, t, &
          (enemed(r)/conta/l, ove(r)/conta, r=1,rep)
!          (enemed(r)/conta/l, ove(r)/conta, ove1(r)/conta, ove2(r)/conta, r=1,rep)        
        !$omp end critical (write_unit40)
        enemed = 0.0;  acc = 0.0;  count = 0.0
        ove = 0.0;  ove1 = 0.0;  ove2 = 0.0;  conta = 0
      end if
    end do

  end do  ! step loop

  !$omp critical (write_unit40)
  write(40,*) ' ';  write(40,*) ' '
  !$omp end critical (write_unit40)

  ! --- Measurement phase (disabled and removed to cleanly exclude modules) ---

end subroutine prog
