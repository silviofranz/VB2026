!======================================================================
! mod_params: compile-time parameters shared across all modules.
! No dependencies. Include with: use mod_params
!======================================================================
module mod_params
  implicit none

  integer, parameter :: expl     = 9
  integer, parameter :: l        = 2**expl       ! # spins
  integer, parameter :: rep      = 11             ! parallel tempering replicas
  integer, parameter :: exig     = 2
  integer, parameter :: invgamma = 2**exig       ! = 4  (hierarchical block size)
  integer, parameter :: connett  = 4             ! graph connectivity (degree)
  integer, parameter :: nclaus   = connett*l/2   ! = 256
  integer, parameter :: lnorm    = l/invgamma    ! = 32 (# of blocks)

  ! ---- Simulation and Output Parameters ----
  integer, parameter :: sample   = 128*10           ! total disorder samples
  integer, parameter :: tw       = 2**17            ! equilibration MC sweeps
  integer, parameter :: nt       = 16            ! measurement phase sweeps (2**nt)
  
  ! ---- Random Number Generator Seeds ----
  integer, parameter :: idummy_seed  = -19827645 ! base seed for MC sweeps
  integer, parameter :: idummy1_seed = -1187163  ! base seed for initialisation/graphs
end module mod_params
