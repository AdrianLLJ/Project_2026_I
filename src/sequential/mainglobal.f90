program mainglobal
    ! cooligula: this is a work in progress of the ultimate "global" main program, 
    ! for now it just allows the sequential execution of the
    ! initialisation and the energy modules

    ! Initialization modules (@J-dot-Barrientos)
    use io_module
    use init_config
    use centerPolymer
    use system
    ! Energy modules (@AdrianLLJ)
    use nonBonded
    use bonded
    use energy
    ! MC modules (@J-dot-Barrientos, @AdrianLLJ, @cooligula)
    use mcproposal
    use mcloop

    implicit none
    ! Initialize RNG
    integer :: nseed
    integer, allocatable :: seed(:)
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 83226
    call random_seed(put=seed)
    ! Initialize RNG

    ! System initialization
    call init_io()
    call readInput()
    call allocateSystem()
    call initDihedrals()
    call initPolymer()
    call centerInBox()
    call writeXYZ("systemConfig.xyz", 0)
    ! System initialization

    ! Energies initialization
    call shiftLenJon()
    print*, 'RC, ECUT:', RC, ECUT
    call totEnergy(En, Eb, Enb)
    print*, "Initial Enb, Eb, En:", Enb, Eb, En
    ! Energies initialization

    ! MC evolution
    call runMC(ntry, naccept)
    print*, "Final: ntry, naccept:", ntry, naccept
    ! MC evolution

    ! Test
    call totEnergy(En, Eb, Enb)
    print*, "Final Enb, Eb, En:", Enb, Eb, En

end program mainglobal