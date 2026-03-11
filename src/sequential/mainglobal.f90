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
	character(len=512) :: path_log
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


    ! Open log file to write basic results of the simulation
    path_log = get_filepath("simulation.log")
    open(91, file = trim(path_log), position="append", status="unknown")
    
    ! Energies initialization
    call shiftLenJon()
    write(91, *) 'RC, ECUT:', RC, ECUT
    call totEnergy(En, Eb, Enb)
    write(91, *), "Initial Enb, Eb, En:", Enb, Eb, En
    ! Energies initialization

    ! MC evolution
    call runMC(ntry, naccept)
    write(91, *), "Final: ntry, naccept:", ntry, naccept
    ! MC evolution

    ! Test final energies and acceptance ratio
    call totEnergy(En, Eb, Enb)
    write(91, *) "Final Enb, Eb, En:", Enb, Eb, En
    write(91, *) "Attempts, accepted, ratio(%):", ntry, naccept, 100.d0*naccept/ntry
    ! Test final energies and acceptance ratio

    close(91)

end program mainglobal