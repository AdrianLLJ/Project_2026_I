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
    integer :: nseed
    integer, allocatable :: seed(:)
    double precision :: t1, t2
	character(len=512) :: path_log

    ! Start CPU counter
    call cpu_time(t1)

    ! Initialize RNG
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

    ! Verlet list initialization
    if (isVlist.eq.1) then 
        call allocVerlet()
        call new_vlist()
    end if
    ! Verlet list initialization
    
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
    write(91, *) "Final running energies:", Enb, Eb, En
    call totEnergy(En, Eb, Enb)
    write(91, *) "Final configuration energies Enb, Eb, En:", Enb, Eb, En
    write(91, *) "Attempts, accepted, ratio(%):", ntry, naccept, 100.d0*naccept/ntry
    ! Test final energies and acceptance ratio

    ! Final CPU counter
    call cpu_time(t2)
    write(91, *) "Total CPU time:", t2-t1

    close(91)

end program mainglobal