program mainglobal
    ! cooligula: this is a work in progress of the ultimate "global" main program, 
    ! for now it just allows the sequential execution of the
    ! initialisation and the energy modules

    use mpi
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
    integer :: ierr, rank, nproc
    integer :: nseed, i
    integer, allocatable :: seed(:)
    double precision :: t1, t2
	character(len=512) :: path_log

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    ! Start CPU counter
    call cpu_time(t1)

    ! Different RNG per replica
    call random_seed(size=nseed)
    allocate(seed(nseed))
    do i = 1, nseed
        seed(i) = 12345 + rank*1000 + i
    end do
    call random_seed(put=seed)
    ! Different RNG per replica

    ! IO (each replica writes in its own folder)
    call init_io(rank)

    ! Leer parámetros solo en rank 0
    if (rank == 0) then
        call readInput()
    end if

    ! Enviar parámetros a todos los procesos
    call broadcastInput()

    ! Verificación
    if (rank == 0) then
        print *, "Input broadcasted to", nproc, "replicas."
    end if

    ! Distribute temperatures
    if (nproc > 1) then
        TEMP = TEMP * (MAX_TEMP / TEMP)**(dble(rank) / dble(nproc-1))
    end if

    ! System initialization
    call allocateSystem()
    call initDihedrals()
    call initPolymer()
    call centerInBox()
    call writeXYZ("systemConfig.xyz", 0)
    ! System initialization

    ! Log file for each replica
    write(path_log, '(A,I4.4,A)') "simulation_", rank, ".log"
    path_log = get_filepath(path_log)
    open(91, file = trim(path_log), position="append", status="unknown")

    ! Verlet list initialization
    if (isVlist.eq.1) then 
        call allocVerlet()
        call new_vlist()
    end if
    ! Verlet list initialization
    
    ! Energies initialization
    call shiftLenJon()
    ! write(91, *) 'RC, ECUT:', RC, ECUT
    call totEnergy(En, Eb, Enb)     
    ! write(91, *), "Initial Enb, Eb, En:", Enb, Eb, En
    ! Energies initialization

    ! MC evolution
    call runMC(rank, nproc, ntry, naccept)
    ! write(91, *), "Final: ntry, naccept:", ntry, naccept
    ! MC evolution

    ! Test final energies and acceptance ratio
    call totEnergy(En, Eb, Enb)
    !write(91, *) "Final Enb, Eb, En:", Enb, Eb, En
    !write(91, *) "Attempts, accepted, ratio(%):", ntry, naccept, 100.d0*naccept/ntry
    ! Test final energies and acceptance ratio

    ! Final CPU counter
    call cpu_time(t2)

    write(91, *) "Rank:", rank
    write(91, *) "Attempts, accepted, ratio(%):", ntry, naccept, 100.d0*naccept/ntry
    write(91, *) "Total CPU time:", t2-t1

    close(91)

    call MPI_Finalize(ierr)

end program mainglobal