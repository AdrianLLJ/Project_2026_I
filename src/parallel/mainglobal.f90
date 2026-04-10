program mainglobal

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
    integer :: ierr, rank_world, nproc_world
    integer :: num_replicas, color, rank_rep, REPLICA_COMM
    integer :: nseed, i
    integer, allocatable :: seed(:)
    double precision :: t1, t2
	character(len=512) :: path_log

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_world, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc_world, ierr)

    ! Leer parámetros solo en rank 0
    if (rank_world == 0) then
        call readInput()
    end if

    ! Enviar parámetros a todos los procesos
    call broadcastInput()

    !--Hierarchical MPI setup
    num_replicas = nproc_world / nproc_per_replica
    ! 'color' identifies the replica group (0, 1, 2...)
    color = rank_world / nproc_per_replica
    ! Split the global communicator into replica sub-communicators
    call MPI_Comm_split(MPI_COMM_WORLD, color, rank_world, REPLICA_COMM, ierr)
    ! Get rank within the local replica (0 to nproc_per_replica-1)
    call MPI_Comm_rank(REPLICA_COMM, rank_rep, ierr)
    ! ---Hierachical MPI setup

    ! Start CPU counter (only on replica masters)
    if (rank_rep .eq. 0) call cpu_time(t1)

    ! Different RNG per replica
    call random_seed(size=nseed)
    allocate(seed(nseed))
    do i = 1, nseed
        seed(i) = 12345 + rank_world*1000 + i
    end do
    call random_seed(put=seed)
    ! Different RNG per replica

    ! IO (each replica writes in its own folder)
    call init_io(color)

    ! Verificación
    if (rank_world == 0) then
        print *, "Input broadcasted to", num_replicas, "replicas."
    end if

    ! Distribute temperatures
    if (num_replicas > 1) then
        TEMP = TEMP * (MAX_TEMP / TEMP)**(dble(color) / dble(num_replicas-1))
    end if

    ! System initialization
    call allocateSystem()
    call initDihedrals()
    call initPolymer()
    call centerInBox()
    call writeXYZ("systemConfig.xyz", 0, color)
    ! System initialization

    ! Log file for Master
    if (rank_rep == 0) then
        write(path_log, '(A,I4.4,A)') "simulation_", color, ".log"
        path_log = get_filepath(path_log)
        open(91, file = trim(path_log), position="append", status="unknown")
    end if

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
    write(91, *), "Initial Enb, Eb, En:", Enb, Eb, En
    ! Energies initialization

    ! MC evolution
    call runMC(rank_world, nproc_world, num_replicas, ntry, naccept, REPLICA_COMM)
    ! write(91, *), "Final: ntry, naccept:", ntry, naccept
    ! MC evolution

    ! Test final energies and acceptance ratio
    call totEnergy(En, Eb, Enb)
    ! write(91, *) "Final Enb, Eb, En:", Enb, Eb, En
    ! write(91, *) "Attempts, accepted, ratio(%):", ntry, naccept, 100.d0*naccept/ntry
    ! Test final energies and acceptance ratio

    ! Final CPU counter
    if (rank_rep .eq. 0) call cpu_time(t2)

    if (rank_rep == 0) then
        write(91, *) "Rank:", color
        write(91, *) 'RC, ECUT:', RC, ECUT
        write(91, *) "Attempts, accepted, ratio(%):", ntry, naccept, 100.d0*naccept/ntry
        write(91, *) "Total CPU time:", t2-t1

        close(91)
    end if

    ! Cleanup
    call MPI_Comm_free(REPLICA_COMM, ierr)
    call MPI_Finalize(ierr)

end program mainglobal