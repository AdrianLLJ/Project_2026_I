program init_replicas
    use mpi
    use io_module
    use init_config
    use centerPolymer
    use system

    implicit none

    integer :: ierr, rank, nprocs
    integer :: nseed, i
    integer, allocatable :: seed(:)
    character(len=512) :: path_log

    ! Inicializar MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    ! Inicializar IO
    call init_io()

    ! Cada réplica escribe en su carpeta
    write(out_dir, '(A,"/replica_",I4.4)') trim(out_dir), rank
    call execute_command_line("mkdir -p " // trim(out_dir))

    ! RNG distinto por réplica
    call random_seed(size=nseed)
    allocate(seed(nseed))
    do i = 1, nseed
        seed(i) = 12345 + rank*1000 + i
    end do
    call random_seed(put=seed)

    ! Leer parámetros solo en rank 0
    if (rank == 0) then
        call readInput()
    end if

    ! Enviar parámetros a todos los procesos
    call broadcastInput()

    ! Verificación
    if (rank == 0) then
        print *, "Input broadcasted to", nprocs, "replicas."
    end if

    ! Inicialización del sistema
    call allocateSystem()
    call initDihedrals()
    call initPolymer()
    call centerInBox()

    ! Guardar configuración
    call writeXYZ("init.xyz", 0)

    ! Calcular centro de masa
    call computeCOM(rcm)

    ! Log de verificación
    path_log = get_filepath("init_check.log")
    open(91, file=trim(path_log), status="unknown")

    write(91,*) "Replica:", rank
    write(91,*) "N:", N
    write(91,*) "Seed example:", seed(1)
    write(91,*) "Center of mass:", rcm
    write(91,*) "First atom:", R(:,1)
    write(91,*) "Last atom:", R(:,N)

    close(91)

    call MPI_Finalize(ierr)

end program init_replicas
