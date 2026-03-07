program mainglobal

    ! use hello_module

    ! call sayHello()

    ! cooligula: this is a work in progress of the ultimate "global" main program, 
    ! for now it just allows the sequential execution of the
    ! initialisation and the energy modules

    ! Jonathan's modules
    use io_module
    use init_config
    use centerPolymer
    use system

    ! Adrian's modules
    use nonBonded
    use bonded
    use energy

    implicit none

    ! Initialization
    call readInput()
    call allocateSystem()
    call initDihedrals()
    call initPolymer()
    call centerInBox()
    call writeXYZ("initial.xyz")
    ! Initialization

    ! Energy
    call shiftLenJon()
    print*, 'RC, ECUT:', RC, ECUT
    call enerNonBond(Enb)
    call enerBonded(Eb)
    call totEnergy(En)
    print*, "Enb, Eb, En:", Enb, Eb, En
    ! Energy

end program mainglobal