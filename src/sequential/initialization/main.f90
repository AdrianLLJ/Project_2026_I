program create_polymer
        use io_module
        use init_config
        use system
        implicit none

        call readInput()

        call allocateSystem()

        call initDihedrals()

        call initPolymer()
        
        call writeXYZ("initial.xyz")
end program create_polymer
