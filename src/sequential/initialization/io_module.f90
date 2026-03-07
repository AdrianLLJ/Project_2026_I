module io_module
        use constants
        use system
        implicit none

        contains
                subroutine readInput()
                        open(10, file = "input.dat")
                        read(10,*) N
                        read(10,*) TEMP
                        read(10,*) BOX
                        read(10,*) BLEN
                        read(10,*) BANG
                        read(10,*) SIG
                        read(10,*) EPS
                        read(10,*) MASS
                        read(10,*) SHIFT
                        read(10,*) RC
                        close(10)

                        ! Convert degrees to radians
                        BANG = BANG * PI / 180.d0

                        ! Define HBOX as half of the simulation box length
                        HBOX = 0.5d0*BOX

                        ! Define actual cut-off radius
                        RC = min(RC*SIG, HBOX)
                end subroutine readInput

                subroutine writeXYZ(filename)
                        character(len=*), intent(in) :: filename
                        integer :: i

                        open(20, file = filename)
                        write(20, '(I8)') N
                        write(20, *) "Polymer initial configuration"

                        do i = 1, N
                            write(20, '(A, 3F12.6)') "C", R(1, i), R(2,i), R(3,i)
                        end do

                        close(20)
                end subroutine writeXYZ
end module io_module
