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
                        read(10,*) maxDih
                        read(10,*) N_MCEQUI
                        read(10,*) N_MCPROD
                        read(10,*) NATTEMPTS
                        close(10)

                        ! Convert degrees to radians
                        BANG = BANG * PI / 180.d0
                        maxDih = maxDih * PI / 180.d0

                        ! Define HBOX as half of the simulation box length
                        HBOX = 0.5d0*BOX

                        ! Define actual cut-off radius
                        RC = min(RC*SIG, HBOX)

                end subroutine readInput

                subroutine writeXYZ(filename, iframe)
                        ! Writes snapshot of the system in the current state (iframe) to visualize in VMD.
                        ! We work in unwrapped coordinates R, but apply PBC to R for visualization only
                        character(len=*), intent(in) :: filename
                        integer :: i, iframe
                        ! append: to concatenate new data; unknown: if it does not exist, it creates the file 
                        open(20, file = filename, status="unknown", position="append", action="write")
                        write(20, '(I8)') N
                        write(20, *) iframe
                        do i = 1, N
                            ! Write coords with PBC: 0<= R(:,i)<= BOX
                            write(20, '(A, 3F12.6)') "C", mod(R(1, i)+BOX,BOX), &
                                                          mod(R(2, i)+BOX,BOX), &
                                                          mod(R(3, i)+BOX,BOX)
                        end do

                        close(20)
                end subroutine writeXYZ
end module io_module
