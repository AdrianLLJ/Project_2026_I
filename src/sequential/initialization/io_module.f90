module io_module
        use constants
        use system
        implicit none

        character(len=256), save :: out_dir

        contains
                subroutine init_io()
                        ! Read the command line argument
                        call get_command_argument(1, out_dir)
                        
                        ! Fallback if no argument is provided
                        if (trim(out_dir) == "") out_dir = "."
                end subroutine init_io

                function get_filepath(filename) result(filepath)
                        ! Helper function: Takes "name.dat" and returns "build/results/name.dat"
                        character(len=*), intent(in) :: filename
                        character(len=512) :: filepath
                        
                        filepath = trim(out_dir) // "/" // trim(filename)
                end function get_filepath

                subroutine readInput()
                        ! open(10, file = "input.dat")
                        ! read(10,*) N
                        ! read(10,*) TEMP
                        ! read(10,*) BOX
                        ! read(10,*) BLEN
                        ! read(10,*) BANG
                        ! read(10,*) SIG
                        ! read(10,*) EPS
                        ! read(10,*) MASS
                        ! read(10,*) SHIFT
                        ! read(10,*) RC
                        ! read(10,*) maxDih
                        ! read(10,*) N_MCEQUI
                        ! read(10,*) N_MCPROD
                        ! read(10,*) NATTEMPTS
                        ! close(10)

                        ! Read directly from standard input instead of a hardcoded file
                        ! This way, we can run ./exec < input.dat
                        read(*,*) N
                        read(*,*) TEMP
                        read(*,*) BOX
                        read(*,*) BLEN
                        read(*,*) BANG
                        read(*,*) SIG
                        read(*,*) EPS
                        read(*,*) MASS
                        read(*,*) SHIFT
                        read(*,*) RC
                        read(*,*) RV
                        read(*,*) maxDih
                        read(*,*) N_MCEQUI
                        read(*,*) N_MCPROD
                        read(*,*) NATTEMPTS
                        read(*,*) NSAVE
                        read(*,*) initRandom
                        read(*,*) isVlist

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
                        character(len=512) :: fullpath

                        ! Get the full path: results folder + filename
                        fullpath = get_filepath(filename)
                        ! append: to concatenate new data; unknown: if it does not exist, it creates the file 
                        open(20, file = trim(fullpath), status="unknown", position="append", action="write")
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
