module io_module
use constants
use system
implicit none

! Ouput directory per process (replica)
character(len=256), save :: out_dir

contains

subroutine init_io(rank)
        ! Initializes the output directory for each process
        ! (replica) based on the command line argument and
        ! the rank of the process.
        integer, intent(in) :: rank
        character(len=512) :: base_dir
        
        ! Read the command line argument
        call get_command_argument(1, base_dir)
        if (trim(base_dir) == "") base_dir = "results"
        
        ! Each process writes in its own folder
        write(out_dir, '(A,"/replica_",I4.4)') trim(base_dir), rank

        ! Create the directory if it doesn't exist
        call system("mkdir -p " // trim(out_dir))
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
        read(*,*) MAX_NEIGH
        read(*,*) maxDih
        read(*,*) N_MCEQUI
        read(*,*) N_MCPROD
        read(*,*) NATTEMPTS
        read(*,*) NSAVE
        read(*,*) initRandom
        read(*,*) isVlist
        read(*,*) MAX_TEMP
        read(*,*) N_SWAP

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
        integer, intent(in) :: iframe
        integer :: i
        character(len=512) :: fullpath

        ! Get the full path: results folder + filename
        !fullpath = trim(out_dir) // "/" // trim(filename)
        fullpath = get_filepath(filename)
        ! append: to concatenate new data; unknown: if it does not exist, it creates the file 
        open(20, file = trim(fullpath), status="unknown", position="append", action="write")
        ! Write the number of atoms and frame number as a comment line for VMD
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

subroutine broadcastInput()
        use mpi
        use system
        implicit none
        integer :: ierr

        call MPI_Bcast(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(SHIFT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(initRandom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(isVlist, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(MAX_NEIGH, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(N_MCEQUI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(N_MCPROD, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NATTEMPTS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(NSAVE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call MPI_Bcast(TEMP, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(BOX, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(HBOX, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(BLEN, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(BANG, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(EPS, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(SIG, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(MASS, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(RC, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(RV, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(maxDih, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(MAX_TEMP, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(N_SWAP, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine broadcastInput

end module io_module
