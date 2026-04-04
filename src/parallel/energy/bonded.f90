module bonded 
    ! ---
    ! This module defines and initialises all the variables related to the bonded interactions
    ! (i.e. torsion potential) as well as all the functions/subroutines for computing the 
    ! bonded energy contribution. 
    ! The current version only considers a C skeleton.
    ! ---
    use constants
    use io_module
    use system
    implicit none
    double precision :: Eb = 0.d0

    contains

    subroutine enerTorsion(dihedral_index, utors)
        ! Computes the torsion energy contribution of one dihedral, from:
        ! William L. Jorgensen, David S. Maxwell, and Julian Tirado-Rives
        ! Journal of the American Chemical Society 1996 118 (45), 11225-11236
        
        implicit none
        integer, intent(in) :: dihedral_index
        double precision, intent(out) :: utors
        double precision :: V1, V2, V3, dihedral

        dihedral = DANG(dihedral_index)

        ! In kcal/mol for the polyethilene would be (values from paper):
        ! V1 = 1.411d0
        ! V2 = -0.271d0
        ! V3 = 3.145d0

        ! In reduced units (e.g. dividing by epsilon 0.0913 kcal/mol *for our particular system*)
        V1 = 15.46d0
        V2 = -2.97d0
        V3 = 34.45d0

        utors = 0.5d0*V1*(1 + cos(dihedral)) + &
                0.5d0*V2*(1 - cos(2.d0*dihedral)) + &
                0.5d0*V3*(1 + cos(3.d0*dihedral))
    end subroutine enerTorsion

    subroutine enerBonded(eb)
        ! Computes the total bonded energy (sum of the torsion energy of each dihedral).
        
        implicit none
        integer :: i
        double precision, intent(out) :: eb
        double precision :: utors
        character(len=256) :: out_dir
        character(len=512) :: filepath

        ! Instead of hardcoding the output file directory, it is obtained from the first CLI argument
        ! First we get the argument
        ! filepath = get_filepath("dihedrals.dat")

        eb = 0.d0
        ! open(30, file = trim(filepath), status="unknown", position="append", action="write")
        ! write(30, '(A)') "i, DANG(i), utors(i)"
            do i = 1, N - 3
                call enerTorsion(i, utors)
                !call writeDihedrals(i, utors)
                eb = eb + utors
            end do
        close(30)
    end subroutine enerBonded

    subroutine writeDihedrals(i, utors)
        ! Writes i, DANG(i) and utors to a file.
        implicit none
        integer, intent(in) :: i
        double precision, intent(in) :: utors
        write(30, '(I8, F12.6, F12.6)') i, DANG(i), utors
    end subroutine writeDihedrals
end module bonded