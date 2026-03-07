module bonded 
    ! ---
    ! This module defines and initialises all the variables related to the bonded interactions
    ! (i.e. torsion potential) as well as all the functions/subroutines for computing the 
    ! bonded energy contribution. 
    ! The current version only considers a C skeleton.
    ! ---
    use constants
    use system
    implicit none
    double precision :: Eb = 0.d0

    contains

    subroutine enerTorsion(dihedral, utors)
        ! Computes the torsion energy contribution of one dihedral, from:
        ! William L. Jorgensen, David S. Maxwell, and Julian Tirado-Rives
        ! Journal of the American Chemical Society 1996 118 (45), 11225-11236
        
        implicit none
        double precision, intent(in) :: dihedral
        double precision, intent(out) :: utors
        double precision :: V1, V2, V3

        V1 = 1.411d0
        V2 = -0.271d0
        V3 = 3.145d0

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

        eb = 0.d0
        do i = 1, N - 3
            call enerTorsion(DANG(i), utors)
            eb = eb + utors
        end do
    end subroutine enerBonded

end module bonded