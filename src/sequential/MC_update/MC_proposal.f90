module mcproposal
    ! This module proposes a new dihedral angle and applies the
    ! a uniform rotation to the molecule to change the positions
    ! of the particles involved.

    use system
    use init_config
    implicit none

    contains

    subroutine proposeDihedral(k, deltaPhi)
        ! Selects which dihedral to change, proposes a change
        ! of the angle and identifies which 4 particles define the dihedral
        implicit none
        double precision :: r
        double precision, intent(out) :: deltaPhi
        integer, intent(out) :: k
        ! Selects a random number between 1 and N-3
        call random_number(r)
        k = 1 + FLOOR((N-3)*r)
        
        ! Angular displacement in (-maxDih, +maxDih)
        call random_number(r)
        deltaPhi = (2.d0*r - 1.d0) * maxDih
    end subroutine proposeDihedral

    subroutine rotatePart(k, deltaPhi)
        ! Rotates particles k+3 to N by deltaPhi on a unit sphere,
        ! where k corresponds to the k-th dihedral.
        implicit none
        integer, intent(in) :: k
        integer :: i
        double precision, intent(in) :: deltaPhi
        double precision :: axis(3), dist2axis(3), r_rot(3), ax_x_dist(3)

        ! Create axis of rotation given k
        axis = R(:, k+1) - R(:, k+2)
        axis = axis / sqrt(sum(axis**2.d0))

        do i = k+3, N
            ! Distance to axis of rotation
            dist2axis = R(:, i) - R(:, k+2)

            ! Compute rotated vector (Rodrigues' formula)
            call cross(axis, dist2axis, ax_x_dist)
            r_rot = dist2axis*cos(deltaPhi) &
                  + ax_x_dist*sin(deltaPhi) &
                  + axis*dot_product(axis, dist2axis)*(1-cos(deltaPhi))

            ! Update position vector given the position of the axis
            R(:, i) = R(:, k+2) + r_rot
        end do
    end subroutine rotatePart
end module mcproposal