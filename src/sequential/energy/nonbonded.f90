module nonBonded
    ! ---
    ! This module defines and initialises all the variables related to the non-bonded interactions
    ! (i.e. Lennard-Jones potential) as well as all the functions/subroutines for computing the 
    ! non-bonded energy contribution. 
    ! The current version only considers a C skeleton in bulk (implemented using PBC's).
    ! ---
    ! Variables/parameters:
    ! EPS   : the energy of the potential well. 
    ! SIG   : distance at which the particle-particle potential energy is 0.
    ! MASS  : mass of the element.
    ! RC    : cutoff distance at which interactions are truncated.
    ! SHIFT : if true shifts the potential such that at the cutoff the energy is 0.
    ! ECUT  : energy shift value.
    ! ---
    ! Reduced units:
    ! if EPS = SIG = MASS = 1, seting the units of energy, distance and mass, respectively.
    ! Other units can be derived as combinations of this triad (see link).
    ! https://en.wikipedia.org/wiki/Lennard-Jones_potential#Dimensionless_(reduced_units)
    ! ---
    use constants
    implicit none

    logical, parameter :: SHIFT=.TRUE.
    double precision, parameter :: EPS=1.d0, SIG=1.d0, MASS=1.d0, RC=min(2.5d0*SIG, HBOX)
    double precision, parameter :: RC2=RC*RC, EPS4=4.d0*EPS, SIG2=SIG*SIG
    double precision :: ECUT, Enb = 0.d0

    contains

    subroutine shiftLenJon()
        ! Sets the energy shift for the Lennard-Jones potential.
        ! This subroutine should be called only once, prior to the MC simulation.
        implicit none

        if (SHIFT) then
            ECUT = 0.d0 
            call enerLenJon(RC2, ECUT)
        end if
    end subroutine shiftLenJon

    subroutine enerLenJon(r2, enij)
        ! Computes the Lennard-Jones interaction between two non-bonded atoms
        ! separated by a squared distance r2, with a cutoff RC and shift.

        implicit none
        double precision, intent(in) :: r2
        double precision, intent(inout) :: enij
        double precision :: r2i, r6i

        if (r2.le.RC2) then
            r2i = SIG2/r2
            r6i = r2i*r2i*r2i
            if (SHIFT) then
                enij = EPS4*(r6i*r6i-r6i) - ECUT
            else
                enij = EPS4*(r6i*r6i-r6i)
            end if
        else
            enij = 0.d0
        end if
    end subroutine enerLenJon

    subroutine minImgConv(dx, dy, dz)
        ! Applies the minimum image convention to assign Lennard-Jones
        ! interactions between pairs of particles in a bulk system.
        implicit none

        double precision, intent(inout) :: dx, dy, dz

        ! E.g. if dx < 0, then sign(BOX, dx) = -BOX, thus we add +BOX
        if (abs(dx).gt.HBOX) dx = dx - sign(BOX, dx)
        if (abs(dy).gt.HBOX) dy = dy - sign(BOX, dy)
        if (abs(dz).gt.HBOX) dz = dz - sign(BOX, dz)
    end subroutine minImgConv

    subroutine enerPart(Xi, Yi, Zi, I, Jb, enI)
        ! Computes the non-bonded energy contribution of one particle due to all other
        ! interacting particles (i.e. within the cutoff region).
        ! Since we are considering a rigid molecule (bond lengths/angles fixed) with
        ! dihedral angles defined by every 4 consecutive particles, we only consider 
        ! non-bonded interactions for particles separated by 4+ consecutive particles.
        implicit none

        integer, intent(in) :: I, Jb
        integer :: j
        double precision, intent(in) :: Xi, Yi, Zi
        double precision, intent(out) :: enI
        double precision :: dx, dy, dz, r2, enIj

        enI = 0.d0
        ! Shift of 3 to ensure we skip particles related by the same dihedral  
        do j = Jb + 3, NPART
            if (j.ne.I) then
                ! Relative displacement between particles I, j
                dx = Xi - X(j)
                dy = Yi - Y(j)
                dz = Zi - Z(j)

                ! PBC (minimum image convention)
                call minImgConv(dx, dy, dz)

                ! Lennard-Jones interaction between particles I, j
                r2 = dx*dx + dy*dy + dz*dz
                call enerLenJon(r2, enIj)
                enI = enI + enIj
            end if
        end do
    end subroutine enerPart

    subroutine enerNonBond(enb)
        ! Computes the total non-bonded energy of the system due to each interacting particle.
        implicit none

        integer :: i, jb
        double precision, intent(out) :: enb
        double precision :: xi, yi, zi, eni

        enb = 0.d0
        do i = 1, NPART - 4
            xi = X(i)
            yi = Y(i)
            zi = Z(i)
            jb = i + 1
            call enerPart(xi, yi, zi, i, jb, eni)
            enb = enb + eni
        end do
    end subroutine enerNonBond
end module nonBonded