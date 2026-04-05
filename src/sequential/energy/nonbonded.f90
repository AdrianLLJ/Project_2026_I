module nonBonded
    ! ---
    ! This module defines and initialises all the variables related to the non-bonded interactions
    ! (i.e. Lennard-Jones potential) as well as all the functions/subroutines for computing the 
    ! non-bonded energy contribution. 
    ! The current version only considers a C skeleton.
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
    use system
    implicit none

    ! ---FIRST SECTION VARIABLES
    double precision :: RC2, EPS4, SIG2
    double precision :: ECUT, Enb = 0.d0

    ! ---SECOND SECTION VARIABLES
    ! nlist and list are allocated in new_vlist. Make sure to initialize them at some point.
    integer, allocatable :: nlist(:), list(:,:)
    double precision, allocatable :: posv(:,:)

    contains

    ! ---FIRST SECTION

    subroutine shiftLenJon()
        ! Sets the energy shift for the Lennard-Jones potential.
        ! This subroutine should be called only once, prior to the MC simulation.
        implicit none

        RC2 = RC*RC
        EPS4 = 4.d0*EPS
        SIG2 = SIG*SIG

        if (SHIFT.eq.1) then
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
            if (SHIFT.eq.1) then
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
        ! do j = Jb, N
        do j = Jb + 3, N
            ! if (abs(j-I).ge.3) then
            if (j.ne.I) then
                ! Relative displacement between particles I, j
                dx = Xi - R(1, j)
                dy = Yi - R(2, j)
                dz = Zi - R(3, j)

                !PBC (minimum image convention)
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
        ! It allows to make the computation with or without Verlet lists (via the global variable isVlist).
        implicit none

        integer :: i, jb
        double precision, intent(out) :: enb
        double precision :: xi, yi, zi, eni

        enb = 0.d0
        do i = 1, N - 4
            xi = R(1, i)
            yi = R(2, i)
            zi = R(3, i)

            if (isVlist.eq.1) then
                call enerPartVlist(Xi, Yi, Zi, i, 0, eni)
            else if (isVlist.eq.0) then
                jb = i + 1 
                call enerPart(xi, yi, zi, i, jb, eni)
            end if

            enb = enb + eni
        end do
    end subroutine enerNonBond

    ! ! ---SECOND SECTION --------------------------------------------------------------------------------
    ! This section contains all the subroutines that create/update the Verlet list,
    ! considering a Verlet radius RV for a system with PBC to consider
    ! the minimum image as neigbhours (based on Bekker et al.'s algorithm).
    ! The author (@AdrianLLJ) chose to add it to the same module since there is a codependence between
    ! subroutines from the different sections (since they are all related to the non-bonded interactions)
    ! ----------------------------------------------------------------------------------------------------

    subroutine checkUpdateVlist()
        ! Checks if the Verlet lists needs to be recomputed pre/post a dihedral rotation.
        !  The subroutine checks which particle has moved the most after a rotation of 1 diherdral.
        ! Then for that particle checks if this distance is greater than the 'skin depth'. See module 'mcloop'.

        implicit none
        integer :: i
        double precision :: dx, dy, dz, dr2, max_dr2

        max_dr2 = 0.d0
        do i = 1, N
            dx = R(1, i) - posv(1, i)
            dy = R(2, i) - posv(2, i)
            dz = R(3, i) - posv(3, i)
            
            call minImgConv(dx, dy, dz)
            dr2 = dx*dx + dy*dy + dz*dz
            if (dr2 .gt. max_dr2) max_dr2 = dr2
        end do

        if (sqrt(max_dr2) .gt. (RV-RC)/2.d0) then 
            call new_vlist()
        end if

    end subroutine checkUpdateVlist

    !  subroutine checkUpdateVlist(k, poskp3)
    !     ! Checks if the Verlet lists needs to be recomputed pre/post a dihedral rotation.
    !     ! The subrotine recives as input which dihedral was changed (k), and checks if 
    !     ! particle k+3 has moved more that the 'skin depth' given by RC, RV. 
    !     ! Given how the the proposing of new dihedrals is done, we only need to check if the first
    !     ! rotating particle has left its RV sphere. See module 'mcloop'.

    !     implicit none
    !     integer, intent(in) :: k
    !     double precision, intent(in) :: poskp3(3)  ! the position (old/new) of the k+3 particle
    !     double precision :: dxkp3, dykp3, dzkp3, drkp3

    !     ! compute displacement between new position and Verlet reference
    !     dxkp3 = poskp3(1) - posv(1, k+3)
    !     dykp3 = poskp3(2) - posv(2, k+3)
    !     dzkp3 = poskp3(3) - posv(3, k+3)

    !     ! apply minimum image convention for PBC
    !     call minImgConv(dxkp3, dykp3, dzkp3)

    !     drkp3 = sqrt(dxkp3*dxkp3 + dykp3*dykp3 + dzkp3*dzkp3)        

    !     ! Check if need to update
    !     ! TODO : check if RV-RC needs to be (RV-RC)/2
    !     if (drkp3.gt.(RV-RC)/2.d0) then 
    !         call new_vlist()
    !     end if

    ! end subroutine checkUpdateVlist

    subroutine allocVerlet()
        ! Blind allocation: trial and error allocation of the number of neigbours for the
        ! Verlet list; depends on the density of particles given an initial configuration
        ! In general, it is safer to add enough, since the num. of neighbours changes during the simulation
        implicit none

        allocate(posv(3, N))
        allocate(nlist(N))
        allocate(list(N, MAX_NEIGH))
    end subroutine allocVerlet

    subroutine new_vlist()
        !---
        ! Creates/reconstructs the Verlet list of each non-bonding interacting particle with PBC.
        ! The subroutine computes the number of neighbours of each interacting particle, nlist, if
        ! the distance between particle i and neigbhour j is less than the Verlet radius.
        ! The full list is given in matrix format "list(particle i,neighbour number)=particle j"
        ! Example:
        ! if particle i=23 has neigbhour j=15 within the Verlet radius and it is the 4-th neighbour
        ! found for particle i=23, then in the Verlet list of particle i=23, particle j=15
        ! will be assigned the neigbhour number 4. And so "j=15=list(i=23, nlist(i=23)=4)".
        !---

        implicit none
        integer :: i, j
        double precision :: dposDist, dpos(3)

        ! Blind allocation: we create a matrix list with as many columns as if the Verlet list
        ! was not used, namely, N choose 2. Just so that we have enough components (to avoid
        ! 'smarter' but more prone to error/less readable allocations).

        ! Copy the positions of each particle
        do i=1, N
            nlist(i) = 0
            posv(:, i) = R(:, i)
        end do
        ! before i=1, N-1
        do i=1, N-4
            ! Shift of 4 to ensure we skip particles related by the same dihedral
            do j=i+4, N
                dpos(:) = R(:, i) - R(:, j)
                ! Apply PBC (minimum image convention) for the neigbhours
                call minImgConv(dpos(1), dpos(2), dpos(3))
                ! Assign neighbours to Verlet lists
                dposDist = sqrt(dpos(1)**2 + dpos(2)**2 + dpos(3)**2)
                if (dposDist.lt.RV) then
                    nlist(i) = nlist(i) + 1
                    nlist(j) = nlist(j) + 1

                    if ((nlist(i).gt.MAX_NEIGH).or.(nlist(j).gt.MAX_NEIGH)) then 
                        write(91, *) "Verlet list overflow"
                    end if

                    list(i, nlist(i)) = j
                    list(j, nlist(j)) = i
                end if 
            end do
        end do
    end subroutine new_vlist

    subroutine enerPartVlist(Xi, Yi, Zi, I, k, enI)
        ! Computes the non-bonded energy of particle I by the interaction with
        ! all its Verlet neigbouring particles.  It does the equivalent of subroutine 'enerPart'.
        ! For total energy computations pass k=0.
        ! For the MC loop, pass the actual value of the dihedral k.
        implicit none

        integer, intent(in) :: I, k
        integer :: j, neighIj
        double precision, intent(in) :: Xi, Yi, Zi
        double precision, intent(out) :: enI
        double precision :: dx, dy, dz, r2, enIj

        enI = 0.d0 
        do neighIj = 1, nlist(I)
            j = list(I, neighIj)
            ! Avoid double counting and particles related by the same dihedral
            if (j.ge.max(k+2, I+4)) then
                ! Relative displacement between particles I, j
                dx = Xi - R(1, j)
                dy = Yi - R(2, j)
                dz = Zi - R(3, j)

                !PBC (minimum image convention)
                call minImgConv(dx, dy, dz)

                ! Lennard-Jones interaction between particles I, j
                r2 = dx*dx + dy*dy + dz*dz
                call enerLenJon(r2, enIj)
                enI = enI + enIj
            end if
        end do
    end subroutine enerPartVlist
end module nonBonded