module verlet 
    ! ---
    ! This module contains all the subroutines that create/update the Verlet list,
    ! considering a Verlet radius RV for a system with PBC to consider
    ! the minimum image as neigbhours (based on Bekker et al.'s algorithm).
    ! ---
    use constants
    use io_module
    use system
    use nonBonded

    implicit none
    
    ! nlist and list are allocated in new_vlist. Make sure to initialize them at some point.
    integer, allocatable :: nlist(:), list(:,:)
    double precision :: posv(3, N)

    contains

    subroutine checkUpdateVlist(k, poskp3, nlist, posv, list)
        ! Checks if the Verlet lists needs to be recomputed pre/post a dihedral rotation.
        ! The subrotine recives as input which dihedral was changed (k), and checks if 
        ! particle k+3 has moved more that the 'skin depth' given by RC, RV. 
        ! Given how the the proposing of new dihedrals is done, we only need to check if the first
        ! rotating particle has left its RV sphere. See module 'mcloop'.

        implicit none
        integer, intent(in) :: k
        integer, allocatable, intent(inout) :: nlist(:), list(:,:)
        double precision, intent(in) :: poskp3(3)  ! the position (old/new) of the k+3 particle
        double precision, intent(inout) :: posv(3, N)
        double precision :: rposkp3, rvposkp3

        ! Define norm of poskp3 and vposkp3 vectors (the latter is the current position of 
        ! particle k+3 in the Verlet scheme).
        rposkp3 = sqrt(poskp3(1)**2 + poskp3(2)**2 + poskp3(3)**2) 
        rvposkp3 = sqrt(posv(1, k+3,1)**2 + posv(2, k+3)**2 + posv(3, k+3)**2)        

        ! Check if need to update
        ! TODO : Check if (RV-RC) needs to be (RV-RC)/2.d0
        if (abs(rposkp3-rvposkp3).gt.(RV-RC)) call new_vlist(nlist, posv, list)

    end subroutine checkUpdateVlist

    subroutine new_vlist(nlist, posv, list)
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
        integer, allocatable, intent(out) :: nlist(:), list(:,:)
        double precision, intent(out):: posv(3, N)

        ! Blind allocation: we create a matrix list with as many columns as if the Verlet list
        ! was not used, namely, N choose 2. Just so that we have enough components (to avoid
        ! 'smarter' but more prone to error/less readable allocations).
        allocate(nlist(N))
        allocate(list(N, int(N*(N-1)/2.d0)))

        ! Copy the positions of each particle
        do i=1, N
            nlist(i) = 0
            posv(:, i) = R(:, i)
        end do

        do i=1, N-1
            ! Shift of 3 to ensure we skip particles related by the same dihedral
            do j=i+3, N
                dpos(:) = R(:, i) - R(:, j)
                ! Apply PBC (minimum image convention) for the neigbhours
                call minImgConv(dpos(1), dpos(2), dpos(3))
                ! Assign neighbours to Verlet lists
                dposDist = sqrt(dpos(1)**2 + dpos(2)**2 + dpos(3)**2)
                if (dposDist.lt.RV) then
                    nlist(i) = nlist(i) + 1
                    nlist(j) = nlist(j) + 1
                    list(i, nlist(i)) = j
                    list(j, nlist(j)) = i
                end if 
            end do
        end do
    end subroutine new_vlist

    subroutine enerPartVlist(Xi, Yi, Zi, I, enI)
        ! Computes the non-bonded energy of particle I by the interaction with
        ! all its Verlet neigbouring particles. It serves the same purpose as subroutine
        ! 'enerPart' in module nonBonded.
        implicit none

        integer, intent(in) :: I
        integer :: j, neighIj
        double precision, intent(in) :: Xi, Yi, Zi
        double precision, intent(out) :: enI
        double precision :: dx, dy, dz, r2, enIj

        enI = 0.d0 
        do neighIj = 1, nlist(I)
            j = list(I, neighIj)
            ! Relative displacement between particles I, j
            dx = Xi - R(1, j)
            dy = Yi - R(2, j)
            dz = Zi - R(3, j)

            ! Lennard-Jones interaction between particles I, j
            r2 = dx*dx + dy*dy + dz*dz
            call enerLenJon(r2, enIj)
            enI = enI + enIj
        end do
    end subroutine enerPartVlist

   
end module verlet