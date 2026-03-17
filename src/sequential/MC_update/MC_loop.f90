module mcloop
    ! Use all substeps of the MC loop:
    ! 1. PROPOSAL:
        ! 1.1 propose new dihedral
        ! 1.2 compute new positions (rotation matrix)
    ! 2. ENERGY DIFFERENCE:
        ! 2.1. store old energy
        ! 2.2. compute new energy (for new dihedral/positions)
        ! 2.3. compute energy difference:
            ! Strategy: 
                ! --------------------
                !  BONDED INTERACTION
                ! --------------------
                ! Since only one k-th-dihedral  is changed per iteration, 
                ! only the torsion energy difference of this k-dihedral
                ! needs to be computed.
                ! ---------------------
                ! NONBONDED INTERACTION
                ! ---------------------
                ! The non-bonded interaction between the atoms to the left
                ! and the atoms to the right of the dihedral are the only
                ! ones that will be different after the rotation.
                ! Thus, we group the particles between the "left" [1,k+1]
                ! and "right" [k+2, N] planes. The non-bonded energy change will be
                ! due *only* to the interaction of each particle in the left plane
                ! all particles in the right plane.
                ! So we *do not* need to compute non-bonded interactions between particles
                ! in the same plane.
                ! ---------------------
    ! 3. ACCEPTANCE TEST:
        ! 3.1. exp(-dE/T)...
        ! if accepted: update positions, dihedrals
        ! if rejected: maintain old positions, dihedrals


    use mcproposal
    use energy
    use system
    use io_module
    implicit none

    double precision :: E_new, dE
    integer :: ntry, naccept

    contains

    subroutine runMC(ntry, naccept) 
        ! Runs the basic Monte Carlo loops and stages (equilibration, production and sampling)
        integer :: i, j, k, n_mcs
        integer, intent(inout) :: ntry, naccept
        double precision :: R_old(3, N), phi_old, deltaPhi, dE

        ! total number of MC steps
        n_mcs = N_MCEQUI+N_MCPROD
        
        ntry = 0
        naccept = 0
        do i=1, n_mcs
            if (mod(i,100).eq.0) print*, "MC step:", i
            ! Each MC steps corresponds to (on average) trying to change
            ! each dihedral once
            do j=1, NATTEMPTS
                ! generate new trial configuration
                call MC_trial_step(k, R_old, phi_old, deltaPhi, dE)
                ! accept or reject the proposed configuration
                call accept_reject(k, R_old, deltaPhi, phi_old, dE, ntry, naccept, E_new)
            end do
            if (i.ge.N_MCEQUI) then
                ! Only store data during "production" phase
                if (mod(i, NSAVE) == 0) then
                    call sample(i)
                    call writeXYZ("systemConfig.xyz", i)
                end if
            end if
        end do
    end subroutine

    subroutine MC_trial_step(k, R_old, phi_old, deltaPhi, dE)
        ! Implements a single MC trial step, that is,
        ! it proposes to change a single dihedral and computes the
        ! energy difference between the old and new configurations
        implicit none
        integer, intent(out) :: k
        integer :: I
        double precision, intent(out) :: R_old(3,N), phi_old, deltaPhi, dE
        double precision :: utors_old_k, utors_new_k, enI, enb_old, enb_new


        ! propose a dihedral change deltaPhi
        call proposeDihedral(k, deltaPhi)
        ! Store old positions and dihedral involved in the update
        R_old = R
        
        ! Check we need to recompute the Verlet list before the rotation
        if (isVlist.eq.1) call checkUpdateVlist(k, R_old(:,k+3), nlist, posv, list)

        ! compute *old* torsion and non-bonded energies
        call enerTorsion(k, utors_old_k)
        enb_old = 0.d0
        do I = 1, k+1
            if (isVlist.eq.1) then
                call enerPartVlist(R(1, I), R(2, I), R(3, I), k-1, enI)
            else if (isVlist.eq.0) then
                call enerPart(R(1, I), R(2, I), R(3, I), I, k-1, enI)
            end if
            enb_old = enb_old + enI
        end do

        ! rotate particles:[k+3,N] by deltaPhi
        call rotatePart(k, deltaPhi)
        ! Update dihedral
        phi_old = DANG(k)
        DANG(k) = DANG(k) + deltaPhi

        ! Check we need to recompute the Verlet list after the rotation
        if (isVlist.eq.1) call checkUpdateVlist(k, R(:,k+3), nlist, posv, list)

        ! compute *new* torsion and non-bonded energies
        call enerTorsion(k, utors_new_k)
        enb_new = 0.d0
        do I = 1, k+1
            if (isVlist.eq.1) then
                call enerPartVlist(R(1, I), R(2, I), R(3, I), k-1, enI)
            else if (isVlist.eq.0) then
                call enerPart(R(1, I), R(2, I), R(3, I), I, k-1, enI)
            end if
            enb_new = enb_new + enI
        end do

        ! compute change of energy
        dE = (enb_new - enb_old) + (utors_new_k - utors_old_k)
        ! print*, "dE(torsion)", utors_new_k - utors_old_k
    end subroutine

    subroutine accept_reject(k, R_old, deltaPhi, phi_old, dE, ntry, naccept, E_new)
        ! Decides wether or not to accept the proposed change of the system given
        ! by a change in the k-th dihedral by an amount deltaPhi.
        ! Assumes En in use... (total energy)
        implicit none
        integer, intent(in) :: k
        integer, intent(inout) :: ntry, naccept
        double precision, intent(in) :: R_old(3,N), deltaPhi, phi_old, dE
        double precision, intent(out) :: E_new
        double precision :: rnd
        
        ntry = ntry + 1
        E_new = En + dE
        ! print*, "deltaPhi:", deltaPhi

        if (dE.eq.0.d0) then
            ! Never accept moves that leave the system unchanged
            DANG(k) = phi_old
            R = R_old 
        else if (dE.lt.0.d0) then
            ! Always accept moves that reduce the energy
            En = E_new
            naccept = naccept + 1
        else if (dE.gt.0.d0) then
            ! Accept moves that increase the energy with prob. given by the Metropolis rule
            call random_number(rnd)
            if (rnd .le. exp(-dE/TEMP)) then
                En = E_new
                naccept = naccept + 1
            else
                ! Restore changed positions and dihedral
                DANG(k) = phi_old
                R= R_old
            end if
        end if
    end subroutine accept_reject

    subroutine sample(step)
        ! Extracts data during the "production" phase for statistical analysis
        implicit none
        integer, intent(in) :: step
        integer :: i
        double precision :: Ree, Rg, Rcm(3)
        character(len=512) :: path_ener, path_tors, path_struct

        ! Energy evolution
        path_ener = get_filepath("energy.dat")
        open(40, file = trim(path_ener), position="append", status="unknown")
        ! Write Step and Total Energy (En)
        write(40, '(I8, F18.6)') step, En
        close(40)

        ! Torsion angle
        path_tors = get_filepath("torsions.dat")
        open(41, file = trim(path_tors), position="append", status="unknown")
        ! Write Step and all dihedral angles
        do i = 1, N-3
            write(41, '(I8, F14.6)') step, DANG(i)
        end do
        close(41)

        !! End-to-end Distance and Radius of Gyration
        ! End-to-end Distance (Ree)
        Ree = sqrt(sum((R(:, N) - R(:, 1))**2.d0))

        ! Center of Mass (Rcm) for radius of gyration
        Rcm = 0.d0
        do i = 1, N
            Rcm(:) = Rcm(:) + R(:, i)
        end do
        Rcm = Rcm / dble(N)

        ! Radius of Gyration (Rg)
        Rg = 0.d0
        do i = 1, N
            Rg = Rg + sum((R(:, i) - Rcm)**2.d0)
        end do
        Rg = sqrt(Rg / dble(N))

        path_struct = get_filepath("structure.dat")
        open(42, file=trim(path_struct), position="append", status="unknown")
        ! Write Step, Ree, and Rg
        write(42, '(I8, F14.6, F14.6)') step, Ree, Rg
        close(42)

    end subroutine sample
end module mcloop

! Safer version of MC_step:
! call enerTorsion(k, utors_old_k) ! store old energy
! call enerNonBond(enb_old)
! call proposeDihedral(k, deltaPhi) ! propose a dihedral change deltaPhi

! ! Store old positions and dihedral involved in the update
! R_old = R(k+3, N)
! DANG_old = DANG(k)

! ! rotate particles by deltaPhi
! call rotatePart(k, deltaPhi)

! ! Compute change of energy
! call enerTorsion(k, utors_new_k) ! only need to compute the contribution of k 
! call enerNonBond(enb_new) ! total non-bonded energy 
! dE = (enb_new - enb_old) + (utors_new_k - utors_old_k)