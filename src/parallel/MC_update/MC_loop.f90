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
    use mpi
    implicit none

    double precision :: E_new, dE
    integer :: ntry, naccept

    contains

    subroutine runMC(rank_world, nproc_world, ntry, naccept, REPLICA_COMM) 
        ! Runs the basic Monte Carlo loops and stages (equilibration, production and sampling)
        integer, intent(in) :: rank_world, nproc_world, REPLICA_COMM
        integer :: ierr, rank_rep
        integer :: i, j, k, n_mcs
        integer, intent(inout) :: ntry, naccept
        double precision :: R_old(3, N), phi_old, deltaPhi, dE

        ! Get local rank within the replica
        call MPI_Comm_rank(REPLICA_COMM, rank_rep, ierr)

        ! total number of MC steps
        n_mcs = N_MCEQUI+N_MCPROD
        ntry = 0
        naccept = 0

        do i=1, n_mcs
            ! Only the master of each replica prints to avoid flooding
            if (rank_rep == 0 .and. mod(i,100).eq.0) print*, "MC step:", i
            ! Each MC steps corresponds to (on average) trying to change each dihedral once
            do j=1, NATTEMPTS
                ! generate new trial configuration
                call MC_trial_step(k, R_old, phi_old, deltaPhi, dE, REPLICA_COMM)
                ! accept or reject the proposed configuration
                call accept_reject(k, R_old, deltaPhi, phi_old, dE, ntry, naccept, E_new, REPLICA_COMM)
            end do
            ! Sample observables during equilibration & production to track changes
            if (mod(i, NSAVE) == 0) then
                call sample(i, rank_world, nproc_world)
            end if
            if ((i.ge.N_MCEQUI).and.(mod(i, NSAVE) == 0)) then
                ! Only store config. for the master of the replica
                if (rank_rep == 0) call writeXYZ("systemConfig.xyz", i, rank_world)
            end if

            ! Try replica exchange
            if (mod(i, N_SWAP) == 0) then
                call replica_exchange(rank_world, nproc_world, i)
            end if
        end do
    end subroutine

    subroutine MC_trial_step(k, R_old, phi_old, deltaPhi, dE, REPLICA_COMM)
        ! Implements a single MC trial step, that is,
        ! it proposes to change a single dihedral and computes the
        ! energy difference between the old and new configurations
        implicit none
        integer, intent(out) :: k
        integer :: I
        double precision, intent(out) :: R_old(3,N), phi_old, deltaPhi, dE
        double precision :: utors_old_k, utors_new_k, enI, enb_old, enb_new

        ! MPI Variables for Intra-replica parallelism
        integer, intent(in) :: REPLICA_COMM
        integer :: rank_rep, nprocs_rep, ierr
        integer :: n_total, chunk_size, remainder
        integer :: i_start, i_end
        double precision :: local_enb_old, local_enb_new

        ! Get local processor rank and size inside this specific replica
        call MPI_Comm_rank(REPLICA_COMM, rank_rep, ierr)
        call MPI_Comm_size(REPLICA_COMM, nprocs_rep, ierr)

        ! Master process of the replica : propose a dihedral change deltaPhi 
        if (rank_rep == 0) then 
            call proposeDihedral(k, deltaPhi)
        end if
        ! Broadcast the chosen dihedral index and rotation angle to all workers in the replica
        call MPI_Bcast(k, 1, MPI_INTEGER, 0, REPLICA_COMM, ierr)
        call MPI_Bcast(deltaPhi, 1, MPI_DOUBLE_PRECISION, 0, REPLICA_COMM, ierr)

        ! Store old positions and dihedral involved in the update
        R_old = R
        ! Check we need to recompute the Verlet list before the rotation
        if (isVlist.eq.1) then 
            call checkUpdateVlist()
        end if

        ! Calculate local loop bounds for I = 1 to k+1
        n_total = k + 1
        chunk_size = n_total / nprocs_rep
        remainder = mod(n_total, nprocs_rep)
        if (rank_rep < remainder) then
            i_start = 1 + rank_rep * (chunk_size + 1)
            i_end = i_start + chunk_size
        else
            i_start = 1 + remainder * (chunk_size + 1) + (rank_rep - remainder) * chunk_size
            i_end = i_start + chunk_size - 1
        end if

        ! compute *old* torsion and non-bonded energies (in parallel)
        call enerTorsion(k, utors_old_k)
        local_enb_old = 0.d0
        do I = i_start, i_end
            if (isVlist.eq.1) then
                ! before was k-1 as input instead of I
                call enerPartVlist(R(1, I), R(2, I), R(3, I), I, k, enI)
            else if (isVlist.eq.0) then
                call enerPart(R(1, I), R(2, I), R(3, I), I, max(k, I+1), enI)
            end if
            local_enb_old = local_enb_old + enI
        end do
        ! Gather total old energy across all ranks in the replica
        call MPI_Allreduce(local_enb_old, enb_old, 1, MPI_DOUBLE_PRECISION, MPI_SUM, REPLICA_COMM, ierr)

        ! rotate particles:[k+3,N] by deltaPhi
        phi_old = DANG(k)
        call rotatePart(k, deltaPhi)
        ! Update dihedral
        DANG(k) = DANG(k) + deltaPhi

        ! Check we need to recompute the Verlet list after the rotation
        if (isVlist.eq.1) then 
            call checkUpdateVlist()
        end if

        ! compute *new* torsion and non-bonded energies
        call enerTorsion(k, utors_new_k)
        local_enb_new = 0.d0
        do I = i_start, i_end
            if (isVlist.eq.1) then
                call enerPartVlist(R(1, I), R(2, I), R(3, I), I, k, enI)
            else if (isVlist.eq.0) then
                call enerPart(R(1, I), R(2, I), R(3, I), I, max(k, I+1), enI)
            end if
            local_enb_new = local_enb_new + enI
        end do
        ! Gather total new energy across all ranks in the replica
        call MPI_Allreduce(local_enb_new, enb_new, 1, MPI_DOUBLE_PRECISION, MPI_SUM, REPLICA_COMM, ierr)

        ! compute change of energy
        dE = (enb_new - enb_old) + (utors_new_k - utors_old_k)
        ! print*, "dE(torsion)", utors_new_k - utors_old_k
    end subroutine

    subroutine accept_reject(k, R_old, deltaPhi, phi_old, dE, ntry, naccept, E_new, REPLICA_COMM)
        ! Decides wether or not to accept the proposed change of the system given
        ! by a change in the k-th dihedral by an amount deltaPhi.
        ! Assumes En in use... (total energy)
        implicit none
        integer, intent(in) :: k, REPLICA_COMM
        integer, intent(inout) :: ntry, naccept
        double precision, intent(in) :: R_old(3,N), deltaPhi, phi_old, dE
        double precision, intent(out) :: E_new
        double precision :: rnd

        integer :: rank_rep, ierr, accepted_local

        call MPI_Comm_rank(REPLICA_COMM, rank_rep, ierr)
        accepted_local = 0
        E_new = En + dE

        ! Only the Master decides to avoid divergence due to different RNG streams
        if (rank_rep .eq. 0) then
            ntry = ntry + 1
            if (dE.eq.0.d0) then
                ! Never accept moves that leave the system unchanged
                accepted_local = 0
            else if (dE.lt.0.d0) then
                ! Always accept moves that reduce the energy
                accepted_local = 1
            else if (dE.gt.0.d0) then
                ! Accept moves that increase the energy with prob. given by the Metropolis rule
                call random_number(rnd)
                if (rnd .le. exp(-dE/TEMP)) then
                    accepted_local = 1
                end if
            end if
        end if
        ! Broadcast the master's decision to all workers in the replica
        call MPI_Bcast(accepted_local, 1, MPI_INTEGER, 0, REPLICA_COMM, ierr)

        ! Apply the decision synchronously across all processors in the replica
        if (accepted_local .eq. 1) then
            En = E_new
            if (rank_rep .eq. 0) naccept = naccept + 1
        else
            ! Restore changed positions and dihedral
            DANG(k) = phi_old
            R = R_old
        end if
    end subroutine accept_reject

    subroutine sample(step, rank, nproc)
        ! Extracts data during the "production" phase for statistical analysis
        implicit none
        integer, intent(in) :: step, rank, nproc
        integer :: i, ierr
        double precision :: Ree, Rg, Rcm(3)
        character(len=512) :: path_ener, path_tors, path_struct
        
        ! Dynamic buffers for gathered values
        double precision :: all_En(nproc), all_Ree(nproc), all_Rg(nproc)

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
        
        ! Gather physical properties from all replicas into Rank 0
        call MPI_Gather(En, 1, MPI_DOUBLE_PRECISION, all_En, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(Ree, 1, MPI_DOUBLE_PRECISION, all_Ree, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(Rg, 1, MPI_DOUBLE_PRECISION, all_Rg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (rank == 0) then
            ! Write step and ALL energies to single file
            path_ener = get_filepath("energy.dat")
            open(40, file = trim(path_ener), position="append", status="unknown")
            write(40, *) step, all_En
            close(40)

            ! Write step and ALL structural properties to single file
            path_struct = get_filepath("structure.dat")
            open(42, file=trim(path_struct), position="append", status="unknown")
            write(42, *) step, all_Ree, all_Rg
            close(42)

            ! Torsion angle (Rank 0 only, no MPI Gathering necessary)
            path_tors = get_filepath("torsions.dat")
            open(41, file = trim(path_tors), position="append", status="unknown")
            do i = 1, N-3
                write(41, '(I8, F14.6)') step, DANG(i)
            end do
            close(41)
        end if
    end subroutine sample

    subroutine replica_exchange(rank, nproc, step)
        use nonBonded, only: new_vlist
        implicit none
        integer, intent(in) :: rank, nproc, step
        integer :: partner, ierr
        integer :: status(MPI_STATUS_SIZE)
        double precision :: E_me, E_partner
        double precision :: T_me, T_partner
        double precision :: swap_prob, rnd
        integer :: accept_swap, swap_decision
        
        partner = -1
        
        ! Determine partner safely based on step counter
        if (mod(step/N_SWAP, 2) == 0) then
            ! Even pairs: 0-1, 2-3, 4-5
            if (mod(rank, 2) == 0) then
                partner = rank + 1
            else
                partner = rank - 1
            end if
        else
            ! Odd pairs: 1-2, 3-4, 5-6
            if (mod(rank, 2) == 1) then
                partner = rank + 1
            else
                partner = rank - 1
            end if
        end if
        
        ! Ensure boundary partners don't go out of bounds
        if (partner >= 0 .and. partner < nproc) then
            ! We have a valid swap pair
            E_me = En
            T_me = TEMP
            
            ! Exchange energy and temperature
            call MPI_Sendrecv(E_me, 1, MPI_DOUBLE_PRECISION, partner, 0, &
                              E_partner, 1, MPI_DOUBLE_PRECISION, partner, 0, &
                              MPI_COMM_WORLD, status, ierr)
                              
            call MPI_Sendrecv(T_me, 1, MPI_DOUBLE_PRECISION, partner, 1, &
                              T_partner, 1, MPI_DOUBLE_PRECISION, partner, 1, &
                              MPI_COMM_WORLD, status, ierr)
                              
            ! Lower rank evaluates the Metropolis criterion
            if (rank == min(rank, partner)) then
                swap_prob = exp( (1.d0/T_me - 1.d0/T_partner) * (E_me - E_partner) )
                
                if (swap_prob >= 1.d0) then
                    accept_swap = 1
                else
                    call random_number(rnd)
                    if (rnd <= swap_prob) then
                        accept_swap = 1
                    else
                        accept_swap = 0
                    end if
                end if
            end if
            
            ! Broadcast decision to the higher rank
            if (rank == min(rank, partner)) then
                call MPI_Send(accept_swap, 1, MPI_INTEGER, partner, 2, MPI_COMM_WORLD, ierr)
                swap_decision = accept_swap
            else
                call MPI_Recv(swap_decision, 1, MPI_INTEGER, partner, 2, MPI_COMM_WORLD, status, ierr)
            end if
            
            ! If accepted, swap the state
            if (swap_decision == 1) then
                call MPI_Sendrecv_replace(R, 3*N, MPI_DOUBLE_PRECISION, partner, 3, partner, 3, MPI_COMM_WORLD, status, ierr)
                if (N > 3) call MPI_Sendrecv_replace(DANG, N-3, MPI_DOUBLE_PRECISION, partner, &
                                 4, partner, 4, MPI_COMM_WORLD, status, ierr)
                call MPI_Sendrecv_replace(En,  1, MPI_DOUBLE_PRECISION, partner, 5, partner, 5, MPI_COMM_WORLD, status, ierr)
                call MPI_Sendrecv_replace(Eb,  1, MPI_DOUBLE_PRECISION, partner, 6, partner, 6, MPI_COMM_WORLD, status, ierr)
                call MPI_Sendrecv_replace(Enb, 1, MPI_DOUBLE_PRECISION, partner, 7, partner, 7, MPI_COMM_WORLD, status, ierr)
                
                ! Need to reconstruct the Verlet list because positions suddenly changed!
                if (isVlist == 1) then
                    call new_vlist()
                end if
            end if
        end if
    end subroutine replica_exchange
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