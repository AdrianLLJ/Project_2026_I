module system
        implicit none

        integer :: N
        double precision :: BLEN, BANG
        double precision, allocatable :: R(:,:)   ! positions (3, N)
        double precision, allocatable :: DANG(:)  ! dihedrals (N-3)
end module system
