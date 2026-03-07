module system
        implicit none

        integer :: N, SHIFT
        double precision :: TEMP, BOX, HBOX, BLEN, BANG
        double precision :: EPS, SIG, MASS, RC
        double precision, allocatable :: R(:,:)   ! positions (3, N)
        double precision, allocatable :: DANG(:)  ! dihedrals (N-3)
end module system
