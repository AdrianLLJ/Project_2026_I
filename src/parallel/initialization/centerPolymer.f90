module centerPolymer
    ! Module to shift the polymer to the center of the box after initialization
    use system
    implicit none
    double precision :: rcm(3)=0.d0, displ(3)=0.d0
    contains
        subroutine shiftPolymer(displ)
            implicit none
            integer :: i
            double precision, intent(in) :: displ(3)

            do i = 1, N
                R(:,i) = R(:,i) + displ
            end do
        end subroutine shiftPolymer


        subroutine computeCOM(rcm)
            implicit none

            double precision, intent(inout) :: rcm(3)

            rcm = sum(R, dim=2) / dble(N)
        end subroutine computeCOM


        subroutine centerInBox()
            implicit none
            
            call computeCOM(rcm)

            displ = (/BOX/2.d0, BOX/2.d0, BOX/2.d0/) - rcm

            call shiftPolymer(displ)
        end subroutine centerInBox
    end module centerPolymer