module init_config
        use constants
        use system
        implicit none

        contains
                subroutine allocateSystem()
                        allocate(R(3,N))

                        if (N > 3) then
                                allocate(DANG(N-3))
                        end if
                end subroutine allocateSystem

                subroutine initDihedrals()
                        double precision :: r
                        integer :: i

                        if (N > 3) then
                                if (initRandom.eq.1) then
                                        ! Random dihedral initialisation
                                        do i = 1, N-3
                                                call random_number(r)
                                                DANG(i) = PI * (2.d0 * r - 1)
                                        end do
                                else
                                        !Initialises all dihedrals to pi
                                        do i = 1, N-3
                                                DANG(i) = PI
                                        end do
                                end if
                        end if
                end subroutine initDihedrals

                subroutine initPolymer()
                        integer :: i
                        double precision :: b1(3), b2(3)
                        double precision :: e1(3), e2(3), e3(3)
                        double precision :: phi

                        ! first atom
                        R(:,1) = (/0.d0, 0.d0, 0.d0/)
                        
                        ! second atom
                        R(:,2) = (/BLEN, 0.d0, 0.d0/)

                        ! third atom
                        R(1,3) = BLEN - BLEN * cos(BANG)
                        R(2,3) = BLEN * sin(BANG)
                        R(3,3) = 0.d0

                        ! remaining atoms
                        do i = 4, N
                            b1 = R(:,i-2) - R(:,i-3)
                            b2 = R(:,i-1) - R(:,i-2)

                            e1 = b2 / sqrt(sum(b2**2))
                            call cross(b1, b2, e2)
                            e2 = e2 / sqrt(sum(e2**2))
                            call cross(e2, e1, e3)

                            phi = DANG(i-3)

                            R(:,i) = R(:,i-1) + BLEN * ( &
                                     - cos(BANG) * e1 &
                                     + sin(BANG) * sin(phi) * e2 &
                                     + sin(BANG) * cos(phi) * e3 )
                        end do
                end subroutine initPolymer

                subroutine cross(a,b,c)
                        ! Computes the cross product c = a x b
                        double precision, intent(in) :: a(3), b(3)
                        double precision, intent(out) :: c(3)

                        c(1) = a(2) * b(3) - a(3) * b(2)
                        c(2) = a(3) * b(1) - a(1) * b(3)
                        c(3) = a(1) * b(2) - a(2) * b(1)
                end subroutine cross
end module init_config
