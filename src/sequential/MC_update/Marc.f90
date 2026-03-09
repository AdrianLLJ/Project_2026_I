program MC_Update
        implicit none
        integer,parameter::kd=8
        integer,parameter::Natoms=6
        real(kd),dimension(Natoms,3)::initialpos
        real(kd),dimension(Natoms,3)::pos
        integer,parameter::Nsteps=100000
        real(kd),parameter::PI=4*atan(1.0_8)
        real(kd),parameter::d_to_r=2*PI/360
        real(kd),parameter::Temp=0.1
        real(kd),parameter::maxDih=5
        real(kd)::dihedral,NewD,E,NewE,rand,dE,b2norm,x,y,norma
        real(kd),dimension(3)::b1,b2,b3,posa1,posa2,posa3,posa4,n1,n2,ax,normax,posorig,posrot
        integer::i,j,currentIt,r,k,step

        real(8) :: Rot(3,3)
        real(8) :: c, s, t

        !Estructura inicial preestablerta de 6 atoms, s'ha de canviar per la funció de inicialització
        ! initialpos(1,:)=(/-4.1441938029, 1.1522778995, 0.6791332440/)
        ! initialpos(2,:)=(/-3.0600337337, 0.3910845926, 1.4121238072/)
        ! initialpos(3,:)=(/-1.7044639241, 0.9346274965, 1.0133719403/)
        ! initialpos(4,:)=(/-0.6202171934, 0.1735068110, 1.7463127311/)
        ! initialpos(5,:)=(/ 0.7353523678, 0.7170534180, 1.3475650660/)
        ! initialpos(6,:)=(/ 1.8195118466,-0.0441309577, 2.0805657877/)
        ! pos=initialpos

        do step=1,Nsteps
        call random_number(rand)
        r = 2 + FLOOR((Natoms-3)*rand) !random dihedral angle to change
        !Mesurement of dihedral angle between atoms 1 and 4 from the set of 4 atoms
        !consecutive in this case
        posa1=pos(r-1,:)
        posa2=pos(r,:)
        posa3=pos(r+1,:)
        posa4=pos(r+2,:)

        b1=posa2-posa1
        b2=posa3-posa2
        b3=posa4-posa3

        n1=cross(b1,b2)
        n2=cross(b2,b3)

        b2norm = sqrt(dot_product(b2, b2))

        x = dot_product(n1, n2)
        y = b2norm * dot_product(b1, n2)

        dihedral = atan2(y, x)/d_to_r
        !Fins aquí només és una funció per mesurar una angle dihedre
        
        !A partir d'aquí és el procés per canviar el dihedre seleccionat a un agle +NewD
        call random_number(rand)
        rand=rand-0.5
        NewD=2*rand*maxDih*d_to_r
        posorig=pos(r+1,:)
        do i=1,Natoms !Canvia el centre de coordenades al atom r+1
                pos(i,1)=pos(i,1)-posorig(1)
                pos(i,2)=pos(i,2)-posorig(2)
                pos(i,3)=pos(i,3)-posorig(3)
        end do

        ax=-pos(r,:)
        normax=sqrt(dot_product(ax,ax))
        ax=ax/normax

        c = cos(NewD)
        s = sin(NewD)
        t = 1.0d0 - c

        ! Build rotation matrix
        Rot(1,1) = c + ax(1)*ax(1)*t
        Rot(1,2) = ax(1)*ax(2)*t - ax(3)*s
        Rot(1,3) = ax(1)*ax(3)*t + ax(2)*s

        Rot(2,1) = ax(2)*ax(1)*t + ax(3)*s
        Rot(2,2) = c + ax(2)*ax(2)*t
        Rot(2,3) = ax(2)*ax(3)*t - ax(1)*s

        Rot(3,1) = ax(3)*ax(1)*t - ax(2)*s
        Rot(3,2) = ax(3)*ax(2)*t + ax(1)*s
        Rot(3,3) = c + ax(3)*ax(3)*t


        do j = r+1, Natoms
             posrot=pos(j,:)
             pos(j,1) = dot_product(Rot(1,:), posrot)
             pos(j,2) = dot_product(Rot(2,:), posrot)
             pos(j,3) = dot_product(Rot(3,:), posrot)
        end do
        enddo
        do k=1,Natoms !print les coordenades transformades per comprovar que funciona
                print *, 'C',pos(k,1),pos(k,2),pos(k,3)
        enddo
        !dihedral=dihedral_0
        !E=energy(dihedral)
        !currentIt=0
        !open(unit=10,file='E1.txt',action='write')
        !write(10,'(3A12)') 'Iteration','Dihedral','Energy'
        !do i=1,N_MCS
        !        currentIt=currentIt+1
        !        call random_number(rand)
        !        rand=rand-0.5
        !        NewD=dihedral+2*rand*maxDih
        !        NewE=energy(NewD)
        !        if (NewE<=E) then
        !                E=NewE
        !                dihedral=NewD
        !        else
        !                call random_number (rand)
        !                dE=NewE-E
        !                if (rand<=exp(-dE/Temp))then
        !                        E=NewE
        !                        dihedral=NewD
        !                endif
        !        endif
        !        write (10,'(I12,2F12.8)') currentIt,dihedral,E
        !enddo

        contains

        FUNCTION cross(a, b) !Function to make the cross product
                integer,parameter::kd=8
                real(kd), dimension(3) :: cross
                real(kd), dimension(3), intent(in) :: a, b

                cross(1) = a(2) * b(3) - a(3) * b(2)
                cross(2) = a(3) * b(1) - a(1) * b(3)
                cross(3) = a(1) * b(2) - a(2) * b(1)
        END FUNCTION cross
end program MC_Update