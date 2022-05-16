! Fortran code to calculate Order Parameter Q
! Zhipan Liu Group  2020
! Use YLMB FORTRAN code from internet

      subroutine Get_Order_parameter(na,iza,cell,xa,atom,Qglobal, ene, sym)

      implicit none
      integer:: namx !,lmax
      parameter(namx=100)
      integer::na,ntol,i,j,i1,i2,i3,k,bond,l,m,super(3),km,atom,totalbond,control,Iiza
      integer :: iza(na),maxiza,miniza, sym
      double precision,allocatable , save ::xac(:,:),xa0(:,:),r0(:,:)
      double precision xa(3,na),r1,cell0(3,3)
      double precision Q1,co1,co2
      double precision alpha,Pi,Qlocal(4),Qglobal(4)
      double precision dist_all,dist,dtmp,r(3),cell(3,3),celli(3,3)
      double precision rij(namx*namx),theta(namx*namx),phi(namx*namx),radi(namx*namx),rmax
      double precision :: ene
      logical, save :: first=.true.
      COMPLEX*16  Y(49),Q(49),Qall(49)
      logical CONT
      integer:: imode=1

      save maxiza,miniza

      Pi=3.14159265359d0
      alpha=1.0d0
      rmax=10.0d0

      if(first) then
      allocate(xac(3,na),xa0(3,na),r0(na,na))
      first = .false.
      maxiza=0
      miniza=100
      do i=1,na
        do j=1,na
              call fastbond(iza(i),iza(j),r0(i,j))
              if(r0(i,j) < 0.05d0) then
                  call species_radius(iza(i),co1)
                  call species_radius(iza(j),co2)
                  r0(i,j) = co1 + co2
                  if(r0(i,j) < 0.5d0) then
                    r0(i,j)=2d0
                  endif
              end if
        enddo
      enddo
      endif
      if(atom==0) then
       Iiza=0 !maxiza
      else
       Iiza=atom
      endif

      xa0=xa
      cell0=cell
      if(Imode>=0) then
      call cellconstraint(na,cell0,xa0,1)
      endif


      call reclat_loc(cell0, celli, 0)
      do i = 1,na
         xac(:,i) = matmul(transpose(celli),xa0(:,i))
         xac(:,i)= modulo(xac(:,i)+100.0d0,1.0d0)
      enddo

      do k = 1,na
        do j = 1,3
          xa0(j,k) = cell0(j,1) * xac(1,k) + cell0(j,2) * xac(2,k) + cell0(j,3) * xac(3,k)
        enddo
      enddo

      CONT=.false.

      Qall=0.0d0
      Qlocal=0d0
      totalbond=0
      super=0
      if(Imode>=0) super(:)= 6

      do km=1 , na
        bond=1
         if(iza(km)==Iiza .or. Iiza==0) then !goto 90
         j=km
!        write(6,'(A,3f12.4)') 'Si ', xa0(1,j),xa0(2,j),xa0(3,j)
         Q=0d0
         do i=1,na
!             if((iza(i)/=iza(j) .or. CONT) ) then !.and. i/=j) then
              do i1=-super(1),super(1)
              do i2=-super(2),super(2)
              do i3=-super(3),super(3)
                    r(:)=(xac(1,i)+dble(i1))*cell0(:,1)+(xac(2,i)+ &
                          dble(i2))*cell0(:,2)+(xac(3,i)+dble(i3))*cell0(:,3) - xa0(:,j)
                    r1=dsqrt(sum(r*r))
                    if(r1 > 20d0) cycle    ! 20 cutoff radius
                    if(r1>r0(i,j)*0.6d0 .and. r1<r0(i,j)*1.3d0) bond=bond+1
                    if(r1>r0(i,j)*0.6d0 .and. r1<r0(i,j)*1.3d0 ) then
!
                    call YLMB(r,6,Y) ! sphe_harm(l,m,theta(i),phi(i))
                    Q=Q+Y
                    Qall=Qall+Y 
                     endif
              enddo
              enddo
              enddo
!        endif
       enddo

       if(bond>1) bond=bond-1

       do j=0,6,2
       Q1=0d0
       do i=j*j+1,(j+1)*(j+1)
       Q1=Q1+Q(i)*DCONJG(Q(i))
       enddo
       Q1=sqrt(4.0d0*Pi*Q1/(2.0d0*dble(j)+1.0d0))
       Q1=Q1/dble(na)/dble(bond)
       Qlocal(j/2+1)=Q1+Qlocal(j/2+1)
       enddo

      endif ! Iiza
      enddo  ! km

      Qglobal=Qlocal

      End subroutine

      subroutine fastbond(zi,zj,bondlen)

      double precision,dimension(53,53),save :: bondmap
      integer ::  zi,zj,i,j
      double precision :: bondlen
      logical ,save :: init=.true.

      if(init) then
        bondmap=0d0
        bondmap(1,1) = 0.80
        bondmap(1,5) = 1.084
        bondmap(1,6) = 1.084
        bondmap(1,7) = 1.001
        bondmap(1,8) = 0.947
        bondmap(1,9) = 0.92
        bondmap(1,14) = 1.48
        bondmap(1,15) = 1.415
        bondmap(1,16) = 1.326
        bondmap(1,17) = 1.28
        bondmap(1,35) = 1.41
        bondmap(1,53) = 1.60
        bondmap(5,8) = 1.512
        bondmap(6,6) = 1.512
        bondmap(6,7) = 1.439
        bondmap(6,8) = 1.393
        bondmap(6,9) = 1.353
        bondmap(6,14) = 1.86
        bondmap(6,15) = 1.84
        bondmap(6,16) = 1.812
        bondmap(6,17) = 1.781
        bondmap(6,35) = 1.94
        bondmap(6,53) = 2.16
        bondmap(7,7) = 1.283
        bondmap(7,8) = 1.333
        bondmap(7,9) = 1.36
        bondmap(7,14) = 1.74
        bondmap(7,15) = 1.65
        bondmap(7,16) = 1.674
        bondmap(7,17) = 1.75
        bondmap(7,35) = 1.90
        bondmap(7,53) = 2.10
        bondmap(8,8) = 1.45
        bondmap(8,9) = 1.42
        bondmap(8,14) = 1.63
        bondmap(8,15) = 1.66
        bondmap(8,16) = 1.470
        bondmap(8,17) = 1.70
        bondmap(8,22) = 2.1   ! Ti-O
        bondmap(8,35) = 1.85
        bondmap(8,53) = 2.05
        bondmap(9,14) = 1.57
        bondmap(9,15) = 1.54
        bondmap(9,16) = 1.55
        bondmap(14,8) = 2.0
        bondmap(14,14) = 2.32
        bondmap(14,15) = 2.25
        bondmap(14,16) = 2.15
        bondmap(14,17) = 2.02
        bondmap(14,35) = 2.19
        bondmap(14,53) = 2.44
        bondmap(15,15) = 2.21
        bondmap(15,16) = 2.10
        bondmap(15,17) = 2.03
        bondmap(15,35) = 2.21
        bondmap(15,53) = 2.47
        bondmap(16,16) = 2.052
        bondmap(16,17) = 2.04
        bondmap(16,35) = 2.24
        bondmap(16,53) = 2.40
        bondmap(17,17) = 1.99
        bondmap(35,35) = 2.28
        bondmap(53,53) = 2.67
        init = .false.
      end if

      if(zi > 53 .or. zj > 53) then
        bondlen = 0.0d0
      else
        if(zi<zj) then
            bondlen = bondmap(zi,zj)
        else
            bondlen = bondmap(zj,zi)
        end if
      end if
      end subroutine

      subroutine species_radius(atom,radius)

      integer::atom
      double precision::radius

        radius=0.0d0
        if(atom==1)   radius=0.400d0
        if(atom==2)   radius=0.490d0
        if(atom==3)   radius=1.400d0
        if(atom==4)   radius=1.400d0
        if(atom==5)   radius=1.170d0
        if(atom==6)   radius=0.650d0
        if(atom==7)   radius=0.620d0
        if(atom==8)   radius=0.600d0
        if(atom==9)   radius=0.570d0
        if(atom==10)  radius=0.510d0
        if(atom==11)  radius=2.230d0
        if(atom==12)  radius=1.720d0
        if(atom==13)  radius=1.820d0
        if(atom==14)  radius=1.222d0
        if(atom==15)  radius=1.230d0
        if(atom==16)  radius=1.120d0
        if(atom==17)  radius=1.120d0
        if(atom==18)  radius=0.880d0
        if(atom==19)  radius=2.770d0
        if(atom==20)  radius=2.230d0
        if(atom==21)  radius=2.090d0
        if(atom==22)  radius=2.000d0
        if(atom==23)  radius=1.920d0
        if(atom==24)  radius=1.850d0
        if(atom==25)  radius=1.790d0
        if(atom==26)  radius=1.720d0
        if(atom==27)  radius=1.670d0
        if(atom==28)  radius=1.350d0
        if(atom==29)  radius=1.570d0
        if(atom==30)  radius=1.530d0
        if(atom==31)  radius=1.810d0
        if(atom==32)  radius=1.520d0
        if(atom==33)  radius=1.330d0
        if(atom==34)  radius=1.170d0
        if(atom==35)  radius=1.120d0
        if(atom==36)  radius=1.030d0
        if(atom==37)  radius=2.210d0
        if(atom==38)  radius=1.860d0
        if(atom==39)  radius=2.270d0
        if(atom==40)  radius=2.160d0
        if(atom==41)  radius=2.080d0
        if(atom==42)  radius=2.010d0
        if(atom==43)  radius=1.950d0
        if(atom==44)  radius=1.340d0
        if(atom==45)  radius=1.340d0
        if(atom==46)  radius=1.340d0
        if(atom==47)  radius=1.340d0
        if(atom==48)  radius=1.340d0
        if(atom==49)  radius=2.000d0
        if(atom==50)  radius=1.720d0
        if(atom==51)  radius=1.530d0
        if(atom==52)  radius=1.420d0
        if(atom==53)  radius=1.320d0
        if(atom==54)  radius=1.240d0
        if(atom==55)  radius=3.340d0
        if(atom==56)  radius=2.780d0
        if(atom==57)  radius=2.740d0
        if(atom==58)  radius=2.700d0
        if(atom==59)  radius=2.670d0
        if(atom==60)  radius=2.640d0
        if(atom==61)  radius=2.620d0
        if(atom==62)  radius=2.590d0
        if(atom==63)  radius=2.560d0
        if(atom==64)  radius=2.540d0
        if(atom==65)  radius=2.510d0
        if(atom==66)  radius=2.490d0
        if(atom==67)  radius=2.470d0
        if(atom==68)  radius=2.450d0
        if(atom==69)  radius=2.420d0
        if(atom==70)  radius=2.400d0
        if(atom==71)  radius=2.250d0
        if(atom==72)  radius=2.160d0
        if(atom==73)  radius=2.090d0
        if(atom==74)  radius=2.020d0
        if(atom==75)  radius=1.370d0
        if(atom==76)  radius=1.320d0
        if(atom==77)  radius=1.370d0
        if(atom==78)  radius=1.350d0
        if(atom==79)  radius=1.350d0
        if(atom==80)  radius=1.760d0
        if(atom==81)  radius=2.080d0
        if(atom==82)  radius=1.810d0
        if(atom==83)  radius=1.630d0
        if(atom==84)  radius=1.530d0
        if(atom==85)  radius=1.430d0
        if(atom==86)  radius=1.340d0

      end subroutine

      subroutine cellconstraint(na,cell,xa,control)

      implicit none

      double precision  :: cell(3,3),xa(3,na),aplusb,zz


      integer    :: iv, i,ix,ia,j,na,k,i0,j0,k0,m,control
      double precision :: cellm(3), ang(3),celang(3),cell2(3,3),celli(3,3),Pi,cell0(3,3)
      double precision :: latt_abc(3),latt_ang0(3),latt_ang(3),alplp,betlp,gamlp,alp,blp,latt_abc2(3)
      double precision :: clp,xac(3),xxx,ab(3),abl
      logical CONT,SSWoutput
      logical , save :: first = .true.
      integer, save :: change(50,4) ,Lastchange(50,4) !,last_m
!     data last_m /0/
      SSWoutput=.false.
      if(first) then
        Lastchange=0
        first=.false.
      endif

      Pi=3.1415926535897932384626d0

      cell0=cell

      call outcell_output_loc(cell,latt_abc,latt_ang)
!      write(6,'(A8,6F15.10)') 'cell0',latt_abc,latt_ang
     ! change=0
     ! endif


      if(.true.) then

!     last_m=m
        change=0
        m=0
291     continue

        m=m+1

        do iv = 1,3
          cellm(iv) = dot_product(cell(:,iv),cell(:,iv))
          cellm(iv) = dabs(cellm(iv))
        enddo
        i=minloc(cellm,1)
        j=maxloc(cellm,1)
        if(i==j) then
          i=1
          j=2
          k=3
        else
          k=6-i-j
        endif
        zz = dot_product(cell(:,i),cell(:,j))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)

        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=1
            change(m,2)=i
            change(m,3)=j
            goto 291
        endif

        ab(:)=cell(:,j)-floor(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)

        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=3
            change(m,2)=i
            change(m,3)=j
            goto 291
        endif

        celang(:) = cell(:,i)+cell(:,k)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,j),celang(:))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=k
            change(m,4)=j
            goto 291
        endif
        ab(:)=cell(:,j)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=k
            change(m,4)=j
            goto 291
        endif


        celang(:) = cell(:,i)-cell(:,k)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,j),celang(:))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=-k
            change(m,4)=j
            goto 291
        endif
        ab(:)=cell(:,j)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=-k
            change(m,4)=j
            goto 291
        endif

        celang(:) = cell(:,i)+cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,k),celang(:))
        ab(:)=cell(:,k)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14)  then
            cell(:,k)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=j
            change(m,4)=k
            goto 291
        endif
        ab(:)=cell(:,k)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14)  then
            cell(:,k)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=j
            change(m,4)=k
            goto 291
        endif


        celang(:) = cell(:,i)-cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,k),celang(:))
        ab(:)=cell(:,k)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14)  then
            cell(:,k)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=-j
            change(m,4)=k
            goto 291
        endif
        ab(:)=cell(:,k)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14)  then
            cell(:,k)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=-j
            change(m,4)=k
            goto 291
        endif


        celang(:) = cell(:,k)+cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,i),celang(:))
        ab(:)=cell(:,i)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14)  then
            cell(:,i)=ab
            change(m,1)=2
            change(m,2)=k
            change(m,3)=j
            change(m,4)=i
            goto 291
        endif
        ab(:)=cell(:,i)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14)  then
            cell(:,i)=ab
            change(m,1)=4
            change(m,2)=k
            change(m,3)=j
            change(m,4)=i
            goto 291
        endif


        celang(:) = cell(:,k)-cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,i),celang(:))
        ab(:)=cell(:,i)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14)  then
            cell(:,i)=ab
            change(m,1)=2
            change(m,2)=k
            change(m,3)=-j
            change(m,4)=i
            goto 291
        endif
        ab(:)=cell(:,i)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14)  then
            cell(:,i)=ab
            change(m,1)=4
            change(m,2)=k
            change(m,3)=-j
            change(m,4)=i
            goto 291
        endif



        zz = dot_product(cell(:,i),cell(:,k))
        ab(:)=cell(:,k)-ceiling(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14)  then
            cell(:,k)=ab
            change(m,1)=1
            change(m,2)=i
            change(m,3)=k
            goto 291
        endif
        ab(:)=cell(:,k)-floor(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14)  then
            cell(:,k)=ab
            change(m,1)=3
            change(m,2)=i
            change(m,3)=k
            goto 291
        endif



        zz = dot_product(cell(:,k),cell(:,j))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/cellm(k))*sign(1.0d0,zz)*cell(:,k)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=1
            change(m,2)=k
            change(m,3)=j
            goto 291
        endif
        ab(:)=cell(:,j)-floor(dabs(zz)/cellm(k))*sign(1.0d0,zz)*cell(:,k)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14)  then
            cell(:,j)=ab
            change(m,1)=3
            change(m,2)=k
            change(m,3)=j
            goto 291
        endif


         if(SSWoutput) then
           do i=1,m-1
             write(6,*) change(i,1:4)
           enddo
       !   write(6,*) 'lattice change in cellconstraint',m-1
         endif

        if(m>1.or.Lastchange(1,1)/=0) then
          do i=1,50
            do j=1,4
            if(change(i,j)/=Lastchange(i,j)) then
             goto 45
             endif
           enddo
          enddo
        endif
45      continue
        Lastchange=change
      else !if(.not.LallowRot) then

        m=1
        do while (change(m,1)/=0)
         do iv = 1,3
           cellm(iv) = dot_product(cell(:,iv),cell(:,iv))
           cellm(iv) = dabs(cellm(iv))
         enddo

         if(change(m,1)==1) then
            j=change(m,3)
            i=change(m,2)
            zz = dot_product(cell(:,i),cell(:,j))
            cell(:,j)=cell(:,j)-ceiling(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
         elseif(change(m,1)==2) then
           i=change(m,2)
           k=change(m,3)
           j=change(m,4)
           if(k>0) then
             celang(:) = cell(:,i)+cell(:,k)
           else
             celang(:) = cell(:,i)-cell(:,k)
           endif
           aplusb=sum(celang*celang)
           aplusb= dabs(aplusb)
           zz = dot_product(cell(:,j),celang(:))
           cell(:,j)=cell(:,j)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
         elseif(change(m,1)==3) then
           j=change(m,3)
           i=change(m,2)
           zz = dot_product(cell(:,i),cell(:,j))
           cell(:,j)=cell(:,j)-floor(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
         elseif(change(m,1)==4) then
           i=change(m,2)
           k=change(m,3)
           j=change(m,4)
           if(k>0) then
             celang(:) = cell(:,i)+cell(:,k)
           elseif(k<0) then
             celang(:) = cell(:,i)-cell(:,k)
           endif
           aplusb=sum(celang*celang)
           aplusb= dabs(aplusb)
           zz = dot_product(cell(:,j),celang(:))
           cell(:,j)=cell(:,j)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
         endif
         m=m+1
        enddo

      endif   ! LallowRot end


292   continue


       cell2=cell
       call outcell_output_loc(cell,latt_abc,latt_ang)
!       write(6,'(A8,6F15.10)') 'cell1',latt_abc,latt_ang

       alplp=latt_ang(3)* pi/180.d0
       betlp=latt_ang(2)* pi/180.d0
       gamlp=latt_ang(1)* pi/180.d0


       alp=latt_abc(1)
       blp=latt_abc(2)
       clp=latt_abc(3)

        cell(3,3) = clp
        cell(2,3) = 0.d0
        cell(1,3) = 0.d0
        cell(3,2) = blp * dcos(alplp)
        cell(2,2) = blp * dsin(alplp)
        cell(1,2) = 0.d0
        cell(3,1) = alp * dcos(betlp)
        xxx = (dcos(gamlp) - dcos(betlp)*dcos(alplp))/dsin(alplp)
        cell(2,1) = alp * xxx
        cell(1,1) = alp * dsqrt(dsin(betlp)*dsin(betlp) - xxx*xxx)

!     now find the rotation matrix

      call reclat_loc(cell, celli, 0)

      call reclat_loc(cell2, celli, 0)
      do ia = 1,na
         xac(:) = matmul(transpose(celli),xa(:,ia))
         if(control==1) xac(:) =modulo(xac(:)+100.0d0,1.0d0)
         do ix = 1,3
           xa(ix,ia) = cell(ix,1) * xac(1) + cell(ix,2) * xac(2) + cell(ix,3) * xac(3)
         enddo
      enddo

      end subroutine

      SUBROUTINE reclat_loc (A,B,IOPT)


      integer :: iopt, i
      DOUBLE PRECISION A(3,3),B(3,3), pi, c , ci
      PI=ACOS(-1.D0)
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*PI
      DO 20 I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
  20  CONTINUE
      END

      SUBROUTINE YLMB(V,LMAX,Y)
      INTEGER            LMAX
      DOUBLE PRECISION   V(3)
      COMPLEX*16         Y(*)
!
!     ..................................................................
! 1.     PROGRAM UNIT 'YLM'
!           Calculates spherical harmonics
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           The spherical harmonics (Condon and Shortley convention)
!             Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!           for vector V (given in Cartesian coordinates)
!           are calculated. In the Condon Shortley convention the
!           spherical harmonics are defined as
!                             +------+
!                        m    |   1     m              im(Phi)
!           Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!                            \| 2(Pi)   l
!                  m
!           where P (cos(Theta)) is the normalized Associated Legendre
!                  l
!           function. Thus,
!                                          m      *
!                            Y(l,-m) = (-1) Y(l,m)
!
!
! 3.     USAGE
!           DOUBLE PRECISION V(3), Y(5*5)
!           V(1) = ...
!           V(2) = ...
!           V(3) = ...
!           CALL YLM(V,4,Y)
!
!        ARGUMENT-DESCRIPTION
!           V      - DOUBLE PRECISION vector, dimension 3        (input)
!                    Must be given in Cartesian coordinates.
!                    Conversion of V to polar coordinates gives the
!                    angles Theta and Phi necessary for the calculation
!                    of the spherical harmonics.
!           LMAX   - INTEGER value                               (input)
!                    upper bound of L for which spherical harmonics
!                    will be calculated
!                    constraint:
!                       LMAX .GE. 0 (not checked)
!           Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!                    contains the calculated spherical harmonics
!                    Y(1)                   for L .EQ. 0 (M = 0)
!                    Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!                    ...
!                    Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!                                           for L .EQ. LMAX
!                                              (M = -L,...,L)
!                    constraint:
!                       Dimension of Y .GE. (LMAX+1)**2 (not checked)
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           Type COMPLEX*16 is used which does not conform to the
!           FORTRAN 77 standard.
!           Also the non-standard type conversion function DCMPLX()
!           is used which combines two double precision values into
!           one double complex value.
!
! 4.     REMARKS
!           none
!
! 5.     METHOD
!           The basic algorithm used to calculate the spherical
!           harmonics for vector V is as follows:
!
!           Y(0,0)
!           Y(1,0)
!           Y(1,1)
!           Y(1,-1) = -Y(1,1)
!           DO L = 2, LMAX
!              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!              DO M = L-2, 0, -1
!                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!                 Y(L,-M)= (-1)**M*Y(L,M)
!              ENDDO
!           ENDDO
!
!           In the following the necessary recursion formulas and
!           starting values are given:
!
!        Start:
!                        +------+
!                        |   1
!           Y(0,0) =  -+ | -----
!                       \| 4(Pi)
!
!                                   +------+
!                                   |   3
!           Y(1,0) =  cos(Theta) -+ | -----
!                                  \| 4(Pi)
!
!                                     +------+
!                                     |   3    i(Phi)
!           Y(1,1) =  - sin(Theta) -+ | ----- e
!                                    \| 8(Pi)
!
!        Formula 1:
!
!           Y(l,l) =
!                           +--------+
!                           | (2l+1)   i(Phi)
!            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!                          \|   2l
!
!        Formula 2:
!                                  +---------------+
!                                  |  (2l-1)(2l+1)
!           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!                                 \|   (l-m)(l+m)
!
!                                    +--------------------+
!                                    |(l-1+m)(l-1-m)(2l+1)
!                              -  -+ |-------------------- Y(l-2,m)
!                                   \|  (2l-3)(l-m)(l+m)
!
!        Formula 3: (not used in the algorithm because of the division
!                    by sin(Theta) which may be zero)
!
!                                    +--------------+
!                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!                      sin(Theta)   \| (l+m+1)(l-m)
!
!                                    +--------------+
!                                    |(l-m-1)(l+m+2)  -2i(Phi)
!                              -  -+ |-------------- e        Y(l,m+2)
!                                   \| (l-m)(l+m+1)
!
! 6.     DATE
!           26. April 1994                                   Version 1.2
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13, PI
      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!
      PI = (4.0D+0)*ATAN(1.0D+0)
!
!        Y(0,0)
!
      YLLR = 1.0D+0/SQRT(4.0D+0*PI)
      YLLI = 0.0D+0
      Y(1) = DCMPLX(YLLR,YLLI)
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      IF (LMAX .LE. 0) GOTO 999
!
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0D+0) THEN
         A = V(1)/ABMAX
         B = V(2)/ABMAX
         AB = SQRT(A*A+B*B)
         COSPH = A/AB
         SINPH = B/AB
      ELSE
         COSPH = 1.0D+0
         SINPH = 0.0D+0
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         AB = A*A + B*B
         ABC = SQRT(AB + C*C)
         COSTH = C/ABC
         SINTH = SQRT(AB)/ABC
      ELSE
         COSTH = 1.0D+0
         SINTH = 0.0D+0
      ENDIF
!
!        Y(1,0)
!
      Y(3) = DCMPLX(SQRT(3.0D+0)*YLLR*COSTH,0.0D+0)
!
!        Y(1,1) ( = -DCONJG(Y(1,-1)))
!
      TEMP1 = -SQRT(1.5D+0)*YLLR*SINTH
      Y(4) = DCMPLX(TEMP1*COSPH,TEMP1*SINPH)
      Y(2) = -DCONJG(Y(4))
!
      DO 20 L = 2, LMAX
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         YL1L1R = DBLE(Y(INDEX-1))
         YL1L1I = DIMAG(Y(INDEX-1))
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
         YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
         YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
         Y(INDEX2) = DCMPLX(YLLR,YLLI)
         Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = SQRT(DBLE(2*L+1))*COSTH
         YLL1R = TEMP2*YL1L1R
         YLL1I = TEMP2*YL1L1I
         Y(INDEX2) = DCMPLX(YLL1R,YLL1I)
         Y(INDEX)  = -MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX2 - 4*L + 2
         I2L  = INDEX2 - 2*L
         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO 10 M = L - 2, 0, -1
!
!        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
            YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
            Y(INDEX2) = DCMPLX(YLMR,YLMI)
            Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
!
            MSIGN  = -MSIGN
            INDEX2 = INDEX2 - 1
            INDEX  = INDEX  + 1
            I4L2   = I4L2   - 1
            I2L    = I2L    - 1
   10    CONTINUE
   20 CONTINUE
!
  999 RETURN
!
!        End of 'YLM'
!
      END

      subroutine outcell_output_loc(cell,cellm,celang)



      implicit none

      integer,parameter :: dp = selected_real_kind(precision(1.0d0))
      real(dp), intent(in)   :: cell(3,3)


      integer    :: iv, ix
      real(dp)   :: cellm(3), celang(3), volume, volcel
      real(dp),                      save :: Ang = 1.0d0/0.529177d0
      real(dp),                      save :: pi = 3.14159265358d0
      external volcel

      do iv = 1,3
        cellm(iv) = dot_product(cell(:,iv),cell(:,iv))
        cellm(iv) = sqrt(cellm(iv))
      enddo


      celang(1) = dot_product(cell(:,1),cell(:,2))
      celang(1) = acos(celang(1)/(cellm(1)*cellm(2)))*180._dp/pi
      celang(2) = dot_product(cell(:,1),cell(:,3))
      celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180._dp/pi
      celang(3) = dot_product(cell(:,2),cell(:,3))
      celang(3) = acos(celang(3)/(cellm(2)*cellm(3)))*180._dp/pi

!      cellm=cellm/Ang


      end subroutine outcell_output_loc

DOUBLE PRECISION FUNCTION VOLCEL( C )
DOUBLE PRECISION C(3,3)
VOLCEL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) + ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) + ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
VOLCEL = ABS( VOLCEL )
END

subroutine abc2cell(alp,blp,clp,latt_angle,lattice)
double precision::alp,blp,clp,latt_angle(3),lattice(3,3),alplp,betlp,gamlp,pi
pi=3.1415926535897932384626d0
lattice=0.0d0
alplp=latt_angle(3)* pi/180.d0
betlp=latt_angle(2)* pi/180.d0
gamlp=latt_angle(1)* pi/180.d0
lattice(1,1) = alp
lattice(1,2) = blp * dcos(gamlp)
lattice(2,2) = blp * dsin(gamlp)
lattice(1,3) = clp * dcos(betlp)
lattice(2,3) = (clp*dcos(alplp) - clp*dcos(gamlp)*dcos(betlp))/dsin(gamlp)
lattice(3,3) = dsqrt(clp*clp - lattice(1,3)*lattice(1,3) - lattice(2,3)*lattice(2,3))
end subroutine

       integer function elemz(element)
         implicit none
         character(len=2)::element
         if(trim(element)=='H') then
          elemz = 1
          return


         else if(trim(element)=='Li') then
          elemz = 3
          return


         else if(trim(element)=='B') then
          elemz = 5
          return

         else if(trim(element)=='N') then
          elemz = 7
          return

         else if(trim(element)=='F') then
          elemz = 9
          return

         else if(trim(element)=='Na') then
          elemz = 11
          return

         else if(trim(element)=='Al') then
          elemz = 13
          return

         else if(trim(element)=='P') then
          elemz = 15
          return

         else if(trim(element)=='Cl') then
          elemz = 17
          return

         else if(trim(element)=='K') then
          elemz = 19
          return

         else if(trim(element)=='Sc') then
          elemz = 21
          return

         else if(trim(element)=='Sc') then
          elemz = 22
          return

         else if(trim(element)=='V') then
          elemz = 23
          return

         else if(trim(element)=='Mn') then
          elemz = 25
          return

         else if(trim(element)=='Co') then
          elemz = 27
          return

         else if(trim(element)=='Cu') then
          elemz = 29
          return

         else if(trim(element)=='Ga') then
          elemz = 31
          return

         else if(trim(element)=='As') then
          elemz = 33
          return

         else if(trim(element)=='Br') then
          elemz = 35
          return

         else if(trim(element)=='Rb') then
          elemz = 37
          return

         else if(trim(element)=='Y') then
          elemz = 39
          return

         else if(trim(element)=='Nb') then
          elemz = 41
          return

         else if(trim(element)=='Tc') then
          elemz = 43
          return

         else if(trim(element)=='Rh') then
          elemz = 45
          return

         else if(trim(element)=='Ag') then
          elemz = 47
          return

         else if(trim(element)=='In') then
          elemz = 49
          return

         else if(trim(element)=='Sb') then
          elemz = 51
          return

         else if(trim(element)=='I') then
          elemz = 53
          return

         else if(trim(element)=='Cs') then
          elemz = 55
          return

         else if(trim(element)=='La') then
          elemz = 57
          return

         else if(trim(element)=='Pr') then
          elemz = 59
          return

         else if(trim(element)=='Pm') then
          elemz = 61
          return

         else if(trim(element)=='Eu') then
          elemz = 63
          return

         else if(trim(element)=='Tb') then
          elemz = 65
          return

         else if(trim(element)=='Ho') then
          elemz = 67
          return

         else if(trim(element)=='Tm') then
          elemz = 69
          return

         else if(trim(element)=='Lu') then
          elemz = 71
          return

         else if(trim(element)=='Ta') then
          elemz = 73
          return

         else if(trim(element)=='Re') then
          elemz = 75
          return

         else if(trim(element)=='Ir') then
          elemz = 77
          return

         else if(trim(element)=='Au') then
          elemz = 79
          return

         else if(trim(element)=='Tl') then
          elemz = 81
          return

         else if(trim(element)=='Bi') then
          elemz = 83
          return

         else if(trim(element)=='At') then
          elemz = 85
          return

         else if(trim(element)=='Fr') then
          elemz = 87
          return

         else if(trim(element)=='Ac') then
          elemz = 89
          return

         else if(trim(element)=='Pa') then
          elemz = 91
          return

         else if(trim(element)=='Np') then
          elemz = 93
          return

         else if(trim(element)=='Am') then
          elemz = 95
          return

         else if(trim(element)=='Bk') then
          elemz = 97
          return

         else if(trim(element)=='Es') then
          elemz = 99
          return

         else if(trim(element)=='Md') then
          elemz = 101
          return

         else if(trim(element)=='Lw') then
          elemz = 103
          return

         else if(trim(element)=='He') then
          elemz = 2
          return

         else if(trim(element)=='Be') then
          elemz = 4
          return

         else if(trim(element)=='C') then
          elemz = 6
          return

         else if(trim(element)=='O') then
          elemz = 8
          return

         else if(trim(element)=='Ne') then
          elemz = 10
          return

         else if(trim(element)=='Mg') then
          elemz = 12
          return

         else if(trim(element)=='Si') then
          elemz = 14
          return

         else if(trim(element)=='S') then
          elemz = 16
          return

         else if(trim(element)=='Ar') then
          elemz = 18
          return

         else if(trim(element)=='Ca') then
          elemz = 20
          return

         else if(trim(element)=='Ti') then
          elemz = 22
          return

         else if(trim(element)=='Cr') then
          elemz = 24
          return

         else if(trim(element)=='Fe') then
          elemz = 26
          return

         else if(trim(element)=='Ni') then
          elemz = 28
          return

         else if(trim(element)=='Zn') then
          elemz = 30
          return

         else if(trim(element)=='Ge') then
          elemz = 32
          return

         else if(trim(element)=='Se') then
          elemz = 34
          return

         else if(trim(element)=='Kr') then
          elemz = 36
          return

         else if(trim(element)=='Sr') then
          elemz = 38
          return

         else if(trim(element)=='Zr') then
          elemz = 40
          return

         else if(trim(element)=='Mo') then
          elemz = 42
          return

         else if(trim(element)=='Ru') then
          elemz = 44
          return

         else if(trim(element)=='Pd') then
          elemz = 46
          return

         else if(trim(element)=='Cd') then
          elemz = 48
          return

         else if(trim(element)=='Sn') then
          elemz = 50
          return

         else if(trim(element)=='Te') then
          elemz = 52
          return

         else if(trim(element)=='Xe') then
          elemz = 54
          return

         else if(trim(element)=='Ba') then
          elemz = 56
          return

         else if(trim(element)=='Ce') then
          elemz = 58
          return

         else if(trim(element)=='Nd') then
          elemz = 60
          return

         else if(trim(element)=='Sm') then
          elemz = 62
          return

         else if(trim(element)=='Gd') then
          elemz = 64
          return

         else if(trim(element)=='Dy') then
          elemz = 66
          return

         else if(trim(element)=='Er') then
          elemz = 68
          return

         else if(trim(element)=='Yb') then
          elemz = 70
          return

         else if(trim(element)=='Hf') then
          elemz = 72
          return

         else if(trim(element)=='W') then
          elemz = 74
          return

         else if(trim(element)=='Os') then
          elemz = 76
          return

         else if(trim(element)=='Pt') then
          elemz = 78
          return

         else if(trim(element)=='Hg') then
          elemz = 80
          return

         else if(trim(element)=='Pb') then
          elemz = 82
          return

         else if(trim(element)=='Po') then
          elemz = 84
          return

         else if(trim(element)=='Rn') then
          elemz = 86
          return

         else if(trim(element)=='Ra') then
          elemz = 88
          return

         else if(trim(element)=='Th') then
          elemz = 90
          return

         else if(trim(element)=='U') then
          elemz = 92
          return

         else if(trim(element)=='Pu') then
          elemz = 94
          return

         else if(trim(element)=='Cm') then
          elemz = 96
          return

         else if(trim(element)=='Cf') then
          elemz = 98
          return

         else if(trim(element)=='Fm') then
          elemz = 100
          return

         else if(trim(element)=='No') then
          elemz = 102
          return
         else
          elemz = 0
          return
         endif

        end function


      subroutine Get_dis_weight_order_parameter(na,iza,cell,xa,atom,Qglobal, ene, sym)

      implicit none
      integer:: namx !,lmax
      parameter(namx=100)
      integer::na,ntol,i,j,i1,i2,i3,k,bond,l,m,super(3),km,atom,totalbond,control,Iiza
      integer :: iza(na),maxiza,miniza, sym
      double precision,allocatable , save ::xac(:,:),xa0(:,:),r0(:,:)
      double precision xa(3,na),r1,cell0(3,3)
      double precision Q1,co1,co2
      double precision alpha,Pi,Qlocal(4),Qglobal(4)
      double precision dist_all,dist,dtmp,r(3),cell(3,3),celli(3,3)
      double precision rij(namx*namx),theta(namx*namx),phi(namx*namx),radi(namx*namx),rmax
      double precision :: ene
      logical, save :: first=.true.
      COMPLEX*16  Y(49),Q(49),Qall(49)
      logical CONT
      integer:: imode=1

      save maxiza,miniza

      Pi=3.14159265359d0
      alpha=1.0d0
      rmax=10.0d0

      if(first) then
      allocate(xac(3,na),xa0(3,na),r0(na,na))
      first = .false.
      maxiza=0
      miniza=100
      do i=1,na
        do j=1,na
              call fastbond(iza(i),iza(j),r0(i,j))
              if(r0(i,j) < 0.05d0) then
                  call species_radius(iza(i),co1)
                  call species_radius(iza(j),co2)
                  r0(i,j) = co1 + co2
                  if(r0(i,j) < 0.5d0) then
                     !write(*,*) 'missing species',iza(i),iza(j)
                    r0(i,j)=2d0
                  endif
              end if
        enddo
      enddo
      endif
      if(atom==0) then
       Iiza=0 !maxiza
      else
       Iiza=atom
      endif

      xa0=xa
      cell0=cell
      if(Imode>=0) then
      call cellconstraint(na,cell0,xa0,1)
      endif


      call reclat_loc(cell0, celli, 0)
      do i = 1,na
         xac(:,i) = matmul(transpose(celli),xa0(:,i))
         xac(:,i)= modulo(xac(:,i)+100.0d0,1.0d0)
      enddo

      do k = 1,na
        do j = 1,3
          xa0(j,k) = cell0(j,1) * xac(1,k) + cell0(j,2) * xac(2,k) + cell0(j,3) * xac(3,k)
        enddo
      enddo

      CONT=.false.

      Qall=0.0d0
      Qlocal=0d0
      totalbond=0
      super=0
      if(Imode>=0) super(:)= 6   ! supercell size

      do km=1 , na
        bond=1
         if(iza(km)==Iiza .or. Iiza==0) then !goto 90
         j=km
         Q=0d0
         do i=1,na
              do i1=-super(1),super(1)
              do i2=-super(2),super(2)
              do i3=-super(3),super(3)
                    r(:)=(xac(1,i)+dble(i1))*cell0(:,1)+(xac(2,i)+ &
                          dble(i2))*cell0(:,2)+(xac(3,i)+dble(i3))*cell0(:,3) - xa0(:,j)
                    r1=dsqrt(sum(r*r))

                    if(r1 > 20d0) cycle  !cutoff

                    if(r1>r0(i,j)*0.6d0 .and. r1<r0(i,j)*1.3d0) bond=bond+1
                    if(r1>r0(i,j)*0.6d0) then
                    r1=dexp(-0.5d0*(max(0d0,(r1-r0(i,j)*0.6d0)/(r0(i,j)*0.6d0))))
                    if(r1>1d-4) then ! 2.0d0*r0(i,j)+0.5d0 .and. r1 > 0.5d0) then
                    call YLMB(r,6,Y) !sphe_harm(l,m,theta(i),phi(i))
                    Q=Q+Y*r1/0.71d0     ! distance weighted  0.71 is a prefactor
                    Qall=Qall+Y !*r1
                     endif
                     endif
              enddo
              enddo
              enddo
!        endif
       enddo

       if(bond>1) bond=bond-1

       do j=0,6,2
       Q1=0d0
       do i=j*j+1,(j+1)*(j+1)
       Q1=Q1+Q(i)*DCONJG(Q(i))
       enddo
       Q1=sqrt(4.0d0*Pi*Q1/(2.0d0*dble(j)+1.0d0))
       Q1=Q1/dble(na)/dble(bond)
       Qlocal(j/2+1)=Q1+Qlocal(j/2+1)
       enddo

      endif ! Iiza
      enddo  ! km

      Qglobal=Qlocal

      End subroutine
