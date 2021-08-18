c     =============================
c     MODULE M_CONTCAR  -  begin
c     =============================
      module m_contcar
c  MODULE FOR MANIPULATE VASP GEOMETRY FILES 
c                                            
c  written by Matheus Paes Lima (07/05/2017) 
c                                            
c  menu: 
c    0) type(contcar)                (06/04/2017)
c    1) read_contcar(var,fname)      (06/04/2017)
c    2) CONTCAR_to_xyz(var)          (06/04/2017)                           
c    3) write_contcar(var,fname)     (07/04/2017)                     
c    4) direct_to_Cartesian(var)     (07/04/2017)                                    
c    5) Cartesian_to_direct(var)     (07/04/2017)
c    6) translate                    (07/04/2017)                
c    7) put_selective(var)           (07/04/2017)
c    8) put_in_box(var)              (07/04/2017)
c    9) erease_atom(var,i)           (10/04/2017)
c   10) supercell(var,nx,ny,nz)      (10/04/2017)
c   11) locpot skills                (19/06/2017)
c   12) supercellM                   (25/07/2018)
c   13) RotateCell                   (25/07/2018)
c   14) Rotate                       (21/03/2019)
c   15) add_atom                     (01/08/2020)
c   16) purge_contcar                (11/08/2020)
c   17) write_contcarF               (11/08/2020)
c   18) CONTCAR_to_xyzF              (11/08/2020)                           
c   19) translate_list(var,t,list,n) (18/08/2021)
c
c
c
c
c
      implicit none
c     variables of a POSCAR or a CONTCAR file of VASP
c     if necessary, more variables can be easly added.
c     The subroutine read_contcar allocate and charge 
c     these values from a POSCAR file named in "fname"
      type, public :: contcar
          character*64              :: head
          real*8                    :: a0
          real*8                    :: cell(3,3)
          integer                   :: nspecs
          character*2, allocatable  :: atm_name(:)
          integer, allocatable      :: atm_name_n(:)
          logical                   :: scaled
          real*8, allocatable       :: x(:), y(:), z(:)
          logical                   :: selec_dyn
          character*1, allocatable  :: SD_x(:), SD_y(:), SD_z(:)
          logical                   :: grid=.false.
          integer                   :: ngx, ngy, ngz
          real*8, allocatable       :: v(:)
      end type contcar


      contains

c     =============================
c     subroutine read_contcar -  begin
c     =============================
      subroutine read_contcar(var,fname)
      implicit none
      type(contcar) :: var
      character*64  :: fname
      
c     .: internal variables 
      character*128 :: string
      character*1   :: st(128), Cdummy
      integer       :: nwords, ierr, ntot, i
      logical       :: long
      
      long=.true.
      long=.false.

      open(21,file=fname,status="unknown")
      
      read(21,'(a128)') string
      var%head=""
      read(string,'(a64)') var%head
          if (long) write(*,*) var%head
      
      read(21,'(a128)') string
      read(string,*) var%a0
          if (long) write(*,16) var%a0
      
      read(21,'(a128)') string
      read(string,*) var%cell(1,:)
          if (long) write(*,16) var%cell(1,:)

      read(21,'(a128)') string
      read(string,*) var%cell(2,:)
          if (long) write(*,16) var%cell(2,:)
   
      read(21,'(a128)') string
      read(string,*) var%cell(3,:)
          if (long) write(*,*) var%cell(3,:)
      
      read(21,'(a128)') string
      st=""
      read (string,*,iostat=ierr) st
      nwords = count(st /= "")
      var%nspecs=nwords
      allocate(var%atm_name(nwords),var%atm_name_n(nwords))
      read(string,*) var%atm_name(:)
          if (long) write(*,'(100a5)')  var%atm_name(:) 

      read(21,'(a128)') string
      read(string,*) var%atm_name_n(:)
          if (long) write(*,*) var%atm_name_n(:) 
      ntot=sum(var%atm_name_n(:))
      allocate(var%x(ntot),var%y(ntot),var%z(ntot))
      allocate(var%SD_x(ntot),var%SD_y(ntot),var%SD_z(ntot))
      var%SD_x(:)="T";var%SD_y(:)="T"; var%SD_z(:)="T" 
      read(21,'(a128)') string
      read(string,*) cdummy
      var%selec_dyn=.false.

      if ((cdummy=="s").or.(cdummy=="S")) then
          var%selec_dyn=.true.
          if (long) write(*,*)  "selective dynamics was found"
          read(21,'(a128)') string
          read(string,*) cdummy
          if ((cdummy=="C").or.(cdummy=="c")) then
          if (long) write(*,*)  "Cartesian corrdinates was found"
          var%scaled=.false.
          elseif ((cdummy=="D").or.(cdummy=="d")) then
          if (long) write(*,*)  "Direct coordinates was found"
          var%scaled=.true.
          else
          stop "unknown option CONTCAR"
          endif
          do i=1,ntot
              read(21,'(a128)') string    
              read(string,*) var%x(i), var%y(i), var%z(i), 
     &                       var%SD_x(i), var%SD_y(i), var%SD_z(i)
          enddo
      else
          if (long) write(*,*) "selective dynamics was NOT found"
          if ((cdummy=="C").or.(cdummy=="c")) then
          if (long) write(*,*) "Cartesian coordinates was found"
          var%scaled=.false.
          elseif ((cdummy=="D").or.(cdummy=="d")) then
          if (long) write(*,*) "direct coordinates was found"
          var%scaled=.true.
          else
          stop "unknown option CONTCAR"
          endif
          do i=1,ntot
              read(21,'(a128)') string    
              read(string,*) var%x(i), var%y(i), var%z(i)
          enddo
      endif

      if (long) then
      do i=1,ntot
      if (var%selec_dyn) then
      write(*,15) var%x(i),var%y(i),var%z(i),
     &        var%SD_x(i),var%SD_y(i),var%SD_z(i)            
      else
      write(*,15) var%x(i),var%y(i),var%z(i)
      endif
      enddo
      endif
15    format(3f21.12,3a5)
16    format(3f21.12)
      close(21)
      return
      end subroutine read_contcar

      
      
      
      
      
c     =============================
c     subroutine CONTCAR_to_xyz  -  begin
c     =============================      
      
      
      subroutine CONTCAR_to_xyz(var)
      implicit none
      type(contcar) var
      
C     .: internal      
      integer i, j, icount,itot
      itot=sum(var%atm_name_n(:))
      write(*,*) itot
      write(*,*) " generated from a POSCAR file "
      icount=1
      do i=1,var%nspecs
         do j=1,var%atm_name_n(i)
             if (var%scaled) then
             write(*,17) var%atm_name(i),
     &       var%x(icount)*var%cell(1,1)+
     &       var%y(icount)*var%cell(2,1)+
     &       var%z(icount)*var%cell(3,1),
     &       var%x(icount)*var%cell(1,2)+
     &       var%y(icount)*var%cell(2,2)+
     &       var%z(icount)*var%cell(3,2),
     &       var%x(icount)*var%cell(1,3)+
     &       var%y(icount)*var%cell(2,3)+
     &       var%z(icount)*var%cell(3,3)
             else
             write(*,17) var%atm_name(i), var%x(icount), 
     &                 var%y(icount),var%z(icount) 
            endif
         icount=icount+1
         enddo      
      enddo
      
17    format(a2,5x,3f21.12)      
      
      end subroutine CONTCAR_to_xyz

c     =============================
c     subroutine write_contcar -  begin
c     =============================      
      
      subroutine write_contcar(var)
      implicit none
      type(contcar) :: var
      integer       :: ntot, icount, i, j
       write(*,'(a64)') var%head
       write(*,16) var%a0
       write(*,16) var%cell(1,:)
       write(*,16) var%cell(2,:)
       write(*,16) var%cell(3,:)
       write(*,'(100a5)') var%atm_name(:) 
       write(*,'(100i5)') var%atm_name_n(:) 
       if (var%selec_dyn)   write(*,'(a18)') "Selective Dynamics"
       if (.not.var%scaled) write(*,'(a9)') "Cartesian"
       if (var%scaled)      write(*,'(a6)') "Direct"
       ntot=sum(var%atm_name_n(:))
       icount=1
       if (var%selec_dyn) then
        do i=1,var%nspecs
           do j=1,var%atm_name_n(i)
               write(*,17) var%x(icount),var%y(icount),var%z(icount),
     &         var%SD_x(icount),var%SD_z(icount),var%SD_z(icount)
               icount=icount+1
           enddo      
        enddo
        else
        do i=1,var%nspecs
           do j=1,var%atm_name_n(i)
               write(*,17) var%x(icount),var%y(icount),var%z(icount) 
               icount=icount+1
           enddo      
        enddo
        
        endif

15    format(3f21.12,3a5)
16    format(3f21.12)
17    format(3f21.12,3x,a1,3x,a1,3x,a1)     
      end subroutine write_contcar
 
c     =============================
c     subroutine direct_to_Cartesian  -  begin
c     ============================= 
 
      subroutine direct_to_Cartesian(var)
      implicit none
      type(contcar)       :: var
      integer             :: ntot, i
      real*8              :: x1, y1, z1
      
      if (.not.var%scaled) stop 'geometry already in Cartesian coord'
      
      var%scaled=.false.
      ntot=sum(var%atm_name_n(:))
      do i=1,ntot
      x1=var%x(i)*var%cell(1,1)+
     &   var%y(i)*var%cell(2,1)+
     &   var%z(i)*var%cell(3,1)
      y1=var%x(i)*var%cell(1,2)+
     &   var%y(i)*var%cell(2,2)+
     &   var%z(i)*var%cell(3,2)
      z1=var%x(i)*var%cell(1,3)+
     &   var%y(i)*var%cell(2,3)+
     &   var%z(i)*var%cell(3,3)
      var%x(i)=x1 ; var%y(i)=y1 ; var%z(i)=z1
      enddo
      end subroutine direct_to_Cartesian
      
c     =============================
c     subroutine Cartesian_to_direct  -  begin
c     =============================

      subroutine Cartesian_to_direct(var)
      implicit none
      type(contcar)       :: var
      integer             :: ntot, i
      real*8              :: x1(3),xd(3)
      real*8              :: cellT(3,3)
      
      if (var%scaled) stop 'geometry already in Direct coord'
      
      var%scaled=.true.
      ntot=sum(var%atm_name_n(:))
      do i=1,ntot
      
      cellT(1,1)=var%cell(1,1)
      cellT(1,2)=var%cell(2,1)
      cellT(1,3)=var%cell(3,1)
      
      cellT(2,1)=var%cell(1,2)
      cellT(2,2)=var%cell(2,2)
      cellT(2,3)=var%cell(3,2)
      
      cellT(3,1)=var%cell(1,3)
      cellT(3,2)=var%cell(2,3)
      cellT(3,3)=var%cell(3,3)
      
      x1(1)=var%x(i);x1(2)=var%y(i);x1(3)=var%z(i)
      call LU(cellT,x1,xd,3)     
      var%x(i)=xd(1) ; var%y(i)=xd(2) ; var%z(i)=xd(3)
      enddo
      
      end subroutine Cartesian_to_direct
 
c     =============================
c     subroutine translate -  begin
c     ============================= 
 
 
      subroutine translate(var,t)
      implicit none
      type(contcar)       :: var
      real*8              :: t(3), cellT(3,3)
      
c     .: internal variables
      real*8              :: td(3)
      
      if (.not.var%scaled) then
          var%x(:)=var%x(:)+t(1)
          var%y(:)=var%y(:)+t(2)
          var%z(:)=var%z(:)+t(3)
      else
          cellT(1,1)=var%cell(1,1)
          cellT(1,2)=var%cell(2,1)
          cellT(1,3)=var%cell(3,1)
      
          cellT(2,1)=var%cell(1,2)
          cellT(2,2)=var%cell(2,2)
          cellT(2,3)=var%cell(3,2)
      
          cellT(3,1)=var%cell(1,3)
          cellT(3,2)=var%cell(2,3)
          cellT(3,3)=var%cell(3,3)   
          
          call LU(cellT,t,td,3)
          
          var%x(:)=var%x(:)+td(1)
          var%y(:)=var%y(:)+td(2)
          var%z(:)=var%z(:)+td(3)
          
      endif      
      return
      end subroutine translate
      
      
c     =============================
c     subroutine put_Selective -  begin
c     =============================      
      subroutine put_Selective(var)
      implicit none
      type(contcar)       :: var
      
      if (var%selec_dyn) stop 'CONTCAR already with selective'
      var%selec_dyn=.true.
      
      end subroutine put_Selective

      
c     =============================
c     subroutine put_in_box -  begin
c     =============================      
      subroutine put_in_box(var)
      implicit none
      type(contcar) var
      logical isin, cart
      integer ntot, i
      
      ntot=sum(var%atm_name_n(:))
 
      cart=(.not.var%scaled)
      if (cart) call Cartesian_to_direct(var)
      do i=1,ntot
c         for x    
          if (var%x(i).ge.1.d0) var%x(i)=var%x(i)-1.d0 
          if (var%x(i).lt.0.d0) var%x(i)=var%x(i)+1.d0
c         for y      
          if (var%y(i).ge.1.d0) var%y(i)=var%y(i)-1.d0 
          if (var%y(i).lt.0.d0) var%y(i)=var%y(i)+1.d0
c         for z      
          if (var%z(i).ge.1.d0) var%z(i)=var%z(i)-1.d0 
          if (var%z(i).lt.0.d0) var%z(i)=var%z(i)+1.d0
      enddo
      if (cart) call direct_to_Cartesian(var)
      return
      end subroutine put_in_box
      
c     =============================
c     subroutine put_in_box -  begin
c     =============================  
      subroutine erease_atom(var,iatm)
      implicit none
      type(contcar)         ::   var 
      integer               ::   iatm
      
c     .: internal variables
      integer               :: i, icount, ic1, ic2
      integer               :: ntot
      double precision,allocatable:: xN(:),yN(:),zN(:)
      character*1,allocatable:: SDxN(:),SDyN(:),SDzN(:)
      
c     .: encontrando especie 
      ic1=1
      do i=1,var%nspecs
      ic2=ic1+var%atm_name_n(i)
c      print*, i,ic1,ic2,iatm
      if((iatm.ge.ic1).and.(iatm.lt.ic2))then
c         print*, "apagar especie", i, var%atm_name(i)
          exit
      endif
      ic1=ic2
      enddo
      write(*,*) "apagar atomo", iatm," especie", i, var%atm_name(i)
      
c     salvando posicoes 
      ntot=sum(var%atm_name_n(:))
      allocate(xN(ntot-1),yN(ntot-1),zN(ntot-1))
      allocate(SDxN(ntot-1),SDyN(ntot-1),SDzN(ntot-1))
      xN(1:iatm-1)=var%x(1:iatm-1)
      yN(1:iatm-1)=var%y(1:iatm-1)
      zN(1:iatm-1)=var%z(1:iatm-1)
      
      xN(iatm:ntot-1)=var%x(iatm+1:ntot)
      yN(iatm:ntot-1)=var%y(iatm+1:ntot)
      zN(iatm:ntot-1)=var%z(iatm+1:ntot)
      
      SDxN(1:iatm-1)=var%SD_x(1:iatm-1)
      SDyN(1:iatm-1)=var%SD_y(1:iatm-1)
      SDzN(1:iatm-1)=var%SD_z(1:iatm-1)
      
      SDxN(iatm:ntot-1)=var%SD_x(iatm+1:ntot)
      SDyN(iatm:ntot-1)=var%SD_y(iatm+1:ntot)
      SDzN(iatm:ntot-1)=var%SD_z(iatm+1:ntot)
      
      deallocate(var%x,var%y,var%z)
      allocate(var%x(ntot-1),var%y(ntot-1),var%z(ntot-1))
      var%x(1:ntot-1)=xN(1:ntot-1)
      var%y(1:ntot-1)=yN(1:ntot-1)
      var%z(1:ntot-1)=zN(1:ntot-1)
      
      var%SD_x(1:ntot-1)=SDxN(1:ntot-1)
      var%SD_y(1:ntot-1)=SDyN(1:ntot-1)
      var%SD_z(1:ntot-1)=SDzN(1:ntot-1)      
      
      var%atm_name_n(i)=var%atm_name_n(i)-1
      
      deallocate(xN,yN,zN)
      deallocate(SDxN,SDyN,SDzN)
      return
      end subroutine erease_atom
      
c     =============================
c     subroutine supercell -  begin
c     =============================      
      subroutine supercell(var,nx,ny,nz)
      implicit none
      type(contcar)        :: var
      integer              :: nx, ny, nz
      
c     .: internal varibles
      integer              :: nat, natS
      integer              :: ix, iy, iz, icount,iatm
      real*8, allocatable  :: xS(:),yS(:),zS(:)
      character*1,allocatable :: SDx(:),SDy(:),SDz(:)
      
      call direct_to_Cartesian(var)
      
      nat=sum(var%atm_name_n(:))
      natS=nat*nx*ny*nz
      
c  .: 1) modifying geometry      
      allocate(xS(natS),yS(natS),zS(natS))
      allocate(SDx(natS),SDy(natS),SDz(natS))
      
      icount=1
      do iatm=1,nat
          do ix=0,nx-1
          do iy=0,ny-1
          do iz=0,nz-1
          xS(icount)=var%x(iatm)+ix*var%cell(1,1)
     .                          +iy*var%cell(2,1)
     .                          +iz*var%cell(3,1)
          yS(icount)=var%y(iatm)+ix*var%cell(1,2)
     .                          +iy*var%cell(2,2)
     .                          +iz*var%cell(3,2)
          zS(icount)=var%z(iatm)+ix*var%cell(1,3)
     .                          +iy*var%cell(2,3)
     .                          +iz*var%cell(3,3)
          SDx(icount)=var%SD_x(iatm)
          SDy(icount)=var%SD_y(iatm)
          SDz(icount)=var%SD_z(iatm)
          icount=icount+1
          enddo !loop in X
          enddo !loop in Y
          enddo !loop in Z
      enddo !loop on the atoms
      
      deallocate(var%x,var%y,var%z)
      deallocate(var%SD_x,var%SD_y,var%SD_z)
      allocate(var%x(natS),var%y(natS),var%z(natS))
      allocate(var%SD_x(natS),var%SD_y(natS),var%SD_z(natS))
      var%x(:)=xS(:); var%y(:)=yS(:); var%z(:)=zS(:)
      var%SD_x(:)=SDx(:);var%SD_y(:)=SDy(:);var%SD_z(:)=SDz(:)

c     2) modifying the number of atoms
      var%atm_name_n(:)=var%atm_name_n(:)*nx*ny*nz
      
c  .: 3) enlarge box      
      var%cell(1,:)=var%cell(1,:)*nx
      var%cell(2,:)=var%cell(2,:)*ny
      var%cell(3,:)=var%cell(3,:)*nz

      call Cartesian_to_direct(var)
      
      deallocate(xS,yS,zS,SDx,SDy,SDz)
      return
      end subroutine supercell

c     =============================
c     subroutine read_locpot -  begin
c     =============================
      subroutine read_locpot(var,fname)
      implicit none
      type(contcar) :: var
      character*64  :: fname
      
c     .: internal variables 
      character*128 :: string
      character*1   :: st(128), Cdummy
      integer       :: nwords, ierr, ntot, i, itot
      logical       :: long
      
      long=.true.
C      long=.false.

      open(21,file=fname,status="unknown")
      
      read(21,'(a128)') string
      var%head=""
      read(string,'(a64)') var%head
          if (long) write(*,*) var%head
      
      read(21,'(a128)') string
      read(string,*) var%a0
          if (long) write(*,16) var%a0
      
      read(21,'(a128)') string
      read(string,*) var%cell(1,:)
          if (long) write(*,16) var%cell(1,:)

      read(21,'(a128)') string
      read(string,*) var%cell(2,:)
          if (long) write(*,16) var%cell(2,:)
   
      read(21,'(a128)') string
      read(string,*) var%cell(3,:)
          if (long) write(*,*) var%cell(3,:)
      
      read(21,'(a128)') string
      st=""
      read (string,*,iostat=ierr) st
      nwords = count(st /= "")
      var%nspecs=nwords
      allocate(var%atm_name(nwords),var%atm_name_n(nwords))
      read(string,*) var%atm_name(:)
          if (long) write(*,'(100a5)')  var%atm_name(:) 

      read(21,'(a128)') string
      read(string,*) var%atm_name_n(:)
          if (long) write(*,*) var%atm_name_n(:) 
      ntot=sum(var%atm_name_n(:))
      allocate(var%x(ntot),var%y(ntot),var%z(ntot))
      allocate(var%SD_x(ntot),var%SD_y(ntot),var%SD_z(ntot))
      var%SD_x(:)="T";var%SD_y(:)="T"; var%SD_z(:)="T" 
      read(21,'(a128)') string
      read(string,*) cdummy
      var%selec_dyn=.false.

      if ((cdummy=="s").or.(cdummy=="S")) then
          var%selec_dyn=.true.
          if (long) write(*,*)  "selective dynamics was found"
          read(21,'(a128)') string
          read(string,*) cdummy
          if ((cdummy=="C").or.(cdummy=="c")) then
          if (long) write(*,*)  "Cartesian corrdinates was found"
          var%scaled=.false.
          elseif ((cdummy=="D").or.(cdummy=="d")) then
          if (long) write(*,*)  "Direct coordinates was found"
          var%scaled=.true.
          else
          stop "unknown option CONTCAR"
          endif
          do i=1,ntot
              read(21,'(a128)') string    
              read(string,*) var%x(i), var%y(i), var%z(i), 
     &                       var%SD_x(i), var%SD_y(i), var%SD_z(i)
          enddo
      else
          if (long) write(*,*) "selective dynamics was NOT found"
          if ((cdummy=="C").or.(cdummy=="c")) then
          if (long) write(*,*) "Cartesian coordinates was found"
          var%scaled=.false.
          elseif ((cdummy=="D").or.(cdummy=="d")) then
          if (long) write(*,*) "direct coordinates was found"
          var%scaled=.true.
          else
          stop "unknown option CONTCAR"
          endif
          do i=1,ntot
              read(21,'(a128)') string    
              read(string,*) var%x(i), var%y(i), var%z(i)
          enddo
      endif

      if (long) then
      do i=1,ntot
      if (var%selec_dyn) then
      write(*,15) var%x(i),var%y(i),var%z(i),
     &        var%SD_x(i),var%SD_y(i),var%SD_z(i)            
      else
      write(*,15) var%x(i),var%y(i),var%z(i)
      endif
      enddo
      endif

C     reading grid
      var%grid=.true.
      write(*,*) var%grid
      read(21,'(a128)') string    
      read(21,'(a128)') string    
      read(string,*) var%ngx, var%ngy, var%ngz
      itot=var%ngx*var%ngy*var%ngz
      if (long) write(*,*) "var%ngx, var%ngy, var%ngz, itot"
      if (long) write(*,*) var%ngx, var%ngy, var%ngz, itot
      allocate(var%v(var%ngx*var%ngy*var%ngz))
      read(21,'(5e18.11)') var%v(:)

15    format(3f21.12,3a5)
16    format(3f21.12)
      close(21)
      return
      end subroutine read_locpot

            
c     =============================
c     subroutine supercellM -  begin
c     =============================      
      subroutine supercellM(var,M)
      implicit none
      
      type(contcar)        :: var
      real*8               :: M(3,3)
      
c  .: internal varibles
      real*8               :: M2(3,3), Q(3,3), cellSC(3,3)
      integer,parameter    :: maxX=50, maxy=50, maxZ=50
      integer              :: nat
      integer              :: ix, iy, iz, icount,iatm
      real*8               :: xS,yS,zS
      real*8               :: xSl,ySl,zSl
      real*8, parameter    :: delta=1.d-2
      integer              :: Ndim, nsc
      real*8, allocatable  :: xSC(:,:), ySC(:,:), zSC(:,:)
      logical              :: dir

c  .: calculating the number of atoms    
      nat=sum(var%atm_name_n(:))
      
      Ndim=(2*maxX+1)*(2*maxX+1)*(2*maxX+1)
      allocate(xSC(Ndim,nat),ySC(Ndim,nat),zSC(Ndim,nat))
c      print*, "(M) supercell matrix"
c      print*, M(1,:)
c      print*, M(2,:)
c      print*, M(3,:)

c  .: direct coordinated are required
      dir = var%scaled
      if (.not.dir) call Cartesian_to_direct(var)

c  .: calculating Q, the inverse of M
      M2=M
      call inverse(M2,Q,3)
c      print*, "(Q) supercell inverse matrix"
c      print*, Q(1,:)
c      print*, Q(2,:)
c      print*, Q(3,:)
      
c  .: calculating supercell vectors
      cellSC(1,1)=M(1,1)*var%cell(1,1)
     .           +M(2,1)*var%cell(2,1)+M(3,1)*var%cell(3,1)
      cellSC(1,2)=M(1,1)*var%cell(1,2)
     .           +M(2,1)*var%cell(2,2)+M(3,1)*var%cell(3,2)
      cellSC(1,3)=M(1,1)*var%cell(1,3)
     .           +M(2,1)*var%cell(2,3)+M(3,1)*var%cell(3,3)
      
      cellSC(2,1)=M(1,2)*var%cell(1,1)
     .           +M(2,2)*var%cell(2,1)+M(3,2)*var%cell(3,1)
      cellSC(2,2)=M(1,2)*var%cell(1,2)
     .           +M(2,2)*var%cell(2,2)+M(3,2)*var%cell(3,2)
      cellSC(2,3)=M(1,2)*var%cell(1,3)
     .           +M(2,2)*var%cell(2,3)+M(3,2)*var%cell(3,3)
      
      cellSC(3,1)=M(1,3)*var%cell(1,1)
     .           +M(2,3)*var%cell(2,1)+M(3,3)*var%cell(3,1)
      cellSC(3,2)=M(1,3)*var%cell(1,2)
     .           +M(2,3)*var%cell(2,2)+M(3,3)*var%cell(3,2)
      cellSC(3,3)=M(1,3)*var%cell(1,3)
     .           +M(2,3)*var%cell(2,3)+M(3,3)*var%cell(3,3)
      
c      print*,"supercell vectors"
c      print*, cellSC(1,:)
c      print*, cellSC(2,:)
c      print*, cellSC(3,:)
      
 
      
c  .: 1) modifying geometry      
      icount=0
      nsc=0

          do ix=-maxX,maxX
          do iy=-maxY,maxY
          do iz=-maxZ,maxZ
          xS=ix; yS=iy; zS=iz
          
c          xS=var%x(iatm)+ix!*var%cell(1,1)
c          yS=var%y(iatm)+iy!*var%cell(1,2)
c          zS=var%z(iatm)+iz!*var%cell(1,3)

          xSl=Q(1,1)*xS+Q(1,2)*yS+Q(1,3)*zS
          ySl=Q(2,1)*xS+Q(2,2)*yS+Q(2,3)*zS
          zSl=Q(3,1)*xS+Q(3,2)*yS+Q(3,3)*zS
c              print*,"=================="
c              print*, xS, yS, zS,    "cell      - peidei"
c              print'(3f18.8,i7)', xSl, ySl, zSl, iatm
c
          if (xSl.ge.0.d0-delta) then
          if (ySl.ge.0.d0-delta) then 
          if (zSl.ge.0.d0-delta) then
          if (xSl.lt.1.d0-delta) then 
          if (ySl.lt.1.d0-delta) then 
          if (zSl.lt.1.d0-delta) then
           nsc=nsc+1
             do iatm=1,nat
                  icount=icount+1
                  
                  xS=var%x(iatm)+ix !*var%cell(1,1)
                  yS=var%y(iatm)+iy !*var%cell(1,2)
                  zS=var%z(iatm)+iz !*var%cell(1,3)
          
                  xSl=Q(1,1)*xS+Q(1,2)*yS+Q(1,3)*zS
                  ySl=Q(2,1)*xS+Q(2,2)*yS+Q(2,3)*zS
                  zSl=Q(3,1)*xS+Q(3,2)*yS+Q(3,3)*zS
                  
c                  print'(3f18.8,3i6)', xSl, ySl, zSl,ix,iy,iz
                  
                  xSC(nsc,iatm)=xSl
                  ySC(nsc,iatm)=ySl
                  zSC(nsc,iatm)=zSl
                  
              enddo !loop on the atoms
          endif
          endif
          endif
          endif
          endif
          endif
          
          enddo !loop in X
          enddo !loop in Y
          enddo !loop in Z

c          print*, "icount=", icount
       deallocate(var%x,var%y,var%z)
       deallocate(var%SD_x,var%SD_y,var%SD_z)
       allocate(var%x(icount),var%y(icount),var%z(icount))
       allocate(var%SD_x(icount),var%SD_y(icount),var%SD_z(icount))
       iy=1
       do iatm=1,nat
       do ix=1,nsc
           var%x(iy)=xSC(ix,iatm)
           var%y(iy)=ySC(ix,iatm)
           var%z(iy)=zSC(ix,iatm)
           var%SD_x(:)="T"
           var%SD_y(:)="T"
           var%SD_z(:)="T"
           iy=iy+1
       enddo
       enddo

c     2) modifying the number of atoms
      var%atm_name_n(:)=var%atm_name_n(:)*nsc
      
c  .: 3) enlarge box      
      var%cell(:,:)=cellSC(:,:)

      if (.not.dir) call Direct_to_cartesian(var)
!      
      deallocate(xSC,ySC,zSC)
      return
      end subroutine supercellM

c     =============================
c     subroutine rotateCell -  begin
c     =============================      
      
      subroutine rotateCell(var,v,theta)
      implicit none
      type(contcar)       :: var
      real*8              :: v(3), theta
      real*8              :: rotM(3,3), modv, cell_rot(3,3)
      logical             :: scaled, cart

c  .: IN THIS ROTATING METHOD, THE POSCAR SHOULD BE IN 
c     DIRECT COORDINATED, AND ONLY THE CELL VECTORS ARE
c     ROTATED. 
      
c     Normalizing v
      modv=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      v(1)=v(1)/modv;v(2)=v(2)/modv; v(3)=v(3)/modv 
      rotM(1,1)=0.0
      
c  .: transform do Direct Coordinates
      cart=(.not.var%scaled)
      if (cart) call Cartesian_to_direct(var)

c  .: calculating the rotating matrix 
      rotM(1,1)=cos(theta)+(1-cos(theta))*v(1)*v(1)
      rotM(1,2)=(1-cos(theta))*v(1)*v(2)+sin(theta)*v(3)
      rotM(1,3)=(1-cos(theta))*v(1)*v(3)-sin(theta)*v(2)


      rotM(2,1)=(1-cos(theta))*v(1)*v(2)-sin(theta)*v(3)
      rotM(2,2)=cos(theta)+(1-cos(theta))*v(2)*v(2)
      rotM(2,3)=(1-cos(theta))*v(2)*v(3)+sin(theta)*v(1)

      rotM(3,1)=(1-cos(theta))*v(3)*v(1)+sin(theta)*v(2)
      rotM(3,2)=(1-cos(theta))*v(3)*v(2)-sin(theta)*v(1)
      rotM(3,3)=cos(theta)+(1-cos(theta))*v(3)*v(3)    
          
c  .: rotating the cell          
  
      cell_rot=matmul(var%cell,rotM)
      var%cell(:,:)=cell_rot(:,:)

      if (cart) call direct_to_Cartesian(var)    
      return
      end subroutine rotateCell


c     =============================
c     subroutine rotate -  begin
c     =============================      
      
      subroutine rotate(var,o,v,theta)
      implicit none
      type(contcar)       :: var
      real*8              :: o(3), v(3), theta, t(3)
      real*8              :: rotM(3,3), modv, pos(3), pos2(3)
      logical             :: scaled, cart
      integer             :: natms, i

c  .: IN THIS ROTATING METHOD, THE POSCAR SHOULD BE IN 
c     CARTESIAN COORDINATES, AND ONLY THE ATOMIC POSITIONS
c     ARE ROTATED. 
      
c     Normalizing v
      modv=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      v(1)=v(1)/modv;v(2)=v(2)/modv; v(3)=v(3)/modv 
      rotM(1,1)=0.0
      
c  .: transform do Cartesian Coordinates
      cart=(var%scaled)
      if (cart) call Direct_to_cartesian(var)

c  .: calculating the rotating matrix 
      rotM(1,1)=cos(theta)+(1-cos(theta))*v(1)*v(1)
      rotM(1,2)=(1-cos(theta))*v(1)*v(2)+sin(theta)*v(3)
      rotM(1,3)=(1-cos(theta))*v(1)*v(3)-sin(theta)*v(2)


      rotM(2,1)=(1-cos(theta))*v(1)*v(2)-sin(theta)*v(3)
      rotM(2,2)=cos(theta)+(1-cos(theta))*v(2)*v(2)
      rotM(2,3)=(1-cos(theta))*v(2)*v(3)+sin(theta)*v(1)

      rotM(3,1)=(1-cos(theta))*v(3)*v(1)+sin(theta)*v(2)
      rotM(3,2)=(1-cos(theta))*v(3)*v(2)-sin(theta)*v(1)
      rotM(3,3)=cos(theta)+(1-cos(theta))*v(3)*v(3)    

c  .: translating structure to the origin
      t(:)=-o(:)
      call translate(var,t)
          
c  .: rotating the cell          
  
      natms=sum(var%atm_name_n(:))

c     print*, "natms = ", natms

      do i=1,natms
          pos(1)=var%x(i)
          pos(2)=var%y(i)
          pos(3)=var%z(i)
          pos2=matmul(rotM,pos)
          var%x(i)=pos2(1)
          var%y(i)=pos2(2)
          var%z(i)=pos2(3)
      enddo
c  .: translating structure back to the original position
      t(:)=o(:)
      call translate(var,t)

      if (cart) call Cartesian_to_direct(var)    
      return
      end subroutine rotate

c  .: begin add_atom
      subroutine add_atom(var,c,t)
      implicit none
      type(contcar)       :: var
      character*2         :: c
      real*8              :: t(3)
      logical             :: exists
      integer             :: i, ifind, ntot, idummy

c  .: buffer
      character*2, allocatable  :: B_atm_name(:)
      integer, allocatable      :: B_atm_name_n(:)
      real*8, allocatable       :: B_x(:), B_y(:), B_z(:)
      character*1, allocatable  :: B_SD_x(:), B_SD_y(:), B_SD_z(:)


      exists=.false.
c  .: Checking if the added atom exist
      do i=1,var%nspecs
          if (c==var%atm_name(i)) then
              exists=.true.
              ifind=i
          endif
      end do

c  .: calculating the number of atoms
      ntot=sum(var%atm_name_n(:))

c     write(*,*) "entrou"
c  .: buffering the changing data
      allocate(B_atm_name(var%nspecs),B_atm_name_n(var%nspecs)) 
      allocate(B_x(ntot), B_y(ntot), B_z(ntot))
      if (var%selec_dyn)
     . allocate(B_SD_x(ntot), B_SD_y(ntot), B_SD_z(ntot))
      B_atm_name(:)=var%atm_name(:)
      B_atm_name_n(:)=var%atm_name_n(:)
      B_x(:)=var%x(:)
      B_y(:)=var%y(:)
      B_z(:)=var%z(:)
      if (var%selec_dyn) then
          B_SD_x(:)=var%SD_x(:)
          B_SD_y(:)=var%SD_y(:)
          B_SD_z(:)=var%SD_z(:)
      endif

c  .: resizing some variables
      deallocate(var%x,var%y,var%z)
      if (var%selec_dyn) deallocate(var%SD_x,var%SD_y,var%SD_z)
      allocate(var%x(ntot+1),var%y(ntot+1),var%z(ntot+1))
      if (var%selec_dyn) then 
          allocate(var%SD_x(ntot+1),var%SD_y(ntot+1),var%SD_z(ntot+1))
      endif
      

c  .: adding the new atom
c     (if is a new specie)
      if (.not.exists) then
          var%nspecs=var%nspecs+1
          deallocate(var%atm_name,var%atm_name_n)
          allocate(var%atm_name(var%nspecs),var%atm_name_n(var%nspecs))
          do i=1,var%nspecs-1
              var%atm_name(i)=B_atm_name(i)
              var%atm_name_n(i)=B_atm_name_n(i)
          enddo
          var%atm_name(var%nspecs)=c
          var%atm_name_n(var%nspecs)=1

          do i=1,ntot
              var%x(i)=B_x(i)
              var%y(i)=B_y(i)
              var%z(i)=B_z(i)
              if (var%selec_dyn) then
                  var%SD_x(i)=B_SD_x(i)
                  var%SD_y(i)=B_SD_y(i)
                  var%SD_z(i)=B_SD_z(i)
              endif
          enddo
          var%x(ntot+1)=t(1); var%y(ntot+1)=t(2); var%z(ntot+1)=t(3)
          if (var%selec_dyn) then
c         (below, change to "F" if necessary)
          var%SD_x(ntot+1)="T"
          var%SD_y(ntot+1)="T"
          var%SD_z(ntot+1)="T"
          endif

c     (if is an existing specie)
      else
          var%atm_name_n(ifind)=var%atm_name_n(ifind)+1
          idummy=sum(B_atm_name_n(1:ifind))
          do i=1,idummy
              var%x(i)=B_x(i)
              var%y(i)=B_y(i)
              var%z(i)=B_z(i)
              if (var%selec_dyn) then
                  var%SD_x(i)=B_SD_x(i)
                  var%SD_y(i)=B_SD_y(i)
                  var%SD_z(i)=B_SD_z(i)
              endif
          enddo
          var%x(idummy+1)=t(1);var%y(idummy+1)=t(2);var%z(idummy+1)=t(3)
          if (var%selec_dyn) then
              var%SD_x(idummy+1)="F"
              var%SD_y(idummy+1)="F"
              var%SD_z(idummy+1)="F"
          endif 
          do i=idummy+2,ntot+1
              var%x(i)=B_x(i-1)
              var%y(i)=B_y(i-1)
              var%z(i)=B_z(i-1)
              if (var%selec_dyn) then
                  var%SD_x(i)=B_SD_x(i-1)
                  var%SD_y(i)=B_SD_y(i-1)
                  var%SD_z(i)=B_SD_z(i-1)
              endif
          enddo
      endif


c  .: cleaning the trash
      deallocate(B_atm_name, B_atm_name_n,B_x, B_y, B_z) 
      if (var%selec_dyn) deallocate(B_SD_x, B_SD_y, B_SD_z) 

      end subroutine add_atom
c  .: end add_atom


      subroutine purge_contcar(var)
      implicit none
      type(contcar) :: var
          var%head=""
          var%a0=0.0
          var%cell(:,:)=0.0
          var%nspecs=0
          deallocate(var%atm_name,var%atm_name_n,var%x,var%y,var%z)
          var%scaled=.false.
          deallocate(var%SD_x,var%SD_y,var%SD_z)
c         if (var%selec_dyn) deallocate(var%SD_x,var%SD_y,var%SD_z)
          var%selec_dyn=.false.
          var%ngx=0; var%ngy=0; var%ngz=0
          if (var%grid) deallocate(var%v)
          var%grid=.false.
      end subroutine purge_contcar      

c     =============================
c     subroutine write_contcarF -  begin
c     =============================      
      
      subroutine write_contcarF(var,fname)
      implicit none
      type(contcar) :: var
      character*64  :: fname
      integer       :: ntot, icount, i, j

       open(3001,file=fname)
       write(3001,'(a64)') var%head
       write(3001,16) var%a0
       write(3001,16) var%cell(1,:)
       write(3001,16) var%cell(2,:)
       write(3001,16) var%cell(3,:)
       write(3001,'(100a5)') var%atm_name(:) 
       write(3001,'(100i5)') var%atm_name_n(:) 
       if (var%selec_dyn)   write(3001,'(a18)') "Selective Dynamics"
       if (.not.var%scaled) write(3001,'(a9)') "Cartesian"
       if (var%scaled)      write(3001,'(a6)') "Direct"
       ntot=sum(var%atm_name_n(:))
       icount=1
       if (var%selec_dyn) then
        do i=1,var%nspecs
           do j=1,var%atm_name_n(i)
               write(3001,17) var%x(icount),var%y(icount),var%z(icount),
     &         var%SD_x(icount),var%SD_z(icount),var%SD_z(icount)
               icount=icount+1
           enddo      
        enddo
        else
        do i=1,var%nspecs
           do j=1,var%atm_name_n(i)
               write(3001,17) var%x(icount),var%y(icount),var%z(icount) 
               icount=icount+1
           enddo      
        enddo
        
        endif

      close(3001)

15    format(3f21.12,3a5)
16    format(3f21.12)
17    format(3f21.12,3x,a1,3x,a1,3x,a1)     
      end subroutine write_contcarF
      
c     =============================
c     subroutine CONTCAR_to_xyzF  -  begin
c     =============================      
      
      
      subroutine CONTCAR_to_xyzF(var,fname)
      implicit none
      type(contcar) var
      character*64 fname
C     .: internal      
      integer i, j, icount,itot
      
      open(3001,file=fname)

      itot=sum(var%atm_name_n(:))
      write(3001,*) itot
      write(3001,*) " generated from a POSCAR file "
      icount=1
      do i=1,var%nspecs
         do j=1,var%atm_name_n(i)
             if (var%scaled) then
             write(*,17) var%atm_name(i),
     &       var%x(icount)*var%cell(1,1)+
     &       var%y(icount)*var%cell(2,1)+
     &       var%z(icount)*var%cell(3,1),
     &       var%x(icount)*var%cell(1,2)+
     &       var%y(icount)*var%cell(2,2)+
     &       var%z(icount)*var%cell(3,2),
     &       var%x(icount)*var%cell(1,3)+
     &       var%y(icount)*var%cell(2,3)+
     &       var%z(icount)*var%cell(3,3)
             else
             write(3001,17) var%atm_name(i), var%x(icount), 
     &                 var%y(icount),var%z(icount) 
            endif
         icount=icount+1
         enddo      
      enddo
      
      close(3001)
17    format(a2,5x,3f21.12)      
      
      end subroutine CONTCAR_to_xyzF

      

      subroutine translate_list(var,t,list,n)
      implicit none
      type(contcar)       :: var
      real*8              :: t(3), cellT(3,3)
      integer             :: n, list(n)

c     .: internal variables
      real*8              :: td(3)
      integer             :: i, iatm

      do i=1,n
      iatm=list(i)
      if (.not.var%scaled) then
          var%x(iatm)=var%x(iatm)+t(1)
          var%y(iatm)=var%y(iatm)+t(2)
          var%z(iatm)=var%z(iatm)+t(3)
      else
          cellT(1,1)=var%cell(1,1)
          cellT(1,2)=var%cell(2,1)
          cellT(1,3)=var%cell(3,1)

          cellT(2,1)=var%cell(1,2)
          cellT(2,2)=var%cell(2,2)
          cellT(2,3)=var%cell(3,2)

          cellT(3,1)=var%cell(1,3)
          cellT(3,2)=var%cell(2,3)
          cellT(3,3)=var%cell(3,3)

          call LU(cellT,t,td,3)

          var%x(iatm)=var%x(iatm)+td(1)
          var%y(iatm)=var%y(iatm)+td(2)
          var%z(iatm)=var%z(iatm)+td(3)

      endif

      enddo
      return
      end subroutine translate_list

      end module m_contcar
