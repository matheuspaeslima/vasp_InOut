      program main
      use m_contcar
      implicit none
      character*64           :: fname
      real*8, parameter    :: pi=3.1415926535807932384d0
      real*8                 :: o(3), v(3), theta
      integer                :: i, natoms
      integer, allocatable   :: list(:)
      type(contcar)          :: data1

c  .: read geometry
      fname="POSCAR"
      call read_contcar(data1,fname)

      print*, "how many atoms do you like to rotate?"
      read(*,*) natoms
      allocate (list(natoms))
      print*, "which atoms do you like to rotate?"
      read(*,*) list(:)
      print*, "please, provide the origin for the rotation"
      read(*,*) o(:)
      print*, "please, provide the direction for the rotation"
      read(*,*) v(:)
      print*, "please, provide the angle (degree) for the rotation"
      read(*,*) theta
c  .: transforming de angle from degree to radians
      theta=theta*pi/180.d0

      print*, "list(:)=", list(:)
      
      call rotate_list(data1,o,v,theta,list,natoms)
      
      call write_contcar(data1)

      end program main
       
