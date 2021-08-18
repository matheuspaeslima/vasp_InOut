      program main
      use m_contcar
      implicit none
      character*64           :: fname
      real*8                 :: t(3)
      integer                :: i, natoms
      integer, allocatable   :: list(:)
      type(contcar)          :: data1

c  .: read geometry
      fname="POSCAR"
      call read_contcar(data1,fname)

      print*, "how many atoms do you like to translate?"
      read(*,*) natoms
      allocate (list(natoms))
      print*, "which atoms do you like to translate?"
      read(*,*) list(:)
      print*, "please, provide the translation vector"
      read(*,*) t(:)

      print*, "list(:)=", list(:)
      
      call translate_list(data1,t,list,natoms)
      
      call write_contcar(data1)

      end program main
       
