      program main
      use m_contcar
      implicit none
      character*64 fname
      real*8 zmax, zmin, t(3)
      real*8, allocatable :: vpot(:)
      integer             :: i, natoms
      type(contcar) data1

c  .: read geometry
      fname="log2.vasp"
      call read_contcar(data1,fname)
c     call direct_to_cartesian(data1)
c     call put_selective(data1)
      t=0.d0
      t(1)=-12.304000000002+1.787181988918
      t(2)=-7.103717712110
      t(3)=-9.515479999998+7.598811759562+1.9
      call translate(data1,t)
      
      call write_contcar(data1)

c  .: EXAMPE of random numbers
c     open(21,file="seed.dat")
c     read(21,*) iseed1, iseed2
c     write(*,*) "the seeds are :", iseed1, iseed2
C  .: initializing random number generator
c     call fillu(iseed1,iseed2)
c  .: generate random number in [0,1]
c     rng=uni64()


      end program main
       
