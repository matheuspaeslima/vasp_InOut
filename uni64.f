!     ===============================================================
!     first call fillu(iseed1,iseed2), whith iseed1 and iseed2 integer
!     after x=uni64()
!     A semente esta salva no common block
!     referencia: 
!     ===============================================================

      double precision function uni64 ()
      integer i,j
      double precision u(97),r,d,c,x
      parameter (r=9007199254740881.d0/9007199254740992.d0)
      parameter (d=362436069876.d0/9007199254740992.d0)
      data c,i,j/0.,97,33/
      common /unicom/ u
      save
      x=u(i)-u(j)
      if (x.lt.0.0) x=x+1.0d0
      u(i)=x
      i=i-1
      if (i.eq.0) i=97
      j=j-1
      if (j.eq.0) j=97
      c=c-d
      if (c.lt.0.0) c=c+r
      x=x-c
      uni64=x
      if (x.lt.0.) uni64=x+1.d0
      return
      end

      subroutine fillu (seed1,seed2)
      integer seed1,seed2,i,j,x,y
      double precision s,t,u(97)
      common /unicom/ u
      save
      x=seed1
      y=seed2
      do 20 i=1,97
        s=0.d0
        t=0.5d0
        do 10 j=1,53
          x=mod(6969*x,65543)
          y=mod(8888*y,65579)
          if (iand(ieor(x,y),32).gt.0) s=s+t
          t=0.5d0*t
10      continue
        u(i)=s
20    continue
      return
      end
