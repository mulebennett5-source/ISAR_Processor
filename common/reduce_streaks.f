C**
C***********************************************************************
C**
      subroutine reduce_streaks ( ci , nx , ny , dconst , rconst ,
     .                            dbmin )
C**
C***********************************************************************
C**
      implicit none
c
      integer nx , ny , i , j
c
      complex ci(nx,ny)
c
      real    dconst , rconst , dbmin , abslim , total , dfacx(ny) ,
     .        dfacy(nx)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     Limit calculation to dbmin of the top intensity
c
      abslim = 0.0
c
      do j = 1 , ny     ! Loop over range cells
c
         do i = 1 , nx
c
            abslim = amax1( abslim , cabs( ci(i,j) ) )
c
         enddo
c
      enddo
c
      abslim = 10**(-abs(dbmin/20)) * abslim
c
c-----------------------------------------------------------------------
c
c     Reduce cross-range streaks
c
      do j = 1 , ny     ! Loop over range cells
c
         total = 0.0
c
         do i = 1 , nx
c
            total = total + log( amax1( abslim, cabs( ci(i,j) ) ) )
c
         enddo
c
         dfacx(j) = exp( - dconst * total / float(nx) )
c
      enddo             ! End of range loop
c
c-----------------------------------------------------------------------
c
c     Reduce range streaks
c
      do i = 1 , nx     ! Loop over range cells
c
         total = 0.0
c
         do j = 1 , ny
c
            total = total + log( amax1( abslim, cabs( ci(i,j) ) ) )
c
         enddo
c
         dfacy(i) = exp( - rconst * total / float(ny) )
c
      enddo             ! End of range loop
c
c-----------------------------------------------------------------------
c
c     Apply corrections
c
      do j = 1 , ny
c      
         do i = 1 , nx
c         
            ci(i,j) = ci(i,j) * dfacy(i) * dfacx(j)
c         
         enddo
c      
      enddo
c      
      return
      end