C**
C***********************************************************************
C**
      subroutine reduce_sidelobes ( image , noise , nx , ny , taywta )
C**
C***********************************************************************
C**
      implicit none
c
      integer nx , ny , i , j
c
      real    image(nx,ny) , noise , taywta , sl , imax
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      sl = 10.0**(abs(taywta/10.0))
c
      do j = 1 , ny
c
         imax = 0.0
         do i = 1 , nx
            imax = max(imax,image(i,j))
         enddo
c
         do i = 1 , nx
            if ( ( image(i,j) .lt. imax/sl ) .and.
     .           ( image(i,j) .gt. noise ) ) image(i,j) = noise
         enddo
c
      enddo
c
      return
      end