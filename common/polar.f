C**
C***********************************************************************
C**
      subroutine polar ( omega0 , afocus , lambda , nt , nr , nf , dt ,
     .                   dr , df , dotdot )
C**
C***********************************************************************
C**
      implicit none
c
      integer nt , nr , nf , k , j , i
c
      real    dotdot(nf,nr,nt) , omega0 , afocus , lambda , t , r , f ,
     .        dt , dr , df
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Compute analytical representation of acceleration given two
c   parameters - a rotation rate (omega0) and a focus value (afocus).
c
c   Modified June 30, 2000 to center about the time, range, and Doppler
c   coordinates
c
      do k = 1 , nt
c
         t = ( float( k - nt / 2 - 1 ) + 0.5 ) * dt
c
         do j = 1 , nr
c
            r = ( float( j - nr / 2 - 1 ) + 0.5 ) * dr
c
            do i = 1 , nf
c
               f             = ( float( i - nf / 2 - 1 ) + 0.5 ) * df
c
               dotdot(i,j,k) = afocus -
     .                         ( ( 2.0 * omega0 ** 2 ) / lambda ) * r
     .                       - ( omega0 ** 2 ) * t * f
     .                       + ( ( omega0 ** 4 ) / lambda ) * r * t * t
c
            enddo
c
         enddo
c
      enddo
c
      return
      end
