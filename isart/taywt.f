c
                           Program TAYWT
c
      implicit none
c
      real        taywt_rt
c
      integer     nrange , i
c
      parameter ( taywt_rt = - 30.0 , nrange = 684 )
c
      real        wts(nrange)
c
      call taylor ( nrange , taywt_rt , wts )
c
      write ( 6 , '(4(i5,f8.3))' ) ( i , wts(i) , i = 1 , nrange )
c
      stop
      end
C**
C***********************************************************************
C**
      subroutine taylor ( nx , sll_db , weights )
C**
C***********************************************************************
C**
c
c   This subroutine calculates the Taylor illumination function as
c   defined on pg. 249 of Barton, HANDBOOK OF RADAR MEASUREMENT.
c
c   Input parameters:
c
c        nx     = # of data points
c        sll_db = sidelobe level in dB (This should negative)
c
c   Output array:
c
c        weights = Taylor illumination function . This is normalized 
c                  to be one at the center.         
c
      implicit none
c
      integer     nmax , nx , nbar , m , n , ix
c
      parameter ( nmax = 30 )
c
      real        weights(nx) , f(nmax) , pi , sll_db , x , a , dwt ,
     .            sigma2 , fact , df , denom , xnum
c
      data pi / 3.14159265 /
c
      x    = 10.0 ** ( abs( sll_db ) / 20.0 )
      x    = x + sqrt( x ** 2 - 1.0 )
      a    = alog( x ) / pi
      nbar = 2 * a ** 2 + 1.5
c
      if ( nbar .gt. ( nmax + 1 ) ) then
         stop ' nbar is too large '
      endif
c
      sigma2 = float( nbar ) ** 2 /
     .         ( a ** 2 + ( float( nbar ) - 0.5 ) ** 2 )
      fact = 1.0
c
      do m = 1 , nbar - 1
c
         dwt = 1.0
c
         do n = 1 , nbar - 1
c
            df = ( 1.0 - float( m ) ** 2 / sigma2 /
     .           ( a ** 2 + ( float( n ) - 0.5 ) ** 2 ) )
c
            if ( n .ne. m ) df = df / ( 1.0 - float( m ) ** 2 /
     .                                        float( n ) ** 2 )
c
            dwt = dwt * df
c
         enddo
c
         f(m) = fact * dwt / 2.0
         fact = - fact
c
      enddo
c
      do ix = 1 , nx
c
         a     = float( ix - 1 ) / float( nx - 1 ) - 0.5
         xnum  = 0.0
         denom = 0.0
c
         do m = 1 , nbar - 1
c
            xnum  = xnum + f(m) * cos( 2.0 * pi * float( m ) * a )
            denom = denom + f(m)
c
         enddo
c
         weights(ix) = ( 1.0 + 2.0 * xnum ) / ( 1.0 + 2.0 * denom )
c
      enddo
c
      return
      end
