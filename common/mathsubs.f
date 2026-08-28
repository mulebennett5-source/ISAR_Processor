c**
c***********************************************************************
c**
      subroutine caputi ( si_npsi , wta )
c**
c***********************************************************************
c**
c
c   Bill Caputi's sub-image weight function
c
      implicit none
c
      integer si_npsi , ip , si_naskip
c
      real    wta(si_npsi) , pi
c
c-----------------------------------------------------------------------
c
      pi        = atan2( 0.0 , - 1.0 )
c
      si_naskip = si_npsi / 2
c
      do ip = 1 , si_npsi
c
         wta(ip) = 0.75 * sin( pi * float( ip ) / float( 1 + si_npsi ) )
c
         if ( ip .le. si_naskip ) then
c
            wta(ip) = wta(ip) + float( ip ) / float( si_naskip )
c
         else
c
            wta(ip) = wta(ip) + 2.0 - float( ip - 1 ) /
     .                                float( si_naskip)
c
         endif
c
         wta(ip) = wta(ip) / 1.75
c
      enddo
c
      return
      end
c**
c***********************************************************************
c**
      real function atan2p ( ci , cr )
c**
c***********************************************************************
c**
c
c   ATAN2 Function with trap for both components zero
c
c-----------------------------------------------------------------------
c
      implicit none
c
      real ci , cr
c
      if ( ci .eq. 0.0 .and. cr .eq. 0.0 ) then
c
         atan2p = 0.0
c
      else
c
         atan2p = atan2( ci , cr )
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      real function phasec ( c )
c**
c***********************************************************************
c**
c
c   Phase of a complex number in degrees
c
c-----------------------------------------------------------------------
c
      implicit none
c
      complex c
c
      real    raddeg , atan2p
c
      raddeg = 180 / atan2( 0.0 , - 1.0 )
c
      phasec = raddeg * atan2p( aimag( c ) , real( c ) )
c
      return
      end
c**
c***********************************************************************
c**
      subroutine fade ( x_filt , xnew , i , n )
c**
c***********************************************************************
c**
c
c   Implement a simple fading memory filter for a real number
c
c     x_filt  :  Filtered version (input and output)
c
c     xnew    :  Latest version of the un-filtered time series
c
c     i       :  Sequence number of the latest data point
c
c     n       :  Length of the fading memory filter
c
c-----------------------------------------------------------------------
c
      implicit none
c
      real    x_filt , xnew
c
      integer i , n
c
c-----------------------------------------------------------------------
c
      if ( i .le. 0 ) then
c
         x_filt = xnew
c
      else if ( i .lt. n ) then
c
         x_filt = ( x_filt * float( i - 1 ) + xnew ) / float( i )
c
      else
c
         x_filt = ( x_filt * float( n - 1 ) + xnew ) / float( n )
c
      endif
c
      return
      end

C**
C***********************************************************************
C**
      subroutine wtfft ( d , wt , fft , work , isign , ishift ,
     .                   nd , nfft )
C**
C***********************************************************************
C**
      implicit none
c
      integer isign , nn(1) , ishift , jshift , nd , nfft , i , j , jsh
c
      real    d(2,nd) , wt(nd) , fft(2,nfft) , work(2,nfft)
c
c-----------------------------------------------------------------------
c
      do i = 1 , nfft
c
         work(1,i) = 0.0
         work(2,i) = 0.0
c
      enddo
c
      if ( ishift .ge. 0 ) then
c
         jshift = ishift
c
      else
c
         jshift = ishift + nfft * ( 1 + ( iabs( ishift ) / nfft ) )
c
      endif
c
      do j = 1 , nd
c
         jsh = 1 + mod( j + jshift - 1 , nfft )
c
         work(1,jsh) = work(1,jsh) + wt(j) * d(1,j)
         work(2,jsh) = work(2,jsh) + wt(j) * d(2,j)
c
      enddo
c
      nn(1) = nfft
c
      call fourt ( work , nn , 1 , isign , 1 , fft , 2 * nfft )
c
c   Rotate DC to the center
c
      do j = 1 , nfft
c
         jsh = j - nfft / 2
c
         if ( jsh .le. 0 ) jsh = jsh + nfft
c
         fft(1,jsh) = work(1,j)
         fft(2,jsh) = work(2,j)
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      SUBROUTINE STATS ( X , N , XMIN , XMAX , AVGX , SDX )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      INTEGER          N , I
C
      REAL             X(N) , XMIN , XMAX , AVGX , SDX
C
      DOUBLE PRECISION SUMX , SUMXSQ
C
C-----------------------------------------------------------------------
C
      IF ( N .LE. 0 ) THEN
C
         AVGX = 0
         XMAX = 0
         XMIN = 0
         SDX  = 0
C
      ELSE
C
         SUMX   = 0.0
         SUMXSQ = 0.0
C
         XMIN = X(1)
         XMAX = X(1)
C
         DO I = 1 , N
C
            SUMX   = SUMX   + X(I)
            SUMXSQ = SUMXSQ + X(I) ** 2
C
            XMIN = AMIN1( XMIN , X(I) )
            XMAX = AMAX1( XMAX , X(I) )
C
         ENDDO
C
         AVGX   = SUMX   / FLOAT( N )
         SUMXSQ = SUMXSQ / FLOAT( N )
C
         SDX = SUMXSQ - ( SUMX / FLOAT( N ) ) ** 2
C
         IF ( SDX .GT. 0.0 ) THEN
C
            SDX = SQRT( SDX )
C
         ELSE
C
            SDX = 0.0
C
         ENDIF
C
      ENDIF
C
      RETURN
      END
c**
c***********************************************************************
c**
      subroutine getasq ( asq , c , ni , nj , ifirst , ilast , jfirst ,
     .                    jlast )
c**
c***********************************************************************
c**
c   Get the average squared value from a section of a complex array
c
      implicit none
c
      real    asq
c
      integer ni , nj , ifirst , ilast , jfirst , jlast , i , j
c
      complex c(ni,nj)
c
      asq = 0.0
c
      if ( ( ilast .ge. ifirst ) .and. ( jlast .ge. jfirst ) ) then
c
         do j = jfirst , jlast
c
            do i = ifirst , ilast
c
               asq = asq + cabs( c(i,j) ) ** 2
c
            enddo
c
         enddo
c
         asq = asq / float( ( ilast + 1 - ifirst ) *
     .                      ( jlast + 1 - jfirst ) )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      complex function cov ( c1 , c2 , n )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i
c
      complex c1(n) , c2(n)
c
c-----------------------------------------------------------------------
c
      cov = cmplx( 0.0 , 0.0 )
c
      do i = 1 , n
c
         cov = cov + conjg( c1(i) ) * c2(i)
c
      enddo
c
      cov = cov / float( n )
c
      return
      end
C**
C***********************************************************************
C**
      complex function cov2 ( c1 , c2 , mdim , ndim , m , n , ilag ,
     .                        jlag , window )
C**
C***********************************************************************
C**
      implicit none
c
      integer mdim , ndim , m , n , ilag , jlag , i , j , jtot , mhalf ,
     .        nhalf , window , ilim1 , ilim2 , jlim1 , jlim2
c
      complex c1(mdim,ndim) , c2(mdim,ndim)
c
c-----------------------------------------------------------------------
c
      cov2 = cmplx( 0.0 , 0.0 )
c
      if ( mdim .ge. m .and. ndim .ge. n ) then
c
         mhalf = ( mdim - m ) / 2
         nhalf = ( ndim - n ) / 2
c
         ilim1 = max(    1 , 1 + mhalf - ilag / 2 )
         ilim2 = min( mdim , mdim - mhalf - ilag / 2 )
c
         jlim1 = max(    1 , 1 + nhalf - jlag / 2 )
         jlim2 = min( ndim , ndim - nhalf - jlag / 2 )
c
         jtot  = 0
c
         do j = jlim1 , jlim2
c
            if ( j .le. ndim / 2 - window / 2 .or.
     .           j .gt. ndim / 2 + window / 2 ) then
c
               jtot = jtot + 1
c
               do i = ilim1 , ilim2
c
                  cov2  = cov2 + conjg( c1(i+ilag,j+jlag) ) * c2(i,j)
c
               enddo
c
            endif
c
         enddo
c
         if ( jtot * ( ilim2 - ilim1 + 1 ) .gt. 0 ) then
c
            cov2 = cov2 / float( jtot * ( ilim2 - ilim1 + 1 ) )
c
         else
c
            write ( 7 , * ) ' No points in covariance window in COV2'
            stop            ' No points in covariance window in COV2'
c
         endif
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine fine_shift ( c1 , c2 , n )
C**
C***********************************************************************
C**
c   Purpose:  To do a sub-pixel shift of one vector to align it with
c             another.  The shift is done in the Fourier domain.
c
      implicit none
c
      integer n , i , nn(1)
c
      complex c1(n) , c2(n) , covp , cov
c
      real    dphase , phase , atan2p , work(2*n)
c
c-----------------------------------------------------------------------
c
      nn(1) = n
c
      call fourt1c ( c1 , nn , 1 , - 1 , 0 , work , 2 * n )
c
      call fourt1c ( c2 , nn , 1 , - 1 , 0 , work , 2 * n )
c
      covp   = cov( c1 , c2(2) , n - 1 )
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
      do i = 1 , n
c
         phase = 0.5 * float( i - 1 - n / 2 ) * dphase
c
         c1(i) = c1(i) * cmplx( cos( phase ) , - sin( phase ) )
c
         c2(i) = c2(i) * cmplx( cos( phase ) , + sin( phase ) )
c
      enddo
c
      covp   = cov( c1 , c2(2) , n - 1 )
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
      call fourt1c ( c1 , nn , 1 , + 1 , 0 , work , 2 * n )
c
      call fourt1c ( c2 , nn , 1 , + 1 , 0 , work , 2 * n )
c
      c1(:) = c1(:) / float( n )
c
      c2(:) = c2(:) / float( n )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine c1smooth ( c , n , niters )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i , iter , niters
c
      complex c(n)
c
c-----------------------------------------------------------------------
c
c   Smooth complex time series niters times with a three point smoother
c
      if ( niters .gt. 0 ) then
c
         do iter = 1 , niters
c
c   Even points except n
c
            do i = 2 , n - 2 , 2
c
               c(i) = 0.25 * ( c(i+1) + c(i-1) + 2.0 * c(i) )
c
            enddo
c
            c(n) = 0.5 * ( c(n-1) + c(n-2) )
c
c   Odd points except 1
c
            do i = 3 , n - 1 , 2
c
               c(i) = 0.25 * ( c(i+1) + c(i-1) + 2.0 * c(i) )
c
            enddo
c
            c(1) = 0.5 * ( c(2)   + c(3)   )
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine r1smooth ( r , n , niters )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i , iter , niters
c
      real    r(n)
c
c-----------------------------------------------------------------------
c
c   Smooth complex time series niters times with a three point smoother
c
      if ( niters .gt. 0 ) then
c
         do iter = 1 , niters
c
c   Even points except n
c
            do i = 2 , n - 2 , 2
c
               r(i) = 0.25 * ( r(i+1) + r(i-1) + 2.0 * r(i) )
c
            enddo
c
c   Odd points except 1
c
            do i = 3 , n - 1 , 2
c
               r(i) = 0.25 * ( r(i+1) + r(i-1) + 2.0 * r(i) )
c
            enddo
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine c1smooth_pg ( c , n , nsmth , dphase )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i , nsmth
c
      complex c(n) , covp , cov
c
      real    dphase , phase , atan2p
c
c-----------------------------------------------------------------------
c
c   Smooth
c
      if ( nsmth .gt. 0 ) call c1smooth ( c , n , 1 )
c
c   Remove mean phase gradient
c
      covp   = cov( c , c(2) , n - 1 )
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
c     write ( 7 , * ) ' dphase = ' , dphase
c
      do i = 1 , n
c
         phase = float( i - 1 - n / 2 ) * dphase
c
         c(i)  = c(i) * cmplx( cos( phase ) , sin( phase ) )
c
      enddo
c
c   Smooth
c
      if ( nsmth .gt. 1 ) call c1smooth ( c , n , nsmth - 1 )
c
c   Replace phase gradient
c
      do i = 1 , n
c
         phase = float( i - 1 - n / 2 ) * dphase
c
         c(i)  = c(i) * cmplx( cos( phase ) , - sin( phase ) )
c
      enddo
c
      covp   = cov( c , c(2) , n - 1 )
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
c     write ( 7 , * ) ' dphase = ' , dphase
c
      return
      end
c**
c***********************************************************************
c**
      subroutine canwgts ( cwt1 , cwt2 , c11 , c22 , c33 , c21 , c23 ,
     .                     c13 , cerror )
c**
c***********************************************************************
c**
c   Purpose:  To calculate the complex weights for a 2-component
c             cancellation algorithm.  Channel 2 is to be cancelled by
c             a linear combination of channels 1 and 3.
c
      implicit none
c
      complex c13 , c23 , c21 , cwt1 , cwt2
c
      real    c11 , c22 , c33 , det , cerror
c
c-----------------------------------------------------------------------
c
      det     = c11 * c33 - cabs( c13 ) ** 2
c
      if ( det .le. 0.0 ) then  !  Should never happen!
c
         write ( 6 , * ) 'Negative determinant:' , c11 , c33 , c13 , det
c
         read ( 5 , * )
c
      endif
c
      cwt1 = ( conjg( c23 ) * c13 - conjg( c21 ) * c33 ) / det
      cwt2 = ( conjg( c21 * c13 ) - conjg( c23 ) * c11 ) / det
c
      cerror  = c22 + c11 * cabs( cwt1 ) ** 2 +
     .              + c33 * cabs( cwt2 ) ** 2 +
     .              + 2.0 * real( c21 * cwt1 )
     .              + 2.0 * real( c23 * cwt2 )
     .              + 2.0 * real( conjg( cwt1 ) * cwt2 * c13 )
c
      if ( cerror .le. 0.0 ) then  !  Should never happen!
c
         write ( 6 , * ) 'Negative error:' , c22 , c11 , c33 , c13 ,
     .                                       cerror
c
         read ( 5 , * )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      SUBROUTINE SMOOTH ( A , W , NX , NY , NXSMTH , NYSMTH , NLOOPS )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      INTEGER NX , NY , I , J , IXSMTH , NXSMTH , IYSMTH , NYSMTH ,
     .        ILOOPS , NLOOPS
C
      REAL    A(NX,NY) , W(NX,NY)
C
C-----------------------------------------------------------------------
C
      DO ILOOPS = 1 , NLOOPS
C
         IF ( NXSMTH .GT. 0 ) THEN
C
            DO IXSMTH = 1 , NXSMTH
C
               DO J = 1 , NY
C
                  DO I = 2 , NX - 1
C
                     W(I,J) = 0.25 *
     .                        ( A(I-1,J) + A(I+1,J) + 2.0 * A(I,J) )
C
                  ENDDO
C
                  W(1,J)  = A(1,J)
                  W(NX,J) = A(NX,J)
C
               ENDDO
C
               CALL RSWAP ( A , W , NX * NY )
C
            ENDDO
C
         ENDIF
C
         IF ( NYSMTH .GT. 0 ) THEN
C
            DO IYSMTH = 1 , NYSMTH
C
               DO I = 1 , NX
C
                  DO J = 2 , NY - 1
C
                     W(I,J) = 0.25 *
     .                        ( A(I,J-1) + A(I,J+1) + 2.0 * A(I,J) )
C
                  ENDDO
C
                  W(I,1)  = A(I,1)
                  W(I,NY) = A(I,NY)
C
               ENDDO
C
               CALL RSWAP ( A , W , NX * NY )
C
            ENDDO
C
         ENDIF
C
      ENDDO
C
      RETURN
      END
C**
C***********************************************************************
C**
      SUBROUTINE RSWAP ( A , B , N )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      INTEGER N , I
C
      REAL    A(N) , B(N) , TEMP
C
C-----------------------------------------------------------------------
C
      DO I = 1 , N
C
         TEMP = A(I)
         A(I) = B(I)
         B(I) = TEMP
C
      ENDDO
C
      RETURN
      END
C**
C***********************************************************************
C**
      REAL FUNCTION DB( X )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      REAL X
C
      DB = 10.0 * ALOG10( X )
C
      RETURN
      END
C**
C***********************************************************************
C**
      subroutine bswap ( n , b )
C**
C***********************************************************************
C**
c
c   Swaps bytes to convert INTEGER*2 from Sun to IBM-PC format
c
      implicit none
c
      integer     n , i
c
      character   b(2,n)*1 , btmp*1
c
c-----------------------------------------------------------------------
c
      do i = 1 , n
c
         btmp   = b(1,i)
         b(1,i) = b(2,i)
         b(2,i) = btmp
c
      enddo
c
      return
      end
