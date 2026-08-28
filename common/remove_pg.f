C**
C***********************************************************************
C**
      subroutine remove_phase_gradient ( c , n )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i
c
      complex c(n) , covp , cov
c
      real    dphase , phase , atan2p
c
c-----------------------------------------------------------------------
c
c   Remove the mean phase gradient from a complex vector
c
c   Compute single-lag covariance value
c
      covp   = cov( c , c(2) , n - 1 )
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
c   Remove phase gradient
c
      do i = 1 , n
c
         phase = float( i - 1 - n / 2 ) * dphase
c
         c(i)  = c(i) * cmplx( cos( phase ) , sin( phase ) )
c
      enddo
c
c   Re-compute single-lag covariance value for output
c
      covp   = cov( c , c(2) , n - 1 )
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine remove_phase_gradient_2 ( c , n , m )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , m , i , j
c
      complex c(n,m) , r(m) , covp , cov , czero
c
      real    dphase , phase , atan2p
c
c-----------------------------------------------------------------------
c
c   Remove the mean phase gradient from a 2-D complex array
c
      czero = cmplx( 0.0 , 0.0 )
c
c   Remove phase gradient in range
c
      covp  = czero
c
      do i = 1 , n
c
         r(:) = c(i,:)
c
         covp = covp + cov( r , r(2) , m - 1 )
c
      enddo
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
      do j = 1 , m
c
         phase = float( j - 1 - m / 2 ) * dphase
c
         c(:,j) = c(:,j) * cmplx( cos( phase ) , sin( phase ) )
c
      enddo
c
c   Remove phase gradient in cross-range
c
      covp  = czero
c
      do j = 1 , m
c
         covp = covp + cov( c(1,j) , c(2,j) , n - 1 )
c
      enddo
c
      dphase = - atan2p( aimag( covp ) , real( covp ) )
c
      do i = 1 , n
c
         phase = float( i - 1 - n / 2 ) * dphase
c
         c(i,:) = c(i,:) * cmplx( cos( phase ) , sin( phase ) )
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine remove_phase_gradient_2_2 ( c1 , c2 , nx , ny , phlag )
c**
c***********************************************************************
c**
      implicit none
c
      integer     nx , ny , i , j , phlag , phlag_use , ii , jj
c
      complex     c1(nx,ny) , c2(nx,ny) , covp , czero
c
      real        xphdelt , yphdelt , p , raddeg , power(2,2,2) ,
     .            power1avg , power2avg , atan2p
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Estimate and remove a phase gradient between two images
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      raddeg = 180.0 / atan2( 0.0 , - 1.0 )
c
      czero  = cmplx( 0.0 , 0.0 )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   First, calculate the phase gradient using a lag of phlag
c
      covp      = czero
c
      phlag_use = min( phlag , nx / 4 )
c
      do j = 1 , ny
c
         do i = 1 , nx - phlag_use
c
            covp = covp + c1(i,j) * conjg( c1(i+phlag_use,j) )
     .                  * c2(i+phlag_use,j) * conjg( c2(i,j) )
c
         enddo
c
      enddo
c
      xphdelt = atan2p( aimag( covp ) , real( covp ) ) /
     .                  float( phlag_use )
c
      write ( 7 , * ) ' X-Phase diff. per pixel - before (deg): ' ,
     .                  xphdelt * raddeg
c
c-----------------------------------------------------------------------
c
c   Remove the estimated phase gradient from the second image
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            p       = - xphdelt * float( i - nx / 2 )
c
            c2(i,j) = c2(i,j) * cmplx( cos( p ) , sin( p ) )
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Verify by re-calculating the gradient
c
      covp = czero
c
      do j = 1 , ny
c
         do i = 1 , nx - phlag_use
c
            covp = covp + c1(i,j) * conjg( c1(i+phlag_use,j) ) *
     .                    c2(i+phlag_use,j) * conjg( c2(i,j) )
c
         enddo
c
      enddo
c
      xphdelt = atan2p( aimag( covp ) , real( covp ) ) /
     .                  float( phlag_use )
c
      write ( 7 , * ) ' X-Phase diff. per pixel - after (deg):  ' ,
     .                  xphdelt * raddeg
c
c-----------------------------------------------------------------------
c
c   First, calculate the phase gradient using a lag of phlag
c
      covp      = czero
c
      phlag_use = min( phlag , ny / 4 )
c
      do j = 1 , ny - phlag_use
c
         do i = 1 , nx
c
            covp = covp + c1(i,j) * conjg( c1(i,j+phlag_use) )
     .                  * c2(i,j+phlag_use) * conjg( c2(i,j) )
c
         enddo
c
      enddo
c
      yphdelt = atan2p( aimag( covp ) , real( covp ) ) /
     .                  float( phlag_use )
c
      write ( 7 , * ) ' Y-Phase diff. per pixel - before (deg): ' ,
     .                  yphdelt * raddeg
c
c-----------------------------------------------------------------------
c
c   Remove the estimated phase gradient from the second image
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            p       = - yphdelt * float( j - ny / 2 )
c
            c2(i,j) = c2(i,j) * cmplx( cos( p ) , sin( p ) )
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Verify by re-calculating the gradient
c
      covp = czero
c
      do j = 1 , ny - phlag_use
c
         do i = 1 , nx
c
            covp = covp + c1(i,j) * conjg( c1(i,j+phlag_use) )
     .                  * c2(i,j+phlag_use) * conjg( c2(i,j) )
c
         enddo
c
      enddo
c
      yphdelt = atan2p( aimag( covp ) , real( covp ) ) /
     .                  float( phlag_use )
c
      write ( 7 , * ) ' Y-Phase diff. per pixel - after (deg):  ' ,
     .                  yphdelt * raddeg
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Diagnostic calculation - calculate the power in each image for each
c   of the four quadrants
c
      power(:,:,:) = 0.0
c
      do j = 1 , ny / 2
c
         jj = j + ny / 2
c
         do i = 1 , nx / 2
c
            ii = i + nx / 2
c
            power(1,1,1) = power(1,1,1) + c1(i,j)   * conjg( c1(i,j) )
c
            power(1,1,2) = power(1,1,2) + c2(i,j)   * conjg( c2(i,j) )
c
            power(2,1,1) = power(2,1,1) + c1(ii,j)  * conjg( c1(ii,j) )
c
            power(2,1,2) = power(2,1,2) + c2(ii,j)  * conjg( c2(ii,j) )
c
            power(1,2,1) = power(1,2,1) + c1(i,jj)  * conjg( c1(i,jj) )
c
            power(1,2,2) = power(1,2,2) + c2(i,jj)  * conjg( c2(i,jj) )
c
            power(2,2,1) = power(2,2,1) + c1(ii,jj) * conjg( c1(ii,jj) )
c
            power(2,2,2) = power(2,2,2) + c2(ii,jj) * conjg( c2(ii,jj) )
c
         enddo
c
      enddo
c
      power(:,:,:) = power(:,:,:) / float( ( nx / 2 ) * ( ny / 2 ) )
c
      power1avg    = 0.25 * ( power(1,1,1) + power(2,1,1) +
     .                        power(1,2,1) + power(2,2,1) )
c
      power2avg    = 0.25 * ( power(1,1,2) + power(2,1,2) +
     .                        power(1,2,2) + power(2,2,2) )
c
      power(:,:,1) = power(:,:,1) / power1avg
c
      power(:,:,2) = power(:,:,2) / power2avg
c
      write ( 7 , '(1x,4f16.4)' ) power
c
      return
      end
