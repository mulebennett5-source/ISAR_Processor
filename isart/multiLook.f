C**
C***********************************************************************
C**
      subroutine multiLook ( imageout , imagein , rLook , fLook ,
     .                       nf , nr )
C**
C***********************************************************************
C**
      implicit none
c
      integer rlook , flook , nf , nr , totalLooks , i , j , ii , jj
c
      complex imagein(nf,nr)
c
      real    imageout(nf,nr,2)
c
c Output:
c       imageOut = Multilooked and median filtered intensity image
c
c Input:
c         imageIn  = Complex image (Doppler X Range)
c         rlook    = Convolution size in range
c         flook    = Convolution size in Doppler
c
c Default [rlook = 1, flook = 1] yields a (3X3 intensity convolution)
c
c In general the convolution is rectangular, (1+2*rlook) X (1+2*flook)
c
c-----------------------------------------------------------------------
c
c     Divide each Doppler line by its median abs
c     intenseR   = median(abs(imageIn),2);
c
c     do i = 1 , nf
c        imageIn(i,:) = imageIn(i,:)/intenseR(i);
c     enddo
c
c-----------------------------------------------------------------------

c     Image intensity for output
      do i = 1 , nf
         do j = 1 , nr
            imageout(i,j,1) = cabs(imageIn(i,j))**2;
         enddo
      enddo

c     Square filter only if rlook or flook > 0

      if ( rlook .gt. 0 .or. flook .gt. 0 ) then

c        Convolve the data with a window of (1+2*rlook) by (1+2*flook) pixels
   
         imageout(:,:,2) = imageout(:,:,1)

         totalLooks      = (1+2*rlook)*(1+2*flook);

         do i = 1+flook,nf-flook
            do j = 1+rlook,nr-rlook
               imageout(i,j,1) = 0;
               do ii = -flook,flook
                  do jj = -rlook,rlook
                     imageout(i,j,1) = imageout(i,j,1)+
     .                                 imageout(i+ii,j+jj,2)/totalLooks
                  enddo
               enddo
            enddo
         enddo
c
      endif
c    
      return
      end