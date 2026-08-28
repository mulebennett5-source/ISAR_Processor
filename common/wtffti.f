C**
C***********************************************************************
C**
      subroutine wtffti ( wt , fft , work , isign , nbands , nd , nfft )
C**
C***********************************************************************
C**
      implicit none
c
      integer isign , nbands , nn(1) , nd , nfft , nfftb , j , jb , jj ,
     .        jsh
c
      real    wt(nd) , fft(2,nfft*nbands) , work(2,nfft*nbands)
c
c-----------------------------------------------------------------------
c
      nfftb = nfft * nbands
c
c   Rotate center to DC
c
      do j = 1 , nfftb
c
         jsh = j - nfftb / 2
c
         if ( jsh .le. 0 ) jsh = jsh + nfftb
c
         work(1,jsh) = fft(1,j)
         work(2,jsh) = fft(2,j)
c
      enddo
c
c   Fourier transform back to time using the opposite sign used in the
c   original call to wtfft
c
      nn(1) = nfftb
c
      call fourt ( work , nn , 1 , - isign , 1 , fft , 2 * nfftb )
c
c   Rotate DC (in time) to center
c
      do j = 1 , nfftb
c
         jsh = j - nfftb / 2
c
         if ( jsh .le. 0 ) jsh = jsh + nfftb
c
         fft(1,jsh) = work(1,j)
         fft(2,jsh) = work(2,j)
c
      enddo
c
c   Divide out the original weight function in time
c
      do j = 1 , nd
c
         do jb = 1 , nbands
c
            jj  = nbands * ( j - 1 ) + jb
c
            jsh = jj + ( nfftb - nd * nbands ) / 2
c
            fft(1,jsh) = fft(1,jsh) / wt(j)
            fft(2,jsh) = fft(2,jsh) / wt(j)
c
         enddo
c
      enddo
c
      return
      end
