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
