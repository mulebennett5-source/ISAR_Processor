c**
c***********************************************************************
c**
      subroutine subimg_rt ( csbimg , nsa , nsr , ntp , ncr , naskip ,
     .                       nakeep , str , wtac , work , ac , cac ,
     .                       nwork )
c**
c***********************************************************************
c**
c   Purpose:  To return a coarse resolution complex image
c
c-----------------------------------------------------------------------
c
      implicit none
c
      include 'sarprm.h'      !  Standard ISAR-T parameters
c
      include 'updates.h'     !  Updates to parameters
c
      integer  nsa , nsr , ntp , ncr , naskip , nakeep ,
     .         ir , ia , isr , nacent , nalose , nwork
c
      real     wtac(3*naskip) , str(2,ntp,ncr) , ac(2,nsa) , work(nwork)
c
      complex  csbimg(nsr,nakeep) , cac(nsa)
c
c-----------------------------------------------------------------------
c
c   Compress in time to form coarse Doppler (sub-image)
c
      nacent = ( 3 * naskip ) / 2
c
c   For ISAR mode the sub-image is simply the Fourier transform with
c   respect to time at fixed range
c
      do isr = 1 , nsr          !  All range cells in sub-image
c
         ir = isr + ( ncr - nsr ) / 2
c
         call wtfft ( str(1,1+naskip/2,ir) , wtac , ac , work ,
     .                + 1 , - nacent , 3 * naskip , nsa )
c
c   Keep only the middle Doppler cells to save memory.  This expression
c   works for both even and odd values of nakeep.  For an odd number it
c   keeps the DC plus an equal number of negative and positive points.
c   For an even number it keeps one more negative than positive.  In
c   the case of the limit of nakeep = nsa, it the first value is the
c   Nyquist.
c
         nalose = max( 0 , ( nsa - nakeep + 1 ) / 2 )
c
         do ia = 1 , nakeep
c
            csbimg(isr,ia) = cac(ia+nalose)
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
      return
      end
