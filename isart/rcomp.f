C**
C***********************************************************************
C**
      subroutine rcomp ( br , clight , dr0 , rc , crc , pt , ntr , ncr ,
     .                   wtr , dtr , work , nwork )
C**
C***********************************************************************
C**
c
c   Purpose:  Range compresses a single pulse given the desired values
c             of the mo-comp center, dr0, and the mo-comp velocity, vr0.
c
c***********************************************************************
c
      implicit none
c
      integer ntr , ncr , ir , nwork
c
      real    br , clight , dr0 , pt , phase , wtr(ntr) , work(nwork) ,
     .        const , dtr
c
      complex cpfast
c
      real    rc(2,ncr)
      complex crc(ncr)
c
      const = ( 2.0 * 1.0E+12 * br * dr0 ) / clight
c
      do ir = 1 , ntr
c
         phase   = pt + const * dtr * float( ir - 1 - ntr / 2 )
c
         crc(ir) = crc(ir) * cpfast( phase )
c
      enddo
c
c   Compress in range and store in range-compressed array.  Allow for
c   up-chirp (br>0) or down-chirp (br<0)
c
      if ( br .gt. 0.0 ) then
c
         call wtfft ( rc , wtr , rc , work , + 1 , - ntr / 2 ,
     .                ntr , ncr )
c
      else
c
         call wtfft ( rc , wtr , rc , work , - 1 , - ntr / 2 ,
     .                ntr , ncr )
c
      endif
c
      return
      end
