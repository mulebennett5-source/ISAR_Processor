C**
C***********************************************************************
C**
      integer function reclda0 ( nbytes , nreals )
C**
C***********************************************************************
C**
c   The purpose of this routine is to return the system-dependent
c   value in the open statement for a direct access file.
c
      implicit none
c
      integer nbytes    !  Number of bytes desired for the record
c
      integer nreals    !  Number of reals for the record
c
c-----------------------------------------------------------------------
c
c     reclda0 = nbytes + 4 * nreals             !  Sun and PC-MicroSoft
c
      reclda0 = ( nbytes + 3 ) / 4 + nreals     !  Vax, Iris, and PC-DEC
c
c     reclda0 = nbytes + 8 * nreals             !  Cray
c
      return
      end
C**
C***********************************************************************
C**
      integer function reclda ( nbytes , nreals )
C**
C***********************************************************************
C**
c   The purpose of this routine is to return the system-dependent
c   value in the open statement for a direct access file.
c
      implicit none
c
      integer       nbytes            !  Number of bytes per record
c
      integer       nreals            !  Number of reals per record
c
      integer       byte_record_size  !  Unit of direct access record
c
      integer       bytes_per_word    !  Bytes per word for computer
c
      logical       first
c
      character     testbytes(16)*1
c
      integer       testint , j
c
      equivalence ( testbytes , testint )
c
      save          first , byte_record_size , bytes_per_word
c
      data          first / .true. /
c
c-----------------------------------------------------------------------
c
      if ( first ) then
c
         open ( 8 , file = 'c:\isar\foofoo.foo' , form = 'unformatted' ,
     .              status = 'unknown' , access = 'direct' , recl = 4 )
c
         write ( 8 , rec = 1 , err = 10 ) testbytes
c
         byte_record_size = 4
c
         go to 100
c
   10    close ( 8 , status = 'delete' )
c
         open ( 8 , file = 'c:\isar\foofoo.foo' , form = 'unformatted' ,
     .              status = 'unknown' , access = 'direct' , recl = 16 )
c
         write ( 8 , rec = 1 , err = 20 ) testbytes
c
         byte_record_size = 1
c
         go to 100
c
   20    write ( 6 , * ) ' Both tests for byte_record_size failed'
c
  100    close ( 8 , status = 'delete' )
c
         do j = 1 , 16
c
            testbytes(j) = char(255)
c
         enddo
c
         testint = 0
c
         bytes_per_word = 0
c
         do j = 2 , 16
c
            if ( ( bytes_per_word .eq. 0 )      .and.
     .           ( testbytes(j) .ne. char(0) ) ) then
c
               bytes_per_word = j - 1
c
            endif
c
         enddo
c
         first = .false.
c
      endif
c
      reclda = ( nbytes + ( byte_record_size - 1 ) ) /
     .                  byte_record_size
     .         + nreals * ( bytes_per_word / byte_record_size )
c
      return
      end
C**
C***********************************************************************
C**
      real function dtime ( rtime )
C**
C***********************************************************************
C**
      implicit none
c
      real rtime(2) , tnow
c
      real tlast / 0.0 /
c
      save tlast
c
      call cpu_time( tnow )
c
      dtime = tnow - tlast
c
      tlast = tnow
c
      return
      end