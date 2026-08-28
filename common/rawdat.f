C**
C***********************************************************************
C**
      subroutine rawdat ( first , rawfil , nfmt , kpulse , nr , hist ,
     .                    rdata , aux , nerror , work , nwork , nrcent )
C**
C***********************************************************************
C**
      implicit none
c
      logical      first
c
      character    rawfil*80
c
      integer      nfmt , kpulse , nr , reclen , j , jj , nerror ,
     .             jfirst , nrcent , offset , reclda , nwork
c
      byte         rawbyt(2,8320)
c
      integer*2    rawi2(2,8192)
c
      integer      hist(2,65536) , idummy , line , frame , nn(2) , shift
c
      real         rdata(2,nr) , aux(10) , work(2,nwork/2) , scale ,
     .             bfloat
c
      real*8       rcray(2,8192)
c
c  Vars associated with reading NRL data
c
      integer      nrlsamps , nrl6 , nsp3l
      parameter  ( nrlsamps = 736 , nrl6 = 2182 , nsp3l = 256 + 8192 )
c
c  Vars associated with reading NRL data
c
      real         ang , hdgang , hdg0 , alt , ew , ns ,
     .             profile_rwts(nrlsamps)
c
      real*8       t
c
      integer      ispike
c
c-----------------------------------------------------------------------
c
c   Lincoln Lab (exp, i, q)
c
      byte         lle(2052) , lli(2052) , llq(2052)
c
      real         realiq(2,4*8192)
c
c-----------------------------------------------------------------------
c
c   NAWC P-3 L-band:  Chirp filter, chirp rate, A/D rate, chirp duration
c
      integer      lchirp , nfilt , ncells
c
      real         filter(2,8192) , brp3l , dtrp3l , trp3l , tdp3l ,
     .             pi , phase
c
      parameter  ( brp3l = - 30.0E+12 , dtrp3l = 1.0 / 125.0E+6 ,
     .             tdp3l = 4.0E-6 )
c
c-----------------------------------------------------------------------
c
      save         reclen , pi , filter , ncells , nfilt , profile_rwts
c
c-----------------------------------------------------------------------
c
c  Open rawfil each time to handle multi-channel data
c
         if (                    ( nr .gt. 6*8192 ) .or.
     .        ( nfmt .eq. 10 .and. nr .gt.   4096 ) .or.
     .        ( nfmt .eq. 11 .and. nr .gt.   4096 ) .or.
     .        ( nfmt .eq. 12 .and. nr .gt.   4096 ) .or.
     .        ( nfmt .eq. 13 .and. nr .gt.   4096 ) .or.
     .        ( nfmt .eq. 14 .and. nr .gt.   2048 ) .or.
     .        ( nfmt .eq. 15 .and. nr .gt.   2048 ) ) then
c
            write ( 6 , * ) ' Pulse too long in RAWDAT!'
            write ( 7 , * ) ' Pulse too long in RAWDAT!'
c
            go to 2000
c
         endif
c
         pi = atan2( 0.0 , - 1.0 )
c
         if ( nfmt .eq.  0 ) reclen = 2 * nr
         if ( nfmt .eq.  1 ) reclen = 4 * nr
         if ( nfmt .eq.  2 ) reclen = 8 * nr
         if ( nfmt .eq.  4 .or. nfmt .eq. 5 ) reclen = 3 * ( 2052 )
         if ( nfmt .eq.  6 ) reclen = nrl6
         if ( nfmt .eq. 10 ) reclen = nsp3l
         if ( nfmt .eq. 11 ) reclen = nsp3l
         if ( nfmt .eq. 12 ) reclen = 8192
         if ( nfmt .eq. 13 ) reclen = 8192
         if ( nfmt .eq. 14 ) reclen = 4096
         if ( nfmt .eq. 15 ) reclen = 4096
         if ( nfmt .eq. 16 ) reclen = 8 * nr
         if ( nfmt .eq. 17 ) reclen = 16 * nr
         if ( nfmt .eq. 18 ) reclen = 16 * nr
c
         if ( nfmt .ne. 3 ) close ( 1 )
c
         if ( first .and. nfmt .eq. 3 ) then
c
            open ( 99 , file = 'profile_rwts.dat' , status = 'old' ,
     .                  form = 'unformatted' , access = 'direct' ,
     .                  recl = reclda( 0 , nrlsamps ) )
            read ( 99 , rec = 1 ) profile_rwts
c
            close ( 99 )
c
         endif
c
         if ( nfmt .ne. 3 .and. nfmt .ne. 16 )
     .                  open ( 1 , file    = rawfil ,
     .                             form    = 'UNFORMATTED' ,
     .                             access  = 'DIRECT' , err = 2001 ,
     .                             status  = 'OLD' ,
     .                             recl    = reclda( reclen , 0 ) )
c
         if ( nfmt .eq. 16 )
     .                  open ( 1 , file    = rawfil ,
     .                             form    = 'UNFORMATTED' ,
     .                             access  = 'DIRECT' , err = 2001 ,
     .                             status  = 'OLD' ,
     .                             convert = 'big_endian' ,
     .                             recl    = reclda( reclen , 0 ) )
c
         if ( nfmt .eq. 17 )
     .                  open ( 1 , file    = rawfil ,
     .                             form    = 'UNFORMATTED' ,
     .                             access  = 'DIRECT' , err = 2001 ,
     .                             status  = 'OLD' ,
     .                             convert = 'cray' ,
     .                             recl    = reclda( reclen , 0 ) )
c
         if ( nfmt .eq. 18 )
     .                  open ( 1 , file    = rawfil ,
     .                             form    = 'UNFORMATTED' ,
     .                             access  = 'DIRECT' , err = 2001 ,
     .                             status  = 'OLD' ,
     .                             convert = 'big_endian' ,
     .                             recl    = reclda( reclen , 0 ) )
c
         go to 2002
 2001    write ( 6 , * ) ' Cannot open file ' , rawfil
         go to 2000
 2002    continue
c
      if ( first ) then
c
         if ( nfmt .eq. 10 .or. nfmt .eq. 11 ) then
c
c   Compute match filter for NAWC P-3
c
            if ( nfmt .eq. 10 ) then
c
               ncells = 4096       !  Use all the swath, low resolution
c
            else
c
               ncells = nr         !  Use middle pixels, full resolution
c
            endif
c
c   Filter is over-sampled by a factor of 2
c
            nfilt = 2 * ncells
c
            do j = 1 , nfilt
c
               filter(1,j) = 0.0
               filter(2,j) = 0.0
c
            enddo
c
            filter(1,1) = 1.0     !  DC
c
c   Form the de-chirp filter at twice the time sampling rate and twice
c   the number of frequency bins
c
            lchirp = ifix( tdp3l / dtrp3l )
c
            do j = 2 , lchirp
c
               trp3l               = 0.5 * float( j - 1 ) * dtrp3l
               phase               = - brp3l * pi * trp3l ** 2
c
               filter(1,j)         = cos( phase )
               filter(2,j)         = sin( phase )
c
               filter(1,nfilt+2-j) = filter(1,j)
               filter(2,nfilt+2-j) = filter(2,j)
c
            enddo
c
c   Swap signs of even elements to shift filter to middle
c
            do j = 2 , nfilt , 2
c            
               filter(1,j)         = - filter(1,j)
               filter(2,j)         = - filter(2,j)
c
            enddo
c
            nn(1) = nfilt
            call fourt ( filter , nn , 1 , - 1 , 1 , realiq , 8192 )
c
         endif
c
         if ( nfmt .ne. 3 ) first = .false.
c
      endif
c
      do j = 1 , 10
c
         aux(j) = 0.0
c
      enddo
c
      if ( nfmt .eq. 0 ) then
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rawbyt(1,j) , rawbyt(2,j) , j = 1 , nr )
c
         do j = 1 , nr
c
            rdata(1,j) = float( rawbyt(1,j) )
            rdata(2,j) = float( rawbyt(2,j) )
            
c           rdata(1,j) = amax1(-32.0,amin1(32.0,rdata(1,j)))
c           rdata(2,j) = amax1(-32.0,amin1(32.0,rdata(2,j)))
c
         enddo
c
      else if ( nfmt .eq. 1 ) then
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rawi2(1,j) , rawi2(2,j) , j = 1 , nr )
c
         do j = 1 , nr
c
            hist(1,rawi2(1,j)+32768) = hist(1,rawi2(1,j)+32768) + 1
            hist(2,rawi2(2,j)+32768) = hist(2,rawi2(2,j)+32768) + 1
c
            rdata(1,j)               = float( rawi2(1,j) )
            rdata(2,j)               = float( rawi2(2,j) )
c
         enddo
c
      else if ( nfmt .eq. 2 .or. nfmt .eq. 16 ) then
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rdata(1,j) , rdata(2,j) , j = 1 , nr )
c
      else if ( nfmt .eq. 17 ) then
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rcray(1,j) , rcray(2,j) , j = 1 , nr )
c
         do j = 1 , nr
c
            rdata(1,j) = rcray(1,j)
            rdata(2,j) = rcray(2,j)
c
         enddo
c
      else if ( nfmt .eq. 18 ) then
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rcray(1,j) , rcray(2,j) , j = 1 , nr )
c
         do j = 1 , nr
c
            rdata(1,j) = rcray(1,j)
            rdata(2,j) = rcray(2,j)
c
         enddo
c
      else if ( nfmt .eq. 3 ) then
c
         if ( nr .lt. nrlsamps ) go to 2000
c
         call rdpuls ( kpulse , work , ang , hdgang , hdg0 , alt , ew ,
     .                 ns , t , first , rawfil , nrlsamps )
c
c   Zeropad pulse out to nr samples
c
         if ( nr .gt. nrlsamps ) then
c
            do j = nrlsamps + 1, nr
c
               work(1,j) = 0.0
               work(2,j) = 0.0
c
            enddo
c
         endif
c
c   "Un-range compress" the data back to the time domain
c
         nn(1) = nr
         call fourt ( work(1,1) , nn , 1 , -1 , 1 , rdata(1,1) , nr )
c
c   Divide by window function
c
         do j = 1, nr
c
            rdata(1,j) = rdata(1,j) * profile_rwts(j)
c
            rdata(2,j) = rdata(2,j) * profile_rwts(j)
c
         enddo 
c
c   Remove spikes from the data
c
c       DC Spike
c
         rdata(1,1+nr/2) = 0.5 * ( rdata(1,nr/2) + rdata(1,2+nr/2) )
         rdata(2,1+nr/2) = 0.5 * ( rdata(2,nr/2) + rdata(2,2+nr/2) )
c  
c       Nyquist Spike
c  
         rdata(1,1)      = 0.5 * ( rdata(1,2) + rdata(1,nr) )
         rdata(2,1)      = 0.5 * ( rdata(2,2) + rdata(2,nr) )
c  
c       Frequency Spike 1
c  
         ispike          = 1 + nr / 4
c
         rdata(1,ispike) = 0.5 * ( rdata(1,ispike-1) +
     .                             rdata(1,ispike+1) )
         rdata(2,ispike) = 0.5 * ( rdata(2,ispike-1) +
     .                             rdata(2,ispike+1) )
c  
c       Frequency Spike 2
c  
         ispike          = 1 + 3 * nr / 4
c
         rdata(1,ispike) = 0.5 * ( rdata(1,ispike-1) +
     .                             rdata(1,ispike+1) )
         rdata(2,ispike) = 0.5 * ( rdata(2,ispike-1) +
     .                             rdata(2,ispike+1) )
c  
c       Frequency Spike 3
c  
         ispike          = 234
c
         rdata(1,ispike) = 0.5 * ( rdata(1,ispike-1) +
     .                             rdata(1,ispike+1) )
         rdata(2,ispike) = 0.5 * ( rdata(2,ispike-1) +
     .                             rdata(2,ispike+1) )
c  
c       Frequency Spike 4
c  
         ispike          = 504
c
         rdata(1,ispike) = 0.5 * ( rdata(1,ispike-1) +
     .                             rdata(1,ispike+1) )
         rdata(2,ispike) = 0.5 * ( rdata(2,ispike-1) +
     .                             rdata(2,ispike+1) )
c  
      else if ( nfmt .eq. 4 .or. nfmt .eq. 5 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Lincoln Lab 2048 sample records plus 12 byte header in 8-8-4 format
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( lle(j) , lli(j) , llq(j) , j = 1 , 2052 )
c
         shift = ( 2048 - nr ) / 2
c
         j = 0
         do jj = 1 , 2052
c
            if ( ( lle(jj) .and. 64 ) .ne. 0 ) then
c
c   If a certain bit of the exponent is set, then the triplet represents
c   auxiliary information - decode it but don't output any of it at this
c   time
c
c-----------------------------------------------------------------------
c
               if ( ( ( lle(jj) .and. 12 ) / 4 ) .eq. 0 ) then
c
                  line = 256 * int( lli(jj) .and. 31 ) +
     .                   ( llq(jj) .and. 255 )
c
               endif
c
               if ( ( ( lle(jj) .and. 12 ) / 4 ) .eq. 1 ) then
c
                  frame = ( llq(jj) .and. 255 )
c
               endif
c
c-----------------------------------------------------------------------
c
            else
c
c   If bit not set, the triplet is a complex data value stored in 8-8-4
c   format
c
c   First, calculate the exponent, using the four least significant
c   bits of lle
c
               j          = j + 1
               scale      = ( 2.0 ** ( lle(jj) .and. 15 ) ) / 4096.0
c
               realiq(1,j) = float( lli(jj) ) * scale
               realiq(2,j) = float( llq(jj) ) * scale
c
            endif
c
         enddo
c
c   Since the above values are range-compressed, FFT it to convert back
c   to fast time
c
         if ( nfmt .eq. 4 ) then
c
c   Preserve swath
c
            nn(1) = 2048
            call fourt ( realiq , nn , 1 , - 1 , 1 , work , 1 )
c
         else
c
c   Preserve resolution
c
            nn(1) = nr
            shift = min( 2048 - nr , nrcent + ( 2048 - nr ) / 2 )
            call fourt ( realiq(1,1+shift) , nn , 1 , - 1 , 1 ,
     .                   work , 1 )
c
         endif
c
c
c   Change sign of alternate pulses so that the DC bias falls at the
c   near and far range of the swath and so that the target is in the
c   center
c
         do j = 1 + shift , nr + shift , 2
c
            realiq(1,j) = - realiq(1,j)
            realiq(2,j) = - realiq(2,j)
c
         enddo
c
c   Finally, load the output I/Q array with the middle part of the
c   fast time series
c
         do j = 1 , nr
c
            rdata(1,j) = realiq(1,j+shift)
            rdata(2,j) = realiq(2,j+shift)
c
c   As a diagnostic, calculate the histogram of I/Q, assuming 16 bit
c   dynamic range
c
            idummy = 32768 +
     .               min( 32768 , max( -32767 , int( rdata(1,j) ) ) )
            hist(1,idummy) = hist(1,idummy) + 1
            idummy = 32768 +
     .               min( 32768 , max( -32767 , int( rdata(2,j) ) ) )
            hist(2,idummy) = hist(2,idummy) + 1
c
         enddo
c
c                      End of Lincoln Lab format
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   New (1993) NRL Format
c
      else if ( nfmt .eq. 6 ) then
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rawbyt(1,j) , rawbyt(2,j) , j = 1 , nrl6 / 2 )
c
c   Make range cell 301 the center if less than 600 cells are asked for
c
         jfirst = max( 1 , 301 - nr / 2 )
c
         do j = jfirst , jfirst - 1 + nr
c
            jj = j + 1 - jfirst
c
c   Flip sign of alternate pulses and byte swap
c
            if ( mod( kpulse , 2 ) .eq. 1 ) then
c
               rdata(2,jj) = + float( rawbyt(1,j+3) )
               rdata(1,jj) = + float( rawbyt(2,j+3) )
c
            else
c
               rdata(2,jj) = - float( rawbyt(1,j+3) )
               rdata(1,jj) = - float( rawbyt(2,j+3) )
c
            endif
c
         enddo
c
c   Perform shift of the data to put center (300) at DC
c
         if ( nr .lt. 1024 ) then
c
            call rshift ( rdata(1,1) , 2 * nr , nr - 2 , work(1,1) )
c
         else
c
            call rshift ( rdata(1,1) , 2 * nr , 600 , work(1,1) )
c
         endif
c
c   "Un-range-compress" the data back to the time domain
c
         nn(1) = nr
         call fourt ( work(1,1) , nn , 1 , - 1 , 1 , rdata(1,1) , nr )
c
c   Perform final shift of the data to correct time format
c
         call rshift ( work(1,1) , 2 * nr , nr - 2 , rdata(1,1) )
c
      else if ( nfmt .eq. 10 .or. nfmt .eq. 11 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                        NAWC P-3 L-Band SAR
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rawbyt(1,j) , rawbyt(2,j) , j = 1 , 4224 )
c
c   Convert to real, skipping header and put DC at first element
c
         offset = min( nrcent , ( 4096 - ncells ) / 2 )
         offset = max( nrcent , - ( 4096 - ncells ) / 2 )
         shift  = 128 + offset + ( 2048 - ncells / 2 )
c
         do j = 1 , ncells / 2
c
            realiq(1,j)          = float( rawbyt(1,j+shift+ncells/2) )
            realiq(2,j)          = float( rawbyt(2,j+shift+ncells/2) )
c
            realiq(1,j+ncells/2) = float( rawbyt(1,j+shift) )
            realiq(2,j+ncells/2) = float( rawbyt(2,j+shift) )
c
         enddo
c
c   Swap signs of even elements to shift DC to middle in output array
c
         do j = 2 , ncells , 2
c
            realiq(1,j) = - realiq(1,j)
            realiq(2,j) = - realiq(2,j)
c
         enddo
c
         nn(1) = ncells
         call fourt ( realiq , nn , 1 , - 1 , 1 , work(1,1) , nr )
c
c   Multiply by de-chirp filter
c
         shift = ( ncells - nr ) / 2
c
         do j = 1 , nr
c
            rdata(1,j) = + realiq(1,j+shift) * 
     .                     filter(1,j+shift+ncells/2)
            rdata(2,j) = + realiq(2,j+shift) *
     .                     filter(2,j+shift+ncells/2)
c
         enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      else if ( nfmt .eq. 12 .or. nfmt .eq. 13 .or.
     .          nfmt .eq. 14 .or. nfmt .eq. 15 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                        SRI UWB SAR
c
         read ( 1 , rec = kpulse , err = 2000 )
     .        ( rawbyt(1,j) , rawbyt(2,j) , j = 1 , reclen / 2 )
c
         if ( nfmt .eq. 12 .or. nfmt .eq. 14 ) then
c
            ncells = reclen / 2   !  Use all the swath, low resolution
c
         else
c
            ncells = nr           !  Use middle pixels, full resolution
c
         endif
c
c   Convert to real, put DC at first element
c
         offset = min( nrcent , ( reclen / 2 - ncells ) / 2 )
         offset = max( nrcent , - ( reclen / 2 - ncells ) / 2 )
         shift  = offset + ( reclen / 4 - ncells / 2 )
c
         do j = 1 , ncells / 2
c
            jj = j + shift + ncells / 2
c
            if ( jj .le. reclen / 2 ) then
c
               realiq(2,j) = bfloat( rawbyt(1,jj) )
               realiq(1,j) = bfloat( rawbyt(2,jj) )
c
            else
c
               realiq(2,j) = 0.0
               realiq(1,j) = 0.0
c
            endif
c
            jj = j + shift
c
            if ( jj .le. reclen / 2 ) then
c
               realiq(2,j+ncells/2) = bfloat( rawbyt(1,jj) )
               realiq(1,j+ncells/2) = bfloat( rawbyt(2,jj) )
c
            else
c
               realiq(2,j) = 0.0
               realiq(1,j) = 0.0
c
            endif
c
         enddo
c
c
c   Swap signs of even elements to shift DC to middle in output array
c
         do j = 2 , ncells , 2
c
            realiq(1,j) = - realiq(1,j)
            realiq(2,j) = - realiq(2,j)
c
         enddo
c
         nn(1) = ncells
         call fourt ( realiq , nn , 1 , - 1 , 1 , work(1,1) , nr )
c
         shift = ( ncells - nr ) / 2
c
         do j = 1 , nr
c
            rdata(1,j) = realiq(1,j+shift)
            rdata(2,j) = realiq(2,j+shift)
c
         enddo
c
c   Zero DC to remove A-to-D bias
c
         rdata(1,1+nr/2) = 0.0
         rdata(2,1+nr/2) = 0.0
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      endif
c
      nerror = 0
      return
c
 1000 nerror = 1
      return
c
 2000 nerror = 2
c
      return
      end
C**
C***********************************************************************
C**
      real function bfloat ( inbyte )
C**
C***********************************************************************
C**
      implicit none
c
      byte inbyte
c
      real rbyte
c
      rbyte = float( inbyte )
c
      if ( rbyte .lt. 0.0 ) then
c
         bfloat = rbyte + 128.0
c
      else
c
         bfloat = rbyte - 128.0
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      SUBROUTINE RDPULS ( JPULSE , RDATA , ANG , HDGANG , HDG0 , ALT ,
     .                    EW , NS , T , INIT , IFILE , NR )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      CHARACTER        IFILE*80
C
      INTEGER          I ,  JPULSE , OLDREC , OFFSTD , NR , RECLDA
C
      REAL             RDATA(2,NR)
C
      REAL             ANG , HDGANG , ALT , BAM , RINPHS , RQUADR , PI ,
     .                 PRNANG , EW , NS , HDG0 , SIGNMULT
C
      DOUBLE PRECISION T , DNTIME
C
      BYTE             BYTEBUFF(1536) , BANNOT(120)
C
      SAVE             PI , BAM , OLDREC , BYTEBUFF
C
      LOGICAL          INIT
C
      INCLUDE         'aps137hd.h'
C
      OFFSTD = 64
C
      IF ( INIT ) THEN
C
         OPEN ( UNIT = 2 , FILE = IFILE , STATUS = 'OLD' ,
     .          FORM = 'UNFORMATTED' , ACCESS = 'DIRECT' ,
     .          RECL = RECLDA( 1536 , 0 ) )
C
         OLDREC = 0
C
         PI     = ATAN2( 0.0 , - 1.0 )
         BAM    = ( 2.0 * PI ) / 4096.0
C
C        WRITE ( 7 , '(/,1X,A48,A31,/)' )
C    .         '    TIME        ACLAT    ACLON    ACALT    RANGE' ,
C    .         '    NSVEL    EWVEL HDG  ANG QU '
C
      ELSE IF ( JPULSE .LT. 0 ) THEN
C
         CLOSE ( 2 )
         RETURN
C
      ENDIF
C
      READ ( 2 , REC = JPULSE ) BYTEBUFF
C
      DO I = 1 , 38
C
         BANNOT(2*I-1) = BYTEBUFF(I)
C
      ENDDO
C
      CALL HDR137 ( BANNOT )
C
      HDGANG = PRNANG( ( 1024 - ACHDG ) * BAM )
C
      IF ( INIT ) THEN
C
         HDG0 = HDGANG
         INIT = .FALSE.
C
      ENDIF
C
      EW     = EWVEL
      NS     = NSVEL
      T      = DNTIME( TIME )
      ANG    = PRNANG( ( 1024 - ANTAZ ) * BAM ) - HDG0
      ANG    = PRNANG( ANG )
      ALT    = ACALT
C
C        -------------------------------------------
C        FLIP THE SIGN OF THE DATA EVERY OTHER PULSE
C        -------------------------------------------
C
      IF ( MOD( JPULSE , 2 ) .EQ. 0 ) THEN
C
         SIGNMULT = -1 
C
      ELSE
C
         SIGNMULT = 1
C
      ENDIF 
C
      DO I = 1 , NR
C
         RINPHS     = BYTEBUFF(2*I-1+OFFSTD)
         RQUADR     = BYTEBUFF(2*I+OFFSTD)
         RDATA(1,I) = SIGNMULT * RINPHS
         RDATA(2,I) = SIGNMULT * RQUADR
C     
      ENDDO
C
      RETURN
      END
C**
C***********************************************************************
C**
      DOUBLE PRECISION FUNCTION DNTIME ( TIME )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      INTEGER*2 TIME(4)
C
      DNTIME = 0.001 * TIME(4) + TIME(3) + 60.0 * TIME(2) +
     .         3600.0 * TIME(1) 
C
      RETURN
      END
C**
C***********************************************************************
C**
      SUBROUTINE HDR137 ( BANNOT )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      INTEGER     I
C
      BYTE        BANNOT(120) , BUFFER(38)
C
      INTEGER*2   W_TEMP(19) , ITEMP , IVAL(120) , IBITS2
C
      EQUIVALENCE ( BUFFER(1) , W_TEMP(1) )
C
      INCLUDE 'aps137hd.h'
C
C Properly arrange bytes
C
      DO I = 1 , 120
C
         IVAL(I) = BANNOT(I)
C
      ENDDO
C
      DO I = 1 , 38 , 2
C
         BUFFER(I)   = BANNOT(((I-1)*2)+1)
         BUFFER(I+1) = BANNOT(((I-1)*2)+3)
C
      ENDDO
C
      DTYPE   = IBITS2(W_TEMP(1),12,4)
      ITEMP   = W_TEMP(1)
C
C sign extend if necessary
C
      IF ( BTEST(W_TEMP(1),9) ) THEN
C
         ITEMP = IOR(ITEMP,64512)
C
      ELSE
C
         ITEMP = IAND(ITEMP,1023)
C
      ENDIF
C
      ACALT   = ITEMP * 200. * .3048
C
      SPARR   = IBITS2(W_TEMP(2),8,8)
C
      DATE(1) = IBITS2(W_TEMP(2),4,4) * 10 + IBITS2(W_TEMP(2),0,4)
C
      DATE(2) = IBITS2(W_TEMP(3),12,4) * 10 + IBITS2(W_TEMP(3),8,4)
C
      DATE(3) = IBITS2(W_TEMP(3),4,4) * 10 + IBITS2(W_TEMP(3),0,4)
C
      ACLAT   = W_TEMP(4) * 0.0054932
C
      ACLON   = W_TEMP(5) * 0.0054932
C
      HIRES   = BTEST(W_TEMP(6),15)
C
      RANGE   = IBITS2(W_TEMP(6),0,15) * .0078125 * 1852.
C
      ACHDG   = IBITS2(W_TEMP(7),4,12)
C
      RMODE   = IBITS2(W_TEMP(7),0,4)
C
      ANTAZ   = IBITS2(W_TEMP(8),4,12)
C
      DQUAL   = IBITS2(W_TEMP(8),1,3)
C
      ASVAL   = BTEST(W_TEMP(8),0)
C
      TAGNO   = IBITS2(W_TEMP(9),6,10)
C
      SSTAT   = IBITS2(W_TEMP(9),3,3)
C
      FPATH   = IBITS2(W_TEMP(9),0,3)
C
      ITEMP   = W_TEMP(10)
C
C   Sign extend if necessary
C
      IF ( BTEST(ITEMP,14) ) THEN
C
         ITEMP = IBSET(ITEMP,15)
C
      ELSE
C
         ITEMP = IBCLR(ITEMP,15)
C
      ENDIF
C
      NSVEL   = ITEMP * .0390625 * 1852. / 3600.
C
      VERTP   = BTEST(W_TEMP(11),15)
C
      ITEMP   = W_TEMP(11)
C
C sign extend if necessary
C
      IF ( BTEST(ITEMP,14) ) THEN
C
         ITEMP = IBSET(ITEMP,15)
C
      ELSE
C
         ITEMP = IBCLR(ITEMP,15)
C
      ENDIF
C
      EWVEL   = ITEMP * .0390625 * 1852. / 3600.
C
      ASANG   = IBITS2(W_TEMP(12),7,9)
C
      ATILT   = IBITS2(W_TEMP(12),2,5)
C
      APRTR   = IBITS2(W_TEMP(12),0,2)
C
      ASERR   = IBITS2(W_TEMP(13),8,8)
C
      ECODE   = IBITS2(W_TEMP(13),0,3)
C
      TIME(3) = IBITS2(W_TEMP(14),12,4) * 10 + IBITS2(W_TEMP(14),8,4)
C
      TIME(4) = IBITS2(W_TEMP(14),4,4) * 10 + IBITS2(W_TEMP(14),0,4)
C
      TIME(1) = IBITS2(W_TEMP(15),12,4) * 10 + IBITS2(W_TEMP(15),8,4)
C
      TIME(2) = IBITS2(W_TEMP(15),4,4) * 10 + IBITS2(W_TEMP(15),0,4)
C
      RNGAV   = BTEST(W_TEMP(16),15)
C
      ITEST   = BTEST(W_TEMP(16),14)
C
      PULSE   = IBITS2(W_TEMP(16),11,3)
C
      PRFNO   = IBITS2(W_TEMP(16),8,3)
C
      IF ( PRFNO .EQ. 1 ) THEN
C
         PRFNO = 479
C
      ELSE IF ( PRFNO .EQ. 3 ) THEN
C
         PRFNO = 509
C
      ELSE IF ( PRFNO .EQ. 4 ) THEN
C
         PRFNO = 418
C
      ELSE IF ( PRFNO .EQ. 5 ) THEN
C
         PRFNO = 1707
C
      ELSE IF ( PRFNO .EQ. 7 ) THEN
C
         PRFNO = 2026
C
      ENDIF
C
      SWATH   = IBITS2(W_TEMP(16),4,4) * 65536
C
      TIME(4) = TIME(4) * 10 + IBITS2(W_TEMP(16),0,4)
C
      SWATH   = SWATH + W_TEMP(17)
C
C   Correct for misinterpretation of swath lsw as being negative
C   if msb is set
C
      IF ( BTEST(W_TEMP(17),15) ) SWATH = SWATH + 65536
C
      TLOCK   = BTEST(W_TEMP(18),15)
C
      PHASE   = IBITS2(W_TEMP(18),9,6) * 5.625
C
      TKCEL   = IBITS2(W_TEMP(18),0,9)
C
      AGCGN   = IBITS2(W_TEMP(19),8,8) * .25
C
      MANGN   = IBITS2(W_TEMP(19),0,8) * .25
C
      RETURN
      END
C**
C***********************************************************************
C**
      REAL FUNCTION PRNANG ( ANG )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      REAL  ANG , COSANG , SINANG
C
      COSANG = COS( ANG )
      SINANG = SIN( ANG )
C
      PRNANG = ATAN2( SINANG , COSANG )
C
      RETURN
      END
C**
C***********************************************************************
C**
      INTEGER*2 FUNCTION IBITS2 ( WT , M , N )
C**
C***********************************************************************
C**
      IMPLICIT NONE
C
      INTEGER*2 WT
C
      INTEGER   M , N , NN
C
      NN     = MIN( N , 16 - M )
C
      IBITS2 = MOD( ( WT + 65536 ) / ( 2 ** M ) , 2 ** NN )
C
      RETURN
      END
