C**
C***********************************************************************
C**
      subroutine addtgt ( time0 , range0 , freq0 , accel0 , dwdth0 ,
     .                    snr0 , source0 )
C**
C***********************************************************************
C**
c   Purpose:  This routine manages the target list.  It takes the latest
c             target state vector given it and adds that information to
c             a large circular buffer.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c   Use the global target list stored in common
c
      include 'tglist.h'
c
c   Standard parameters
c
      include 'sarprm.h'
c
c   Updated parameters
c
      include 'updates.h'
c
c-----------------------------------------------------------------------
c              
c   Passed Variables
c
      real         time0 , range0 , freq0 , accel0 , dwdth0 , snr0
c
      character    source0*1
c
c   Need a new set of variables for reading from disk because the above
c   call arguments are often set as hard constants - trying to use them
c   as internal variables causes an 'access violation'.
c
      real         time1 , range1 , freq1 , accel1 , dwdth1 , dbsnr1
c
      character    source1*1
c
c   time0 (real,in)                  Time (s) for target 
c   range0 (real,in)                 Range (m) for target 
c   freq0 (real,in)                  Frequency (Hertz) for target.
c   accel0 (real,in)                 Acceleration (Hertz/s) for target 
c   dwdth0 (real,in)                 Doppler width (Hertz) for target
c   snr0 (real,in)                   SNR for target 
c   source0 (real,in)                Letter code for target.  The code
c                                    for a particular target indicates
c                                    the origin of state vector for that
c                                    target:
c
c                                     't' = Came from Target subroutine
c                                     'p' = Came from PGA subroutine
c                                     's' = Came from sub-image array
c 
c
c   (time0,range0,freq0,accel0,dwdth0,snr0,source0) is the state vector
c   for one target.  This information is added to the target list.
c
      integer      iloop , ndot , lastdot
c
      real         db
c
      character    lastfile*80 , filenm*80 , tgtline*80
c
      data lastfile / ' ' /
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      ndot   = lastdot( ofile ) - 1
c
      filenm = ofile(1:ndot) // '.tgt'
c
c   Bail here if we have already read the targets from the file
c
      if ( ( rd_tgt .ne. 0 ) .and. ( filenm .eq. lastfile ) ) then
c
c   If we have already read the targets from disk then write out new
c   targets to the '.tgtout' file to analyze them without using them
c   offline.
c
         write ( 97 , '(1x,6f9.2,3x,a)' ) db( snr0 ) , accel0 , dwdth0 ,
     .                                time0 , range0 , freq0 , source0
c
         return
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Initialize time array with all minus one's  if the routine is called
c   with a negative time
c
      if ( time0 .lt. 0 ) then 
c
         do iloop = 1 , nlist_d
c
            time(iloop)    = - 1.0
            range(iloop)   = 0.0
            freq(iloop)    = 0.0
            accel(iloop)   = 0.0
            dwdth(iloop)   = 0.0
            snr(iloop)     = 0.0
            source(iloop)  = 'n'
            xtgt(iloop)    = 0.0
            ytgt(iloop)    = 0.0
            ddot(iloop)    = 0.0
            ipixtg(iloop)  = 0
            jpixtg(iloop)  = 0
c
         enddo
c
         ndead          = 0
c
         a0             = 0.0
         at             = 0.0
         ar             = 0.0
         af             = 0.0
         art            = 0.0
         aft            = 0.0
c
         iptr           = 0
c
         nlist          = 1
c
         if ( rd_tgt .eq. 0 ) then
c
            lastfile = filenm
c
            return
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   If the flag to read targets from the target file is set and the
c   file name is new then read all the targets from the file.
c
      if ( ( rd_tgt .ne. 0 ) .and. ( filenm .ne. lastfile ) ) then
c
c               Unit 88 - Target information
c
         do iptr = 1 , nlist_d
c
            read ( 88 , '(a)' , end = 100 ) tgtline
c
            if ( index( tgtline , 'pass' ) .eq. 0 ) then
c
               read ( tgtline , * ) dbsnr1 , accel1 , dwdth1 , time1 ,
     .                              range1 , freq1
c
               if      ( index( tgtline , 't' ) .ne. 0 ) then
c
                  source1 = 't'
c
               else if ( index( tgtline , 'p' ) .ne. 0 ) then
c
                  source1 = 'p'
c
               else
c
                  source1 = 's'
c
               endif
c
               time(iptr)   = time1
               range(iptr)  = range1
               freq(iptr)   = freq1
               accel(iptr)  = accel1
               dwdth(iptr)  = dwdth1
               snr(iptr)    = 10.0 ** ( 0.1 * dbsnr1 )
               source(iptr) = source1
c
            endif
c
         enddo
c
  100    continue  !  Come here on end-of-file for target reads
c
         lastfile = filenm
c
         nlist    = max( iptr , 1 )
c
         write ( 6 , * ) nlist , ' Targets read'
c
         write ( 7 , * ) nlist , ' Targets read'
c
         if ( nlist .ge. nlist_d ) then
c
            write ( 6 , * ) ' End of target memory while reading'
c
            write ( 7 , * ) ' End of target memory while reading'
c
         endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      else
c
c   Write target to the sub-image and PGA files and add to the internal
c   list
c
         if ( source0 .eq. 's' ) write ( 88 , '(1x,6f11.3)' )
     .   db( snr0 ) , accel0 , dwdth0 , time0 , range0 , freq0
c
         if ( source0 .eq. 'p' ) write ( 99 , '(1x,6f11.3)' )
     .   db( snr0 ) , accel0 , dwdth0 , time0 , range0 , freq0
c
         if ( ( source0 .eq. 'p' ) .or. ( source0 .eq. 's' ) .or.
     .        ( mod( curve , 3 ) .eq. 0 ) ) then
c
c   Bump pointer ahead to next target record
c
            iptr         = iptr + 1
c
            nlist        = min( nlist + 1 , nlist_d )
c
c   Make sure the pointer wraps around when it is at the end of the
c   circular buffer
c
            if ( iptr .gt. nlist_d ) then
c
               write ( 6 , * ) ' Target list over-run' , iptr , nlist_d
c
               write ( 7 , * ) ' Target list over-run' , iptr , nlist_d
c
               iptr = 1
c
            endif
c
c   Add this target to the list   
c
            time(iptr)   = time0
            range(iptr)  = range0
            freq(iptr)   = freq0
            accel(iptr)  = accel0
            dwdth(iptr)  = dwdth0
            snr(iptr)    = snr0
            source(iptr) = source0
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      return
      end
