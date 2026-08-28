C**
C***********************************************************************
C**
      subroutine getraw ( ntr , ntp , isub0 , isub , nt1 , naskip ,
     .                    spulse , npulse , ncr , rc_raw , crc_raw ,
     .                    cstr , nerror , work , nwork , first ,
     .                    nchuse , nr1 , nre , curtime )
C**
C***********************************************************************
C**
c
c   Purpose:  Fills a portion of the raw data array with a burst of
c             pulses and rearranges the old pulses.  The pulses are
c             obtained either from a disk file or from the internal
c             simulator.
c
c***********************************************************************
c
      implicit none
c
      include     'sarprm.h'
c
      include     'updates.h'
c
      integer      ntr , ntp , ncr , isub0 , isub , nt1 , naskip , ich ,
     .             spulse , npulse , nerror , nedit , npdone , ip , ir ,
     .             nwork , nchuse , off0 , nchtot , nchtotrc , nr1 ,
     .             nre , offset , offset_max , ncancel , nleft ,
     .             ncenter , nright , toffset , toffset_max , ntrwin ,
     .             ntpwin
c
c  For multi-channel mode, output channel is stored as the first layer
c
      real         rc_raw(2,ntr,ntp,nchuse+1)  !  Raw data array -
      complex      crc_raw(ntr,ntp,nchuse+1)   !  possibly multi-channel
c
      complex      cstr(ntp,ncr,nchuse+1)      !  Range-compressed data
c
      real         auxil(10) , work(nwork)
c
      character    rawfil*80
c
      complex      c21 , c13 , c23 , cw1 , cw2 , cov2 , cwt1 , cwt2
c
      real         c11 , c22 , c33 , c21max , cerror , cerrmin ,
     .             curtime , db , dbcancel , dbexpect
c
      logical      first , cfirst
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c            Stage 1 - rearrange the pulses in memory
c
      if ( isub .le. 1 ) then
c
c   For first sub-image load all pulses
c
         nt1 = 1
c
      else
c
         nt1 = 1 + ntp - naskip
c
c   After first sub-image, move over the raw data for the pulses
c   already compressed
c
         if ( nchuse .eq. 1 ) then
c
            nchtot = 1              !  1 channel, 1 level
c
         else
c
            nchtot = nchuse + 1     !  Multiple channels plus room
c                                   !  for cancel channel
c
         endif
c
c   Move over raw data
c
         do ich = 1 , nchtot
c
            do ip = 1 , ntp - naskip
c
               do ir = 1 , ntr
c
                  crc_raw(ir,ip,ich) = crc_raw(ir,ip+naskip,ich)
c
               enddo
c
            enddo
c
         enddo
c
c   Move the compressed pulses also if there are any - this allows
c   this routine to be called by raw data analysis programs as well
c   as the ISAR processor
c
         if ( ncr .gt. 0 ) then
c
            if ( nchuse .gt. 1 .and. cmodel .gt. 0 ) then
c
               nchtotrc = nchuse + 1
c
            else
c
               nchtotrc = 1
c
            endif
c
            do ich = 1 , nchtotrc
c
               do ir = 1 , ncr
c
                  do ip = 1 , ntp - naskip
c
                     cstr(ip,ir,ich) = cstr(ip+naskip,ir,ich)
c
                  enddo
c
               enddo
c
            enddo
c
         endif     !  ncr .gt. 0
c
      endif        !  isub .le. 1
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c          Stage 2 - Initialization of multi-channel alignment
c
      if ( first ) then
c
c   Adaptively find the best offset to the cancellation channel
c
         do ich = 1 , maxchan
c
            pulse_offset(1,ich) = 0    !  Fast-time pulse offset
            pulse_offset(2,ich) = 0    !  Slow-time pulse offset
c
         enddo
c
         if ( nchuse .eq. 2 .and. cmodel .gt. 1 ) then
c
c   Fill the buffer with data channels
c
            cfirst = first
c
            call fillbuff ( 1 , ntp , spulse , npulse , isub , ntr ,
     .                      nchuse , nerror , naskip , nwork , rawfil ,
     .                      crc_raw , rc_raw , cfirst , auxil , work )
c
            if ( nerror .ne. 0 ) return
c
c-----------------------------------------------------------------------
c
c   High-pass the data to remove the effect of the target.  Use the
c   first level as work space for now (later it is used for the cancel
c   channel).  This is done only for the purpose of calculating the
c   optimum lags - the whole raw data buffer is re-read from disk below.
c
            do ich = 2 , nchuse + 1     !  Just the data
c
               do ip = 1 , ntp          !  All pulses in the buffer
c
                  do ir = 2 , ntr - 1   !  All but endpoints
c
                     crc_raw(ir,ip,1)   = 0.25 * crc_raw(ir+1,ip,ich)
     .                                  + 0.25 * crc_raw(ir-1,ip,ich)
     .                                  +  0.5 * crc_raw(ir,ip,ich)
c
                  enddo                 !  RF
c
                  crc_raw(1,ip,1)       = cmplx( 0.0 , 0.0 )
c
                  crc_raw(ntr,ip,1)     = cmplx( 0.0 , 0.0 )
c
c   Subtract off low-pass field
c
                  do ir = 1 , ntr
c
                     crc_raw(ir,ip,ich) = crc_raw(ir,ip,ich)
     .                                  - crc_raw(ir,ip,1)
c
                  enddo
c
               enddo                    !  Pulses
c
            enddo                       !  Channels
c
c-----------------------------------------------------------------------
c
c                   Find the optimum offset
c
c   The original data is stored starting at the second level so that
c   the cancelled channel is always in the first level.
c
            ncancel  = 1   !  Cancelled data used in the ISAR processor
            nleft    = 2   !  'Left' channel
            ncenter  = 3   !  'Center' channel
            nright   = 4   !  'Right' channel
c
            off0     = 4 + naskip / 2
c
            ntpwin   = ntp - off0 - 2
c
            ntrwin   = ntr - off0 - 2
c
            c21max   = 0.0
c
            dbexpect = 0.0
c
            dbcancel = 0.0
c
c   Variance of the cancel channel
c
            c11 = cabs( cov2( crc_raw(1,1,nleft) ,
     .                        crc_raw(1,1,nleft) ,
     .                        ntr , ntp , ntrwin ,
     .                        ntpwin , 0 , 0 , 0 ) )
c
            c33 = c11
c
c   Variance of the primary channel
c
            c22 = cabs( cov2( crc_raw(1,1,ncenter) ,
     .                        crc_raw(1,1,ncenter) ,
     .                        ntr , ntp , ntrwin ,
     .                        ntpwin , 0 , 0 , 0 ) )
c
c   Covariance of the above two cancellation channels
c
            c13 = cov2( crc_raw(1,1,nleft) ,
     .                  crc_raw(1,1,nleft) ,
     .                  ntr , ntp , ntrwin ,
     .                  ntpwin , 0 , - 1 , 0 )
c
            do offset = - off0 , off0
c
               do toffset = - off0 , off0
c
c   Covariance and variance of left channel at previous pulse
c
                  c21 = cov2( crc_raw(1,1,ncenter) ,
     .                        crc_raw(1,1,nleft) ,
     .                        ntr , ntp , ntrwin ,
     .                        ntpwin , toffset , offset , 0 )
c
c   Covariance and variance of right channel at next pulse
c
                  c23 = cov2( crc_raw(1,1,ncenter) ,
     .                        crc_raw(1,1,nleft) ,
     .                        ntr , ntp , ntrwin ,
     .                        ntpwin , toffset , offset - 1 , 0 )
c
c   Compute the weights
c
                  if ( cmodel .eq. 2 ) then
c
                     cw1      = - conjg( c21 ) / c11
c
                     dbcancel = db( 1.0 - cabs( cw1 ) ** 2 )
c
                     if ( cabs( c21 / c11 ) .gt. c21max ) then
c
                        c21max      = cabs( c21 / c11 )
                        offset_max  = offset
                        toffset_max = toffset
c
                        dbexpect    = dbcancel
c
                        cwt1        = cw1
                        cwt2        = cmplx( 0.0 , 0.0 )
c
                     endif
c
                  else if ( cmodel .eq. 3 ) then 
c
                     call canwgts ( cw1 , cw2 , c11 , c22 , c33 , c21 ,
     .                              c23 , c13 , cerror )
c
                     dbcancel = db( cerror / c22 )
c
                     if ( ( toffset .eq. - off0 .and.
     .                      offset .eq. - off0 )
     .                     .or. ( cerror .lt. cerrmin ) ) then
c
                        cerrmin     = cerror
                        offset_max  = offset
                        toffset_max = toffset
c
                        dbexpect    = dbcancel
c
                        cwt1        = cw1
                        cwt2        = cw2
c
                     endif
c
                  endif
c
                  write ( 83 , '(1x,a30,2i6,6f10.3)' )
     .                  ' Offset, Cancellation (dB):  ' ,
     .                    toffset , offset , dbcancel
c
                  write (  6 , '(1x,a30,2i6,6f10.3)' )
     .                  ' Offset, Cancellation (dB):  ' ,
     .                    toffset , offset , dbcancel
c
c                 write ( 83 , 1 )
c                 write ( 83 , 1 ) ' Cross correlation:          ' , c21
c
c                 write ( 83 , 1 ) ' Secondary variance:         ' , c11
c
c                 write ( 83 , 1 ) ' Cross correlation:          ' , c23
c
c                 write ( 83 , 1 ) ' Covar. of cancel channels:  ' , c13
c
c                 write ( 83 , 1 ) ' Complex weight:             ' , cw1
c
c                 write ( 83 , 1 ) ' Complex weight:             ' , cw2
c                 write ( 83 , 1 )
c
    1             format ( a33 , 2f18.4 )
c
               enddo
c
            enddo     !  Search over lags
c
c   In fast time, always shift the cancel channel
c
            pulse_offset(1,1) = - toffset_max
            pulse_offset(1,2) = 0
c
c   In slow time, shift whichever channel will give a positive
c   shift value
c
            if ( offset_max .gt. 0 ) then
c
c   Shift primary channel
c
               pulse_offset(2,2) = + offset_max
c
            else
c
c   Shift cancel channel
c
               pulse_offset(2,1) = - offset_max
c
            endif
c
            write ( 83 , * )
            write ( 83 , * ) ' Lags of max. correlation: ' ,
     .                         toffset_max , offset_max
c
            write ( 83 , * )
            write ( 83 , 1 ) ' Complex weight:             ' , cwt1
c
            write ( 83 , 1 ) ' Complex weight:             ' , cwt2
            write ( 83 , 1 )
c
            write ( 83 , * ) ' Expected Cancellation:      ' , dbexpect
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   For isub <= 1, the raw data is the same - skip this section except
c   for the first time
c
      if ( isub .eq. isub0 .or. isub .gt. 1 ) then
c
c   Read in all raw data for a target buffer and sub-image
c
         cfirst = first
c
         call fillbuff ( nt1 , ntp , spulse , npulse , isub , ntr ,
     .                   nchuse , nerror , naskip , nwork , rawfil ,
     .                   crc_raw , rc_raw , cfirst , auxil , work )
c
         if ( nerror .ne. 0 ) return
c
c-----------------------------------------------------------------------
c
c   Remove anomalous points, i and q offsets, rotate i and q to
c   remove correlation, and rescale i and q to have the same 
c   variances.  Apply only for single channel cases.
c
         if ( fixiq .ne. 0 .and. nchuse .eq. 1 ) then
c
            npdone = nt1 - 1
c
            call fixem ( rc_raw , ntr , ntp , work(1) , work(ntr+1) ,
     .                   work(2*ntr+1) , work(3*ntr+1) , work(4*ntr+1) ,
     .                   nedit , npdone )
c
         endif
c
      endif            !  if ( isub .eq. isub0 .or. isub .gt. 1 )
c
c-----------------------------------------------------------------------
c
      if ( nchuse .gt. 1 ) then
c
c   For cancellation in the un-compressed domain, cancel needs to be
c   called only when the raw data changes.  For compressed data, it
c   needs to be called every time.
c
         if ( isub .eq. isub0 .or.
     .        isub .gt. 1     .or.
     .        ( nchuse .eq. 2 .and. cmodel .ge. 1 ) ) then
c
c              Perform multi-channel clutter cancellation
c
            call cancel ( crc_raw , ntr , ntp , cstr , ncr , nt1 ,
     .                    nchuse , cmodel , nr1 , nre , finemc ,
     .                    curtime , isub0 , isub ) 
c
         endif
c
      endif
c
c   Finally, clear initialization flag after it has been passed to
c   fillbuff (and downward to rawdat and iqgen).
c
      first = .false.  
c
      return
      end
