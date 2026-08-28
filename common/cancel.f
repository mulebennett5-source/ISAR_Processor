C**
C***********************************************************************
C**
      subroutine cancel ( crc_raw , ntr , ntp , cstr , ncr , nt1 ,
     .                    nchuse , cmodel , nr1 , nre , finemc ,
     .                    curtime , isub0 , isub )
C**
C***********************************************************************
C**
c      Purpose:  To create clutter-cancelled radar data for input to
c                an Inverse-SAR processor
c
      implicit none
c
      integer      ntr , ntp , ncr , nt1 , ip , ir , nchuse , ntpwin ,
     .             cmodel , ncancel , nleft , ncenter , nright , nr1 ,
     .             nre , finemc , nrfm , ncrwin , isub0 , isub
c
      complex      crc_raw(ntr,ntp,nchuse+1)   !  Uncompressed raw data
c
      complex      cstr(ntp,ncr,nchuse+1)      !  Range-compressed data
c
      complex      cwt1 , c21 , c23 , c13 , cwt2 , cov2 , c21t ,
     .             c23t , c13t
c
      real         oldwt , newwt , c11 , c22 , c33 , c00 , db ,
     .             cerror , pi , c11t , c22t , c33t , curtime
c
      save               !  Save all estimates so they can be updated
c
      if ( isub .le. 1 ) then                      !  First sub-image
c
         pi      = atan2( 0.0 , - 1.0 )
c
c   Zero running totals
c
         c21t    = cmplx( 0.0 , 0.0 )
         c11t    = 0.0
         c23t    = cmplx( 0.0 , 0.0 )
         c33t    = 0.0
         c22t    = 0.0
         c13t    = cmplx( 0.0 , 0.0 )
c
c   The original data is stored starting at the second level so that
c   the cancelled channel is in the first level.
c
         ncancel = 1        !  Cancelled data used in the ISAR processor
         nleft   = 2        !  'Left' channel
         ncenter = 3        !  'Center' channel
         nright  = 4        !  'Right' channel
c
         newwt   = 1.0
         oldwt   = 0.0
c
         ntpwin  = 3 * ( ntp / 4 )
c
      else
c
         newwt   = 0.25
         oldwt   = 1.0 - newwt
c
         ntpwin  = ntp / 4
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Fundamental choice of number of channels to process
c
      if ( nchuse .eq. 3 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c                      Cancellation Method 1
c
c
c   Perform three-channel cancellation with weights based on the raw
c   data before range compression.  Cancel center channel with left and
c   right channels.
c
c   This method works for the Norden APG-76 Ku-Band radar.
c
c   Covariance and variance of left channel at previous pulse
c
         c21     = cov2( crc_raw(1,1,ncenter) ,
     .                   crc_raw(1,1,nleft) ,
     .                   ntr , ntp , ntr , ntpwin ,
     .                   0 , 1 , 0 )
c
         c11     = cabs( cov2( crc_raw(1,1,nleft) ,
     .                         crc_raw(1,1,nleft) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Covariance and variance of right channel at next pulse
c
         c23     = cov2( crc_raw(1,1,ncenter) ,
     .                   crc_raw(1,1,nright) ,
     .                   ntr , ntp , ntr , ntpwin ,
     .                   0 , - 1 , 0 )
c
         c33     = cabs( cov2( crc_raw(1,1,nright) ,
     .                         crc_raw(1,1,nright) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Variance of main channel
c
         c22     = cabs( cov2( crc_raw(1,1,ncenter) ,
     .                         crc_raw(1,1,ncenter) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Covariance of the above two cancellation channels
c
         c13     = cov2( crc_raw(1,1,nleft) ,
     .                   crc_raw(1,1,nright) ,
     .                   ntr , ntp , ntr , ntpwin ,
     .                   0 , - 2 , 0 )
c
c   Variance of cancel channel (diagnostic)
c
         if ( isub .eq. isub0 ) then
c
            c00     = c22
c
         else
c
            c00     = cabs( cov2( crc_raw(1,1,ncancel) ,
     .                            crc_raw(1,1,ncancel) ,
     .                            ntr , ntp , ntr , ntpwin ,
     .                            0 , 0 , 0 ) )
c
         endif
c
         c21t    = oldwt * c21t + newwt * c21
         c11t    = oldwt * c11t + newwt * c11
         c23t    = oldwt * c23t + newwt * c23
         c33t    = oldwt * c33t + newwt * c33
         c22t    = oldwt * c22t + newwt * c22
         c13t    = oldwt * c13t + newwt * c13
c
c   Update the weights
c
         call canwgts ( cwt1 , cwt2 , c11t , c22t , c33t , c21t , c23t ,
     .                  c13t , cerror )
c
         write ( 83 , 1 )
         write ( 83 , * ) ' 3-Channel cancellation: Uncompressed data'
         write ( 83 , 1 )
         write ( 83 , 1 ) ' Cross correlation:             ' , c21t
c
         write ( 83 , 1 ) ' Secondary variance:            ' , c11t
c
         write ( 83 , 1 ) ' Cross correlation:             ' , c23t
c
         write ( 83 , 1 ) ' Secondary variance:            ' , c33t
c
         write ( 83 , 1 ) ' Covar. of cancel channels:     ' , c13t
c
         write ( 83 , 1 ) ' Complex weight:                ' , cwt1
c
         write ( 83 , 1 ) ' Complex weight:                ' , cwt2
         write ( 83 , 1 )
c
    1    format ( a33 , 4f18.4 )
c
         cerror = amax1( cerror , 1.0e-6 * c22t )
c
c-----------------------------------------------------------------------
c
c   Leave the first two pulses as zeros
c
         do ip = max( 3 , nt1 ) , ntp
c
            do ir = 1 , ntr
c
               crc_raw(ir,ip,ncancel) = crc_raw(ir,ip-1,ncenter)
     .                                + crc_raw(ir,ip-2,nleft)  * cwt1
     .                                + crc_raw(ir,ip,  nright) * cwt2
c
            enddo
c
         enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      else if ( nchuse .eq. 2 .and. cmodel .eq. 0 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c                      Cancellation Method 2
c
c
c   Perform two channel cancellation using weights calculated from raw
c   data before range compression.  Cancel center channel with left
c   channel.
c
c   This method also works for the Norden APG-76, but not as well as the
c   three-channel cancellation above.
c
c   Covariance of left channel at previous pulse
c
         c21     = cov2( crc_raw(1,1,ncenter) ,
     .                   crc_raw(1,1,nleft) ,
     .                   ntr , ntp , ntr , ntpwin ,
     .                   0 , 1 , 0 )
c
c   Variance of left channel
c
         c11     = cabs( cov2( crc_raw(1,1,nleft) ,
     .                         crc_raw(1,1,nleft) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Variance of primary channel (diagnostic)
c
         c22     = cabs( cov2( crc_raw(1,1,ncenter) ,
     .                         crc_raw(1,1,ncenter) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Variance of cancel channel (diagnostic)
c
         if ( isub .eq. isub0 ) then
c
            c00     = c22
c
         else
c
            c00     = cabs( cov2( crc_raw(1,1,ncancel) ,
     .                            crc_raw(1,1,ncancel) ,
     .                            ntr , ntp , ntr , ntpwin ,
     .                            0 , 0 , 0 ) )
c
         endif
c
         c21t    = oldwt * c21t + newwt * c21
         c11t    = oldwt * c11t + newwt * c11
         c22t    = oldwt * c22t + newwt * c22
c
c   Update the weights
c
         cwt1    = - conjg( c21t ) / c11t
c
         cerror  = c22t - c11t * cabs( cwt1 ) ** 2
c
         write ( 83 , 1 )
         write ( 83 , * ) ' 2-Channel cancellation: Uncompressed data'
         write ( 83 , 1 )
         write ( 83 , 1 ) ' Cross correlation (lagged):    ' , c21t
c
         write ( 83 , 1 ) ' Secondary, primary variances:  ' , c11t ,
     .                                                         c22t
c
         write ( 83 , 1 ) ' Complex weight:                ' , cwt1
         write ( 83 , 1 )
c
c-----------------------------------------------------------------------
c
c   Leave the first pulse as zero
c
         do ip = max( 2 , nt1 ) , ntp
c
            do ir = 1 , ntr
c
c   For multiple-channel modes, the output channel is stored in level 1
c
               crc_raw(ir,ip,ncancel) = crc_raw(ir,ip,ncenter)
     .                                + crc_raw(ir,ip-1,nleft) * cwt1
c
            enddo
c
         enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      else if ( nchuse .eq. 2 .and. cmodel .eq. - 1 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c               Cancellation Method 2 - Sum Channel
c
c
c   Perform two channel cancellation using weights calculated from raw
c   data before range compression.
c
c   Covariance of left channel at same pulse
c
         c21     = cov2( crc_raw(1,1,ncenter) ,
     .                   crc_raw(1,1,nleft) ,
     .                   ntr , ntp , ntr , ntpwin ,
     .                   0 , 0 , 0 )
c
c   Variance of left channel
c
         c11     = cabs( cov2( crc_raw(1,1,nleft) ,
     .                         crc_raw(1,1,nleft) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Variance of primary channel (diagnostic)
c
         c22     = cabs( cov2( crc_raw(1,1,ncenter) ,
     .                         crc_raw(1,1,ncenter) ,
     .                         ntr , ntp , ntr , ntpwin ,
     .                         0 , 0 , 0 ) )
c
c   Variance of cancel channel (diagnostic)
c
         if ( isub .eq. isub0 ) then
c
            c00     = c22
c
         else
c
            c00     = cabs( cov2( crc_raw(1,1,ncancel) ,
     .                            crc_raw(1,1,ncancel) ,
     .                            ntr , ntp , ntr , ntpwin ,
     .                            0 , 0 , 0 ) )
c
         endif
c
         c21t    = oldwt * c21t + newwt * c21
         c11t    = oldwt * c11t + newwt * c11
         c22t    = oldwt * c22t + newwt * c22
c
c   Update the weights
c
         cwt1    = - conjg( c21t ) / c11t
c
         cerror  = c22t - c11t * cabs( cwt1 ) ** 2
c
         write ( 83 , 1 )
         write ( 83 , * ) ' 2-Channel cancellation: Uncompressed data'
         write ( 83 , 1 )
         write ( 83 , 1 ) ' Cross correlation (lagged):    ' , c21t
c
         write ( 83 , 1 ) ' Secondary, primary variances:  ' , c11t ,
     .                                                         c22t
c
         write ( 83 , 1 ) ' Complex weight:                ' , cwt1
         write ( 83 , 1 )
c
c-----------------------------------------------------------------------
c
c   Leave the first pulse as zero
c
         do ip = max( 2 , nt1 ) , ntp
c
            do ir = 1 , ntr
c
c   For multiple-channel modes, the output channel is stored in level 1
c
c   To create sum-channel simply change sign from cancellation case.
c   Use no shift in time.
c
               crc_raw(ir,ip,ncancel) = crc_raw(ir,ip,ncenter)
     .                                - crc_raw(ir,ip,nleft) * cwt1
c
            enddo
c
         enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      else if ( nchuse .eq. 2 .and. cmodel .gt. 0 ) then
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c                      Cancellation Method 3
c
c
c   Perform two channel cancellation.  Cancel center channel with left
c   channel in the range-compressed domain.
c
         if ( isub .gt. isub0 )  then
c
c   No range-compressed data on first sub-image yet
c
            if ( finemc .ne. 0 ) then
c
               nrfm = ( nre + 1 - nr1 )   !   Mo-comp window
c
            else
c
               nrfm = ncr / 4
c
            endif
c
            ncrwin = 3 * nrfm
c
c   Covariance between main and cancel channels
c
            c21  = cov2( cstr(1,1,ncenter) , cstr(1,1,nleft) ,
     .                   ntp , ncr , ntpwin , ncrwin , 0 , 0 , nrfm )
c
c   Variance of left channel
c
            c11  = cabs( cov2( cstr(1,1,nleft) , cstr(1,1,nleft) ,
     .                         ntp , ncr , ntpwin , ncrwin , 0 , 0 ,
     .                         nrfm ) )
c 
            if ( cmodel .gt. 2 ) then
c
               c23 = cov2( cstr(1,1,ncenter) , cstr(1,1,nleft) ,
     .                       ntp , ncr , ntpwin , ncrwin , - 1 ,
     .                       0 , nrfm )
c
c   Variance of right channel
c
               c33 = c11
c
c   Covariance of the above two cancellation channels
c
               c13 = cov2( cstr(1,1,nleft) , cstr(1,1,nleft) ,
     .                     ntp , ncr , ntpwin , ncrwin , - 1 , 0 ,
     .                     nrfm )
c
            endif
c
c   Variance of primary channel (diagnostic)
c
            c22  = cabs( cov2( cstr(1,1,ncenter) , cstr(1,1,ncenter) ,
     .                         ntp , ncr , ntpwin , ncrwin , 0 , 0 ,
     .                         nrfm ) )
c
c   Variance of cancel channel to measure effectiveness
c
            c00  = cabs( cov2( cstr(1,1,ncancel) , cstr(1,1,ncancel) ,
     .                         ntp , ncr , ntpwin , ncrwin , 0 , 0 ,
     .                         nrfm ) )
c
            c21t = oldwt * c21t + newwt * c21
            c11t = oldwt * c11t + newwt * c11
            c23t = oldwt * c23t + newwt * c23
            c33t = oldwt * c33t + newwt * c33
            c22t = oldwt * c22t + newwt * c22
            c13t = oldwt * c13t + newwt * c13
c
            write ( 83 , 1 )
            write ( 83 , * )
     .            ' 2-Channel cancellation: Range-compressed data'
            write ( 83 , 1 )
            write ( 83 , 1 ) ' Cross correlation (lagged):   ' , c21t
c
            write ( 83 , 1 ) ' Secondary, primary variances: ' , c11t ,
     .                                               c22t , c11t / c22t
c
            if ( cmodel .gt. 2 ) then
c
               write ( 83 , 1 ) ' Cross correlation:            ' , c23t
c
               write ( 83 , 1 ) ' Covar. of cancel channels:    ' , c13t
c
            endif
c
            if ( cmodel .eq. 1 .or. cmodel .eq. 2 ) then
c
               cwt1   = - conjg( c21t ) / c11t
c
               cerror = c22t - c11t * cabs( cwt1 ) ** 2
c
            else if ( cmodel .eq. 3 ) then
c
               call canwgts ( cwt1 , cwt2 , c11t , c22t , c33t , c21t ,
     .                        c23t , c13t , cerror )
c
            endif
c
            cerror = amax1( cerror , 1.0e-6 * c22t )
c
            write ( 83 , 1 ) ' Complex weight:  #1#            ' ,
     .                         curtime , cwt1
c
            write ( 83 , 1 ) ' Complex weight:  #2#            ' ,
     .                         curtime , cwt2
c
            write ( 83 , * )
c
         endif
c                
c-----------------------------------------------------------------------
c
c   Data correction phase: use information from the compressed pulses
c   to correct the uncompressed data.
c
         do ip = nt1 , ntp
c
            do ir = 1 , ntr
c
c   For multiple-channels mode, the active channel is stored in level 1
c
c   Definition of cmodel:
c
c      cmodel = 0   :  Cancellation coefficients are calculated in the
c                      uncompressed domain in the previous method
c
c      cmodel = 1   :  Coefficients are calculated from range-compressed
c                      data with no auto-shift
c
c      cmodel = 2   :  Same as cmodel = 1, except that in routine
c                      'getraw' the files are lined up to maximize the
c                      correlation coefficient
c
c      cmodel = 3   :  Same as cmodel = 2, except that two adjacent
c                      lags are used
c
               if ( isub .eq. isub0 ) then
c
c   First time just put the center channel into the active position
c
                  crc_raw(ir,ip,ncancel) = crc_raw(ir,ip,ncenter)
c
               else if ( cmodel .eq. 1 .or. cmodel .eq. 2 ) then
c
c   Cancel center channel with left channel lagged one pulse
c
                  crc_raw(ir,ip,ncancel) = crc_raw(ir,ip,ncenter)
     .                                   + crc_raw(ir,ip,nleft) * cwt1
c
               else if ( cmodel .eq. 3 ) then
c
c   Cancel center channel with left channel lagged and at same pulse
c
                  crc_raw(ir,ip,ncancel) = crc_raw(ir,ip-1,ncenter)
     .                                   + crc_raw(ir,ip-1,nleft) * cwt1
     .                                   + crc_raw(ir,ip,nleft)   * cwt2
c
               endif
c
            enddo
c
         enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      endif                      !  Cancellation algorithm options
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      if ( isub .gt. 0 ) then
c
         write ( 82 , '(1x,a5,3f12.4)' ) ' #c# ' , curtime ,
     .                 db( cerror / c22t ) , db( c00 / c22 )
c
         write ( 83 , '(1x,a5,3f12.4)' ) ' #c# ' , curtime ,
     .                 db( cerror / c22t ) , db( c00 / c22 )
c
         write (  6 , '(1x,a33,3f12.4)' )
     .               ' Cancellation efficiency:        ' , curtime ,
     .                 db( cerror / c22t ) , db( c00 / c22 )
c
      endif
c
      return
      end
