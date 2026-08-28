C**
C***********************************************************************
C**
      subroutine movie0 ( ntr , ncr , nsr , ntp , nsa , naskip ,
     .                    nakeep , spulse , npulse , nspf , nspt ,
     .                    nafill , nfr , nfa , rc_raw , crc_raw , rc ,
     .                    crc , str , cstr , ac , cac , csbimg , wtr ,
     .                    wtac , wtaf , ascan , cframe , bytes , ctgt ,
     .                    pgatgt , dotdot , work , iwork , nwork ,
     .                    nabuff , mrrat , mprat , nchuse , mpass ,
     .                    isub , iframe , nafill_a , nfocus )
C**
C***********************************************************************
C**                                                                   **
C**               Special Version for Step-Chirp Data                 **
C**                                                                   **
C***********************************************************************
c
c   This is the main routine for the RDRTec Inverse Synthetic Aperture
c   Radar processor, ISAR-T.  In a normal run of ISAR-T this subroutine
c   is called once; in the multi-pass mode it is called twice.
c
c             mode = 101 : Normal ISAR mode - simply FFT the data at
c                          fixed range and produce a movie following
c                          a moving target.  No adaptive motion
c                          compensation
c
c             mode = 102 : Same as mode = 101 but with adaptive motion
c                          compensation based on sub-image centroid
c                          estimated from the most recent sub-image
c
c   In the ISAR mode, the definition of the three types of calls to the
c   routine is:
c
c             mpass = 0  : Single pass method - don't buffer data
c
c             mpass = 1  : First pass of multi-pass method.  Calculate
c                          the sub-image and write it to disk.
c
c             mpass = 2  : Second pass of multi-pass method.  Do not
c                          calculate sub-image but read it from disk.
c
c   This routine operates on a fundamental cycle of NASKIP pulses of the
c   radar.  Every NASKIP pulses a new coarse doppler sub-image is formed
c   and the decision is made as to whether to form an full resolution
c   image frame.  Usually a full resolution frame is formed every other
c   sub-image until enough data has been collected for the full desired
c   coherent integration time.  At this point a new frame is formed
c   every NFSKIP sub-images.  An exception to this rule is for 'exact
c   time mode' for which NFSKIP=0.  In this mode new image frames are
c   formed at fixed time intervals; the decision to form one at the end
c   of a given cycle is made based on whether a frame is due or not.
c
c***********************************************************************
c
      implicit none
c
      include     'sarprm.h'      !  Standard ISAR-T parameters
c
      include     'updates.h'     !  Updates to parameters from the
c                                    first major release of the code
c
      include     'realtime.h'    !  Real-time parameters
c
      include     'tglist.h'      !  Target parameters and motion
c                                    estimates
c
      include     'motime.h'      !  CPU timing variables - diagnostics
c
c***********************************************************************
c
c        Array sizes and raw data length passed as arguments
c
      integer      ntr , ncr , nsr , ntp , nsa , naskip , nakeep ,
     .             npulse , nspf , nspt , nabuff , nafill , nfr , nfa ,
     .             mrrat , mprat , nchuse , nafill_a
c
c        ntr       :  Range samples in fast time
c        ncr       :  Range cells compressed
c        nsr       :  Range cells per sub-image
c        ntp       :  Pulses in the sub-image buffer
c        nsa       :  Azimuth cells in a sub-image
c        naskip    :  Pulses between sub-images
c        nakeep    :  Sub-image azimuth cells saved
c        npulse    :  Total pulses to process for movie
c        nspf      :  Sub-images per frame
c        nspt      :  Sub-images per target
c        nabuff    :  Number of sub-images saved
c        nafill    :  Total sub-images per frame including zero-fill
c        nafill_a  :  Total sub-images per frame including zero-fill and
c                     augmented by extra sub-images required to compute
c                     two frames at once in the real-time mode
c        nfr       :  Range cells in a frame
c        nfa       :  Azimuth cells in a frame
c        mrrat     :  Number of range cells in the 'depth of focus'
c        mprat     :  Number of doppler cells between 'master particles'
c        nchuse    :  Number of data channels used (usually 1, but can
c                     be 2 or 3 for clutter cancellation modes) 
c        nfocus    :  Number of alternate focus settings
c
c***********************************************************************
c
c                  Real arrays which are equivalenced in the main
c                  program to the following complex arrays:
c
      real      rc_raw(2,ntr,ntp,nchuse+1) !  Raw data buffer
      complex   crc_raw(ntr,ntp,nchuse+1)  !  (Equivalenced complex)
c
      real      rc(2,ncr)                  !  Work array for range
      complex   crc(ncr)                   !  compression
c
      real      str(2,ntp,ncr,nchuse+1)    !  Range-compressed signal
      complex   cstr(ntp,ncr,nchuse+1)     !  history at present
c
c                                          !  Work array for Doppler
      real      ac(2,(nakeep+nsa)*ntr+8*nfa+nafill_a*mrrat)
      complex   cac((nakeep+nsa)*ntr+8*nfa+nafill_a*mrrat) 
c
      complex   ctgt(ncr,13)
c
c***********************************************************************
c
c                  Other variables passed as arguments
c
      integer   spulse , nwork , mpass , maxitu , nfocus
c
      real      wtr(ntr) , wtac(3*naskip) , wtaf(nspf) ,
     .          ascan(nfr+nfa) , work(nwork) , pgatgt(nsr,5) ,
     .          dotdot(1+nfa/mprat,nfr/mrrat,nabuff) , vwt(ncr)
c
      complex   cframe(nfa,nfr,2+2*nfocus)   ! Fine resolution complex
c                                            ! image plus work space for
c                                            ! the PGA
c
      complex   csbimg(nakeep,nsr,nabuff)    ! 3-D buffer of sub-images
c
      integer   iwork(nwork)
c
      character bytes(nfa,nfr)*1
c
c***********************************************************************
c
c                  Local variables
c
      real      tr , bw , rhor , dr0 , vr0 , pt , vadd , dfmp , f0 ,
     .          vrbar , dts , dtf , pi , vrsave , drfbig , ctrast ,
     .          fcmax , taywtc , athrsh , rwindo , vwindo , sigsq ,
     .          tframe , pnoise , tdummy , rcen , fcen ,
     .          astemp , dopcen , dopcen_s , drfreq , lambdap ,
     .          ptb(nbands) , rfc(nbands) , efhz , doff(nbands) ,
     .          ptbuse , dr0use , tinteg_use , sum_frsq
c
      integer   nsubim , ir , ip , ia , isub , iisub , ipulse , ntg ,
     .          nerror , nt1 , iframe , nr1 , nre , nr1fm , isub0 ,
     .          nrefm , actual , pframe , ich , nchtot , ipm1 , nbdiv ,
     .          ntrcom , nrfsmth , iband , irp , ntrband , nftotal ,
     .          rflose , iib , ndop , nspf_use , kk
c
      logical   first , frame , qfirst
c
c                  Variables for Dave Brown's pga routine
c
      integer   icount , nptgt
c
      real      dtpga , db , corg , cnew , c_reduce
c
c                  Variables for getacc Routine
c
      real      ti0 , tv0 , curtime
c
c   Control parameters for virtual sub-images needed for the exact
c   time implementation of the frame generation process
c
      integer   nspfv          !  No. of virtual sub-images for the
                               !  present frame
c
      real      dtv            !  Time separation for virtual
                               !  sub-images
c
      real      tint           !  Local integration time (<= tinteg)
c
      character csfile*80
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c***********************************************************************
c
c                  Start Executable Statements
c
c***********************************************************************
c
c   The raw data file is opened in subroutine RAWDAT on the first call
c   using the input file name and extension '.raw'
c
      first  = .true.                    !  Initialize data files
c
      qfirst = .true.                    !  Initialize quantizer
c
      rcen    = 0.0
      fcen    = 0.0
c
c-----------------------------------------------------------------------
c
c   Pre-compute constants and Taylor weight arrays
c
      dts    = dtp * float( naskip )     !  Time between sub-images
c
      pnoise = 25.0                      !  Noise floor percentile
c
      athrsh = 0.5 / ( dts ** 2 )        !  Acceleration limit
c
      pi     = atan2( 0.0 , - 1.0 )      !  Pi
c
      tr     = dtr * float( ntr )        !  Range integration time
c
      bw     = 1.0E+12 * abs( br ) * tr  !  Bandwidth (Hz)
c
      rhor   = clight / ( 2.0 * bw )     !  Range resolution (m)
c
      drc    = rhor * float( ntr ) /     !  Range sampling ( m )
     .                float( ncr )       !  of range compressed data
c
      efhz   = efghz * 1.0E+9            !  Center frequency in Hz
c
c-----------------------------------------------------------------------
c
c   In ISAR mode the sub-image range cell size is the same as for
c   the range-compressed data
c
      drs    = drc
c
c-----------------------------------------------------------------------
c
c   Coarse Doppler sampling rate (Hz) and minus Nyquist
c
      dfc    = 1.0 / ( float( nsa ) * dtp )
c
      fc0    = - float( nakeep / 2 ) * dfc
c
      fcmax  = fc0 + dfc * float( nakeep - 1 )
c
c-----------------------------------------------------------------------
c
c   Sampling rates in range and Doppler of fine resolution images
c
      drf    = rhor / overrg
      drfbig = drf * float( mrrat )
c
      if ( nfskip .gt. 0 ) then
c
c   In the simplest mode, the frames output time and the integration
c   time are an integer multiple of the sub-image separation time
c
         dff    = 1.0 / ( float( nafill ) * dts )
         dtv    = dts
c
         dtf    = float( nfskip ) * dts
c
      else
c
c   In 'exact-time' mode, the frame rate and the integration time are
c   independent numbers with no restrictions.
c
         dff    = ( 1.0 / tinteg ) * float( nspf ) / float( nafill )
         dtv    = tinteg / float( nspf )
c
         dtf    = 1.0 / frrate
         tframe = 2.0 * dts            !  Fine resolution frame time
c
      endif
c
c-----------------------------------------------------------------------
c
c   A master particle is an image point for which the integration is
c   done very accurately by summing the interpolated sub-images.  In
c   between these points an FFT of the master particle time series is
c   used to interpolate the image values.
c
      dfmp   = dff * float( mprat )   !  Master particle frequency
c                                        separation
c-----------------------------------------------------------------------
c
c   Sampling rate for to pga routine.  This is the time step for the
c   effective signal history obtained by Fourier transforming the
c   complex image.
c
      dtpga  = 1.0 / ( float( nfa ) * dff )
c
c-----------------------------------------------------------------------
c
c   Min and max range of sub-images and fine resolution images
c
      rsmin  = - float( nsr / 2 ) * drs
      rsmax  = rsmin + float( nsr - 1 ) * drs
      rfmin  = - float( nfr / 2 ) * drf
      rfmax  = rfmin + float( nfr - 1 ) * drf
c
c   Fine resolution window in range and velocity - used to thin the
c   targets in the front-end coherent detector (subroutine target)
c
      rwindo = float( nsr ) * drs
      vwindo = float( nakeep ) * dfc * 0.5 * lambda
c
c   Display some of this information to the user and log file
c
      call movie2 ( drs , rsmin , rsmax , drf , rfmin , rfmax ,
     .              mprat , dfc , dff , dfmp )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Compute Taylor weight array for range
c
      taywtr  = - abs( taywtr )
c
      wtr(:)  = 0.0
c
c   Number of fast-time samples per RF band
c
      ntrband = ntr / nbands
c
      if ( mod( ntr , nbands ) .ne. 0 ) stop 'ntr, nbands inconsistent'
c
      if ( mod( ntrband , 2 ) .ne. 0  ) stop 'ntr/nbands not even'
c
c   Size of each RF pixel (Hz)
c
      drfreq  = br * 1.0E+12 * dtr
c
c   Find the total number of bands between the centers of the first and
c   last bands after they are interpolated to an integer number of bands
c
      call ntrcomb ( nftotal , nbdiv , nbands , rfdelt , drfreq )
c
c   Total number of active RF pixels in the overlapping bands
c
      ntrcom = nftotal + ntrband
c
c   Trim the RF bands to allow for a few to be set to zero during the
c   interpolation process
c
      rflose = nbands
c
c   Use a Taylor weight over the expected full-amplitude region
c
      call taylor ( ntrcom - 2 * rflose , taywtr , wtr(1+rflose) )
c
c-----------------------------------------------------------------------
c
c   The bands are assumed to be fired in the order they are stored in
c   the data file and are equally spaced in time over the full PRI
c
      do iband = 1 , nbands
c
c   Offset to center frequencies for the bands in the original data
c
         rfc(iband)  = rfdelt * ( float( iband - 1 - nbands / 2 )
     .                    - 0.5 * float( 1 - mod( nbands , 2 ) ) )
c
c   Time offset to each band
c
         doff(iband) = ( dtp / float( nbands ) ) *
     .                 ( float( iband - 1 - nbands / 2 )
     .                   - 0.5 * float( 1 - mod( nbands , 2 ) ) )
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Compute weight array for coarse Doppler
c
c   For modes which are not real-time, hardwire the Taylor weight used
c   in sub-image formation since it is matched to the overlap factor of
c   2/3
c
      taywtc = - 40.0
c
      call taylor ( 3 * naskip , taywtc , wtac )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   There are nsubim possible sub-images; the index 0 is used for a
c   first pass at range compression for the purpose of determining
c   the average velocity.  Afterwards the previous velocity is used
c   to approximately correct for the average velocity.
c
      nsubim = ( iabs( npulse ) - ntp ) / naskip
c
      iframe = 0                       !  Fine resolution frame counter
c
      vr0    = 0.0                     !  Mean radial velocity
c
      vrsave = 0.0                     !  vr0 before target detection
c
      vrbar  = 0.0                     !  Target mean velocity
c
      dr0    = 0.0                     !  Offset of object from MCP
c
      ntg    = 0                       !  Targets being tracked
c
      if ( notch .eq. 1 .or. notch .eq. 3 ) vr0 = vnotch
c
      if ( notch .eq. 2 .or. notch .eq. 3 ) dr0 = rnotch
c
      dopcen   = 0.0
c
      dopcen_s = 0.0
c
c   Slip velocity correction passed to IMGEN
c
      vadd     = 0.0
c
c   Compute number of channels to compress.  Only compress the first
c   channel unless there are multiple channels and we are doing an
c   advanced algorithm which requires range-compressed data.
c
      if ( nchuse .gt. 1 .and. cmodel .gt. 0 ) then
c
         nchtot = nchuse + 1
c
      else
c
         nchtot = 1
c
      endif
c
c                       END OF INITIALIZATION
c
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
c
c                       START OF CALCULATIONS
c
c   The fundamental cycle of the routine is the sub-image time interval
c   which is equal to NASKIP radar pulses.  However, there is a ramp-up
c   number of cycles allowed for the motion compensation and target
c   acquisition algorithms to work before any image frames are output.
c
c   Start the sub-image counter at a negative value to allow cycles to
c   resolve the aliases.  The negative and zero sub-images physically
c   cover the same time period as the first one but the data may be
c   re-compensated based on an iterative estimate of the target's
c   location and radial velocity.
c
      isub0    = - 1
c
      do 1000 isub = isub0 , nsubim
c
c   Treat sub-image array as a circular buffer with pointer = iisub
c
         if ( isub .le. 1 ) then
c
            iisub = 1             !  Multiple passes for first sub-image
c
         else
c
            iisub = iisub + 1     !  Afterwards, increment pointer
c
         endif
c
c   Make pointer circular - repeating ( 1 --> nabuff )
c
         if ( iisub .gt. nabuff ) iisub = iisub - nabuff
c
c***********************************************************************
c
         if ( mpass .lt. 2 ) then
c
c                      Range compress raw data
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
c   Compute mid-point time of current target buffer
c
            curtime = float( 2 + max( isub - 1 , 0 ) ) * dts - 0.5 * dts
c
            ohtime = ohtime + dtime( timer )
c
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c----------------------  Read in raw data  -----------------------------
c
            call getraw ( ntr , ntp , isub0 , isub , nt1 , naskip ,
     .                    spulse , npulse , ncr , rc_raw , crc_raw ,
     .                    cstr , nerror , work , nwork , first ,
     .                    nchuse , nr1 , nre , curtime )
c
            if ( nerror .ne. 0 ) go to 5000
c
            iqtime = iqtime + dtime( timer )
c
c----------------------  Read in raw data  -----------------------------
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c--------------  Range Compression and Fine Mo-Comp  -------------------
c
            do ip = nt1 , ntp
c
c   Number the pulses according to the order they were read in, ignoring
c   the fact that the time order could have been reversed above.
c
               ipulse = naskip * max( isub - 1 , 0 ) + ip +
     .                  ( spulse - 1 )
c
c-----------------------------------------------------------------------
c
c                  First Pulse Initialization
c
               if ( ipulse .eq. spulse ) then
c
c   Initialize the first range-compressed pulse to zero to start the
c   pre-summer
c
                  do ir = 1 , ncr
c
                     do ich = 1 , nchtot
c
                        cstr(1,ir,ich) = cmplx( 0.0 , 0.0 )
c
                     enddo
c
                  enddo
c
c   Set the phase estimate and number of targets to zero
c
                  pt     = 0.0
                  ptb(:) = 0.0
                  ntg    = 0
c
               endif
c 
c-----------------------------------------------------------------------
c
c   Update the mo-comp position taking into account both the normal
c   mo-comp velocity, vr0.
c
               dr0 = dr0 + ( vr0 + v_corr ) * dtp
c
               if ( alias .eq. - 5 ) then
c
                  read  ( 87 , '(3f16.6)' )
     .                    tdummy , dr0 , vr0
c
               else
c
                  write ( 87 , '(3f16.6)' )
     .                    float( ipulse - spulse ) * dtp , dr0 , vr0
c
               endif
c
               if ( nbands .eq. 1 ) then
c
                  f0 = vr0 * ( 2.0 / lambda )  !  Mean Doppler frequency
c
                  pt = pt + f0 * dtp
c
c   Compress the pulse for all channels
c
                  do ich = 1 , nchtot
c
                     crc(:) = cmplx( 0.0 , 0.0 )
c
                     do ir = 1 , ntr
c
                        crc(ir) = crc_raw(ir,ip,ich)
c
                     enddo
c
                     call rcomp ( br , clight , dr0 , rc , crc , pt ,
     .                            ntr , ncr , wtr , dtr , mode , work ,
     .                            nwork )
c
c   Finally, load the compressed pulse into the range-compressed buffer
c   array.
c
                     ipm1 = max( ip - 1 , 1 )
c
                     do ir = 1 , ncr
c
                        cstr(ip,ir,ich) = presum * cstr(ipm1,ir,ich) +
     .                                    ( 1.0 - presum ) * crc(ir)
c
                     enddo
c
                  enddo  !  Loop over channels
c
               else
c
                  do iband = 1 , nbands
c
                     iib        = ntrband * ( iband - 1 )
c
                     lambdap    = clight / ( efhz + rfc( iband ) )
c
                     f0         = vr0 * ( 2.0 / lambdap )
c
                     ptb(iband) = ptb(iband) + f0 * dtp
c
                     ptbuse     = ptb(iband) + f0 * doff(iband)
c
                     dr0use     = dr0 + vr0 * doff(iband)
c
c   Compress the pulse for all channels
c
                     do ich = 1 , nchtot
c
                        crc(:) = cmplx( 0.0 , 0.0 )
c
                        do ir = 1 , ntrband
c
                           irp     = ir + iib
c
                           crc(ir) = crc_raw(irp,ip,ich)
c
                        enddo
c
                        call rcomp ( br , clight , dr0use , rc , crc ,
     .                               ptbuse , ntrband , ncr , wtr ,
     .                               dtr , mode , work , nwork )
c
                        if ( qantiz .ne. 0 )
     .                     call quant ( qfirst , iqlsb , rc , 2 * ncr )
c
                        do ir = 1 , ntrband
c
                           irp              = ir + iib
c
                           cstr(ip,irp,ich) = crc(ir)
c
                        enddo
c
                     enddo    !  Loop over channels
c
                  enddo       !  Loop over bands
c
               endif
c
c-----------------------------------------------------------------------
c
c                   Adaptive Motion Compensation
c
               if ( mode .gt. 102 .and. finemc .ne. 0 ) then
c
                  call mocomp ( vr0 , ctgt , ncr , ip , isub , isub0 ,
     .                          naskip , nfa , nr1 , nre , nr1fm ,
     .                          nrefm , dr0 , vrsave , ntp , sigsq ,
     .                          ipulse , spulse , cstr , vwt )
c
               endif
c
            enddo           !  do ip = nt1 , ntp    ( Loop over pulses )
c
            rctime = rctime + dtime( timer )
c
         endif     !  ( mpass .lt. 2 ) - Range-compress raw data
c
c--------------  Range Compression and Fine Mo-Comp  -------------------
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c-------------------  Sub-image Compression  ---------------------------
c
c   For zero-th sub-image, skip to the end to repeat range compression
c   using the velocity and range offset estimates from the dominant
c   targets
c
c   Non-Real-time version of sub-image compression
c
         if ( rt_img .eq. 0 .and. mode .ne. 0 .and. isub .ge. 1 ) then
c
            nrfsmth = 0
c
            ndop    = max( nsa , nakeep )
c
            call subimg0 ( mpass , isub , csbimg(1,1,iisub) , nsa ,
     .                     nsr , ntp , ncr , naskip , nakeep , str ,
     .                     cstr , wtac , work , ac , cac , nwork ,
     .                     nerror , wtr , ntr , rc , crc , curtime ,
     .                     dopcen , ntrcom , nrfsmth , ndop )
c
            if ( mode .eq. 102 .and. finemc .ne. 0 )
     .           vr0 = vr0 + 0.5 * alpha0 * dopcen * 0.5 * lambda
c
            if ( nerror .ne. 0 ) go to 5000
c
            sitime = sitime + dtime( timer )
c
         endif
c
c   Diagnostic calculation - rarely used
c
c        if ( mode .eq. 9 ) then
c
c   Make a movie of the sub-images
c
c           csfile = 'subimgs'
c
c           call oframe ( csbimg(1,1,iisub) , work , bytes , - isub ,
c    .                    csfile , nsr , nakeep , 0 , 0 , dbinc ,
c    .                    nlocal , pixbar , color , 0 , 0.0 , 0 , 0 )
c
c        endif
c
c-------------------  Sub-image Compression  ---------------------------
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c--------------  Frame Integration and Autofocus  ----------------------
c
c   If in standard mode (nfskip>0), form a fine resolution frame every
c   other sub-image until the sub-image buffer is full.  Afterwards,
c   form one every nfskip subimages.
c
c   If in 'exact time' mode (nfskip=0), then check every sub-image to
c   see if enough time has elapsed to form another frame.
c
c   Compute the number of sub-images to process for a fine resolution
c   frame
c
c   To summarize the integration time variables:
c
c   Total number of sub-images in the buffer
c
         actual  = min( isub , nabuff )
c
c   Times for last and first sub-images in the buffer
c
         curtime = ( 2.0 + float( isub - 1 ) ) * dts
         ti0     = curtime - float( actual - 1 ) * dts
c
c-----------------------------------------------------------------------
c
         if ( nfskip .gt. 0 ) then
c
c   Standard mode - create a fine resolution frame every other sub-image
c   until there are nspf sub-images.  Afterwards, create a frame every
c   nfskip sub-images
c
            frame = ( isub .ge. nspf .and. mod( isub , nfskip ) .eq. 0 )
     .         .or. ( isub .le. nspf. and. mod( isub , 2 ) .eq. 0 .and.
     .                isub .ge. 2 )
c
            if ( frame ) then
c
               tframe     = 0.5 * ( curtime + ti0 )
               nspfv      = min( actual , nspf )
               tint       = float( nspfv ) * dts
               tinteg_use = tint
               nspf_use   = nspfv
c
            endif
c
         else
c
c-----------------------------------------------------------------------
c
c   In exact time mode, form another frame when the center time of the
c   sub-image stack is within half of its time increment of the desired
c   time
c
            frame = ( ( 0.5 * ( curtime + ti0 ) .ge.
     .                tframe + dtf - 0.5 * dts ) )
     .                 .and. ( isub .ge. 2 )
c
c   If nfskip<=0, then form a frame at frrate per second.  If nfskip=0
c   then the integration time is fixed at tinteg seconds.  However, if
c   nfskip<0, choose the integration time based on the target's rotation
c   rate.
c
            if ( frame ) then
c

               dff    = ( 1.0 / tinteg ) * float( nspf ) /
     .                                     float( nafill )
c
               dtpga  = 1.0 / ( float( nfa ) * dff )
c
               dtv    = tinteg / float( nspf )
c
               dfmp   = dff * float( mprat )
c
c   Update the frame time by 1/frrate
c
               tframe = tframe + dtf
c
c   Estimate approximate integration time, then convert it to an even
c   number of virtual sub-images and then recalculate it
c
               if ( mode .eq. 101 ) then
c
                  tint = amin1( 2.0 * ( tframe - 2.0 * dts ) ,
     .                          tinteg )
c
               else
c
                  tint = tinteg
c
               endif
c
               nspfv  = min( nspf , nint( tint / dtv ) )
               nspfv  = max( 2 , nspfv - mod( nspfv , 2 ) )
c
               tint   = float( nspfv ) * dtv
c
            endif
c
         endif
c
         if ( mode .ne. 0 .and. frame ) then
c
c   Form fine resolution ISAR movie frame from sub-images
c
            iframe = iframe + 1
c
c***********************************************************************
c-----------------------------------------------------------------------
c----------------  Coherent Frame Integration  -------------------------
c
            IF ( ( quiet              .gt. 1 ) .or.
     .           ( mod( iframe , 20 ) .eq. 0 ) ) THEN
c
            write ( 6 , * )
            write ( 6 , * ) '------------------------------------------'
            write ( 6 , '(a24,3i6)' ) ' Pass, Sub-image, Frame:' ,
     .                                  mpass , isub , iframe
            write ( 6 , * ) '------------------------------------------'
c
            ENDIF
c
            write ( 7 , * )
            write ( 7 , * ) '------------------------------------------'
            write ( 7 , '(a24,3i6)' ) ' Pass, Sub-image, Frame:' ,
     .                                  mpass , isub , iframe
            write ( 7 , * ) '------------------------------------------'
c
c-----------------------------------------------------------------------
c
c   Don't do centering for first frames.  Afterwards, correct slowly.
c   The parameter 'center' should be greater than or equal to zero.
c
            if ( ( iframe .lt. nspf / 2 ) .or. ( rt_img .ne. 0 ) ) then
c
               rcenuse = 0.0
               fcenuse = 0.0
c
            else
c
               rcenuse = rcenuse + center * ( rcen - rcenuse )
               fcenuse = fcenuse + center * ( fcen - fcenuse )
c
            endif
c
c***********************************************************************
c***********************************************************************
c
c   Get estimate of the acceleration field
c
            tv0 = tframe - 0.5 * dtv * float( nspfv - 1 )
c
c-----------------------------------------------------------------------
c
c   Estimate target reports from the sub-image field
c
            if ( tgt_si .gt. 0 .and. mode .ne. 10 )               
     .         call gettgt_s ( csbimg , nakeep , nsr , nabuff ,
     .                         pnoise , pgatgt , iwork , dts ,
     .                         iisub , isub , 0.0 , work , nptgt ,
     .                         work , curtime , 0.0 , 0.0 )
c
c   Estimate the acceleration array
c
            call getacc   ( 1 + nfa / mprat , nfr / mrrat ,
     .                      nspfv , nspt , dotdot ,
     .                      fcenuse - 0.5 * float( nfa ) * dff ,
     .                      rcenuse + rfmin + 0.5 * drf * mrrat ,
     .                      tv0 , dfmp , drfbig , dtv , athrsh ,
     .                      rcen , fcen , iframe )
c
c   Done with acceleration field
c
c***********************************************************************
c***********************************************************************
c
c   Generate Taylor weights for fine resolution image
c
            taywta = - abs( taywta )
c
            call taylor ( nspfv , taywta , wtaf )
c
c   Generate complex image
c
            call imgen    ( cframe , nfa , nfr , dff , drf ,
     .                      csbimg , nakeep , nsa , nsr , dfc ,
     .                      drs , fc0 , actual , nabuff , nafill ,
     .                      iisub , dts , ac , cac , mprat ,
     .                      wtaf , work , nwork , br , lambda ,
     .                      clight , dotdot , mrrat , rgwalk ,
     .                      vadd , nspfv , dtv , ti0 , tframe ,
     .                      iwork , strtch , rcenuse , fcenuse ,
     .                      reduce )
c
            frtime = frtime + dtime( timer )
c
c----------------  Coherent Frame Integration  -------------------------
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c-----------------  Phase Gradient Autofocus ---------------------------
c
            if ( mpass .eq. 0 .or. mpass .eq. 2 ) then
c
c   This is the last pass - autofocus the frame before output if
c   required
c
               maxitu = abs( maxit )
c
            else
c
c   This is not the last pass before output - don't bother to autofocus
c
               maxitu = 0
c
            endif
c
c   Do range autofocus
c
            if ( ( maxit .lt. 0 ) .and. ( nfa .eq. nfr ) ) then
c
               cframe(:,:,2) = cframe(:,:,1)
c
               do ir = 1 , nfr
c
                  do ia = 1 , nfa
c
                     cframe(ia,ir,1) = cframe(ir,ia,2)
c
                  enddo
c
               enddo
c
               call pga ( cframe , nfa , nfr , iwork , iwork(nfr+1) ,
     .                    work , work(nfa+1) , cac , ac , cac(nfa+1) ,
     .                    work(2*nfa+1) , work(4*nfa+1) ,
     .                    work(5*nfa+1) , work(6*nfa+1) , icount ,
     .                    dtpga , corg , cnew , work(7*nfa+1) , nptgt ,
     .                    pnoise , pgatgt , work(8*nfa+1) , dotdot ,
     .                    mprat , mrrat , nafill , nspfv , dfmp ,
     .                    tframe , rcenuse , fcenuse , 0 )
c
               do ir = 1 , nfr
c
                  do ia = 1 , nfa
c
                     cframe(ia,ir,2) = cframe(ir,ia,1)
c
                  enddo
c
               enddo
c
               cframe(:,:,1) = cframe(:,:,2)
c
               pgtime = pgtime + dtime( timer )
c
            endif
c
            if ( ( maxitu .ne. 0 .or. curve .ne. 0 ) .and.
     .           ( rt_img .eq. 0 .and. rt_pga .eq. 0 ) ) then
c
c   Call PGA if needed for autofocus for this frame or for target
c   reports for GETACC to use.  Not called in real-time mode.
c
               call pga ( cframe , nfa , nfr , iwork , iwork(nfr+1) ,
     .                    work , work(nfa+1) , cac , ac , cac(nfa+1) ,
     .                    work(2*nfa+1) , work(4*nfa+1) ,
     .                    work(5*nfa+1) , work(6*nfa+1) ,icount ,
     .                    dtpga , corg , cnew , work(7*nfa+1) , nptgt ,
     .                    pnoise , pgatgt , work(8*nfa+1) , dotdot ,
     .                    mprat , mrrat , nafill , nspfv , dfmp ,
     .                    tframe , rcenuse , fcenuse , nfocus )
c
               pgtime = pgtime + dtime( timer )
c
            else
c
               cnew   = ctrast( cframe(1,1,1+nfocus) , nfa , nfr )
c
               corg   = cnew
c
            endif 
c
            if ( mod(mode,100) .eq. 1 .and. reduce .gt. 1 ) then
c
               call reduce_streaks ( cframe , nfa , nfr , 1.0 , 0.25  )
c
               c_reduce = ctrast( cframe , nfa , nfr )
c
            else
c
               c_reduce = cnew
c
            endif
c            
            write ( 82 , '(1x,5f12.4)' )  tframe , db( corg ) ,
     .                                    db( cnew ) , db( c_reduce ) ,
     .                                    vadd
c
            if ( quiet .gt. 1 ) write (  6 , '(1x,a26,5f12.4)' )
     .                          ' Autofocus efficiency:    ' , tframe ,
     .                            db( corg ) , db( cnew ) ,
     .                            db( c_reduce ) , vadd
c
c-----------------  Phase Gradient Autofocus ---------------------------
c-----------------------------------------------------------------------
c***********************************************************************
c
c   Output complex image to disk and also write an eight-bit dB
c   scaled version of it for display
c
c-----------------------------------------------------------------------
c
            if ( mpass .eq. 0 .or. mpass .eq. 2 ) then
c
               if ( otype .eq. 0 ) then
c
                  pframe = + iframe
c
               else
c
                  pframe = - iframe
c
               endif
c
               call oframe ( cframe , work , bytes , pframe , ofile ,
     .                       nfa , nfr * ( 1 + 2 * nfocus ) , 0 , 0 ,
     .                       dbinc , nlocal , pixbar , color , imgslk ,
     .                       taywta , rlook , flook )
c
               write ( 7 , '(/,a,i6,2(a,f10.3),/)' )
     .                 ' Contrast for frame:' , iframe , ' = ' , cnew ,
     .                 ' at time: ' , tframe
c
c   Compute and output the A-Scan vector
c
               do ir = 1 , nfr + nfa
c
                  ascan(ir) = 0.0
c
               enddo
c
               do ir = 1 , nfr
c
                  do ia = 1 , nfa
c
                     astemp =        cframe(ia,ir,1+nfocus) *
     .                        conjg( cframe(ia,ir,1+nfocus) )
c
                     ascan(ir)     = ascan(ir)     + astemp
c
                     ascan(nfr+ia) = ascan(nfr+ia) + astemp
c
                  enddo
c
               enddo
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
               write ( 89 , rec = iframe ) ascan
c
               write ( 92 , '(4f14.2,4i9,4f14.2)' )
     .               rcenuse  , fcenuse , tframe , tinteg_use , nspf ,
     .               nspf_use , iframe , ais(1) , ais(2) , ais(3) , cnew
c
            endif       !  End of frame generation process
c
         endif
c
c--------------  Frame Integration and Autofocus  ----------------------
c-----------------------------------------------------------------------
c***********************************************************************
c***********************************************************************
c
         ohtime = ohtime + dtime( timer )
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
 1000 continue       !  End of loop over sub-images
c
c
c                       END OF CALCULATIONS
c
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
 5000 continue       !  File I/O error processing
c
c   Check to see if there was a data error before a frame was computed.
c   If so, then set pass to - 1 so that the second pass is not done.
c
      if ( nerror .ne. 0 .and. iframe .eq. 0 ) mpass = - 1
c
      return
      end
