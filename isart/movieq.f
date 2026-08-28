C**
C***********************************************************************
C**
      subroutine movie  ( ntr , ncr , nsr , ntp , nsa , naskip ,
     .                    nakeep , spulse , npulse , nspf , nafill ,
     .                    nfr , nfa , rc_raw , crc_raw , rc , crc ,
     .                    str , cstr , ac , cac , csbimg , wtr , wtac ,
     .                    wtaf , ascan , cframe , bytes , tgt , ctgt ,
     .                    pgatgt , dotdot , work , iwork , nwork ,
     .                    nabuff , mrrat , mprat , nchuse , pass ,
     .                    isub , iframe , nafill_a , nfocus )
C**
C***********************************************************************
C**
c   This is the main routine for the SAIC Inverse Synthetic Aperture
c   Radar processor, ISAR-T.  In a normal run of ISAR-T this subroutine
c   is called once; in the multi-pass mode it is called twice.
c
c   There are two main types of 'movies' which can be produced with this
c   routine.  Both modes use a two-stage time integration to form the
c   SAR image frames.  In the first stage, a series of coarse doppler
c   images is formed by the following algorithms:
c
c             mode = 1  :  Normal ISAR mode - simply FFT the data at
c                          fixed range and produce a movie following
c                          a moving target
c
c             mode = 2  :  Ultra-Wide Band strip-map mode - do an
c                          accurate direct integration on a linear path
c                          in range versus time and then produce a
c                          strip-map image as a disk file of indefinite
c                          length
c
c             mode = 3  :  Telephonics Stripmap - no real-time mods
c
c             mode = 4  :  Telephonics Stripmap - using Norden real-time 
c                          5-second mode.
c
c             mode = 5  :  Telephonics Stripmap - no real-time mods with
c                          variable platform speed, mo-comp point speed,
c                          and range.  Intended for data analysis.
c
c             mode = 6  :  Generalized stripmap - compute a series of
c                          spotlight patches and lay them down at the
c                          speed of the beam footprint, focusing the
c                          pulses based on the speeds of the platform
c                          and the beam.  This algorithm assumes that
c                          the aim point is always exactly at broadside.
c
c             mode = 7  :  Spotlight
c
c             mode = 8  :  Spotlight with speed-up
c
c   The primary use of this routine is for the ISAR mode.  However, the
c   two modes are very similar in operation.  The only differences for
c   the UWB strip-map mode are:
c
c             1.  The sub-image array is calculated at a range sampling
c                 of twice that of ISAR mode
c
c             2.  The sub-image array uses a different algorithm which
c                 takes into account the movement through range cells
c
c             3.  The image is formed by a variation on the main
c                 image generation routine, IMGEN, called IMGENU.
c
c             4.  Instead of producing a range-doppler image at each
c                 cycle the UWB algorithm produces a section of a
c                 strip-map image which is as wide in cross-range as
c                 the platform displacement during the time between
c                 frames.
c 
c   In the ISAR mode, the definition of the three types of calls to the
c   routine is:
c
c             pass = 0  :  Single pass method - don't buffer data
c
c             pass = 1  :  First pass of multi-pass method.  Calculate
c                          the sub-image and write it to disk.
c
c             pass = 2  :  Second pass of multi-pass method.  Do not
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
      include     'kalman.h'      !  Kalman Filtering Variables
c
c***********************************************************************
c
c        Array sizes and raw data length passed as arguments
c
      integer      ntr , ncr , nsr , ntp , nsa , naskip , nakeep ,
     .             npulse , nspf , nabuff , nafill , nfr , nfa , mrrat ,
     .             mprat , nchuse , nafill_a
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
      real        rc_raw(2,ntr,ntp,nchuse+1) !  Raw data buffer
      complex     crc_raw(ntr,ntp,nchuse+1)  !  (Equivalenced complex)
c
      real        rc(2,ncr)                  !  Work array for range
      complex     crc(ncr)                   !  compression
c
      real        str(2,ntp,ncr,nchuse+1)    !  Range-compressed signal
      complex     cstr(ntp,ncr,nchuse+1)     !  history at present
c
c                                            !  Work array for Doppler
      real        ac(2,nsa*ntr+8*nfa+nafill_a*mrrat)
      complex     cac(nsa*ntr+8*nfa+nafill_a*mrrat) 
c
      real        tgt(ncr,26)                !  Work for TARGET
      complex     ctgt(ncr,13)
c
c***********************************************************************
c
c                  Other variables passed as arguments
c
      integer     spulse , nwork , pass , maxitu , nfocus
c
      real        wtr(ntr) , wtac(3*naskip) , wtaf(nspf) ,
     .            ascan(nfr+nfa) , work(nwork) , pgatgt(nsr,5) ,
     .            dotdot(1+nfa/mprat,nfr/mrrat,nabuff) 
c
      complex     cframe(nfa,nfr,2+2*nfocus) ! Fine resolution complex
c                                            ! image plus work space for
c                                            ! the PGA
c
      complex     csbimg(nakeep,nsr,nabuff) ! 3-D buffer of sub-images
c
      integer     iwork(nwork)
c
      character   bytes(nfa,nfr)*1
c
      character   csfile*80
c
c***********************************************************************
c
c                  Local variables
c
      real        tr , bw , rhor , dr0 , vr0 , pt , vadd , dfmp , f0 ,
     .            vrbar , r0bar , r0min , r0max , dts , dtf , pi ,
     .            vrsave , rfmmin , rfmmax , drfbig , ctrast , fcmax ,
     .            taywtc , athrsh , rwindo , vwindo , sigsq , vslip ,
     .            eslip , tframe , pnoise , snoise , noise0 , afocus ,
     .            tgrbar , tgrmin , tgrmax , tdummy , rcen , fcen ,
     .            omega0 , x , df_subimg , const , phase , asmin ,
     .            asmax , asbar , assdv , astemp , flight_path ,
     .            ground_path , image_path , mtrs2hz , vspot0 , dop ,
     .            vspot1 , vplat0 , vplat1 , omeg_0 , omeg_1 , dopcen
c
      integer     nsubim , ir , ip , ia , isub , iisub , ipulse ,
     .            ntg , nerror , nt1 , iframe , nr1 , nre , nalias ,
     .            nr1fm , isub0 , nrefm , actual , notchp , iline ,
     .            ne , ounit , pframe , ich , nchtot , ipm1 , mcd_idx ,
     .            nspf_use , itemp , rt_frames , nn(1) , nxmid ,
     .            nxtotal , gpulse , gpulse0 , gpulse1
c
      complex     cpfast
c
      logical     first , frame , qfirst
c
c                 Variables for Dave Brown's pga routine
c
      integer     icount , nptgt
c
      real        dtpga , db , corg , cnew 
c
c                 Variables for getacc Routine
c
      real        ti0 , tv0 , curtime , tinteg_use
c
c   Set maximum number of aliases to search in the acquisition mode -
c   should be an odd number so that the same number of plus and minus
c   aliases are considered
c
      integer     maxal
      parameter ( maxal = 5 )    !  5 = One main alias, 2 plus, 2 minus
c
c   Control parameters for virtual sub-images needed for the exact
c   time implementation of the frame generation process
c
      integer     nspfv          !  No. of virtual sub-images for the
                                 !  present frame
c
      real        dtv            !  Time separation for virtual
                                 !  sub-images
c
      real        tint           !  Local integration time (<= tinteg)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Program parameters not important enough to export to the user
c
      logical     smooth_az_ascan / .true. /
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
      first   = .true.                    !  Initialize data files
c
      qfirst  = .true.                    !  Initialize quantizer
c
      mcd_idx = 0                         !  Index for mo-comp file
c
c-----------------------------------------------------------------------
c
c   Pre-compute constants and Taylor weight arrays
c
      dts     = dtp * float( naskip )     !  Time between sub-images
c
      pnoise  = 25.0                      !  Noise floor percentile
c
      athrsh  = 0.5 / ( dts ** 2 )        !  Acceleration limit
c
      pi      = atan2( 0.0 , - 1.0 )      !  Pi
c
      tr      = dtr * float( ntr )        !  Range integration time
c
      bw      = 1.0E+12 * abs( br ) * tr  !  Bandwidth (Hz)
c
      rhor    = clight / ( 2.0 * bw )     !  Range resolution (m)
c
      drc     = rhor * float( ntr ) /     !  Range sampling ( m )
     .                 float( ncr )       !  of range compressed data
c
c-----------------------------------------------------------------------
c
      if ( mode .eq. 1 .or. mode .eq. 4 .or. mode .eq. 6 .or.
     .     mode .eq. 7 .or. mode .eq. 8 .or. mode .eq. 10 ) then
c
c   In ISAR mode the sub-image range cell size is the same as for
c   the range-compressed data
c
         drs = drc
c
         flight_path = 0.0
c
         ground_path = 0.0
c
         image_path  = 0.0
c
         nxtotal     = 0
c
      else
c
c   In Ultra Wide Band strip-map mode the sub-image range cell size is
c   twice that of the range-compressed data
c
         drs = 2.0 * drc
c
      endif
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
c   Set noise values for Kalman information filter
c
      pnois1 = 0.001 * athrsh ** 2
c
      pnoisr = 0.001 * ( athrsh / ( float( nfr ) * drf ) ) ** 2
c
      pnoisf = 0.025 * ( athrsh / ( float( nfa ) * dff ) ) ** 2
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
      taywtr = - abs( taywtr )
c
      call taylor ( ntr , taywtr , wtr )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Compute weight array for coarse Doppler
c
      if ( rt_img .gt. 0 ) then
c
c   In Telephonics real-time mode fill the array with Bill Caputi's
c   weight function
c
         wtac(:) = 0.0
c
         if ( mod( rt_img , 2 ) .eq. 0 ) then
c
            call caputi ( 2 * naskip , wtac(1+naskip/2) )
c
         else
c
c   In SAIC real-time mode fill the array with a cosine-on-a-pedestal
c   weight function
c
            do ip = 1 + naskip / 2 , 3 * naskip - naskip / 2
c
               x        =  ( float( ip - 3 * ( naskip / 2 ) ) - 0.5 ) /
     .                              float( 2 * naskip )
c
               wtac(ip) = 0.54 + 0.46 * cos( 2.0 * pi * x )
c
            enddo
c
            df_subimg = 1.0 / ( float( 2 * naskip ) * dtp )
c
c   Determine correction for roll-off of filter response versus Doppler
c   frequency - this correction is incorporated into the interpolation
c   stage in 'imgen_rt'
c
            roll_off  = 0.6 * ( dfmp / df_subimg ) ** 2
c
            write ( 7 , * ) ' roll_off = ' , roll_off
c
         endif
c
c-----------------------------------------------------------------------
c
c             Diagnostic for Sub-Image Weight and IPR
c
         write ( 7 , * )
c
         write ( 7 , * ) '---------------------------------------'
c
         write ( 7 , * )
c
         do ip = 1 + naskip / 2 , 3 * naskip - naskip / 2
c
            write ( 7 , '(i10,2f10.4)' ) ip , wtac(ip)
c
         enddo
c
         cac(:) = cmplx( 0.0 , 0.0 )
c
         do ip = 1 , naskip
c
            cac(ip)          = cmplx( wtac(ip+3*(naskip/2)) , 0.0 )
c
            if ( mod(ip,2) .eq. 0 ) cac(ip) = - cac(ip)
c
            cac(4*naskip+1-ip) = - cac(ip)
c
         enddo
c
         nn(1) = 4 * naskip
c
         call fourt( ac , nn , 1 , 1 , 1 , work , nwork )
c
         const = cabs( cac(1+2*naskip) )  !  Divide by DC
c
         write ( 7 , * )
c
         write ( 7 , * ) '---------------------------------------'
c
         write ( 7 , * )
c
         write ( 7 , * ) '  Harmonic  Doppler     Ampl.    Phase '
c
         write ( 7 , * )
c
         do ip = 1 , 4 * naskip
c
            dop = float( ip - 1 - 2 * naskip) /
     .            ( float( 4 * naskip ) * dtp )
c
            write ( 7 , '(i10,f10.2,2f10.4)' ) ip , dop ,
     .      db( amax1( 1.0e-6, cabs( cac(ip) ) / const ) ) ,
     .      ( 180.0 / pi ) * atan2( aimag( cac(ip) ) , real( cac(ip) ) )
c
         enddo
c
         write ( 7 , * )
c
         write ( 7 , * ) '---------------------------------------'
c
         write ( 7 , * )
c
      else
c
c   For modes which are not real-time, hardwire the Taylor weight used
c   in sub-image formation since it is matched to the overlap factor of
c   2/3
c
         taywtc = - 40.0
c
         call taylor ( 3 * naskip , taywtc , wtac )
c
      endif
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
      iline  = 0                       !  Strip-map line counter
c
      vr0    = 0.0                     !  Mean radial velocity
c
      vrsave = 0.0                     !  vr0 before target detection
c
      vrbar  = 0.0                     !  Target mean velocity
c
      vslip  = 0.0                     !  Corr. for incorrect range walk
c
      dr0    = 0.0                     !  Offset of object from MCP
c
      ntg    = 0                       !  Targets being tracked
c
c   Slip velocity correction passed to IMGEN
c
      vadd   = 0.0
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
c-----------------------------------------------------------------------
c
c   New initialization - KAM 7-3-98
c
      theta_t = pi * angle / 180.0
c
c-----------------------------------------------------------------------
c
c                       END OF INITIALIZATION
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
      isub0 = 1 - ( maxal / 2 )
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
         if ( pass .lt. 2 ) then
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
c   Aliases to search: start at maxal and decrease by two each time
c   until 3
c
            nalias  = max( 3 , 1 + 2 * ( 1 - isub ) )
c
c   Report information to the radar control system so that it can
c   modify the way the data is being collected and so it can adjust
c   the range center of the image
c
            call adjust ( isub , dr0 , vr0 , ntp , ncr , nsr , nspf ,
     .                    nakeep , curtime , dts , r0min , r0max ,
     .                    r0bar , rfmmin , rfmmax , cstr , csbimg ,
     .                    omega0 , vslip , eslip , vadd , vrbar ,
     .                    vrsave ,  pt , ntg , isub0 , afocus , nabuff ,
     .                    tgrbar , tgrmin , tgrmax , nr1 , nre , nr1fm ,
     .                    nrefm )
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
                  pt  = 0.0
                  ntg = 0
c
               endif
c 
c-----------------------------------------------------------------------
c
c   Update the mo-comp position taking into account both the normal
c   mo-comp velocity, vr0, and the slip correction velocity, vslip.
c
               dr0 = dr0 + ( vr0 + vslip + v_corr ) * dtp
c
               if ( alias .eq. - 5 ) then
c
                  read  ( 87 , '(3f16.6)' )
     .                    tdummy , dr0 , vr0
c
               else
c
                  write ( 87 , '(4f16.6)' )
     .                    float( ipulse - spulse ) * dtp , dr0 , vr0 ,
     .                    vslip
c
               endif
c
               f0 = vr0 * ( 2.0 / lambda )  !  Mean Doppler frequency
c
               pt = pt + f0 * dtp
c
c   Compress the pulse for all channels
c
               do ich = 1 , nchtot
c
c*********** Write out Mocomped Data ******** - KAM 10/5/98 ************
c
                  if ( wrtmcd .eq. 1 ) then
c
                     const = ( 2.0 * 1.0E+12 * br * dr0 ) / clight
c
                     do ir = 1 , ntr
c
                        phase   = pt + const * dtr *
     .                                 float( ir - 1 - ntr / 2 )
c
                        crc(ir) = crc_raw(ir,ip,ich) * cpfast( phase )
c
                     enddo
c
                     mcd_idx = mcd_idx + 1
c
                     write ( 93, rec = mcd_idx )
     .                     ( crc(ir) , ir = 1 , ntr )
c
                  endif
c
c***********************************************************************
c
                  do ir = 1 , ntr
c
                     crc(ir) = crc_raw(ir,ip,ich)
c
                  enddo
c
                  call rcomp ( br , clight , dr0 , rc , crc , pt , ntr ,
     .                         ncr , wtr , dtr , mode , work , nwork )
c
                  call quant ( qfirst , 0.2 , rc , ncr )
c
c   Finally, load the compressed pulse into the range-compressed buffer
c   array.
c
                  ipm1 = max( ip - 1 , 1 )
c
                  do ir = 1 , ncr
c
                     cstr(ip,ir,ich) = presum * cstr(ipm1,ir,ich) +
     .                                 ( 1.0 - presum ) * crc(ir)
c
                  enddo
c
               enddo  !  Loop over channels
c
c-----------------------------------------------------------------------
c
c                   Adaptive Motion Compensation
c
               if ( finemc .ne. 0 ) then
c
                  if ( finemc .eq. 1 ) then
c
                     nr1fm = nr1
c
                     nrefm = nre
c
                  endif
c
                  call mocomp ( vr0 , vslip , ctgt , ncr , ip , isub ,
     .                          isub0 , naskip , nfa , nr1 , nre ,
     .                          nr1fm , nrefm , dr0 , vrsave , ntp ,
     .                          sigsq , ipulse , spulse , cstr )
c
               endif
c
            enddo           !  do ip = nt1 , ntp    ( Loop over pulses )
c
            rctime = rctime + dtime( timer )
c
c--------------  Range Compression and Fine Mo-Comp  -------------------
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c------------------  Front-end Target Detector  ------------------------
c
c   Run algorithm to detect and estimate the motion of dominant
c   scatterers.  To save computer time, call every other time.
c
            if ( isub .le. 2 .or. mod( isub , 2 ) .eq. 0 ) then
c
c   Apply the notching logic only to the acquisition phase
c
               if ( isub .le. 2 ) then
                  notchp = notch
               else
                  notchp = 0
               endif
c
c   Calculate the noise level to scale the SNR values in Target
c
               ip = 1 + ntp / 2
c
               do ir = 1 , nre - nr1 + 1
c
                  work(ir)  = cabs( cstr(ip,ir+nr1-1,1) ) ** 2
                  iwork(ir) = ir
c
               enddo
c
               snoise = noise0( work , iwork , nre-nr1+1 , pnoise )
c
               if ( snoise .eq. 0.0 ) snoise = 1.0
c
               if ( isub .lt. 1 ) then
c
c   Reset the target information by calling addtgt with a negative time
c
                  call addtgt ( - 1.0 , 0.0 , 0.0 , 0.0 , 0.0  , 0.0 ,
     .                          't' )
c
               endif
c
               if ( rt_rco .ne. 0 ) then
c
                  call target_rt ( cstr , ncr , ntp , nr1 , nre , 1 ,
     .                             ntp , drc , dtp , cac , lambda ,
     .                             nalias , ntg , tgt(1,1) , tgt(1,2) ,
     .                             tgt(1,3) , tgt(1,4) , tgt(1,5) ,
     .                             tgt(1,6) , tgt(1,7) , tgt(1,8) ,
     .                             tgt(1,9) , tgt(1,10) , tgt(1,11) ,
     .                             tgt(1,12) , tgt(1,13) , tgt(1,14) ,
     .                             vrbar , tgrmin , tgrmax , tgrbar ,
     .                             curtime , snoise , vr0 , dr0 ,
     .                             snrmin , athrsh , notchp , vnotch ,
     .                             dvntch , rnotch , drntch , rwindo ,
     .                             vwindo , eslip )
c
               else
c  
                  call target    ( cstr , ncr , ntp , nr1 , nre , 1 ,
     .                             ntp , drc , dtp , cac , lambda ,
     .                             nalias , ntg , tgt(1,1) , tgt(1,2) ,
     .                             tgt(1,3) , tgt(1,4) , tgt(1,5) ,
     .                             tgt(1,6) , tgt(1,7) , tgt(1,8) ,
     .                             tgt(1,9) , tgt(1,10) , tgt(1,11) ,
     .                             tgt(1,12) , tgt(1,13) , tgt(1,14) ,
     .                             vrbar , tgrmin , tgrmax , tgrbar ,
     .                             curtime , snoise , vr0 , dr0 ,
     .                             snrmin , athrsh , notchp , vnotch ,
     .                             dvntch , rnotch , drntch , rwindo ,
     .                             vwindo , eslip )
c  
               endif
c
               tgtime = tgtime + dtime( timer )
c
c   Assign range-compression and target detection time to a separate
c   CPU time account for acquisition
c
               if ( alias .ne. 0 .and. isub .le. 0 ) then
c
                  actime = actime + rctime + tgtime
                  rctime = 0.0
                  tgtime = 0.0
c
               endif
c
            endif
c
         endif     !  ( pass .lt. 2 ) - Range-compress raw data
c
c------------------  Front-end Target Detector  ------------------------
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
            call subimg ( pass , isub , csbimg(1,1,iisub) , nsa , nsr ,
     .                    ntp , ncr , naskip , nakeep , str , cstr ,
     .                    wtac , work , ac , cac , nwork , nerror )
c
            if ( nerror .ne. 0 ) go to 5000
c
            sitime = sitime + dtime( timer )
c
         endif
c
c   Real-time version of sub-image compression
c
         if ( rt_img .ne. 0 .and. mode .ne. 0 .and. isub .ge. 1 ) then
c
            call subimg_rt ( csbimg(1,1,iisub) , nsa , nsr , ntp , ncr ,
     .                       naskip , nakeep , str , wtac , work , ac ,
     .                       cac , nwork )
c
            sitime = sitime + dtime( timer )
c
         endif
c
c   Diagnostic calculation - rarely used
c
         if ( mode .eq. 10 ) then
c
c   Make a movie of the sub-images
c
            csfile = 'subimgs'
c
            call oframe ( csbimg(1,1,iisub) , work , bytes , - isub ,
     .                    csfile , nsr , nakeep , 0 , 0 , dbinc ,
     .                    nlocal , pixbar , color )
c
         endif
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
            if ( mode .eq. 7 .or. mode .eq. 8 ) frame =
     .              ( isub .ge. nspf .and. mod( isub , nfskip ) .eq. 0 )
c
            if ( frame ) then
c
               tframe = 0.5 * ( curtime + ti0 )
               nspfv  = min( actual , nspf )
               tint   = float( nspfv ) * dts
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
c   In variable integration time mode determine the present Doppler
c   pixel size, the master particle frequency separation and the time
c   step between the virtual sub-images
c
               if ( ( nfskip .lt. 0 ) .and.
     .              omega_valid .and. ( angdot .gt. 0.0 ) ) then
c
                  if ( nfskip .eq. - 1 ) then
c
                     tinteg_use = tinteg *
     .                            amin1( 1.5   ,
     .                                  amax1( 0.666 ,
     .                                         angdot / abs( omega ) ) )
c
                     nspf_use   = nspf
c
                  else if ( nfskip .eq. - 2 ) then
c
                     tinteg_use = tinteg
c
                     nspf_use   = nint( float( nspf ) *
     .                            amin1( 1.3333   ,
     .                                  amax1( 0.75 ,
     .                                       abs( omega ) / angdot ) ) )
c
c
                  else
c
                     tinteg_use = tinteg
c
                     nspf_use   = nspf
c
                  endif
c
               else
c
                  tinteg_use = tinteg
c
                  nspf_use   = nspf
c
               endif
c
               dff    = ( 1.0 / tinteg_use ) * float( nspf_use ) /
     .                                         float( nafill )
c
               dtpga  = 1.0 / ( float( nfa ) * dff )
c
               dtv    = tinteg_use / float( nspf_use )
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
               if ( mode .eq. 1 .or. mode .eq. 6 .or. mode .eq. 7 .or.
     .              mode .eq. 8 ) then
c
                  tint = amin1( 2.0 * ( tframe - 2.0 * dts ) ,
     .                          tinteg_use )
c
               else
c
                  tint = tinteg_use
c
               endif
c
               nspfv  = min( nspf_use , nint( tint / dtv ) )
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
     .                                  pass , isub , iframe
            write ( 6 , * ) '------------------------------------------'
c
            ENDIF
c
            write ( 7 , * )
            write ( 7 , * ) '------------------------------------------'
            write ( 7 , '(a24,3i6)' ) ' Pass, Sub-image, Frame:' ,
     .                                  pass , isub , iframe
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
c   Kalman Filter State Initialization/Maintenance
c   Added by KAM 7-2-98
c
            if ( dokalm .gt. 0 ) then
c
               if ( nfskip .gt. 0 ) then
c
c                 Normal Mode
c
                  if ( isub .ge. nspf ) then
c
                     call kalman_maintain ( dtf , acoefs )    
c
                  else          
c
                     if ( isub .eq. 2 ) then
c
c                       First Time
c
                        call kalman_initialize ( dts , acoefs )              
c
                     else
c
c                       Afterwards
c
                        call kalman_maintain ( dts , acoefs )
c
                     endif              
c
                  endif
c
               else
c
c                 Exact Time Mode
c
                  if ( isub .eq. 2 ) then
c
c                    First Time
c
                     call kalman_initialize ( dtf , acoefs )              
c
                  else
c
c                    Afterwards
c
                     call kalman_maintain ( dtf , acoefs )
c
                  endif
c
               endif
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
            if ( mode .eq. 6 ) then
c
c   Read strip map parameters
c
               gpulse0 = max( spulse + 1 , spulse +
     .                              ifix( ( tv0 - 0.5 * tint ) / dtp ) )
c
               read ( 94 , rec = gpulse0 ) strip
c
               vplat0  = strip(2) * sin( strip(5) )
c
               vspot0  = strip(3)
c
               gpulse1 = max( spulse , spulse +
     .                              ifix( ( tv0 + 0.5 * tint ) / dtp ) )
c
               read ( 94 , rec = gpulse1 ) strip
c
               vplat1  = strip(2) * sin( strip(5) )
c
               vspot1  = strip(3)
c
               if ( gpulse1 .gt. gpulse0 ) then
c
                  ddt_vplat = ( vplat1 - vplat0 ) /
     .                        ( float( gpulse1 - gpulse0 ) * dtp )
c
                  ddt_vspot = ( vspot1 - vspot0 ) /
     .                        ( float( gpulse1 - gpulse0 ) * dtp )
c
               endif
c
               gpulse  = max( spulse , spulse + ifix( tv0 / dtp ) )
c
               read ( 94 , rec = gpulse ) strip
c
               write ( 6  , '(1x,4i6,3f7.2,4f8.2)' ) iframe , gpulse ,
     .                 gpulse0 , gpulse1 , ddt_vplat , ddt_vspot ,
     .                 strip(5) , vplat0 , vplat1 , vspot0 , vspot1
c
               write ( 66 , '(1x,4i6,3f7.2,4f8.2)' ) iframe , gpulse ,
     .                 gpulse0 , gpulse1 , ddt_vplat , ddt_vspot ,
     .                 strip(5) , vplat0 , vplat1 , vspot0 , vspot1
c
            endif
c
            if ( mode .eq. 7 ) then
c
c   Read spotlight parameters
c
               gpulse0 = max( spulse + 1 , spulse +
     .                             ifix( ( tv0 - 0.25 * tint ) / dtp ) )
c
               read ( 94 , rec = gpulse0 ) strip
c
               vplat0  = strip(2) * sin( strip(5) )
c
               vspot0  = strip(3)
c
               omeg_0  = abs( vplat0 / strip(1) )
c
               gpulse1 = max( spulse , spulse +
     .                             ifix( ( tv0 + 0.25 * tint ) / dtp ) )
c
               read ( 94 , rec = gpulse1 ) strip
c
               vplat1  = strip(2) * sin( strip(5) )
c
               vspot1  = strip(3)
c
               omeg_1  = abs( vplat1 / strip(1) )
c
               if ( gpulse1 .gt. gpulse0 ) then
c
                  ddt_vplat      = ( vplat1 - vplat0 ) /
     .                             ( float( gpulse1 - gpulse0 ) * dtp )
c
                  ddt_vspot      = ( vspot1 - vspot0 ) /
     .                             ( float( gpulse1 - gpulse0 ) * dtp )
c
                  omegadot_omega = 2.0 * ( ( omeg_1 - omeg_0 ) /
     .                                     ( omeg_1 + omeg_0 ) ) /
     .                             ( float( gpulse1 - gpulse0 ) * dtp )
c
               else
c
                  ddt_vplat      = 0.0
c
                  ddt_vspot      = 0.0
c
                  omegadot_omega = 0.0
c
               endif
c
               gpulse  = max( spulse , spulse + ifix( tv0 / dtp ) )
c
               read ( 94 , rec = gpulse ) strip
c
               write ( 6  , '(1x,4i6,3f7.2,4f8.2)' ) iframe , gpulse ,
     .                 gpulse0 , gpulse1 , ddt_vplat , ddt_vspot ,
     .                 strip(5) , vplat0 , vplat1 , vspot0 , vspot1
               write ( 6 , '(1x,f12.5)' ) omegadot_omega
c
               write ( 66 , '(1x,4i6,3f7.2,4f8.2)' ) iframe , gpulse ,
     .                 gpulse0 , gpulse1 , ddt_vplat , ddt_vspot ,
     .                 strip(5) , vplat0 , vplat1 , vspot0 , vspot1
               write ( 66, '(1x,f12.5)' ) omegadot_omega
c
            endif
c
c-----------------------------------------------------------------------
c
            if ( alias .ne. - 4 .and. alias .ne. - 6 ) then
c
c   Estimate target reports from the sub-image field
c
               if ( tgt_si .gt. 0 )
     .         call gettgt_s ( csbimg , nakeep , nsr , nabuff , pnoise ,
     .                         pgatgt , iwork , dts , iisub , isub ,
     .                         work , nptgt , work , dopcen ,
     .                         0.5 * ( curtime + ti0 +
     .                         dts * float( ( 2 * nabuff ) / 3 ) ) )
c
               if ( rt_img .eq. 0 .or. color .ne. 0 ) then
c
                  call getacc    ( 1 + nfa / mprat , nfr / mrrat , 
     .                             nspfv , dotdot ,
     .                             fcenuse - 0.5 * float( nfa ) * dff ,
     .                             rcenuse + rfmin + 0.5 * drf * mrrat ,
     .                             tv0 , dfmp , drfbig , dtv , athrsh ,
     .                             rcen , fcen )
c
               else
c
c                 Real-time (and not color mode)
c
                  call getacc_rt ( 1 + nfa / mprat , nfr / mrrat ,
     .                             nspfv + nafill_a - nafill , dotdot ,
     .                             fcenuse - 0.5 * float( nfa ) * dff ,
     .                             rcenuse + rfmin + 0.5 * drf * mrrat ,
     .                             tv0 , dfmp , drfbig , dtv , athrsh ,
     .                             rcen , fcen )
c
               endif
c
            else
c
c   Constant rotation case (alias=-4) - omega was input manually.  The
c   acceleration is proportional to range_squared.  When combined with
c   the range walk correction in imgen, this is equivalent to the polar
c   format approximation for constant rotation rate.
c
c   Constant focus correction (alias=-6) - afocus was input manually
c
               call polar ( omega0 , afocus , lambda , nspfv ,
     .                      nfr / mrrat , 1 + nfa / mprat , dtv ,
     .                      drfbig , dff * float( mprat ) , dotdot )
c
            endif
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
            if ( mode .eq. 1 .or. mode .eq. 4 .or. mode .eq. 6 .or.
     .           mode .eq. 7 .or. mode .eq. 8 ) then
c
c              Inverse-SAR Mode and ISAR-like Strip Map Processor
c
               if ( rt_img .ne. 0 ) then
c
c-----------------------------------------------------------------------
c
c                 Real-time Mode
c
c   Flip sign of odd weights to center the image after the FFT
c
                  do ip = 1 , nspfv , 2
c
                     wtaf(ip) = - wtaf(ip)
c
                  enddo
c
c   Real-time image formation processor with in-line PGA
c
                  call imgen_rt ( cframe , nfa , nfr , dff , drf ,
     .                            csbimg , nakeep , nsa , nsr , dfc ,
     .                            drs , fc0 , actual , nabuff , nafill ,
     .                            iisub , cac , mprat , wtaf , work ,
     .                            nwork , lambda , dotdot , mrrat ,
     .                            nspfv , dtv , tframe , rt_nmp ,
     .                            cframe(1,1+nfr/2,3) ,
     .                            cframe(1,1+nfr/2+nfr/4,3) ,
     .                            pgatgt , iwork , iwork(nfr+1) ,
     .                            cframe(1,1,3) , work(2*nafill+2) ,
     .                            nafill_a , rc , crc )
c
               else
c
c-----------------------------------------------------------------------
c
c                 Normal mode
c
                  call imgen    ( cframe , nfa , nfr , dff , drf ,
     .                            csbimg , nakeep , nsa , nsr , dfc ,
     .                            drs , fc0 , actual , nabuff , nafill ,
     .                            iisub , dts , ac , cac , mprat ,
     .                            wtaf , work , nwork , br , lambda ,
     .                            clight , dotdot , mrrat , rgwalk ,
     .                            vadd , nspfv , dtv , ti0 , tframe ,
     .                            iwork , strtch , rcenuse , fcenuse )
c
               endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
            else if ( mode .eq. 2 .or. mode .eq. 3 ) then
c
c              Ultra-Wide-Band Strip-map
c
               call imgenu ( cframe , nfa , nfr , dff , drf , csbimg ,
     .                       nakeep , nsr , dfc , drs , fc0 , actual ,
     .                       nabuff , nafill , iisub , dts , ac , cac ,
     .                       mprat , wtaf , work , nwork , br , lambda ,
     .                       clight , dotdot , mrrat , rgwalk , vadd ,
     .                       nspfv , dtv , ti0 , tframe , iwork ,
     .                       strtch , vfocus , slant0 )
c
            endif
c
            frtime = frtime + dtime( timer )
c
c----------------  Coherent Frame Integration  -------------------------
c-----------------------------------------------------------------------
c***********************************************************************
c-----------------------------------------------------------------------
c-----------------  Phase Gradient Autofocus ---------------------------
c
            if ( pass .eq. 0 .or. pass .eq. 2 ) then
c
c   This is the last pass - autofocus the frame before output if
c   required
c
               maxitu = maxit
c
            else
c
c   This is not the last pass before output - don't bother to autofocus
c
               maxitu = 0
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
     .                    work(5*nfa+1) , work(6*nfa+1) , icount ,
     .                    dtpga , corg , cnew , work(7*nfa+1) , nptgt ,
     .                    pnoise , pgatgt , work(8*nfa+1) , dotdot ,
     .                    mprat , mrrat , nafill , nspfv , dfmp ,
     .                    tframe , rcenuse , fcenuse , nfocus )
c
               pgtime = pgtime + dtime( timer )
c
            else
c
               cnew   = ctrast( cframe , nfa , nfr )
c
               corg   = cnew
c
            endif 
c
            write ( 82 , '(1x,a3,3f12.4)' ) '#p#' , tframe ,
     .                                       db( corg ) , db( cnew )
c
            if ( quiet .gt. 1 ) write (  6 , '(1x,a26,3f12.4)' )
     .                          ' Autofocus efficiency:    ' , tframe ,
     .                            db( corg ) , db( cnew )
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
            if ( mode .eq. 6 ) then
c
c   Strip-map Mode
c
               flight_path = flight_path + strip(2) * dtf
c
               ground_path = ground_path + strip(3) * dtf
c
               mtrs2hz     = ( 2.0 / lambda ) * strip(2) / strip(1)
c
               nxmid       = 2 * ifix ( mtrs2hz *
     .                                  ( ground_path - image_path )
     .                                 / ( 2.0 * dff ) )
c
               nxtotal     = nxtotal + nxmid
c
               image_path  = image_path  + dff * float( nxmid ) /
     .                                           mtrs2hz
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
               write ( 6  , '(1x,4i8,7f11.2)' )
     .               iframe , gpulse , nxmid  , nxtotal    ,
     .               tv0         , strip(1)   , strip(2)   , strip(3) ,
     .               flight_path , image_path , ground_path - image_path
c
               write ( 66 , '(1x,4i8,7f11.2)' )
     .               iframe , gpulse , nxmid  , nxtotal    ,
     .               tv0         , strip(1)   , strip(2)   , strip(3) ,
     .               flight_path , image_path , ground_path - image_path
c
c   Only write out pixels if the integration time is up to steady state
c
               if ( ( gpulse - spulse ) .gt.
     .              ifix( 0.5 * tinteg / dtp ) )
     .            call oframe ( cframe , work , bytes , pframe , ofile ,
     .                          nfa , nfr , nxmid , 0 , dbinc , nlocal ,
     .                          pixbar , color )
c
            else if ( pass .eq. 0 .or. pass .eq. 2 ) then
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
     .                       dbinc , nlocal , pixbar , color )
c
               write ( 7 , '(/,a,i6,2(a,f10.3),/)' )
     .                 ' Contrast for frame:' , iframe , ' = ' , cnew ,
     .                 ' at time: ' , tframe
c
               if ( ( rt_img .ne. 0 ) ) then
c
c   In real-time mode there are two frames computed at once - write out
c   the second one now
c
                  iframe = iframe + 1
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
                  call oframe ( cframe(1,1,2) , work , bytes , pframe ,
     .                          ofile , nfa , nfr , 0 , 0 , dbinc ,
     .                          nlocal , pixbar , color )
c
                  cnew = ctrast( cframe(1,1,2) , nfa , nfr )
c
                  corg = cnew
c
                  write ( 7 , '(/,a,i6,2(a,f10.3),/)' )
     .                 ' Contrast for frame:' , iframe , ' = ' , cnew ,
     .                 ' at time: ' , tframe
c
               endif
c
c   Compute and output the A-Scan vector
c
               if ( rt_img .eq. 0 ) then
c
                  rt_frames = 1
c
               else
c
                  rt_frames = 2
c
               endif
c
               do itemp = 1 , rt_frames
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
     .                           conjg( cframe(ia,ir,1+nfocus) )
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
c   This smoothing routine should only be used for cases for which the
c   amplitude modulations are expected to be large scale
c
c   This is strictly a diagnostic procedure to make the display of the
c   ascan plot look better - it does not affect the image
c
                  if ( smooth_az_ascan ) then
c
                     call stats ( ascan(1+nfr) , nfa , asmin , asmax ,
     .                            asbar , assdv )
c
                     do ia = 1 , nfa
c
                        ascan(ia+nfr) = min( ascan(ia+nfr) ,
     .                                       2.0 * asbar )
c
                     enddo
c
                     call smooth ( ascan(1+nfr) , work , nfa , 1 , 3 ,
     .                             0 , 1 )
c
                  endif
c
c-----------------------------------------------------------------------
c
                  write ( 89 , rec = iframe ) ascan
c
               enddo
c
               if ( ( mode .gt. 1 ) .and. ( mode .lt. 6 ) ) then
c
c   Write out frame to strip-map file
c
                  call iocimg ( ounit , iline , nfa , nfr , 1 , cframe ,
     .                          ne )
c
                  iline = iline + nfa  !  Total lines written to file
c    
               endif
c
               write ( 92 , '(6f12.2,3i9,2f12.4,2f14.2)' )
     .               rcenuse    , fcenuse    , length_est , width_est  ,
     .               tframe     , tinteg_use , nspf       , nspf_use   ,
     .               iframe     , ais(1)     , ais(2)     , ais(3)     ,
     .               cnew
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
         if ( pass .le. 1 ) then
c
            if ( alias .eq. - 5 ) then
c
               read  ( 84 , '(6f12.4)' ) tgrbar , tgrmin , tgrmax ,
     .                                   r0min  , r0max  , r0bar
c
            else
c
               write ( 84 , '(6f12.4)' ) tgrbar , tgrmin , tgrmax ,
     .                                   r0min  , r0max  , r0bar
c
            endif
c
         endif
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
      if ( nerror .ne. 0 .and. iframe .eq. 0 ) pass = - 1
c
      return
      end
