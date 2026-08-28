C**
C***********************************************************************
C**
      subroutine defaults ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                      npulse , nspf , nspt , nafill , nfr , nfa ,
     .                      mrrat , mprat , nchuse , nfocus )
C**
C***********************************************************************
C**
c   This routine sets the default parameters.  It is called at the start
c   of program ISAR-T and when the command 'reset' is given.  It is also
c   called at the start subroutine 'radars' which is invoked by the
c   command "radar = Radar_Type".
c
      implicit none
c
c***********************************************************************
c
      integer      ntr , ncr , naskip , spulse , npulse , nspf , nspt ,
     .             nafill , nfr , nfa
c
      real         sampr , prf
c
c***********************************************************************
c
      include     'sarprm.h'       !  Standard ISAR-T parameters
c
      include     'updates.h'      !  Updates to parameters
c
      include     'realtime.h'     !  Real-Time parameters
c
      include     'kalman.h'       !  Kalman Filter Parameters
c
      include     'tglist.h'       !  Target list definitions
c
      integer      nchuse , mrrat , mprat , nfocus
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c        Set default parameters for ISART
c
c-----------------------------------------------------------------------
c
c        File names
c
      ifile  = 'test'    !  Input and output file names for internal
      ofile  = 'test'    !  test case
c
c-----------------------------------------------------------------------
c
c   Parameters which control the reading and quality control of
c   raw data.
c
      nfmt   = - 1       !  Internally generated I/Q
      mode   = 1         !  ISAR Mode
      fixiq  = 0         !  No data repair functions
      alias  = 1         !  Alias correction option
      nrcent = 0         !  Shift in range center for some data types
c
      efghz  = 10.0      !  X-Band
c
      br     = 15.0      !  15 MHz/Microsecond
      sampr  = 5.12E+7   !  51.2 MHz Sampling in fast time
      strtch = 1         !  Stretch processing in range compression
      ntr    = 512       !  Fast time samples
      prf    = 512.0     !  512 Hz Pulse Repetition Frequency
      ncr    = 1536      !  Compressed samples
      overrg = 1.0       !  No over-sampling in range in final image
      naskip = 32        !  Pulses between sub-images
      spulse = 1         !  Starting pulse
      npulse = 672       !  Number of pulses to be processed
      nspf   = 16        !  Sub-Images per frame
      nspt   = 16        !  Sub-images per target
      nafill = nspf      !  No over-sampling in Doppler
      taywtr = - 40.0    !  Taylor weight for range compression
      taywta = - 40.0    !  Taylor weight for Doppler
      nfr    = 256       !  Range cells in fine resolution image frame
      nfa    = 256       !  Doppler cells in fine resolution image frame
      wrtmcd = 0         !  Don't write out mocomped raw data
      v_corr = 0.0       !  Additional range-walk velocity
      a_corr = 0.0       !  Additional acceleration for strip-map focus
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c        Parameters used to determine 'large' scales used to simplify
c        the calculation of the fine resolution image
c
      mprat  = 4         !  Master particle ratio
      mrrat  = 8         !  'dotdot' array interpolation ratio
c
c        Normally, the frame rate is determined by an integer number of
c        sub-images
c
      nfskip = 1         !  No. of sub-images to skip between fine
c                           resolution frames
c
c        Parameters used for 'exact time' mode in which the frames are
c        calculated at exactly the desired frame rate and integration
c        time independent of the PRF and sub-image rates.
c
      frrate = 10.0      !  Frame rate ( if nfskip=0 )
      tinteg = 1.0       !  Integration time ( if nfskip=0 )
c
c-----------------------------------------------------------------------
c
c        Fine motion compensation parameters used in subroutine 'mocomp'
c
      finemc = 4         !  Fine motion compensation option
      beta   = 0.2       !  Smoothing coefficient for raw data stats.
      alpha0 = 0.5       !  Fine motion compensation relaxation constant
      nlag   = 1         !  Lag number for motion compensation
c
c-----------------------------------------------------------------------
c
c        Parameters for computing eight-bit intensity image frames
c        Details in subroutine 'oframe'
c
      dbinc  = 0.25      !  dB increment for output (>99 = auto-scale)
      nlocal = 0         !  Automatic choice of size of local region
      pixbar = 64        !  Average grey level
      color  = 0         !  Do not do color coding of targets
c
c-----------------------------------------------------------------------
c
c        Acquisition mode parameters - used in subroutine 'target'
c
      snrmin = 100.0     !  Min. SNR for target declaration
c
      notch  = 0         !  Flag for absolute velocity notching (land
                         !  clutter rejection)
c
c   If notch < 0, the notch region is ignored; if notch > 0, the region
c   outside the notch region is ignored
c
      vnotch = 0.0       !  Absolute velocity to notch
      dvntch = 0.0       !  Velocity width of notch
      rnotch = 0.0       !  Absolute range to notch
      drntch = 0.0       !  Range width of notch
c
c-----------------------------------------------------------------------
c
c        Phase Gradient Autofocus (PGA) algorithm parameters - used in
c        subroutine 'pga'
c
      nrkeep = 16        !  Num of range cells to form phase error est.
      npass  = 61        !  Num of passband bins for PGA
      ntaper = 10        !  Num of taper bins for PGA
      maxit  = 20        !  Num max iterations for stop of PGA
      slbkil = 0         !  No Sidelobe Target Editing KAM 7/30/98
c
      pgtype = 0         !  Traditional PGA - pgtype = 2 for Jakowatz
c                        !  new algorithm
c
c-----------------------------------------------------------------------
c
c        Image formation parameter - the parameter rgwalk is used
c        only for testing the importance of the range walk effects.
c        However, it can also be used to enhance the detection of
c        vibration echoes.
c
      rgwalk = 1.0       !  Set to 0.0 to disable range walk
c
c-----------------------------------------------------------------------
c
c        Parameters for fitting acceleration estimates - used in
c        subroutine 'getacc'
c
      curve  = 4         !  Flag for getacc: 0---> fill acc array with 0
c                        !                   1---> estimate acceleration
c                        !                         from front-end
c                        !                   2---> Use only focus info
c                        !                         from the PGA
c                        !                   3---> Use both sources
c                        !
c                        !  Note: Higher values of curve are used for
c                        !  motion models which apply physical
c                        !  constraints.  The above numbers are then
c                        !  interpreted as the mod(4) of curve.
c
      outlie = 3.0       !  No. of std. devs. for outlier control
      acoefs = 6         !  No. of coeficients in acceleration fit
      aghost = 1.0       !  Weight for acceleration ghost particles
      center = 0.2       !  Weight for auto-centering algorithm
      dokalm = 0         !  Don't do Kalman filter
c
      tgt_si = 1         !  Use sub-image targets for GMM
c
      rd_tgt = 0         !  Don't read target list from disk
c
      firsts = 2         !  First subimage for frames
      
      uwb    = 0         !  Do not use Ultra-Wide-Band sub-images
c
      rlook  = 0         !  Multilooks in range
      flook  = 0         !  Multilooks in Doppler
c
      addvib = 0         !  Flag for phase vibration
      addamp = 0         !  Flag for amplitude vibration
      fixamp = 0         !  Flag for amplitude correction
c
      rcentr = 0.0       !  Range offset for Mode 10
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c        Variables used only for internal I/Q generation for
c        testing algorithms - used in subroutine 'iqgen'
c
c-----------------------------------------------------------------------
c
c        Distribution of test particles for nfmt = - 1
c
      ntest  = 17        !  No. of test targets
      drtest = 12.0      !  Range difference between test targets
      dvtest = 0.1       !  Velocity difference between test targets
      rbar   = 0.0       !  Offset added to initial range
c
c-----------------------------------------------------------------------
c
c        Velocity parameters for all simulation types
c
      vbar   = 0.0       !  Average velocity of test targets
      dvdt   = 0.0       !  Acceleration rate of test targets
      vamp   = 0.0       !  Amplitude of velocity vibration (m/s)
      vfreq  = 0.0       !  Frequency of velocity vibration (Hz)
      vamp2  = 0.0       !  Amplitude of vibration (m/s) (secondary)
      vfreq2 = 0.0       !  Frequency of vibration (Hz) (secondary)
c
c-----------------------------------------------------------------------
c
c        Amplitude parameters for all simulation types
c
      aamp   = 0.0       !  Amplitude of amplitude vibration
      dadt   = 0.0       !  Increase of amplitude with time
      afreq  = 0.0       !  Frequency of amplitude vibration (Hz)
c
      noise  = 1.0       !  Gaussian noise level ( RMS i or q )
c
c-----------------------------------------------------------------------
c
c        Velocity error sources
c
      vntime = 1.0 / prf !  Time constant for smoothing velocity noise
      vnoise = 0.0       !  Magnitude of random velocity noise
      slipv  = 0.0       !  Slip velocity error
c
c-----------------------------------------------------------------------
c
c        Error in acceleration versus frequency, range, and time
c
      daccdf = 0.0       !  Derivative of acceleration (Hz/s) wrt freq.
      daccdr = 0.0       !  Derivative of acceleration (Hz/s) wrt range
      daccdt = 0.0       !  Derivative of acceleration (Hz/s) wrt time
c
c-----------------------------------------------------------------------
c
c        Variables used in generation of clutter for simulated data.
c
c               vclut  = clutter velocity (m/s)
c               aclut  = clutter amplitude (m^2)
c               vdclut = clutter velocity variablity (radians)
c               pclut  = clutter distribution (percent) 
c               dclut  = duration of each clutter object (seconds)
c
      vclut  = 2.0
      aclut  = 1.0
      vdclut = 0.5
      pclut  = 0.0
      dclut  = 1.0
c
c-----------------------------------------------------------------------
c
c        Strip-map parameters
c
      vplat  = 100.0     !  Platform velocity used in simulator
      vfocus = 100.0     !  Platform velocity used in focus
      slant0 = 1200.0    !  Slant range to scene center
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c        Parameters added after last major release
c
c-----------------------------------------------------------------------
c
c        Pre-sum parameters
c
      prfrat = 1.0       !  Assume that no pre-summing has been done
c                           before input to this program - used in
c                           subroutine 'target' and 'target_rt'
c
      presum = 0.0       !  No pre-sum filter in this program
c
c-----------------------------------------------------------------------
c
c        Multiple channels for clutter cancellation
c
      nchuse = 1         !  Use only one data channel
c
      cmodel = 0         !  Simplest cancellation model
c
c-----------------------------------------------------------------------
c
c        Corrections to account for definitions of other data
c
      moddat = 0         !  No data corrections required - details in
c                           subroutine fillbuff.f
c
c-----------------------------------------------------------------------
c
c        Movie data type: 0 = each frame is a different file
c
      otype  = 2         !  Single large file for output frames, with
c                           target identification
c
c-----------------------------------------------------------------------
c
c        Multi-pass option - allows global motion model to use future
c        data as well as past data
c
      multip = 0         !  Multi-pass control parameter
c
c-----------------------------------------------------------------------
c
      editac = 2         !  Level of editing for acceleration estimates
c                           used in gettacc.f and getacc_rt.f
c
c-----------------------------------------------------------------------
c
      quiet  = 1         !  Set to 1 to reduce screen output; also makes
c                           the code run faster
c
c-----------------------------------------------------------------------
c
c        A Priori estimates of target parameters
c
      length = - 1.0     !  Expected length of target: negative = N/A
c
      width  = - 1.0     !  Expected width of target:  negative = N/A
c
      angle  =   90.0    !  Expected target angle w.r.t. range (degrees)
c
      ti_con =   0.0     !  Parameter for variable integration time
c
c-----------------------------------------------------------------------
c
c        Step-frequency parameters
c
      nbands = 1         !  Number of bands
c
      rfdelt = 0.0       !  Difference in RF between bands
c
      alt_km = 0.0       !  Altitude (km)
c
c-----------------------------------------------------------------------
c
c        Alternative focus settings
c
      nfocus = 0
c
      accorr = 0.0
c
      arcorr = 0.0
c
      n_vadd = 0
c
      d_vadd = 0.0
c
c-----------------------------------------------------------------------
c
c        Quantizer settings
c
      qantiz = 0
c
      iqlsb  = 0.2
c
c-----------------------------------------------------------------------
c
c        Clutter reduction
c
      reduce = 0
c
      imgslk = 0
c
c-----------------------------------------------------------------------
c
c        Real-time parameters
c
      real_t = 0         !  No general real-time simplifications
c
      rt_rco = 0         !  No real-time simplifications - range
c                           compression
c
      rt_img = 0         !  No real-time simplifications - image
c                           formation
c
      rt_pga = 0         !  No real-time simplifications in PGA
c
      rt_io  = 0         !  No real-time simplifications in I/O
c
c-----------------------------------------------------------------------
c
      return
      end
