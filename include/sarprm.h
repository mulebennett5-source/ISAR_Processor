c
c***********************************************************************
c
c                          File: sarprm.h
c
c        Purpose:  To store parameters for the SAIC/Telephonics ISAR
c                  processor which are needed by a variety of routines
c
c***********************************************************************
c
c             Basic data and Motion Compensation parameters
c
c     Parameter       Type     FUNCTION
c
c-----------------------------------------------------------------------
c
c     nfmt           Integer   Data format: < 0 --> Internally simulated
c
c                                        - 1 : 1-D pattern of scatterers
c
c                                        - 2 : 2-D pattern of scatterers
c
c                                        - 3 : General 3-D object
c
c                                        - 4 : Idealized strip-map
c
c                                        - 5 : Exact strip-map
c
c     [Details of simulated data types are given in subroutine iqgen]
c 
c                                       >= 0 --> Data read from disk
c
c                                          0 :  8-Bit I and Q
c
c                                          1 : 16-Bit I and Q
c
c                                          2 : 32-Bit I and Q (Float)
c
c                                          3 : NRL Format - 768 samples
c
c                                          4 : Lincoln Lab - swath
c
c                                          5 : Lincoln Lab - resolution
c
c                                          Etc.
c
c     [Details of real data types are given in subroutine rawdat]
c
c     mode            Integer  Flag for processing mode
c
c                           Examples:
c
c                              Mode = 0 : Survey
c                              Mode = 1 : ISAR
c                              Mode = 2 : Strip-Map
c
c                           Details of other modes in subroutine movie.f
c
c     strtch          Integer  Flag for stretch processing
c
c                              Strtch = 0 : Raw data was fully range-
c                                           compressed
c
c                              Strtch = 1 : Raw data was fast-time for
c                                           a linear-swept-frequency
c                                           radar system
c
c     fixiq           Integer  Flag for data 'cleanup' (used only for
c                              extremely bad data)
c
c     finemc          Integer  Flag for fine mo-comp (details in routine
c                              mocomp.f)
c
c     alias           Integer  Flag for alias detection (details in
c                              routine adjust.f)
c
c-----------------------------------------------------------------------
c
c              Velocity and range notch parameters - 
c              details in routine target.f
c
c     notch           Integer  Flag for notching out clutter in velocity
c
c     vnotch          Real     Center of velocity range to be used (m/s)
c
c     dvntch          Real     Width of velocity range to be used (m/s)
c
c     rnotch          Real     Center of range to be used (m)
c
c     dvntch          Real     Width of range to be used (m)
c
c-----------------------------------------------------------------------
c
c     nlag            Integer  Lag number used in fine mo-comp
c
c     nfskip          Integer  No. of sub-images to skip between frames
c
c     efghz           Real     Center frequency (MHz) [User choice]
c
c     clight          Real     Speed of light (m/s)   [Internal]
c
c     lambda          Real     Wavelength (m)         [From efghz]
c
c     br              Real     Chirp rate (MHz/Microsecond)
c
c     dtr             Real     Fast time sample time (sec)
c
c     dtp             Real     Slow time sample time (sec)
c
c     dbinc           Real     dB increment for output intensity images
c
c     nlocal          Integer  Size of local region for auto-scale
c
c     pixbar          Integer  Average grey level for auto-scale
c
c     color           Integer  Flag for color coding of targets in the
c                              8-bit output image
c
c     acoefs          Integer  Number of coefficients uses in routine
c                              GETACC
c
c     aghost          Real     Scale factor for weighting 'ghost'
c                              particles in routine GETACC
c
c     beta            Real     Relaxation constant for coherent
c                              smoothing of time series before computing
c                              mo-comp statistics
c
c     alpha0          Real     Relaxation constant for fine mo-comp
c
c     snrmin          Real     Minimum Signal-to-Noise-Ratio for a
c                              target range cell to be accepted
c
c-----------------------------------------------------------------------
c
      integer           nfmt   , fixiq  , finemc , alias  , notch  ,
     .                  nlag   , nfskip , nlocal , pixbar , acoefs ,
     .                  nrcent , mode   , strtch , color
c
      real              efghz  , clight , lambda , br     , dtr    ,
     .                  dtp    , vnotch , dvntch , rnotch , drntch ,
     .                  dbinc  , aghost , beta   , alpha0
c
c-----------------------------------------------------------------------
c
c             Parameter which control frame generation
c
c     curve           Integer  Flag for Global Motion Model options
c
c     taywtr          Real     Taylor weight for range compression
c
c     taywta          Real     Taylor weight for Doppler compression
c
c     overrg          Real     Over-sampling ratio for range
c
c     outlie          Real     No. of standard deviations for target
c                              rejection (used in getacc.f)
c
c     rgwalk          Real     Magnitude of range walk (usually 1.0)
c                              Used for advanced testing and vibration
c                              detection
c
c     frrate          Real     Frame rate (used only if nfskip=0)
c
c     tinteg          Real     Integration time (used only if nfskip=0)
c
      integer           curve
c
      real              taywtr , taywta , overrg , outlie , rgwalk ,
     .                  frrate , tinteg , snrmin
c
c-----------------------------------------------------------------------
c
c             Common block for passing standard parameters
c
      common / sarprm / nfmt   , fixiq  , finemc , nlag   , alias  ,
     .                  efghz  , clight , lambda , br     , dtr    ,
     .                  dtp    , curve  , nfskip , dbinc  , nlocal ,
     .                  pixbar , taywtr , taywta , overrg , notch  ,
     .                  snrmin , vnotch , dvntch , rnotch , drntch ,
     .                  outlie , acoefs , aghost , beta   , alpha0 ,
     .                  rgwalk , nrcent , frrate , tinteg , mode   ,
     .                  strtch , color
c
c-----------------------------------------------------------------------
c
c         Parameters for the Phase Gradient Autofocus algorithm
c
c     nrkeep          Integer  Maximum number of range cells for phase
c                              error iteration
c
c     npass           Integer  Number of Doppler pixels in blur function
c
c     ntaper          Integer  Number of Doppler cells in taper from
c                              the blur function to the surrounding
c
c     maxit           Integer  Maximum number of iterations
c
c     slbkil          Integer  Option to reduce the influence of
c                              sidelobes
c
c     percon          Real     Percent of range cell histogram used to
c                              adaptively set the noise floor
c
      integer           nrkeep , npass , ntaper , maxit , slbkil
c
      real              percon
c
c-----------------------------------------------------------------------
c
c	     Common block for passing standard parameters
c
      common / pgaprm / nrkeep , npass , ntaper , maxit , slbkil ,
     .                  percon

c
c-----------------------------------------------------------------------
c
c***********************************************************************
c***********************************************************************
c
c       Common blocks for file names and the raw data histogram
c
      character         ifile*80 , ofile*80
c
      common / files  / ifile , ofile
c
      integer           hist(2,65536)
c
      common / histgm / hist
c
c***********************************************************************
c***********************************************************************
c
c        Common block for passing parameters to subroutine 'iqgen',
c        the built-in radar data simulator
c
c     ntest           Integer  Number of test targets used
c
c     ntestm          Integer  Maximum number of targets allowed
c                              (parameter)
c
c     dvdt            Real     Acceleration added to all targets
c
c     vamp            Real     Periodic velocity fluctuation added to
c                              all targets
c
c     vfreq           Real     Frequency (Hz) of velocity fluctuation
c
c     vamp2           Real     Periodic velocity fluctuation added to
c                              all targets (secondary vibration)
c
c     vfreq2          Real     Frequency (Hz) of velocity fluctuation
c                              (secondary vibration)
c
c     drtest          Real     Range separation of test pattern
c                              (nfmt = - 1)
c
c     dvtest          Real     Velocity separation of test pattern
c                              (nfmt = - 1)
c
c     v0(ntestm)      Real     Initial velocity of targets (nfmt = - 2)
c
c     v(ntestm)       Real     Velocity of targets (nfmt = - 2)
c
c     r(ntestm)       Real     Range of targets (nfmt = - 2)
c
c     p(ntestm)       Real     Work array for computing target signal
c
c     c(ntestm)       Real     Work array for computing target signal
c
c     noise           Real     RMS I or Q value of receiver noise
c
c     dadt            Real     Rate of change of amplitude applied to
c                              all targets
c
c     aamp            Real     Periodic amplitude fluctuation added to
c                              all targets
c
c     afreq           Real     Frequency (Hz) of amplitude fluctuation
c
c     daccdf          Real     Derivative of acceleration with respect
c                              to doppler frequency
c
c     daccdr          Real     Derivative of acceleration with respect
c                              to range
c
c     daccdt          Real     Derivative of acceleration with respect
c                              to time
c
c     tgtamp(ntestm)  Real     Target amplitude work array
c
c     tgtamp0(ntestm) Real     Target amplitude reference values
c
c     tgdvdt(ntestm)  Real     Target acceleration work array
c
c     vnoise          Real     Platform velocity RMS noise (m/s)
c
c     vntime          Real     De-correlation time for platform velocity
c
c     slipv           Real     Doppler velocity inconsistency error
c
c-----------------------------------------------------------------------
c
c     vclut           Real     Typical velocity of clutter spikes
c
c     aclut           Real     Typical amplitude of clutter spikes
c
c     vdclut          Real     Standard deviation of the velocity of
c                              clutter spikes
c
c     pclut           Real     Percentage of range cells with clutter
c                              spikes
c
c     dclut           Real     Typical duration of clutter spikes
c
c-----------------------------------------------------------------------
c
c     vplat           Real     Platform speed for strip-map mode
c
c     vspot           Real     Velocity of illumination pattern
c
c     slant0          Real     Broadside range for strip-map
c
c-----------------------------------------------------------------------
c
      integer           ntestm , ntest
c
      parameter       ( ntestm = 10000 )
c
      real              dvdt , vbar , vamp , vfreq , drtest , dvtest ,
     .                  v0(ntestm) , r(ntestm) , v(ntestm) , p(ntestm) ,
     .                  c(ntestm) , noise , dadt , aamp , afreq ,
     .                  daccdf , daccdr , daccdt , tgtamp(ntestm),  
     .                  tgtamp0(ntestm) , vnoise , vntime ,
     .                  tgdvdt(ntestm) , slipv , rbar , vclut , aclut ,
     .                  vdclut , pclut , dclut , cr(ntestm) ,
     .                  cv(ntestm) , cp(ntestm) , vamp2 , vfreq2 ,
     .                  cc(ntestm) , vplat , vspot , slant0
c
      common / iqgprm / dvdt   , vbar   , r      , v0      , v      ,
     .                  p      , c      , tgtamp , tgtamp0 , vamp   ,
     .                  vfreq  , drtest , dvtest , ntest   , noise  ,
     .                  dadt   , aamp   , afreq  , daccdf  , daccdr ,
     .                  daccdt , vnoise , vntime , tgdvdt  , slipv  ,
     .                  rbar   , vclut  , aclut  , vdclut  , pclut  ,
     .                  dclut  , cr     , cv     , cc      , cp     ,
     .                  vplat  , vspot  , slant0 , vamp2   , vfreq2
c
c***********************************************************************
c***********************************************************************
c
c           Derived parameters stored in common to ensure
c           that all routine use the same values
c
      real              rfmin  , rfmax  , rsmin  , rsmax  , fc0    ,
     .                  dfc    , dff    , drf    , drc    , drs
c
      common / derprm / rfmin  , rfmax  , rsmin  , rsmax  , fc0    ,
     .                  dfc    , dff    , drf    , drc    , drs
c
c***********************************************************************
c***********************************************************************
c
