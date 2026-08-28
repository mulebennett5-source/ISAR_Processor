C**
C***********************************************************************
C**
      subroutine norden ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                    npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                    mprat , nchuse , nfocus , radar )
C**
C***********************************************************************
C**
c   The purpose of this routine is to facilitate the development of
c   logical groups of parameters for the RDRTec I-SAR processor.
c
c   This routine is invoked by the command: radar=radar_type
c   where radar_type is a character string.
c
      implicit none
c
c***********************************************************************
c
      integer   ntr , ncr , naskip , spulse , npulse , nspf , nfr ,
     .          nfa , nafill , mprat , mrrat , nchuse , nfocus
c
      real      sampr , prf
c
      character radar*80
c
c***********************************************************************
c
      include  'sarprm.h'
c
      include  'updates.h'
c
      include  'realtime.h'
c
c***********************************************************************
c
c   Set default parameters
c
      call defaults ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                mprat , nchuse , nfocus )
c
c   Maritime ISAR: General parameters
c
      efghz  = 9.70               !  Center frequency (GHz)
c
      lambda = clight /
     .         ( 1.0E+9 * efghz )  !  Wavelength (m)
c
c   This is 4 times the real chirp rate since it is designed to create
c   0.5 meter pixels from only 1024 points.  This is the core number of
c   range cells after the de-chirp and downsample stages.
c
      br     = 43.94531            !  Chirp rate (MHz/Microsecond)
c
c   This A/D rate has very little flexibility
c
      sampr  = 200.0E+6            !  A/D sampling (MHz)
c
      strtch = 0                   !  Disable stretch corrections
c
      nfmt   = 0                   !  8-bit data
c
      ntr    = 1024                !  Samples per pulse
c
c   The value of ncr=2048 gives 0.5 meter pixels at the sub-image.
c
      ncr    = 2048                !  Range cells per pulse
c
      taywtr = - 40.0              !  Range Taylor weight
c
      taywta = - 40.0              !  Cross-range Taylor weight
c
      npass  = 64                  !  PGA Bandwidth
c
      ntaper = 10                  !  PGA taper size
c
      acoefs = 4                   !  GETACC Coefficients
c
      curve  = 2                   !  Use only PGA targets
c
      aghost = 0.25                !  Stabilizer for global motion model
c
      nfr    = 256                 !  Range cells in fine image
c
      nfa    = 256                 !  Doppler cells in fine image
c
      mrrat  = 32                  !  Group of range cells
c
      nfskip = 0                   !  Exact time mode
c
      frrate = 1.5625              !  Frame rate 1.5 Hz but 2 frames
c                                  !  are computed at once to yield
c                                  !  3 frames per second
c
      prf    = 400.0
c
      naskip = 8
c
c   Norden Apy-6 ISAR: Parameters for Modes
c
c   There are 5 integration time modes - 0.32, 0.64, 1.28, 2.56, 5.12
c   seconds.  The Doppler pixel spacing for the 3 larger integration
c   time modes is 1/T.  However, since the PRF to the image processor is
c   limited to 400 Hz, the shorter integration time modes use a .78 Hz
c   spacing independent of the integration time.
c
c   In the image processor the master particle spacing is fixed at 12.5 Hz
c   - the same as the sub-image spacing.  Thus, for the 3 shorter
c   integration time modes the number of master particles is 31.  There
c   are 15 for the two second mode and 7 for the four second mode.  This
c   forces the boundary master particles to be zero causing a slight
c   fade to black at the cross-range edges.
c
      if      ( radar(1:8) .eq. 'norden_1' ) then
c
c   Quarter second mode
c
         rt_mode = 1
c
         nspf    = 16
c
         nafill  = 64
c
      else if ( radar(1:8) .eq. 'norden_2' ) then
c
c   Half second mode
c
         rt_mode = 2
c
         nspf    = 32
c
         nafill  = 64
c
      else if ( radar(1:8).eq. 'norden_3' ) then
c
c   One second mode
c
         rt_mode = 3
c
         nspf    = 64
c
         nafill  = 64
c
      else if ( radar(1:8) .eq. 'norden_4' ) then
c
c   Two second mode
c
         rt_mode = 4
c
         nspf    = 128
c
         nafill  = 128
c
      else if ( radar(1:8) .eq. 'norden_5' ) then
c
c   Four second mode
c
         rt_mode = 5
c
         nspf    = 256
c
         nafill  = 256
c
      endif
c
      tinteg  = float( naskip * nspf ) / prf
c
c   Real-time versions - set all the real-time flags
c
      if ( radar(9:9) .eq. 'r' ) then
c
         mprat   = nafill / 4
c
         rt_nmp  = ( nfa / mprat ) - 1
c
         real_t  = 32767
c
         rt_rco  = 32767
c
         rt_img  = 32767
c
         rt_pga  = 32767
c
         rt_io   = 32767
c
         maxit   = 0        !  PGA Iterations
c
c   I/Q Imbalance parameters
c
         ibar    = 0.0
c
         qbar    = 0.0
c
         iqratio = 1.0
c
         iqmix   = 0.0
c
      else
c
c   Do the case without real-time simplifications
c
         real_t  = 0        !  No general real-time simplifications
c
         rt_rco  = 0        !  No simplifications - range compression
c
         rt_img  = 0        !  No simplifications - image formation
c
         rt_pga  = 0        !  No simplifications in PGA
c
         rt_io   = 0        !  No simplifications in I/O
c
         maxit   = 10       !  PGA Iterations
c
         mprat   = 4        !  Small value of master particle ratio
c
c   Double the frame rate since the real-time version does two frames at
c   once
c
         frrate  = 2.0 * frrate
c
c   I/Q Imbalance parameters
c
         ibar    = 0.0
c
         qbar    = 0.0
c
         iqratio = 1.0
c
         iqmix   = 0.0
c
      endif
c
c   Define survey mode
c
      if ( radar(9:9) .eq. 's' ) then
c
         finemc = 0
c
         alias  = 0
c
         maxit  = 0
c
         curve  = 0
c
         rgwalk = 0.0
c
      endif
c
      return
      end
