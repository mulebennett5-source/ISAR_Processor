C**
C***********************************************************************
C**
      subroutine maritime ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                      npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                      mprat , nchuse , nfocus , radar )
C**
C***********************************************************************
C**
c   The purpose of this routine is to facilitate the development of
c   logical groups of parameters for the RDRTec I-SAR processor.
c
c   This routine is invoked by the command: radar=radar_type
c   where radar_type is a character string.
c
c   This version sets the parameters for emulation of the Telephonics
c   Maritime ISAR Processor.
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
      efghz  = 9.470               !  Center frequency (GHz)
c
      lambda = clight /
     .         ( 1.0E+9 * efghz )  !  Wavelength (m)
c
c   This is 6 times the real chirp rate since it is designed to create
c   0.666 meter pixels from only 1024 points.  This is the core number
c   of range cells after the de-chirp and downsample stages.
c
      br     = 43.94531            !  Chirp rate (MHz/Microsecond)
c
c   This A/D rate has very little flexibility
c
      sampr  = 200.0E+6            !  A/D sampling ( 200 MHz)
c
      strtch = 0                   !  Disable stretch corrections to
c                                  !  save copmutations and reduce the
c                                  !  complexity of the frame generation
c                                  !  routine.  These are negligible
c                                  !  for the Maritime Radar waveform.
c
      nfmt   = 0                   !  8-bit data
c
      ntr    = 684                 !  Samples per pulse after 6-to-1
c                                  !  down-sample
c
c   The value of ncr=2048 gives 0.333 meter pixels at the sub-image.
c   In the actual processor this is accomplished via a 1024-point FFT
c   to get 0.666 meter pixels and an interpolation to 0.333 meters.
c
      ncr    = 2048                !  Range cells per pulse
c
      taywtr = - 30.0              !  Range Taylor weight
c
      taywta = - 30.0              !  Cross-range Taylor weight
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
      dbinc  = - 0.25
c
      nfr    = 512                 !  Range cells in fine image
c
      nfa    = 512                 !  Doppler cells in fine image
c
      mrrat  = 32                  !  Group of range cells
c
      nfskip = 0                   !  Exact time mode
c
      frrate = 2.0                 !  Frame rate 2 Hz but 2 frames
c                                  !  are computed at once to yield
c                                  !  4 frames per second
c
      prf    = 512.0
c
      naskip = 8
c
c   Maritime ISAR: Parameters for Modes
c
c   There are 5 integration time modes - 0.25, 0.5, 1.0, 2.0, 4.0
c   seconds.  The Doppler pixel spacing for the 3 larger integration
c   time modes is 1/T.  However, since the PRF to the image processor is
c   limited to 512 Hz, the shorter integration time modes use a 1 Hz
c   spacing independent of the integration time.
c
c   In the image processor the master particle spacing is fixed at 16 Hz
c   - the same as the sub-image spacing.  Thus, for the 3 shorter
c   integration time modes the number of master particles is 31.  There
c   are 15 for the two second mode and 7 for the four second mode.  This
c   forces the boundary master particles to be zero causing a slight
c   fade to black at the cross-range edges.
c
      if      ( radar(1:10) .eq. 'maritime_1' ) then
c
c   Quarter second mode
c
         rt_mode = 1
c
         tinteg  = 0.25
c
         nspf    = 16
c
         nafill  = 64
c
      else if ( radar(1:10) .eq. 'maritime_2' ) then
c
c   Half second mode
c
         rt_mode = 2
c
         tinteg  = 0.5
c
         nspf    = 32
c
         nafill  = 64
c
      else if ( radar(1:10) .eq. 'maritime_3' ) then
c
c   One second mode
c
         rt_mode = 3
c
         tinteg  = 1.0
c
         nspf    = 64
c
         nafill  = 64
c
      else if ( radar(1:10) .eq. 'maritime_4' ) then
c
c   Two second mode
c
         rt_mode = 4
c
         tinteg  = 2.0
c
         nspf    = 128
c
         nafill  = 128
c
      else if ( radar(1:10) .eq. 'maritime_5' ) then
c
c   Four second mode
c
         rt_mode = 5
c
         tinteg  = 4.0
c
         nspf    = 256
c
         nafill  = 256
c
      endif
c
c   Real-time versions - set all the real-time flags
c
      if ( radar(11:11) .eq. 'r' ) then
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
      if ( radar(11:11) .eq. 's' ) then
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
