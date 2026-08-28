C**
C***********************************************************************
C**
      subroutine radars ( sampr , ntr , prf , ncr , naskip , spulse ,
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
      if ( radar .eq. 'survey' ) then
c
         nfr    = ntr / 2
c
         overrg = 0.5
c
         naskip = 16
c
         nspf   = 16
c
         nafill = 16
c
         nfa    = 256
c
         nfskip = 8
c
         finemc = 0
c
         alias  = 0
c
         maxit  = 0
c
         curve  = 0
c
      else if ( radar .eq. 'erim_dcs_mte' ) then
c
c   ERIM DCS Radar:  An X-Band system configured for the DARPA MTE
c                    program
c
         efghz  = 9.66                !  Center frequency (GHz)
c
         lambda = clight /
     .            ( 1.0E+9 * efghz )  !  Wavelength (m)
c
         br     = 13.18359            !  Chirp rate (MHz/Microsecond)
c
         sampr  = 45.0E+6             !  A/D sampling (MHz)
c
         prf    = 976.5625            !  Pulse repetition rate (Hz)
c
         nfmt   = 0                   !  8-bit data
c
         moddat = 1                   !  Conjugate input data
c
         ntr    = 2048                !  Samples per pulse
c
      else if ( radar .eq. 'mte_survey' ) then
c
c   Survey mode for the MTE program - very large scale images with all
c   adaptive algorithms disabled
c
         ncr    = 3072                !  Compressed rabge cells
c
         nfa    = 256                 !  Doppler cells in image
c
         nspf   = 16                  !  Sub-images per frame
c
         naskip = 16                  !  Pulses between sub-images
c
         nafill = 16                  !  nspf plus zero-fill
c
         nfr    = 256                 !  Range cells in image
c
         overrg = 0.2                 !  Range oversampling factor
c                                        ( < 1.0 ==> pixel > IPR )
c
         alias  = 0                   !  No target acquisition
c
         finemc = 0                   !  No motion compensation
c
         curve  = 0                   !  No global motion model
c
         maxit  = 0                   !  No auto-focus
c
         acoefs = 0                   !  No accel. coefficients
c
         editac = 0                   !  Do no target editing
c
         outlie = 2.0                 !  Let in almost all targets
c
         multip = 0                   !  No multi-pass
c
         presum = 0.0                 !  Pre-sum filter weight
c
         nfskip = 4                   !  Output every 4th frame
c
         spulse = 1                   !  Start pulse in data file
c
         npulse = 1000000             !  Run until all pulses used
c
         nchuse = 1                   !  Single channel
c
         cmodel = 0                   !  No cancel model needed
c
      else if ( radar .eq. 'mte_target_slow' ) then
c
c   A long integration time mode appropriate for slow MTE targets
c
         ncr    = 4096                !  Compressed range cells
c
         taywtr = - 30.0              !  Taylor weight in range
c
         taywta = - 40.0              !  Taylor weight in cross-range
c
         nfa    = 256                 !  Doppler cells in images
c
         nspf   = 24
c
         naskip = 48
c
         nafill = 32
c
         nfr    = 256
c
         overrg = 3.0
c
         alias  = 1
c
         finemc = 4
c
         alpha0 = 0.5
c
         beta   = 0.1
c
         curve  = 8                   !  Use only PGA targets in global
c                                        motion model and use physical
c                                        constraints
c
         acoefs = 3                   !  First three coefficients -
c                                        up to rotation
c
         editac = 2                   !  Don't do target screening
c                                        for rotation
c
         outlie = 2.5                 !  Let in almost all targets
c
         multip = 3                   !  Multi-pass method
c
         presum = 0.6                 !  Pre-sum filter weight
c
         maxit  = 20                  !  Max. iterations of PGA
c
         nfskip = 1                   !  Output every frame
c
         spulse = 1                   !  Start at first pulse in file
c
         npulse = 1000000             !  Use all pulses in file
c
         nchuse = 1                   !  Single channel
c
         cmodel = 0                   !  Cancel model (0 = no cancel)
c
      else if ( radar .eq. 'mte_target_fast' ) then
c
c   A shorter integration time mode appropriate for fast MTE targets
c
         ncr    = 4096                !  Compressed range cells
c
         taywtr = - 30.0              !  Taylor weight in range
c
         taywta = - 40.0              !  Taylor weight in cross-range
c
         nfa    = 256                 !  Doppler cells in images
c
         nspf   = 24
c
         naskip = 24
c
         nafill = 32
c
         nfr    = 256
c
         overrg = 3.0
c
         alias  = 1
c
         finemc = 4
c
         alpha0 = 0.5
c
         beta   = 0.2
c
         curve  = 8                   !  Use only PGA targets in global
c                                        motion model and use physical
c                                        constraints
c
         acoefs = 3                   !  First three coefficients -
c                                        up to rotation
c
         editac = 2                   !  Don't do target screening
c                                        for rotation
c
         outlie = 2.5                 !  Let in almost all targets
c
         multip = 3                   !  Multi-pass method
c
         presum = 0.6                 !  Pre-sum filter weight
c
         maxit  = 20
c
         nfskip = 2
c
         spulse = 1
c
         npulse = 1000000
c
         nchuse = 1
c
         cmodel = 0
c
      else if ( radar .eq. 'mte_target_straight' ) then
c
c   A very long integration time mode appropriate for MTE targets on
c   straight roads
c
         ncr    = 4096                !  Compressed range cells
c
         taywtr = - 30.0              !  Taylor weight in range
c
         taywta = - 40.0              !  Taylor weight in cross-range
c
         nfa    = 256                 !  Doppler cells in images
c
         nspf   = 48
c
         naskip = 48
c
         nafill = 64
c
         nfr    = 256
c
         overrg = 3.0
c
         alias  = 1
c
         finemc = 4
c
         alpha0 = 0.5
c
         beta   = 0.1
c
         curve  = 8                   !  Use only PGA targets in global
c                                        motion model and use physical
c                                        constraints
c
         acoefs = 3                   !  First three coefficients -
c                                        up to rotation
c
         editac = 2                   !  Don't do target screening
c                                        for rotation
c
         outlie = 2.5                 !  Let in almost all targets
c
         multip = 3                   !  Multi-pass method
c
         presum = 0.6                 !  Pre-sum filter weight
c
         maxit  = 20                  !  Max. iterations of PGA
c
         nfskip = 1                   !  Output every frame
c
         spulse = 1                   !  Start at first pulse in file
c
         npulse = 1000000             !  Use all pulses in file
c
         nchuse = 1                   !  Single channel
c
         cmodel = 0                   !  Cancel model (0 = no cancel)
c
      else if ( radar(1:5) .eq. 'scatr' ) then
c
c   Set default parameters
c
         call defaults ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                   npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                   mprat , nchuse , nfocus )
c
c   NRaD Radar:  An X-Band system configured for the ONR Small Craft
c                Automatic Target Recognition (SCATR) Program
c
         efghz  = 9.25                !  Center frequency (GHz)
c
         lambda = clight /
     .            ( 1.0E+9 * efghz )  !  Wavelength (m)
c
         br     = 50.0                !  Chirp rate (MHz/Microsecond)
c
         sampr  = 5.12E+7             !  A/D sampling (MHz)
c
         prf    = 200.0               !  Pulse repetition rate (Hz)
c
         nfmt   = 0                   !  8-bit data
c
         moddat = 0                   !  Don't conjugate input data
c
         ntr    = 512                 !  Samples per pulse
c
         taywtr = 0.0                 !  Taylor weight applied in READJJ
c
         strtch = 0                   !  Fully compressed, then FFT in r
c
         curve  = 8                   !  Image targets only, advanced
c
         acoefs = 4                   !  a0, at, ar, af
c
         npass  = 129
c
         dbinc  = - 0.2               !  Count up from the noise level
c
         overrg = 1.5
c
         naskip = 16
c
         nspf   = 10
c
         nafill = 16
c
         center = 0.05
c
         taywta = - 30.0
c
         quiet  = 1                   !  Turn off most screen output
c
         if ( radar(7:12) .eq. 'survey' ) then
c
            ncr    = 512              !  Compressed range cells
c
            nfa    = 256              !  Doppler cells in image
c
            nspf   = 16               !  Sub-images per frame
c
            naskip = 16               !  Pulses between sub-images
c
            nafill = 16               !  nspf plus zero-fill
c
            nfr    = 256              !  Range cells in image
c
            overrg = 0.5              !  Range oversampling factor
c                                        ( < 1.0 ==> pixel > IPR )
c
            alias  = 0                !  No target acquisition
c
            finemc = 0                !  No motion compensation
c
            curve  = 0                !  No global motion model
c
            maxit  = 0                !  No auto-focus
c
            rgwalk = 0.0              !  No range walk
c
            acoefs = 0                !  No accel. coefficients
c
            editac = 0                !  Do no target editing
c
            outlie = 2.0              !  Let in almost all targets
c
            multip = 0                !  No multi-pass
c
            presum = 0.0              !  Pre-sum filter weight
c
            nfskip = 4                !  Output every 4th frame
c
            spulse = 1                !  Start pulse in data file
c
            npulse = 1000000          !  Run until all pulses used
c
            nchuse = 1                !  Single channel
c
            cmodel = 0                !  No cancel model needed
c
         endif
c
      else if ( radar(1:6) .eq. 'norden' ) then
c
         call norden ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                 npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                 mprat , nchuse , nfocus , radar )
c

      else if ( radar(1:8) .eq. 'maritime' ) then
c
         call maritime ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                   npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                   mprat , nchuse , nfocus , radar )
c
      else if ( radar(1:5) .eq. 'lamps' ) then
c
         call lamps ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                mprat , nchuse , nfocus , radar )
c
      else
c
         call rcderr ( ' Unknown radar system.' )
c
      endif
c
      return
      end
