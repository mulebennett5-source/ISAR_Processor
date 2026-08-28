c
c***********************************************************************
c***********************************************************************
c
c                          File updates.h
c
c***********************************************************************
c***********************************************************************
c
c     pfrat      Real      Ratio of PRF to original PRF before pre-sum
c
c     vfocus     Real      Platform velocity used in focusing image
c
c     moddat     Integer   Flag to modify raw data
c
c     cmodel     Integer   Level of cancellation model
c
c     ifile2     Char*80   First cancel channel file name
c
c     ifile3     Char*80   Second cancel channel file name
c
c     otype      Integer   Flag for output data type
c
c     multip     Integer   Multi-pass mode control parameter
c
c     editac     Integer   Edit level for acceleration estimates
c
c     presum     Real      Pre-sum exponential weight: 0.0 for none
c
c     length     Real      Expected target length (meters)
c
c     width      Real      Expected target width (meters)
c
c     angle      Real      Expected target angle w.r.t. range
c                          (degrees from broadside orientation)
c
c     ti_con     Real      Parameter for variable integration time
c
c     wrtmcd     Integer   Flag for writing the mocomped raw data
c
c     center     Real      Weight for auto-centering algorithm
c
c     v_corr     Real      Slip velocity correction due to previous
c                          processing errors
c
c     quiet      Integer   Control parameter for reducing screen output
c
c     pgtype     Integer   Control parameter Phase Gradient Autofocus
c
c     tgt_si     Integer   Control parameter for sub-image targets
c
c     a_corr     Real      Additional acceleration
c
c     strip(5)   Real      Strip map parameters ( Range, V_Plat, V_Spot,
c                                                 Altitude, Squint )
c
c     rd_tgt     Integer   Flag for reading target list from disk
c
c     rfdelt     Real      Difference in RF between bands
c
c     nbands     Integer   Number of RF bands for step-chirp systems
c
c     alt_km     Real      Altitude (km)
c
c     n_vadd     Integer   Number of additional values of range-walk to
c                          search each image frame
c
c     d_vadd     Real      Velocity change for range-walk search (m/s)
c
c     qantiz     Integer   Flag for quantizing the range-compressed data
c
c     iqrms      Real      If qantiz non-zero then use this factor times
c                          the observed RMS I/Q to quantize
c
c     reduce     Integer   Level of clutter reduction
c
c     imgslk     Integer   Flag for sidelobe reduction
c
c     accorr     Real      Increment for 'nfocus' settings
c
c     arcorr     Real      Increment for 'nfocus' settings - derivative
c                          with respect to range
c
c     radar      Integer   Radar number
c
c     firsts     Integer   First subimage to consider for frames
c
c     uwb        Integer   Type of sub-image [0=standard, 1=UWB]
c
c     rlook      Integer   Multilooks in range
c
c     flook      Integer   Multilooks in Doppler
c
c     addvib     Integer   Flag to add phase vibration
c
c     addamp     Integer   Flag to add amplitude vibration
c
c     fixamp     Integer   Flag to correct amplitude vibration
c
c     rcentr     Real      Range offset for Mode 10
c
c-----------------------------------------------------------------------
c
      real               prfrat , vfocus , presum , length , width  ,
     .                   angle  , ti_con , center , v_corr , a_corr ,
     .                   strip(5)        , ddt_vspot       , ddt_vplat ,
     .                   omegadot_omega  , accorr , arcorr , rcentr
c
      integer            cmodel , moddat , otype  , multip , editac ,
     .                   wrtmcd , quiet  , pgtype , tgt_si , rd_tgt
c
      integer            maxchan
      parameter        ( maxchan = 3 )
c
      integer            pulse_offset(2,maxchan)
c
c   Cancel channel file names
c
      character          ifile2*80 , ifile3*80
c
      real               rfdelt , alt_km , d_vadd , iqlsb
c
      integer            nbands , n_vadd , qantiz , reduce , imgslk ,
     .                   firsts , uwb , rlook , flook , addvib ,
     .                   addamp , fixamp
c
      common / updates / prfrat , vfocus , moddat , cmodel , ifile2 ,
     .                   ifile3 , otype  , multip , editac , presum ,
     .                   pulse_offset    , length , width  , angle  ,
     .                   ti_con , wrtmcd , center , v_corr , quiet  ,
     .                   pgtype , tgt_si , a_corr , strip  , ddt_vspot ,
     .                   ddt_vplat       , omegadot_omega  , rd_tgt ,
     .                   rfdelt , nbands , alt_km , n_vadd , d_vadd ,
     .                   qantiz , iqlsb  , reduce , imgslk , accorr ,
     .                   arcorr , firsts , uwb    , rlook  , flook  ,
     .                   addvib , addamp , fixamp , rcentr
c
c***********************************************************************
c***********************************************************************
c
