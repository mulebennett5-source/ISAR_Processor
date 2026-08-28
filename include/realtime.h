c
c                         File realtime.h
c
c***********************************************************************
c***********************************************************************
c
c          Control Parameters for Real-Time Approximations
c
c     real_t     Integer   Real-time approximation level - general
c
c     rt_rco     Integer   Real-time approximation - range compression
c
c     rt_img     Integer   Real-time approximation - frame generation
c
c     rt_pga     Integer   Real-time approximation - PGA
c
c     rt_io      Integer   Real-time approximation - I/O routines
c
c-----------------------------------------------------------------------
c
      integer              real_t , rt_rco , rt_img , rt_pga , rt_io
c
c***********************************************************************
c***********************************************************************
c
c          Variables used for Real-Time Approximations
c
c     roll_off   Real      Correction for roll-off of sub-image with
c                          frequency used in IMGEN
c
c     rt_mode    Integer   Real-time mode number
c
c                             1  :  0.25 second
c
c                             2  :  0.5  second
c
c                             3  :  1.0  second
c
c                             4  :  2.0  second
c
c                             5  :  4.0  second
c
c     rt_nmp     Integer   Real-time number of master particles
c
c     rt_max_in  Integer   Real-time limit on targets for GETACC_RT
c
c-----------------------------------------------------------------------
c
      real                 roll_off
c
      integer              rt_mode , rt_nmp , rt_max_in
c
      parameter          ( rt_max_in = 10000 )
c
c
c***********************************************************************
c***********************************************************************
c
c          Variables used for Real-Time Approximations
c
c
c     ibar       Real      Mean I
c
c     qbar       Real      Mean Q
c
c     iqratio    Real      Ratio of I to Q
c
c     iqmix      Real      Fraction of I to mix with Q
c
c-----------------------------------------------------------------------
c
      real                 ibar , qbar , iqratio , iqmix
c
c***********************************************************************
c***********************************************************************
c
c          Common Block to store all Real-Time Variables
c
      common / realtime / real_t , rt_rco , rt_img , rt_pga , rt_io ,
     .                    roll_off , rt_mode , rt_nmp , ibar , qbar ,
     .                    iqratio , iqmix
c
c***********************************************************************
c***********************************************************************
c
