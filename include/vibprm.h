c
c***********************************************************************
c
c                          File: vibprm.h
c
c        Purpose:  To store parameters for the SAIC vibration detection
c                  processor which are needed by a variety of routines
c
c***********************************************************************
c
c             Basic data and parameters
c
c     Parameter       Type     FUNCTION
c
c-----------------------------------------------------------------------
c
c     xoffset        Integer   X-offset to region center
c
c     yoffset        Integer   Y-offset to region center
c
c     thresh         Real      Threshhold for echo detection
c
c     topsg0         Real      Top value of local cross section allowed
c
c     tgtsg0         Real      Minimum target strength
c
c     oneside        Real      Minimum ratio of threshold needed for
c                              either right or left echoes
c
c     platdir        Real      Platform direction of travel
c                              (right look = -1.0, left = +1.0)
c
c     fmin           Real      Minimum frequency to be searched (Hz)
c
c     fmax           Real      Maximum frequency to be searched (Hz)
c
c     drslant        Real      Slant range pixel size (ft)
c
c     ft_to_m        Real      Conversion from feet to meters
c
c     overaz         Real      Azimuth oversampling factor
c
c     ftype          Char*80   Data type ('F' = Float, 'I' = IQ4)
c
c-----------------------------------------------------------------------
c
      integer           xoffset , yoffset
c
      real              thresh , topsg0 , tgtsg0  , oneside , platdir ,
     .                  fmin   , fmax   , drslant , ft_to_m , overaz
c
      character         ftype*80
c
c***********************************************************************
c***********************************************************************
c
c   'Top-Ten' list is really up to 100000 long
c
      integer           ten
c
      parameter       ( ten = 100000 )
c
c   Information about the potential vibrator
c
      real              topten(ten,16) , topten_c(ten,16)
c
c   Pixel location of the potential vibrator
c
      integer           ipix(ten) , jpix(ten)
c
c***********************************************************************
c***********************************************************************
c
      common / vibprm / xoffset , yoffset  , thresh  , topsg0 ,
     .                  tgtsg0  , oneside  , platdir , fmin   ,
     .                  fmax    , drslant  , ft_to_m , overaz ,
     .                  topten  , topten_c , ipix    , jpix   , ftype
c
c***********************************************************************
c***********************************************************************
c
