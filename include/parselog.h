c
c                           parselog.h
c
c     ISART log file parse structures
c
c                                 123456789012345678901234567890
c
      character   rminstr*29     / 'Final image minimum range (m)' /
      character   rmaxstr*29     / 'Final image maximum range (m)' /
      character   rincstr*24     / 'Final range sampling (m)' /
      character   dfinestr*26    / 'Fine Doppler sampling (Hz)' /
      character   mrratstr*8     / 'Mrrat  =' /
      character   mpratstr*8     / 'Mprat  =' /
      character   efghzstr*8     / 'Efghz  =' /
      character   prfstr*8       / 'Prf    =' /
      character   samprstr*8     / 'Sampr  =' /
      character   brstr*8        / 'Br     =' /
      character   ntrstr*8       / 'Ntr    =' /
      character   nfrstr*8       / 'Nfr    =' /
      character   nfastr*8       / 'Nfa    =' /
      character   nspfstr*8      / 'Nspf   =' /
      character   dbincstr*8     / 'Dbinc  =' /
      character   curvestr*8     / 'Curve  =' /
      character   overrgstr*8    / 'Overrg =' /
      character   nafillstr*8    / 'Nafill =' /
      character   naskipstr*8    / 'Naskip =' /
      character   nfskipstr*8    / 'Nfskip =' /
      character   contraststr*27 / 'Contrast for frame:       =' /
c
      integer     maxframes
c
      parameter ( maxframes = 9999 )
c
c   Structure to hold information gleaned from the log file
c
      structure / log_dat /
c
         real      rmin 
         real      rmax
         real      rinc
         real      dfine
         integer   nfr
         integer   nfa
         integer   nspf
         real      dbinc
         integer   mrrat
         integer   mprat
         real      efghz
         real      prf
         real      sampr
         real      br
         integer   ntr
         real      overrg
         integer   curve
         integer   nafill
         integer   naskip
         integer   nfskip
         real      contrast(maxframes)
         real      time(maxframes)
         character strcontrast(maxframes)*80
c
      end structure 
c      
