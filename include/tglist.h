c 
c***********************************************************************
c
c                          File: tglist.h
c
c***********************************************************************
c
c        Purpose:  To store the composite list of targets detected by
c                  the front-end target detector or the Phase Gradient
c                  Autofocus routine.  These values are used to focus
c                  fine resolution image frames and to estimate target
c                  motion parameters.
c
c***********************************************************************
c
c   Allowed size of the scatterer list
c
      integer           nlist_d
c
      parameter       ( nlist_d = 100000 )
c
c   Number to search through
c
      integer           nlist
c
c-----------------------------------------------------------------------
c
c               Fundamental measured parameters
c
c   time         Real       Time from start of processor of report
c
c   range        Real       Reported range (meter)
c
c   freq         Real       Reported Doppler (Hz)
c
c   accel        Real       Reported acceleration (Hz/second)
c
c   dwdth        Real       Reported Doppler width (Hz)
c
c   snr          Real       Reported Signal-to-Noise Ratio
c
c   source       Character  Reporting source: 't' for front-end target
c                           detector, 'p' for PGA, 'r' from Real-Time
c                           frame generation routine
c
      real              time(nlist_d) , range(nlist_d) , freq(nlist_d)
c
      real              accel(nlist_d) , dwdth(nlist_d) , snr(nlist_d)
c
      character         source(nlist_d)*1
c
c-----------------------------------------------------------------------
c
c               Derived parameters
c
c   iflag           Integer  Flag for edit process
c
c   ipixtg          Integer  Cross-range pixel location of target
c
c   jpixtg          Integer  Range pixel location of target
c
c   omega_valid     Integer  Flag for validity of the estimated rotation
c                            rate
c
c   ndead           Integer  Number of targets which have been
c                            permanently eliminated
c
c   iptr            Integer  Pointer into the circular target buffer
c
c   xtgt, ytgt      Real     Target location (floating point)
c
c   ddot            Real     Rate of change of distance from the cluster
c                            center with time
c
c   tbar            Real     Time for the center of the frame
c
c   fcenuse         Real     Doppler of the scene center for this frame
c
c   rcenuse         Real     Range of the scene center for this frame
c
c   corner          Real     Work array for box containing the target
c
c   a0              Real     Acceleration model coefficient: constant
c
c   at              Real     Acceleration model coefficient: linear in
c                            time
c
c   ar              Real     Acceleration model coefficient: linear in
c                            range
c
c   af              Real     Acceleration model coefficient: linear in
c                            Doppler frequency
c
c   art             Real     Acceleration model coefficient:
c                            proportional to range and time
c
c   aft             Real     Acceleration model coefficient:
c                            proportional to Doppler and time
c
c   ais             Real     3 AIS score values
c
      integer           iflag(nlist_d) , ipixtg(nlist_d) ,
     .                  jpixtg(nlist_d) , omega_valid , ndead , iptr
c
      real              xtgt(nlist_d) , ytgt(nlist_d) , ddot(nlist_d) ,
     .                  tbar , fcenuse , rcenuse , corner(4,2,2) ,
     .                  a0 , at , ar , af , art , aft , ais(3)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c               Put all arrays in common
c
      common / tglist / time , range , freq , accel , dwdth , snr ,
     .                  iptr , source , xtgt , ytgt , iflag , ddot ,
     .                  ipixtg , jpixtg , fcenuse , rcenuse , corner ,
     .                  tbar , ar , af , aft , art , at , ais , nlist
c
c***********************************************************************
c
