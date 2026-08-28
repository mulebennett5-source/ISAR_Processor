c
c***********************************************************************
c
c                              Program ISAR-T
c
c                       An Inverse-SAR Movie Processor
c
c                              Version 3.95
c
c                            January 17, 2018
c
c
c                    Copyright 1994-2017 John R. Bennett
c
c
c    Authors:  John R. Bennett, Kenneth A. Melendez, David S. Brown,
c              Christopher J. Campo, Barton P. Schade, Julie M. Jauregui
c
c              John R. Bennett
c              1945 Briargate Place
c              Escondido, California 92029
c
c              Phone:  858-922-9732
c
c              E-Mail: mulebennett5@gmail.com
c
c***********************************************************************
c
c    Purpose:  To process data from an Inverse Synthetic Aperture
c              Radar of a moving target to create a series of over-
c              lapping movie frames of radar images
c
c    Input:    The input data is raw I/Q signal history which has been
c              corrected in phase to a moving point, the Mo-Comp
c              Point (MCP), by the radar pointing and front-end digital
c              hardware.  The data consists of complex samples
c              spaced uniformly in two time coordinates, fast
c              (range) and slow (azimuth) time.  The fast time is
c              the time within a single pulse of the radar; it is
c              zero at the center of each pulse and its increment is
c              the fundamental A/D sampling interval.  Slow time is
c              zero only at the start of the first pulse; its increment
c              is the pulse repetition interval (1/PRF). 
c
c    Output:   The movie is a series of complex images separated in
c              time by a short interval of time (approx. 0.1 sec) plus
c              certain intermediate files and diagnostic information.
c
c    Assumptions:
c
c              1.  Each radar pulse is assumed to consist of a
c                  constant amplitude signal with a linear change of
c                  frequency with time.  Range compression thus is
c                  accomplished by Fourier transformation of the
c                  data with respect to fast time.
c
c              2.  Each image point is the correlation of the raw
c                  data with the phase function of a unit amplitude
c                  scatterer which moves along an estimate of the path
c                  of scatterers on the object to be imaged.
c
c              3.  If we use the term 'target' for a real scatterer
c                  detected on the object and 'particle' to mean a
c                  hypothetical scatterer used to form an image
c                  point, then the algorithm assumes that the paths
c                  of the particles can be estimated from the paths of
c                  the targets.
c
c              4.  Some terms in the phase function proportional to
c                  the ratio of particle velocity and the speed of
c                  light are neglected.
c
c              5.  The image coordinates are nominal range and
c                  doppler frequency of the particles whose paths
c                  are used to form the image point.
c
c    Algorithm Summary:
c
c              The algorithm consists of three primary processes,
c              Range Compression, Sub-Image Formation, and Frame
c              Generation which are done sequentially and a fourth
c              process, Focus Estimation, which is distributed in
c              nature.  The Focus Estimation process collects
c              information at the output of the other processes and
c              uses it to compute a global model of the object shape,
c              location, and motion which is used to focus the radar
c              images.
c
c              In the Range Compression process a phase correction
c              is applied to cause the target to be approximately
c              centered in the image range and Doppler coordinates
c              and the raw I/Q data is Fourier transformed from
c              fast time to range.
c
c              In the Sub-Image Formation Process a short time record
c              of pulses is used to form a coarse Doppler image.
c
c              In the Focus Estimation process the range-compressed
c              signal history is used to identify point targets on
c              the object to be imaged.  The range, range velocity,
c              and acceleration rate of the targets are estimated with
c              a covariance algorithm.  The average values of range
c              and velocity are fed back to the range compression
c              stage to keep the object centered in range and Doppler.
c              The deviations from these average values are used to
c              estimate the paths of the hypothetical particles used in
c              the Frame Generation Process.
c                  
c              In the Frame Generation process the array of coarse
c              resolution sub-images and the information computed by the
c              the Focus Estimation process is used to compute a fine
c              resolution image frame.
c
c    Program Structure:
c
c              The main program consists of two logical sections - the
c              formal Fortran MAIN program plus a dedicated subroutine,
c              MOVIE.  MAIN takes care of user input and memory
c              allocation.  It then passes control to MOVIE which reads
c              the data files, implements the algorithms and produces
c              the desired output files.  The comments in MAIN are
c              limited to a general discription of the algorithm and the
c              user inputs.  Detailed documentation of the algorithm and
c              the names and uses of the variables are in MOVIE.
c
c***********************************************************************
c
      implicit none
c
c***********************************************************************
c
c   Define all array memory by pointing into a single large array
c
c   Memory sizes:
c
c         nbig   :  The number of words of memory used for real,
c                   complex or integer arrays.
c
c         nbigc  :  Number of complex elements, = nbig / 2
c
c         nbigb  :  Number of bytes used for 8-bit variables
c
      integer       nbig , nbigc , nbigb
c
      parameter   ( nbig  = 320000 * 1024 , nbigc = nbig / 2 ,
     .              nbigb = 4 * nbig )
c
      real          big(nbig)
c
      complex       bigc(nbigc)
c
      integer       bigi(nbig)
c
      character     bigb(nbigb)*1
c
      equivalence ( big , bigc ) , ( big , bigi ) , ( big , bigb )
c
c   Pointers for arrays
c
      integer       imem , irc_raw , icrc_raw , irc , icrc , itr ,
     .              ictr , iac , icac , idd , isbimg , icsbim , iwtr ,
     .              iwtac , iwtaf , iascan , iframe , itgt , ictgt ,
     .              ipgatg , iwork , iiwork , nwork , ibytes
c
c   Make the sub-image array and the image array allocatable
c
      complex, allocatable :: csbimg(:) , cframe(:)
c
      real,    allocatable :: dd(:) , tendat(:,:)
c
      include      'sarprm.h'
c
c***********************************************************************
c
c   Local variables
c
      character     filenm*80
c
      integer       ndot , ntr , ncr , nsr , ntp , nsa , naskip ,
     .              nakeep , nspf , nafill , nfr , nfa , npulse ,
     .              nsubim , frames , ipfile , spulse , nabuff , mrrat ,
     .              mprat , reclda , nchuse , nchmem , mpass , pass1 ,
     .              passn , kframe , ksub , twopwr , nspf_s , nafill_s ,
     .              nafill_a , nfocus , k , ndoti , nspt
c
      real          sampr , prf , prf0 , totdop
c
      logical       prmerr
c
      integer*2     w , h
c
      double precision clight_avg  !  Returns average down to surface
c
c***********************************************************************
c
c   Updates to standard ISART-T processor
c
      include      'updates.h'
c
      include      'motime.h'
c
      include      'realtime.h'
c
      include      'kalman.h'          !  KAM 7-30-98
c
c***********************************************************************
c***********************************************************************
c
      w = 80
c
      h = 500
c
      call wininit ( w , h )  !  Initialize Windows (MS-Windows only)
c
c***********************************************************************
c***********************************************************************
c
c   There are two types of log files.  The first, 'sarnnn.log' is opened
c   as FORTRAN Unit 66 every time the program runs and remains open
c   until the program ends.  All commands are echoed to this file and 
c   some incidental output is also sent to it.  The second type of log
c   file is opened as Unit 7 for each movie created to accept all output
c   related to the case.
c
      call open_sarlog
c
c   Set default parameters
c
      call defaults ( sampr , ntr , prf0 , ncr , naskip , spulse ,
     .                npulse , nspf , nspt , nafill , nfr , nfa ,
     .                mrrat , mprat , nchuse , nfocus )
c
c   Write opening credits to permanent log file and to screen
c
      call isarlogo ( 66 )
      call isarlogo (  6 )
c
      read ( 5 , '(a)' )
c
c   Read parameters and make movies until exit command
c
      do 999 ipfile = 1 , 1000000  !  Large number of cases
c
         prmerr = .false.
c
         call insar ( sampr , ntr , prf0 , ncr , naskip , spulse ,
     .                npulse , nspf , nspt , nafill , nfr , nfa ,
     .                mrrat , mprat , nchuse , nfocus )
c
c   Parameters derived from fundamental parameters
c
         prf    = prf0 / float( nbands )
c
         dtr    = 1.0 / sampr  !  Range sampling interval ( sec. )
         dtp    = 1.0 / prf    !  Pulse repetition interval ( sec. )
c
         clight = clight_avg ( alt_km )
c
         lambda = clight / ( 1.0E+9 * efghz )  !  Wavelength (m)
c
         if ( mod(mode,100) .eq. 1 ) tinteg = float(naskip*nspf) / prf
c
c***********************************************************************
c***********************************************************************
c
c   At this point all parameters have been set from the parameter file
c   or the command parser.  The program now will automatically try to
c   use these parameters to form an ISAR movie.
c
c   Open log file (unit 7)
c
         ndot   = index( ofile , ' ' ) - 1
         filenm = ofile(1:ndot) // '.log'
c
         open ( 7 , file = filenm , form = 'formatted' ,
     .              status = 'unknown' )
c
         call isarlogo ( 7 )
c
         ndoti  = index( ifile , ' ' ) - 1
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Ensure that the fill size is at least as large as the number of
c   sub-images per frame
c
         if ( nafill .lt. nspf ) then
c
            nafill = nspf
c
            write ( 6 , '(/,7x,a40,i8,/)' )
     .         'Parameter Override: nafill set to       ' , nspf
c
            write ( 7 , '(/,7x,a40,i8,/)' )
     .         'Parameter Override: nafill set to       ' , nspf
c
         endif
c
c   Set the Doppler over-sampling ratio to make sure the image does not
c   wrap around Nyquist and that there is no black space
c
c   In 'exact time' mode (nfskip<=0), the two integers NSPF (Number of
c   Sub-images Per Frame) and NAFILL are the number of virtual
c   (interpolated) sub-images used and the number used plus zero-fill.
c   These ratio of these numbers is the rate of over-sampling of the
c   Doppler coordinate.
c
         if ( ( mode .eq. 1 .or. mode .eq. 4 .or. mode .eq. 6
     .          .or. mode .gt. 100 ) .and. nfskip .le. 0 ) then
c
c   Ensure that the frame rate is smaller than the rate of sub-image
c   formation so that the process doesn't run out of data
c
            if ( frrate .gt. prf / float( naskip ) ) then
c
               frrate = prf / float( naskip )
c
               write ( 6 , * ) ' Frame rate reduced to: ' , frrate
               write ( 7 , * ) ' Frame rate reduced to: ' , frrate
c
            endif
c
            nafill = max( nafill , twopwr( nspf ) )
c
c   Doppler pixel size accounting for the integration time and the
c   over-sampling factor
c
            dff    = ( 1.0 / tinteg ) * float( nspf ) /
     .                                  float( nafill )
c
c   Total Doppler width of image
c
            totdop = float( nfa ) * dff
c
c   If image is wider than PRF then modify parameters
c
            if ( totdop .gt. 3 * prf ) then
c
               write ( 6 , * ) ' NSPF and NAFILL before changes: ' ,
     .                           nspf , nafill
               write ( 7 , * ) ' NSPF and NAFILL before changes: ' ,
     .                           nspf , nafill
c
c   Reduce NSPF to expand the image to fill the pixels with only
c   one alias
c
               nspf = ifix( float( nspf ) * prf / totdop )
c
            endif
c
c   If this makes the time step greater than the sub-image time step,
c   then double NSPF and NAFILL until this is not true
c
            do while ( ( tinteg / float( nspf ) ) .gt.
     .                 ( 1.01 * naskip / prf )     )
c
               nspf      = 2 * nspf
c
               nafill    = 2 * nafill
c
            enddo
c
c   Use only even values of NSPF and NAFILL
c
            if ( mod( nspf , 2 )   .eq. 1 ) nspf   = nspf   - 1
c
            if ( mod( nafill , 2 ) .eq. 1 ) nafill = nafill - 1
c
            write ( 6 , * ) ' NSPF and NAFILL changed to:     ' ,
     .                        nspf , nafill
            write ( 7 , * ) ' NSPF and NAFILL changed to:     ' ,
     .                        nspf , nafill
c
c   For the purpose of allocating memory, calculate the number of
c   sub-images per integration time and the equivalent zero-fill value
c
            nspf_s   = nint( tinteg / ( float( naskip ) / prf ) )
c
            nafill_s = ( nafill * nspf_s ) / nspf
c
            if ( rt_img .ne. 0 ) then
c
c   In real-time mode, two frames are computed at once and the sub-image
c   buffer needs to be increased to support this.
c
               nspf_s   = nspf + nint( ( 0.5 / frrate ) /
     .                                 ( float( naskip ) / prf ) )
c
               nafill_a = nafill + nspf_s - nspf
c
               write ( 6 , * ) ' NSPF_S changed to:     ' , nspf_s
               write ( 7 , * ) ' NSPF_S changed to:     ' , nspf_s
c
            endif
c
         else
c
            nspf_s   = nspf
c
            nafill_s = nafill
c
         endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
         call memory ( ntr , ncr , naskip , nspf_s , nafill_s , nfr ,
     .                 nfa , nsr , nsa , nakeep , mrrat , ntp , uwb ,
     .                 npulse , overrg , mode , nabuff , nchuse ,
     .                 nfskip , cmodel , tgt_si , nfocus , nspt , 6 )
c
         call memory ( ntr , ncr , naskip , nspf_s , nafill_s , nfr ,
     .                 nfa , nsr , nsa , nakeep , mrrat , ntp , uwb ,
     .                 npulse , overrg , mode , nabuff , nchuse ,
     .                 nfskip , cmodel , tgt_si , nfocus , nspt , 7 )
c
         nsubim = nint( float( ( iabs( npulse ) - ntp ) / naskip ) )
c
         if ( nfskip .gt. 0 ) then
c
            firsts = min( firsts , nspf )
c 
            frames = max( ( nsubim - nspf ) / nfskip +
     .                    ( ( nspf - firsts ) / 2 ) , 0 )
c
            if ( mode .eq. 7 .or. mode .eq. 8 ) then
c
               if ( nsubim .ge. nspf ) then
c
                  frames = max( 1 + ( nsubim - nspf ) / nfskip , 0 )
c
               else
c
                  frames = 0
c
               endif
c
            endif
c
         else
c
            frames = frrate * ( dtp * float( nsubim * naskip )
     .                          - 0.5 * tinteg )
c
         endif
c
         if ( rt_img .gt. 0 ) frames = 2 * frames ! 2 frames at once
c
         if ( frames .lt. 1 ) then
c
            write ( 6 , * ) ' Not enough data to make one frame'
            write ( 7 , * ) ' Not enough data to make one frame'
c
            prmerr = .true.
c
         else
c
            write ( 6 , '(//,i9,a25)' ) frames ,
     .                             ' Frames will be generated'
            write ( 7 , '(//,i9,a25)' ) frames ,
     .                             ' Frames will be generated'
c
         endif
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Done reading and checking parameters - now write out parameters to
c   the log file and then compute the array memory pointers
c

         call sarst ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                npulse , nspf , nspt , nafill , nfr , nfa ,
     .                mrrat , mprat , nchuse , nfocus , 7 )
c
         imem = 0       !  Set memory pointer to 0
c
c   Work array for raw data - possibly multi-channel
c
         if ( nchuse .gt. 1 ) then
c
            nchmem = nchuse + 1
c
         else
c
            nchmem = 1
c
         endif
c
         icrc_raw = 1 + ( imem + 1 ) / 2
         irc_raw  = ( 2 * icrc_raw ) - 1
         imem     = imem + 2 * ntr * ntp * nchmem
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Work array for range compression
c
         icrc   = 1 + ( imem + 1 ) / 2
         irc    = ( 2 * icrc ) - 1
         imem   = imem + 2 * ncr
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Work array for range compressed signal history - possibly
c   multi-channel
c
         if ( nchuse .gt. 1 .and. cmodel .gt. 0 ) then
c
            nchmem = nchuse + 1
c
         else
c
            nchmem = 1
c
         endif
c
         ictr   = 1 + ( imem + 1 ) / 2
         itr    = ( 2 * ictr ) - 1
         imem   = imem + 2 * ncr * ntp * nchmem
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Work array for Doppler compression
c
         icac   = 1 + ( imem + 1 ) / 2
         iac    = ( 2 * icac ) - 1
c
         if ( rt_img .ne. 0 ) then
c
            imem = imem + 2 * rt_nmp * nafill_a * nfr
c
         else
c
            imem = imem + 2 * ( ( nakeep + nsa ) * ntr +
     .                          8 * nfa + nafill * mrrat )
c
         endif
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Taylor weight array for range compression
c
         iwtr   = imem + 1
         imem   = imem + ntr                   !  Real, ntr x 1
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Taylor weight array for coarse resolution azimuth compression
c
         iwtac  = imem + 1
         imem   = imem + 3 * naskip            !  Real, nsa x 1
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Taylor weight array for fine resolution azimuth compression
c
         iwtaf  = imem + 1
         imem   = imem + nspf                  !  Real, nspf x 1
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   A-Scan array
c
         iascan = imem + 1
         imem   = imem + nfr + nfa             !  Real, nfr + nfa
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Target work array ( Real, 26 * ncr )
c
         ictgt  = 1 + ( imem + 1 ) / 2
         itgt   = ( 2 * ictgt ) - 1
         imem   = imem + 26 * ncr
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   PGA target work array ( Real, 5 * ( nfr + nsr ) * 2 )
c
         ipgatg = imem + 1
         imem   = imem + 5 * ( nfr + nsr ) * 2
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Real work array
c
         iwork  = imem + 1
         nwork  = nfr * nfa * ( 2 + 2 * nfocus )
     .          + 5 * ntr + 2 * ( ncr + nsa )
         nwork  = max( nwork , nabuff * nakeep )
         imem   = imem + nwork              !  Real
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Integer work array
c
         iiwork = imem + 1
         imem   = imem + nwork              !  Integer, 2 x (ncr+nsa)
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c   Byte work array
c
         ibytes = 4 * imem + 1
         imem   = imem + ( nfr * nfa * ( 1 + 2 * nfocus ) + 3 ) / 4
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Work array for acceleration
c
c        idd    = imem + 1
c        imem   = imem + ( 1 + nfa / mprat ) * ( nfr / mrrat )
c    .                                       * nabuff
c
c
c   3-D circular buffer of sub-images
c
c        icsbim = 1 + ( imem + 1 ) / 2
c        isbimg = ( 2 * icsbim ) - 1
c        imem   = imem + 2 * nakeep * nsr * nabuff
c                                      !  Complex, nakeep x nsr x nabuff
c
c   Fine resolution image ( Complex, 2 * nfa * nfr )
c
c        iframe = ( imem / 2 ) + 1
c        imem   = imem + ( 2 + 2 * nfocus ) * 2 * nfr * nfa
c
         if ( allocated(dd) ) deallocate(dd)
         allocate(dd(1:(1+nfa/mprat)*(nfr/mrrat)*nabuff))
c
         if ( allocated(csbimg) ) deallocate(csbimg)
         allocate(csbimg(1:nakeep*nsr*nabuff))
c
         if ( allocated(cframe) ) deallocate(cframe)
         allocate(cframe(1:(2+2*nfocus)*nfr*nfa))
c
         if ( allocated(tendat) ) deallocate(tendat)
         allocate(tendat(npulse,11))
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
         write ( 6 , '(1x,i9,a22)' ) imem , ' Memory words required'
         write ( 7 , '(1x,i9,a22)' ) imem , ' Memory words required'
c
         if ( imem .gt. nbig ) then
c
            write ( 6 , * ) ' Not enough memory.' , imem , nbig
            write ( 7 , * ) ' Not enough memory.' , imem , nbig
c
            prmerr = .true.
c
         endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   If no parameter or memory errors, call the movie processor to
c   compute the ISAR movie
c
         if ( .not. prmerr ) then
c
c-----------------------------------------------------------------------
c
c            Unit 77 - Other Diagnostic Information for the MCA program
c
            filenm = ofile(1:ndot) // '.mca2'
c
            open ( 77 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c-----------------------------------------------------------------------
c
c           Unit 82 - PGA Efficiency Information
c
            filenm = ofile(1:ndot) // '.atr'
c
            open ( 82 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c-----------------------------------------------------------------------
c
c           Unit 83 - Clutter Cancellation Information
c
            if ( nchuse .gt. 1 ) then
c
               filenm = ofile(1:ndot) // '.can'
c
               open ( 83 , file = filenm , form = 'formatted' ,
     .                     status = 'unknown' )
c
            endif
c
c-----------------------------------------------------------------------
c
c           Unit 84 - Motion Compensation Limits
c
            filenm = ofile(1:ndot) // '.mcl'
            open ( 84 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c-----------------------------------------------------------------------
c
c           Unit 85 - Acceleration information
c
            filenm = ofile(1:ndot) // '.gmm'
            open ( 85 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c-----------------------------------------------------------------------
c
c           Unit 86 - Sub-Image array
c
            if ( multip .gt. 0 ) then
c
               filenm = ofile(1:ndot) // '.sbi' 
c
               open ( 86 , file = filenm , form = 'unformatted' ,
     .                     status = 'unknown' , access = 'direct' ,
     .                     recl = reclda( 0 , 2 * nakeep * nsr ) )
c
c              filenm = ofile(1:ndot) // '.sbo' 
c
c              open ( 36 , file = filenm , form = 'unformatted' ,
c    .                     status = 'unknown' , access = 'direct' ,
c    .                     recl = reclda( 0 , 2 * nakeep * nsr ) )
c
            endif
c
c-----------------------------------------------------------------------
c
c           Unit 87 - Motion compensation information
c
            filenm = ofile(1:ndot) // '.moc' 
c
            open ( 87 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown'  )
c
c-----------------------------------------------------------------------
c
            if ( mode .eq. 10 ) then
c
c           Unit 98 - Motion compensation information
c
               filenm = ifile(1:ndoti) // '.ten' 
c
               open ( 98 , file = filenm , form = 'unformatted' ,
     .                     status = 'old' , access = 'direct' ,
     .                     recl = reclda( 0 , 11 ) )
c
               do k = 1 , npulse
                  read( 98 , rec = k+spulse-1 ) tendat(k,1:11)              
               enddo  
c
c              Define relative time
c
               tendat(:,2) = tendat(:,2)-tendat(1,2)
c
               close ( 98 )
c
            endif
c
c-----------------------------------------------------------------------
c
c           Units 88, 99 - Target information
c
            filenm = ofile(1:ndot) // '.tgt'  ! Sub-image targets
c
            open ( 88 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
            filenm = ofile(1:ndot) // '.tgp'  ! PGA targets
c
            open ( 99 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c   Reset the target information by calling addtgt with a negative time
c
            call addtgt ( - 1.0 , 0.0 , 0.0 , 0.0 , 0.0  , 0.0 , 't' )
c
            filenm = ofile(1:ndot) // '.cov'  ! Target Covariances
c
            open ( 100 , file = filenm , form = 'formatted' ,
     .                   status = 'unknown' )
c
            filenm = ofile(1:ndot) // '.gta'  ! Target Covariances
c
            open ( 101 , file = filenm , form = 'unformatted' ,
     .                   status = 'unknown' , access = 'direct' ,
     .                   recl = reclda( 0 , 14 ) )
c
c-----------------------------------------------------------------------
c
c           Unit 89 - A-Scan information
c
            filenm = ofile(1:ndot) // '.as'
c
            open ( 89 , file = filenm , form = 'unformatted' ,
     .                  status = 'unknown' , access = 'direct' ,
     .                  recl = reclda( 0 , nfr + nfa ) )
c
c-----------------------------------------------------------------------
c
c           Unit 90 - Complex strip-map image
c
            if ( mode .gt. 1 .and. mode .lt. 6 ) then
c
               filenm = ofile(1:ndot) // '.ci8'
c
               open ( 90 , file = filenm , form = 'unformatted' ,
     .                status = 'unknown' , access = 'direct' ,
     .                recl = reclda( 0 , 2 * nfr ) )
c
            endif
c
c-----------------------------------------------------------------------
c
c           Unit 92 - File for Scott Musman, IMSI
c           contains fcenuse, rcenuse, length_est, width_est
c
            filenm = ofile(1:ndot) // '.sct'
c
            open ( 92 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c-----------------------------------------------------------------------
c
c           Unit 93 - Write out Mocomped Raw Data
c
            if ( wrtmcd .gt. 0 ) then
c
               filenm = ofile(1:ndot) // '.mcd'
c
               open ( 93 , file = filenm , form = 'unformatted' ,
     .                     status = 'unknown', access = 'direct' ,
     .                     recl = reclda( 0 , 2 * ntr ) )
c
            endif
c
c-----------------------------------------------------------------------
c
c           Unit 94 - Read in Strip Map Data
c
            if ( mode .eq. 6 .or. mode .eq. 7 ) then
c
               filenm = ifile(1:index(ifile,' ')-1) // '.smd'
c
               open ( 94 , file = filenm , form = 'unformatted' ,
     .                     status = 'unknown', access = 'direct' ,
     .                     recl = reclda( 0 , 5 ) )
c
            endif  
c
c-----------------------------------------------------------------------
c
c           Unit 95 - File for Rotation Rate Estimates
c
            if ( quiet .gt. 1 ) then
c
               filenm = ofile(1:ndot) // '.omg'
c
               open ( 95 , file = filenm , form = 'formatted' ,
     .                     status = 'unknown' )                      
c
            endif
c
c-----------------------------------------------------------------------
c
c           Unit 96 - Plot File for Adaptive Phase Estimates
c
            if ( nbands .gt. 1 ) then
c
               filenm = ofile(1:ndot) // '.plt'
c
               open ( 96 , file = filenm , form = 'formatted' ,
     .                     status = 'unknown' )                      
c
               write ( 96 , * )
c
               write ( 96 , * )
c
               if ( nbands .eq. 2 ) then
c
                  write ( 96 , * ) '1,2,2'
c
               else
c
                  write ( 96 , * ) '1,2,3'
c
               endif
c
               write ( 96 , * ) 'F2F Phase Correction'
c
               write ( 96 , * ) 'Sub-image Number'
c
               write ( 96 , * ) 'F1 to F2'
c
               if ( nbands .eq. 2 ) then
c
                  write ( 96 , * ) 'F1 to F2'
c
               else
c
                  write ( 96 , * ) 'F2 to F3'
c
               endif
c
            endif
c
c-----------------------------------------------------------------------
c
c           Unit 97 - Target information
c
            filenm = ofile(1:ndot) // '.tgtout'
c
            open ( 97 , file = filenm , form = 'formatted' ,
     .                  status = 'unknown' )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Finally, generate an ISAR movie
c
            if ( multip .eq. 0 ) then
c
               pass1 = 0
               passn = 0
c
            else if ( multip .eq. 1 ) then
c
               pass1 = 1
               passn = 1
c
            else if ( multip .eq. 2 ) then
c
               pass1 = 2
               passn = 2
c
            else if ( multip .eq. 3 ) then
c
               pass1 = 1
               passn = 2
c
            endif
c
            call timing ( ksub , kframe , big(iwork) , nwork , frrate ,
     .                    0 )
c
            do mpass = pass1 , passn
c
               if ( mode .lt. 100 ) then
c
c   Original version
c
                  call movie  ( ntr , ncr , nsr , ntp , nsa , naskip ,
     .                          nakeep , spulse , npulse , nspf , nspt ,
     .                          nafill , nfr , nfa , big(irc_raw) ,
     .                          bigc(icrc_raw) , big(irc) , bigc(icrc) ,
     .                          big(itr) , bigc(ictr) , big(iac) ,
     .                          bigc(icac) , csbimg , big(iwtr) ,
     .                          big(iwtac) , big(iwtaf) , big(iascan) ,
     .                          cframe , bigb(ibytes) ,
     .                          big(itgt) , bigc(ictgt) , big(ipgatg) ,
     .                          dd , big(iwork) , bigi(iiwork) ,
     .                          nwork , nabuff , mrrat , mprat ,
     .                          nchuse , mpass , ksub , kframe ,
     .                          nafill_a , nfocus , tendat )
c
               else
c
c   Special version for step-chirp data
c
                  call movie0 ( ntr , ncr , nsr , ntp , nsa , naskip ,
     .                          nakeep , spulse , npulse , nspf , nspt ,
     .                          nafill , nfr , nfa , big(irc_raw) ,
     .                          bigc(icrc_raw) , big(irc) , bigc(icrc) ,
     .                          big(itr) , bigc(ictr) , big(iac) ,
     .                          bigc(icac) , csbimg , big(iwtr) ,
     .                          big(iwtac) , big(iwtaf) , big(iascan) ,
     .                          cframe , bigb(ibytes) ,
     .                          bigc(ictgt) , big(ipgatg) , dd ,
     .                          big(iwork) , bigi(iiwork) , nwork ,
     .                          nabuff , mrrat , mprat , nchuse ,
     .                          mpass , ksub , kframe , nafill_a ,
     .                          nfocus )
c
               endif
c
               if ( mpass .lt. 0 ) go to 111
c
c   Mo-comp is only done on first pass
c
               if ( mpass .eq. pass1 ) close ( 87 )   !  Mo-comp file
c
            enddo
c
  111       continue

            call timing ( ksub , kframe , big(iwork) , nwork , frrate ,
     .                    1 )
c
            close ( 77 )   !  Other diagnostic information file
            close ( 82 )   !  PGA Efficiency file
            close ( 83 )   !  Clutter cancellation file
            close ( 85 )   !  Acceleration file
            close ( 86 )   !  Sub-image file
c           close ( 36 )   !  Sub-image file
            close ( 88 )   !  Target file
            close ( 89 )   !  A-Scan file
            close ( 90 )   !  Strip-map image file
            close ( 92 )   !  Scott Musman ATR info
            close ( 93 )   !  Mo-comped data file
            close ( 94 )   !  Strip map data file
            close ( 95 )   !  Rotation Rate Information
            close ( 96 )   !  Plot File - Adaptive Phase Estimate
            close ( 97 )   !  Output targets if rd_tgt .ne. 0
            close ( 99 )   !  Target file for PGA
            close ( 100 )  !  Covariance text file
            close ( 101 )  !  Covariance binary file
c
            close (  7 )   !  Log file for this case
c
         endif
c
  999 continue             !  End of loop over cases
c
 1000 close ( 66 )         !  Global log file
c
      stop                 !  End of program
      end
C**
C***********************************************************************
C**
      subroutine open_sarlog
C**
C***********************************************************************
C**
c
c   Open a log file for ISAR-T incrementing the number until a new file
c   name can be opened.
c
      implicit none
c
      character logfil*10
c
      character num*3
c
      integer   i
c
      do i = 1 , 100
c
         num = '000'
c
         if ( i .lt. 10 )
     .        write ( num(3:3) , '(i1)' ) i
c
         if ( i .ge. 10 .and. i .lt. 100 )
     .        write ( num(2:3) , '(i2)' ) i
c
         if ( i .ge. 100 .and. i .lt. 1000 )
     .        write ( num(1:3) , '(i3)' ) i
c
         logfil = 'sar' // num // '.log'
c
         open ( 66 , file = logfil , form = 'formatted' ,
     .               status = 'NEW' , shared , err = 10 )
c
         return
c
   10    continue
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine insar ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                   npulse , nspf , nspt , nafill , nfr , nfa ,
     .                   mrrat , mprat , nchuse , nfocus )
C**
C***********************************************************************
C**
      implicit none
c
c***********************************************************************
c
c   Program parameters which can be set by any user
c
      integer   ntr , ncr , naskip , spulse , npulse , nspf , nspt ,
     .          nfr , nfa , nafill , nfocus
c
      real      sampr , prf
c
c***********************************************************************
c
c   Program parameters which can only be set by command parser
c
      include  'sarprm.h'
c
      integer   mrrat , mprat
c
c********************** Local Variables ********************************
c
c
      integer   ifirst , iline , ncmd
c
      character tmpstr*80 , lhs*80 , rhs*80
c
      character prompt*3 , ext*3
c
c******************* Program control variables *************************
c
      logical  eflag
c
c   These variables determine the allowable commands and control the
c   program's response to errors detected in the user input or from the
c   system.  Their usage is as follows:
c
c      EFLAG   :  True if the command string has an equal sign in it;
c                 thus, it is an attempt to set a parameter rather than
c                 to execute a word command.
c
c*********************** SAR Commands **********************************
c
c   The user enters commands consisting of strings of characters
c   separated by blanks.  If the string contains an equal sign, it is
c   assumed that a parameter is being set.  If there is not an equal
c   sign, the string is assumed to be a 'word' command.
c
c   Set up the arrays of legal 'word' commands and parameter names
c
      integer     comnum , maxkey
c
      parameter ( comnum = 11 , maxkey = 131 )
c
      character   coms(comnum)*80 , keys(maxkey)*80
c
c   Work variables and names of functions used
c
      integer     nkey , comchk , lenrhs , len , geti4 , getr4 ,
     .            i4val , nchuse
c
      real        r4val
c
      character   hlpfil*80 , hlpcmd*80 , radar*80
c
c***********************************************************************
c
      include    'updates.h'   !  Updates after main release (Ver. 1.5)
c
      include    'realtime.h'  !  Real-time parameters
c
      include    'kalman.h'    !  Kalman Filter Parameters KAM 7-30-98
c
c***********************************************************************
c
      save keys , coms , ifirst , prompt , ext
c
      data ifirst / 0 /
c
      data keys / 'ifile'  , 'ofile'  , 'nfmt'   , 'fixiq'  , 'beta'   ,
     .            'alias'  , 'ntr'    , 'ncr'    , 'naskip' , 'spulse' ,
     .            'npulse' , 'nspf'   , 'nfr'    , 'nfa'    , 'efghz'  ,
     .            'br'     , 'sampr'  , 'prf'    , 'taywtr' , 'taywta' ,
     .            'mprat'  , 'mrrat'  , 'finemc' , 'maxit'  , 'nrkeep' ,
     .            'npass'  , 'ntaper' , 'nrcent' , 'curve'  , 'nlag'   ,
     .            'nfskip' , 'dbinc'  , 'nlocal' , 'pixbar' , 'outlie' ,
     .            'notch'  , 'vnotch' , 'dvntch' , 'rnotch' , 'drntch' ,
     .            'alpha0' , 'overrg' , 'nafill' , 'snrmin' , 'acoefs' ,
     .            'aghost' , 'rgwalk' , 'frrate' , 'tinteg' , 'mode'   ,
     .            'strtch' , 'prfrat' , 'otype'  , 'cmodel' , 'nchuse' ,
     .            'ifile2' , 'ifile3' , 'moddat' , 'radar'  , 'multip' ,
     .            'editac' , 'presum' , 'length' , 'width'  , 'angle'  ,
     .            'ti_con' , 'color'  , 'real_t' , 'rt_rco' , 'rt_img' ,
     .            'rt_pga' , 'rt_io'  , 'dokalm' , 'slbkil' , 'v_corr' ,
     .            'quiet'  , 'pgtype' , 'tgt_si' , 'a_corr' , 'rd_tgt' ,
     .            'nbands' , 'rfdelt' , 'alt_km' , 'imgslk' , 'reduce' ,
     .            'vfocus' , 'wrtmcd' , 'nfocus' , 'accorr' , 'arcorr' ,
     .            'center' , 'n_vadd' , 'd_vadd' , 'qantiz' , 'iqlsb'  ,
     .            'firsts' , 'uwb'    , 'rlook'  , 'flook'  , 'addvib' ,
     .            'addamp' , 'fixamp' , 'rcentr' , 'nspt'   ,
     .            'ntest'  , 'rbar'   , 'vbar'   , 'dvdt'   , 'dadt'   ,
     .            'aamp'   , 'afreq'  , 'daccdf' , 'daccdr' , 'daccdt' ,
     .            'drtest' , 'dvtest' , 'vamp'   , 'vfreq'  , 'vamp2'  ,
     .            'vfreq2' , 'noise'  , 'vnoise' , 'vntime' , 'slipv'  ,
     .            'vclut'  , 'vdclut' , 'aclut'  , 'pclut'  , 'dclut'  ,
     .            'vplat'  , 'slant0' /
c
      data coms / ' ' , 'go' , 'bye' , 'exit' , 'stat' , 'reset' ,
     .            'help' , 'man' , 'quit' , 'realtime' , 'norealtime' /
c
c   These commands are:
c
c     ' '           :  Blank - no effect
c     'go'          :  Create movie
c     'exit'        :  Exits program
c     'quit'        :  Exits program
c     'bye'         :  Synonym for 'exit'
c     'stat'        :  Displays the values of parameters
c     'reset'       :  Resets parameters to their default values
c     'realtime'    :  Use all real-time approximations
c     'norealtime'  :  Use no real-time approximations
c     'help'        :  Accesses help file if available
c     'man'         :  Synonym for 'help'
c
c***********************************************************************
c
      if ( ifirst .eq. 0 ) then
c
         call sarst ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                npulse , nspf , nspt , nafill , nfr , nfa ,
     .                mrrat , mprat , nchuse , nfocus , 6 )
c
         prompt = 'SAR'
         ext    = 'sar'
c
      endif
c
c   Execute commands forever
c
      do 1000 iline = 1 , 64000
c
c   Given a prompt, a file extension and a counter ( PROMPT , EXT ,
c   IFIRST ), return the next command from the user ( TMPSTR ).  If the
c   command has an equal sign in it, set EFLAG=.TRUE. and break the
c   command into the parts on the left and right sides of the '='
c   ( LHS , RHS ).
c
      call getcmd ( prompt , ext , ifirst , tmpstr , lhs , rhs , eflag ,
     .              0 )
c
      if ( tmpstr .ne. ' ' )
     .             write ( 66 , * ) prompt , ' > ' , tmpstr(1:40)
c
c   Handle the commands strings with equal signs
c
      if ( eflag ) then
c
c   Parse command string to obtain input variables
c
         lenrhs = len( rhs ) 
c
         if ( comchk( lhs , keys , maxkey , nkey ) .ne. 0 ) then
c
            call rcderr ( 'Bad keyword.' )
c
c   Choose the default parameters
c
         else if ( lhs .eq. 'radar' ) then
c
            radar = rhs
c
            call radars ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                    npulse , nspf , nafill , nfr , nfa , mrrat ,
     .                    mprat , nchuse , nfocus , radar )
c
         else if ( lhs .eq. 'ifile' ) then
c
            ifile = rhs
c
         else if ( lhs .eq. 'ifile2' ) then
c
            ifile2 = rhs
c
         else if ( lhs .eq. 'ifile3' ) then
c
            ifile3 = rhs
c
         else if ( lhs .eq. 'ofile' ) then
c
            ofile = rhs
c
         else if ( lhs .eq. 'nfmt' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nfmt = i4val
            else
               call rcderr ( 'Bad nfmt.' )
            endif
c
         else if ( lhs .eq. 'strtch' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               strtch = i4val
            else
               call rcderr ( 'Bad strtch.' )
            endif
c
         else if ( lhs .eq. 'mode' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               mode = i4val
c
               if ( mode .eq. 2 ) then
c
                  rgwalk = 0.0
                  write (  6 , * ) ' Strip-map mode: Rgwalk = 0'
                  write ( 66 , * ) ' Strip-map mode: Rgwalk = 0'
c
               endif
c
            else
               call rcderr ( 'Bad mode.' )
            endif
c
         else if ( lhs .eq. 'fixiq' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               fixiq = i4val
            else
               call rcderr ( 'Bad fixiq.' )
            endif
c
         else if ( lhs .eq. 'beta' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               beta = r4val
            else
               call rcderr ( 'Bad beta.' )
            endif
c
         else if ( lhs .eq. 'alias' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               alias = i4val
            else
               call rcderr ( 'Bad alias.' )
            endif
c
         else if ( lhs .eq. 'ntr' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               ntr = i4val
            else
               call rcderr ( 'Bad ntr.' )
            endif
c
         else if ( lhs .eq. 'ncr' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               ncr = i4val
            else
               call rcderr ( 'Bad ncr.' )
            endif
c
         else if ( lhs .eq. 'naskip' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               naskip = i4val
            else
               call rcderr ( 'Bad naskip.' )
            endif
c
         else if ( lhs .eq. 'spulse' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               spulse = i4val
            else
               call rcderr ( 'Bad spulse.' )
            endif
c
         else if ( lhs .eq. 'npulse' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               npulse = i4val
            else
               call rcderr ( 'Bad npulse.' )
            endif
c
         else if ( lhs .eq. 'nspf' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nspf = i4val
            else
               call rcderr ( 'Bad nspf.' )
            endif
c
         else if ( lhs .eq. 'nspt' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nspt = i4val
            else
               call rcderr ( 'Bad nspt.' )
            endif
c
         else if ( lhs .eq. 'nfr' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nfr = i4val
            else
               call rcderr ( 'Bad nfr.' )
            endif
c
         else if ( lhs .eq. 'nfa' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nfa = i4val
            else
               call rcderr ( 'Bad nfa.' )
            endif
c
         else if ( lhs .eq. 'efghz' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               efghz  = r4val
            else
               call rcderr ( 'Bad efghz.' )
            endif
c
         else if ( lhs .eq. 'br' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               br = r4val
            else
               call rcderr ( 'Bad br.' )
            endif
c
         else if ( lhs .eq. 'sampr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               sampr = r4val
            else
               call rcderr ( 'Bad sampr.' )
            endif
c
         else if ( lhs .eq. 'prf' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               prf = r4val
            else
               call rcderr ( 'Bad prf.' )
            endif
c
         else if ( lhs .eq. 'prfrat' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               prfrat = r4val
               write (  6 , * ) ' Prfrat changed to: ' , prfrat
               write ( 66 , * ) ' Prfrat changed to: ' , prfrat
            else
               call rcderr ( 'Bad prfrat.' )
            endif
c
         else if ( lhs .eq. 'presum' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               presum = r4val
            else
               call rcderr ( 'Bad presum.' )
            endif
c
         else if ( lhs .eq. 'taywtr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               taywtr = r4val
            else
               call rcderr ( 'Bad taywtr.' )
            endif
c
         else if ( lhs .eq. 'taywta' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               taywta = r4val
            else
               call rcderr ( 'Bad taywta.' )
            endif
c
         else if ( lhs .eq. 'mprat' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               mprat = i4val
            else
               call rcderr ( 'Bad mprat.' )
            endif
c
         else if ( lhs .eq. 'mrrat' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               mrrat = i4val
            else
               call rcderr ( 'Bad mrrat.' )
            endif
c
         else if ( lhs .eq. 'finemc' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               if ( abs(i4val) .le. 8 ) then
                  finemc = i4val
               else
                  call rcderr ( 'Bad finemc.' )
               endif
            endif
c
         else if ( lhs .eq. 'alpha0' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               alpha0 = r4val
            else
               call rcderr ( 'Bad alpha0.' )
            endif
c
         else if ( lhs .eq. 'notch' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               notch = i4val
            else
               call rcderr ( 'Bad notch.' )
            endif
c
         else if ( lhs .eq. 'vnotch' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vnotch = r4val
            else
               call rcderr ( 'Bad vnotch.' )
            endif
c
         else if ( lhs .eq. 'dvntch' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               dvntch = r4val
            else
               call rcderr ( 'Bad dvntch.' )
            endif
c
         else if ( lhs .eq. 'rnotch' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               rnotch = r4val
            else
               call rcderr ( 'Bad rnotch.' )
            endif
c
         else if ( lhs .eq. 'drntch' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               drntch = r4val
            else
               call rcderr ( 'Bad drntch.' )
            endif
c
         else if ( lhs .eq. 'nrkeep' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nrkeep = i4val
            else
               call rcderr ( 'Bad nrkeep.' )
            endif
c
         else if ( lhs .eq. 'npass' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               npass = i4val
            else
               call rcderr ( 'Bad npass.' )
            endif
c
         else if ( lhs .eq. 'ntaper' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               ntaper = i4val
            else
               call rcderr ( 'Bad ntaper.' )
            endif
c
         else if ( lhs .eq. 'maxit' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               maxit = i4val
            else
               call rcderr ( 'Bad maxit.' )
            endif
c
         else if ( lhs .eq. 'slbkil' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               slbkil = i4val
            else
               call rcderr ( 'Bad slbkil.' )
            endif
c
         else if ( lhs .eq. 'pgtype' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               pgtype = i4val
            else
               call rcderr ( 'Bad pgtype.' )
            endif
c
         else if ( lhs .eq. 'tgt_si' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               tgt_si = i4val
            else
               call rcderr ( 'Bad tgt_si.' )
            endif
c
         else if ( lhs .eq. 'curve' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               curve = i4val
            else
               call rcderr ( 'Bad curve.' )
            endif
c
         else if ( lhs .eq. 'center' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               center = r4val
            else
               call rcderr ( 'Bad center.' )
            endif
c
         else if ( lhs .eq. 'length' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               length = r4val
            else
               call rcderr ( 'Bad length.' )
            endif
c
         else if ( lhs .eq. 'width' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               width = r4val
            else
               call rcderr ( 'Bad width.' )
            endif
c
         else if ( lhs .eq. 'angle' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               angle = r4val
            else
               call rcderr ( 'Bad angle.' )
            endif
c
         else if ( lhs .eq. 'ti_con' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               ti_con = r4val
            else
               call rcderr ( 'Bad ti_con.' )
            endif
c
         else if ( lhs .eq. 'overrg' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               overrg = r4val
            else
               call rcderr ( 'Bad overrg.' )
            endif
c
         else if ( lhs .eq. 'nafill' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nafill = i4val
            else
               call rcderr ( 'Bad nafill.' )
            endif
c
         else if ( lhs .eq. 'nlag' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nlag = i4val
            else
               call rcderr ( 'Bad nlag.' )
            endif
c
         else if ( lhs .eq. 'nrcent' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nrcent = i4val
            else
               call rcderr ( 'Bad nrcent.' )
            endif
c
         else if ( lhs .eq. 'nfskip' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nfskip = i4val
            else
               call rcderr ( 'Bad nfskip.' )
            endif
c
         else if ( lhs .eq. 'dbinc' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               dbinc = r4val
            else
               call rcderr ( 'Bad dbinc.' )
            endif
c
         else if ( lhs .eq. 'snrmin' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               snrmin = r4val
            else
               call rcderr ( 'Bad snrmin.' )
            endif
c
         else if ( lhs .eq. 'outlie' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               outlie = r4val
            else
               call rcderr ( 'Bad outlie.' )
            endif
c
         else if ( lhs .eq. 'aghost' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               aghost = r4val
            else
               call rcderr ( 'Bad aghost.' )
            endif
c
         else if ( lhs .eq. 'rgwalk' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               rgwalk = r4val
            else
               call rcderr ( 'Bad rgwalk.' )
            endif
c
         else if ( lhs .eq. 'acoefs' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               acoefs = i4val
            else
               call rcderr ( 'Bad acoefs.' )
            endif
c
         else if ( lhs .eq. 'editac' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               editac = i4val
            else
               call rcderr ( 'Bad editac.' )
            endif
c
         else if ( lhs .eq. 'nlocal' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nlocal = i4val
            else
               call rcderr ( 'Bad nlocal.' )
            endif
c
         else if ( lhs .eq. 'pixbar' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               pixbar = i4val
            else
               call rcderr ( 'Bad pixbar.' )
            endif
c
         else if ( lhs .eq. 'frrate' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               frrate = r4val
            else
               call rcderr ( 'Bad frrate.' )
            endif
c
         else if ( lhs .eq. 'tinteg' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               tinteg = r4val
            else
               call rcderr ( 'Bad tinteg.' )
            endif
c
         else if ( lhs .eq. 'nchuse' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nchuse = i4val
            else
               call rcderr ( 'Bad nchuse.' )
            endif
c
         else if ( lhs .eq. 'cmodel' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               cmodel = i4val
            else
               call rcderr ( 'Bad cmodel.' )
            endif
c
         else if ( lhs .eq. 'moddat' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               moddat = i4val
            else
               call rcderr ( 'Bad moddat.' )
            endif
c
         else if ( lhs .eq. 'otype' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               otype = i4val
            else
               call rcderr ( 'Bad otype.' )
            endif
c
         else if ( lhs .eq. 'color' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               color = i4val
            else
               call rcderr ( 'Bad color.' )
            endif
c
         else if ( lhs .eq. 'multip' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               multip = i4val
            else
               call rcderr ( 'Bad multip.' )
            endif
c
         else if ( lhs .eq. 'v_corr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               v_corr = r4val
            else
               call rcderr ( 'Bad v_corr.' )
            endif
c
         else if ( lhs .eq. 'a_corr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               a_corr = r4val
            else
               call rcderr ( 'Bad a_corr.' )
            endif
c
         else if ( lhs .eq. 'quiet' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               quiet = i4val
            else
               call rcderr ( 'Bad quiet.' )
            endif
c
         else if ( lhs .eq. 'rd_tgt' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               rd_tgt = i4val
            else
               call rcderr ( 'Bad rd_tgt.' )
            endif
c
         else if ( lhs .eq. 'nbands' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nbands = i4val
            else
               call rcderr ( 'Bad nbands.' )
            endif
c
         else if ( lhs .eq. 'rfdelt' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               rfdelt = r4val * 1.0E+6  !  MHz
            else
               call rcderr ( 'Bad rfdelt.' )
            endif
c
         else if ( lhs .eq. 'alt_km' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               alt_km = r4val
            else
               call rcderr ( 'Bad alt_km.' )
            endif
c
         else if ( lhs .eq. 'nfocus' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               nfocus = i4val
            else
               call rcderr ( 'Bad nfocus.' )
            endif
c
         else if ( lhs .eq. 'accorr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               accorr = r4val
            else
               call rcderr ( 'Bad accorr.' )
            endif
c
         else if ( lhs .eq. 'arcorr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               arcorr = r4val
            else
               call rcderr ( 'Bad arcorr.' )
            endif
c
         else if ( lhs .eq. 'n_vadd' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               n_vadd = i4val
            else
               call rcderr ( 'Bad n_vadd.' )
            endif
c
         else if ( lhs .eq. 'd_vadd' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               d_vadd = r4val
            else
               call rcderr ( 'Bad d_vadd.' )
            endif
c
         else if ( lhs .eq. 'qantiz' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               qantiz = i4val
            else
               call rcderr ( 'Bad qantiz.' )
            endif
c
         else if ( lhs .eq. 'iqlsb' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               iqlsb = r4val
            else
               call rcderr ( 'Bad iqlsb.' )
            endif
c
         else if ( lhs .eq. 'imgslk' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               imgslk = i4val
            else
               call rcderr ( 'Bad imgslk.' )
            endif
c
         else if ( lhs .eq. 'reduce' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               reduce = i4val
            else
               call rcderr ( 'Bad reduce.' )
            endif
c
         else if ( lhs .eq. 'firsts' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               firsts = i4val
            else
               call rcderr ( 'Bad firsts.' )
            endif
c
         else if ( lhs .eq. 'uwb' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               uwb    = i4val
            else
               call rcderr ( 'Bad uwb.' )
            endif
c
         else if ( lhs .eq. 'rlook' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               rlook   = i4val
            else
               call rcderr ( 'Bad rlook.' )
            endif
c
         else if ( lhs .eq. 'flook' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               flook   = i4val
            else
               call rcderr ( 'Bad flook.' )
            endif
c
         else if ( lhs .eq. 'addvib' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               addvib  = i4val
            else
               call rcderr ( 'Bad addvib.' )
            endif
c
         else if ( lhs .eq. 'addamp' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               addamp  = i4val
            else
               call rcderr ( 'Bad addamp.' )
            endif
c
         else if ( lhs .eq. 'fixamp' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               fixamp  = i4val
            else
               call rcderr ( 'Bad fixamp.' )
            endif
c
         else if ( lhs .eq. 'rcentr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               rcentr  = r4val
            else
               call rcderr ( 'Bad rcentr.' )
            endif
c
c-----------------------------------------------------------------------
c------------------------ Real-time Parameters -------------------------
c
         else if ( lhs .eq. 'real_t' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               real_t = i4val
            else
               call rcderr ( 'Bad real_t.' )
            endif
c
         else if ( lhs .eq. 'rt_rco' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               rt_rco = i4val
            else
               call rcderr ( 'Bad rt_rco.' )
            endif
c
         else if ( lhs .eq. 'rt_img' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               rt_img = i4val
            else
               call rcderr ( 'Bad rt_img.' )
            endif
c
         else if ( lhs .eq. 'rt_pga' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               rt_pga = i4val
            else
               call rcderr ( 'Bad rt_pga.' )
            endif
c
         else if ( lhs .eq. 'rt_io' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               rt_io = i4val
            else
               call rcderr ( 'Bad rt_io.' )
            endif
c
c-----------------------------------------------------------------------
c-------------------------- Kalman Filter Parameters -------------------
c
         else if ( lhs .eq. 'dokalm' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               dokalm = i4val
            else
               call rcderr ( 'Bad dokalm.' )
            endif
c
c-------------------------- IQGEN Parameters ---------------------------
c
         else if ( lhs .eq. 'ntest' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               ntest = i4val
            else
               call rcderr ( 'Bad ntest.' )
            endif
c
         else if ( lhs .eq. 'rbar' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               rbar = r4val
            else
               call rcderr ( 'Bad rbar.' )
            endif
c
         else if ( lhs .eq. 'vbar' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vbar = r4val
            else
               call rcderr ( 'Bad vbar.' )
            endif
c
         else if ( lhs .eq. 'dvdt' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               dvdt = r4val
            else
               call rcderr ( 'Bad dvdt.' )
            endif
c
         else if ( lhs .eq. 'dadt' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               dadt = r4val
            else
               call rcderr ( 'Bad dadt.' )
            endif
c
         else if ( lhs .eq. 'aamp' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               aamp = r4val
            else
               call rcderr ( 'Bad aamp.' )
            endif
c
         else if ( lhs .eq. 'afreq' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               afreq = r4val
            else
               call rcderr ( 'Bad afreq.' )
            endif
c
         else if ( lhs .eq. 'daccdf' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               daccdf = r4val
            else
               call rcderr ( 'Bad daccdf.' )
            endif
c
         else if ( lhs .eq. 'daccdr' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               daccdr = r4val
            else
               call rcderr ( 'Bad daccdr.' )
            endif
c
         else if ( lhs .eq. 'daccdt' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               daccdt = r4val
            else
               call rcderr ( 'Bad daccdt.' )
            endif
c
         else if ( lhs .eq. 'drtest' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               drtest = r4val
            else
               call rcderr ( 'Bad drtest.' )
            endif
c
         else if ( lhs .eq. 'dvtest' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               dvtest = r4val
            else
               call rcderr ( 'Bad dvtest.' )
            endif
c
         else if ( lhs .eq. 'vamp' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vamp = r4val
            else
               call rcderr ( 'Bad vamp.' )
            endif
c
         else if ( lhs .eq. 'vfreq' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vfreq = r4val
            else
               call rcderr ( 'Bad vfreq.' )
            endif
c
         else if ( lhs .eq. 'vamp2' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vamp2 = r4val
            else
               call rcderr ( 'Bad vamp2.' )
            endif
c
         else if ( lhs .eq. 'vfreq2' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vfreq2 = r4val
            else
               call rcderr ( 'Bad vfreq2.' )
            endif
c
         else if ( lhs .eq. 'noise' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               noise = r4val
            else
               call rcderr ( 'Bad noise.' )
            endif
c
         else if ( lhs .eq. 'vnoise' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vnoise = r4val
            else
               call rcderr ( 'Bad vnoise.' )
            endif
c
         else if ( lhs .eq. 'vntime' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vntime = r4val
            else
               call rcderr ( 'Bad vntime.' )
            endif
c
         else if ( lhs .eq. 'slipv' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               slipv = r4val
            else
               call rcderr ( 'Bad slipv.' )
            endif
c
         else if ( lhs .eq. 'vclut' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vclut = r4val
            else
               call rcderr ( 'Bad vclut.' )
            endif
c
         else if ( lhs .eq. 'vdclut' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vdclut = r4val
            else
               call rcderr ( 'Bad vdclut.' )
            endif
c
         else if ( lhs .eq. 'aclut' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               aclut = r4val
            else
               call rcderr ( 'Bad aclut.' )
            endif
c
         else if ( lhs .eq. 'pclut' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               pclut = r4val
            else
               call rcderr ( 'Bad pclut.' )
            endif
c
         else if ( lhs .eq. 'dclut' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               dclut = r4val
            else
               call rcderr ( 'Bad dclut.' )
            endif
c
c-----------------------------------------------------------------------
c
c   Strip-map parameters
c
         else if ( lhs .eq. 'vplat' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vplat  = r4val
               vfocus = vplat
               write ( 6  , * ) 'Vfocus has also been changed to Vplat'
               write ( 66 , * ) 'Vfocus has also been changed to Vplat'
            else
               call rcderr ( 'Bad vplat.' )
            endif
c
         else if ( lhs .eq. 'vfocus' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               vfocus = r4val
            else
               call rcderr ( 'Bad vfocus.' )
            endif
c
         else if ( lhs .eq. 'slant0' ) then
c
            if ( getr4( rhs , lenrhs , r4val ) .eq. 0 ) then
               slant0 = r4val
            else
               call rcderr ( 'Bad slant0.' )
            endif
c
         else if ( lhs .eq. 'wrtmcd' ) then
c
            if ( geti4( rhs , lenrhs , i4val ) .eq. 0 ) then
               wrtmcd = i4val
            else
               call rcderr ( 'Bad wrtmcd.' )
            endif
c
c-------------------------- IQGEN Parameters ---------------------------
c
c   End of all parameters
c
         endif
c
      else
c
c   If there isn't an equal sign, check the command against the list of
c   allowable commands ( COMS ).  If it is legal, process it; if not,
c   give an error message.
c
         if ( comchk( tmpstr , coms , comnum , ncmd ) .ne. 0 ) then
c
            call rcderr ( 'Unrecognized expression.' )
c
         else
c
c***********************************************************************
c
c   Display values of the parameters for the user
c
            if ( tmpstr .eq. 'stat' ) then
c
               call sarst ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                      npulse , nspf , nspt , nafill , nfr , nfa ,
     .                      mrrat , mprat , nchuse , nfocus , 6 )
c
            endif
c
            if ( tmpstr .eq. 'reset' ) then
c
               call defaults ( sampr , ntr , prf , ncr , naskip ,
     .                         spulse , npulse , nspf , nspt , nafill ,
     .                         nfr , nfa , mrrat , mprat , nchuse ,
     .                         nfocus )
c
            endif
c
            if ( tmpstr .eq. 'go' ) return
c
            if ( tmpstr .eq. 'help' .or. tmpstr .eq. 'man' ) then
c
               hlpcmd = 'help'
               hlpfil = 'isart.hlp'
c
               call helper ( hlpcmd , hlpfil )
c
            endif
c
            if ( tmpstr .eq. 'bye' .or. tmpstr .eq. 'exit' .or.
     .           tmpstr .eq. 'quit' ) then
c
               close ( 66 )  !  Close log file
c
               stop ' User quit program ISART'
c
            endif
c
            if ( tmpstr .eq. 'realtime' ) then
c
               real_t = 32767
               rt_rco = 32767
               rt_img = 32767
               rt_pga = 32767
               rt_io  = 32767
c
            else if ( tmpstr .eq. 'norealtime' ) then
c
               real_t = 0
               rt_rco = 0
               rt_img = 0
               rt_pga = 0
               rt_io  = 0
c
            endif
c
c   End of loop over word commands
c
         endif
c
c   End of processing command strings with and without equal signs
c
      endif
c
 1000 continue
c
c   Come here on exit
c
 2000 continue
C
      close ( 66 )  !  Close log file
c
      stop ' User quit program ISART'
c
      end
C**
C***********************************************************************
C**
      subroutine sarst ( sampr , ntr , prf , ncr , naskip , spulse ,
     .                   npulse , nspf , nspt , nafill , nfr , nfa ,
     .                   mrrat , mprat , nchuse , nfocus , lun )
C**
C***********************************************************************
C**
      implicit none
c
c***********************************************************************
c
c   Program parameters which can be set by any user
c
      integer      ntr , ncr , naskip , spulse , npulse , nspf , nspt ,
     .             nfr , nfa , nafill , nchuse , nfocus
c
      real         sampr , prf
c
c***********************************************************************
c
c   Program parameters which can only be set by command parser
c
      include     'sarprm.h'
c
      include     'updates.h'
c
      include     'realtime.h'    !  Real-time parameters
c
      include     'kalman.h'      !  Kalman Filter Parameters
c
      integer      mrrat , mprat
c
      integer      lun
c
      write ( lun , '(/)' )
      write ( lun , * ) '                          Standard Parameters'
c
      write ( lun , '(6x,a9,3x,a30,a9,i9)' ) 'Ifile  =' , ifile  ,
     .                                       'Nfmt   =' , nfmt
c
      if ( nchuse .gt. 1 ) then
c
         write ( lun , '(6x,a9,3x,a30,a9,i9)' )
     .                                       'Ifile2 =' , ifile2 ,
     .                                       'Nchuse =' , nchuse
c
         write ( lun , '(6x,a9,3x,a30,a9,i9)' )
     .                                       'Ifile3 =' , ifile3 ,
     .                                       'Cmodel =' , cmodel
c
      endif
c
      write ( lun , '(6x,a9,3x,a30,a9,i9)' ) 'Ofile  =' , ofile  ,
     .                                       'Strtch =' , strtch
c
      write ( lun , * )
c
      write ( lun , '(6x,2(a9,i9,3x),a9,f9.2)' )
     .                                       'Otype  =' , otype ,
     .                                       'Moddat =' , moddat ,
     .                                       'Prfrat =' , prfrat   
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Finemc =' , finemc ,
     .                                       'Fixiq  =' , fixiq  ,
     .                                       'Alias  =' , alias 
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Efghz  =' , efghz  ,
     .                                       'Br     =' , br     ,
     .                                       'Sampr  =' , sampr / 1.0E+6
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Prf    =' , prf    ,
     .                                       'Taywtr =' , taywtr ,
     .                                       'Taywta =' , taywta
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Ntr    =' , ntr    ,
     .                                       'Ncr    =' , ncr    ,
     .                                       'Naskip =' , naskip
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Spulse =' , spulse ,
     .                                       'Npulse =' , npulse ,
     .                                       'Nspf   =' , nspf
c
      write ( lun , '(6x,a9,f9.2,3x,a9,f9.2,3x,a9,i9)' )
     .                                       'Beta   =' , beta   ,
     .                                       'Alpha0 =' , alpha0 ,
     .                                       'Nafill =' , nafill
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Nfr    =' , nfr    ,
     .                                       'Nfa    =' , nfa    ,
     .                                       'Mprat  =' , mprat
c
      write ( lun , '(6x,2(a9,f9.2,3x),a9,i9)' )
     .                                       'Overrg =' , overrg ,
     .                                       'Rgwalk =' , rgwalk ,
     .                                       'Nrcent =' , nrcent
c
      write ( lun , '(6x,2(a9,i9,3x),a9,f9.2,3x)' )
     .                                       'Mrrat  =' , mrrat  ,
     .                                       'UWB    =' , 0      ,
     .                                       'Snrmin =' , snrmin
c
      write ( lun , '(6x,a9,i9,3x,a9,f9.2,3x,a9,i9)' )   
     .                                       'Maxit  =' , maxit  ,
     .                                       'Center =' , center ,
     .                                       'Curve  =' , curve 
c
      write ( lun , '(6x,a9,f9.2,3x,a9,i9,3x,a9,f9.2)' )
     .                                       'Outlie =' , outlie ,
     .                                       'Acoefs =' , acoefs ,
     .                                       'Aghost =' , aghost 
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'Editac =' , editac ,
     .                                       'Length =' , length ,
     .                                       'Width  =' , width
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Presum =' , presum ,
     .                                       'Angle  =' , angle  ,
     .                                       'Ti_con =' , ti_con
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'MultiP =' , multip ,
     .                                       'Nlocal =' , nlocal ,
     .                                       'Pixbar =' , pixbar 
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Nrkeep =' , nrkeep ,
     .                                       'Npass  =' , npass  ,
     .                                       'Ntaper =' , ntaper
c
      write ( lun , '(6x,2(a9,i9,3x),a9,f9.2)' )
     .                                       'Nlag   =' , nlag   ,
     .                                       'Nfskip =' , nfskip ,
     .                                       'Dbinc  =' , dbinc    
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'Notch  =' , notch  ,
     .                                       'Vnotch =' , vnotch ,
     .                                       'Dvntch =' , dvntch
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'Notch  =' , notch  ,
     .                                       'Rnotch =' , rnotch ,
     .                                       'Drntch =' , drntch
c
      write ( lun , '(6x,2(a9,f9.2,3x),a9,i9)' )
     .                                       'Frrate =' , frrate ,
     .                                       'Tinteg =' , tinteg ,
     .                                       'Mode   =' , mode
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Color  =' , color  ,
     .                                       'Real_t =' , real_t ,
     .                                       'Rt_rco =' , rt_rco
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Rt_img =' , rt_img ,
     .                                       'Rt_pga =' , rt_pga ,
     .                                       'Rt_io  =' , rt_io
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Slbkil =' , slbkil , 
     .                                       'dokalm =' , dokalm ,
     .                                       'wrtmcd =' , wrtmcd
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Pgtype =' , pgtype ,
     .                                       'Tgt_si =' , tgt_si ,
     .                                       'Rd_tgt =' , rd_tgt
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'quiet  =' , quiet  ,
     .                                       'v_corr =' , v_corr ,
     .                                       'a_corr =' , a_corr 
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'nbands =' , nbands ,
     .                                       'alt_km =' , alt_km ,
     .                                       'rfdelt =' , 1.e-6 * rfdelt
c
      write ( lun , '(6x,2(a9,i9,3x),a9,f9.2)' )
     .                                       'Reduce =' , reduce ,
     .                                       'N_vadd =' , n_vadd ,
     .                                       'd_vadd =' , d_vadd
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'Nfocus =' , nfocus ,
     .                                       'accorr =' , accorr ,
     .                                       'arcorr =' , arcorr
c
      write ( lun , '(6x,2(a9,i9,3x),a9,f9.2)' )
     .                                       'Nspt   =' , nspt   ,
     .                                       'Imgslk =' , imgslk ,
     .                                       'Vfocus =' , vfocus   
c
      write ( lun , '(6x,2(a9,i9,3x),a9,f9.2)' )
     .                                       'Rlook  =' , rlook  ,
     .                                       'Flook  =' , flook  ,
     .                                       'Rcentr =' , rcentr
c
      write ( lun , '(6x,3(a9,i9,3x))' )     'Addvib =' , addvib ,
     .                                       'Addamp =' , addamp ,
     .                                       'Fixamp =' , fixamp 
c
      write ( lun , '(6x,a9,i9,3x,a9,f9.2,3x,a9,i9)' )
     .                                       'Qantiz =' , qantiz ,
     .                                       'iqlsb  =' , iqlsb  ,
     .                                       'firsts =' , firsts
c
      write ( lun , * )
c
      if ( nfmt .ge. 0 ) return
c
c   nfmt < 0 ==> Internally generated data from IQGEN
c
      write ( lun , '(/,a43)' )
     .   '                           IQGEN Parameters'
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Vbar   =' , vbar   ,
     .                                       'Dvdt   =' , dvdt   ,
     .                                       'Vamp   =' , vamp
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Dadt   =' , dadt   ,
     .                                       'Aamp   =' , aamp   ,
     .                                       'Afreq  =' , afreq
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Daccdf =' , daccdf ,
     .                                       'Daccdr =' , daccdr ,
     .                                       'Daccdt =' , daccdt
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Drtest =' , drtest ,
     .                                       'Dvtest =' , dvtest ,
     .                                       'Vfreq  =' , vfreq
c
      write ( lun , '(6x,a9,i9,3x,2(a9,f9.2,3x))' )
     .                                       'Ntest  =' , ntest  ,
     .                                       'Noise  =' , noise  ,
     .                                       'Vnoise =' , vnoise
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Vntime =' , vntime ,
     .                                       'Slipv  =' , slipv  ,
     .                                       'Vamp2  =' , vamp2
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )
     .                                       'Vclut  =' , vclut  ,
     .                                       'Vdclut =' , vdclut ,
     .                                       'Aclut  =' , aclut
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Pclut  =' , pclut  ,
     .                                       'Dclut  =' , dclut  ,
     .                                       'Vfreq2 =' , vfreq2
c
      write ( lun , '(6x,3(a9,f9.2,3x))' )   'Vplat  =' , vplat  ,
     .                                       'Slant0 =' , slant0 ,
     .                                       'Rbar   =' , rbar 
c
      write ( lun , * )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine memory ( ntr , ncr , naskip , nspf_s , nafill_s , nfr ,
     .                    nfa , nsr , nsa , nakeep , mrrat , ntp , uwb ,
     .                    npulse , overrg , mode , nabuff , nchuse ,
     .                    nfskip , cmodel , tgt_si , nfocus , nspt ,
     .                    lun )
C**
C***********************************************************************
C**
      implicit none
c
      integer  ntr , ncr , naskip , nspf_s , nafill_s , nfr , nfa ,
     .         nsr , nsa , nakeep , mrrat , nafine , ntp , npulse ,
     .         nabuff , twopwr , lun , mode , nchuse , cmodel , nchmem ,
     .         nfskip , tgt_si , nfocus , uwb , nspt
c
      real     mw , overrg
c
      include 'realtime.h'
c
      write ( lun , * )
c
c   Ensure that output image size is an integer times mrrat
c
      if ( mod( nfr , mrrat ) .ne. 0 ) then
c
         nfr = mrrat * ( 1 + ( nfr / mrrat ) )
c
         write ( lun , '(/,7x,a40,i8,/)' )
     .      ' Parameter Override: nfr rounded up to  ' , nfr
c
      endif
c
c   Set the pulse buffer to four times the skip length
c
      ntp = 4 * naskip
c
c   Ensure that there is enough over-sampling in range
c
      ncr = max( ncr , 3 * ntr )
c
c   Ensure that the number of pulses is at least as large as the
c   raw data buffer
c
      if ( npulse .lt. ntp ) then
c
         npulse = ntp
c
         write ( lun , '(/,7x,a40,i8,/)' )
     .      'Parameter Override: npulse set equal to ' , npulse
c
      endif
c
c   Choose enough sub-image range cells so that the final range extent
c   is covered plus a buffer zone of 256 at each end
c
      nsr = 48 +
     .      ifix( ( float( nfr ) * float( ncr ) ) /
     .            ( float( ntr ) * overrg ) )
c
      nsr = nsr - mod( nsr , 2 )
c
      nsr = min( nsr , ncr )
c
c   Ensure that the sub-image size is a power of two at least as large
c   as three times the skip length
c
      if ( rt_img .ne. 0 ) then
c
         nsa    = twopwr( 3 * naskip )
c
         nakeep = nsa
c
         nabuff = nspf_s
c
      else if ( uwb .gt. 0 ) then
c
c   For ultra wide band mode space the sub-image points twice as far apart
c
         nsr    = nsr / 2
c
c   In strip-map mode, keep half the frequency bins plus 4
c
         nsa    = 3 * naskip
c
         nakeep = 4 + ( nsa / 2 ) + mod( nsa / 2 , 2 )
c
c   Enlarge time storage
c
         nabuff = nspf_s + ( nfa / naskip )
c
      else
c
c   In ISAR mode, keep only those frequencies which contribute to the
c   fine resolution image
c
c      nafine = no. of coarse Doppler cells to cover fine res. image
c
c      Add a buffer region to nafine which corresponds to 4 coarse
c      Doppler cells (2 on each side)
c
         nsa    = twopwr( 3 * naskip )
c
         nafine = ( ( nfa * nsa ) / ( nafill_s * naskip ) ) 
c
         nakeep = min( nsa , 4 + nafine )
c
         if ( mod( nakeep , 2 ) .ne. 0 ) nakeep = nakeep + 1
c
         nabuff = nspf_s
c
      endif
c
c   For variable integration time mode over-ride the above logic to
c   keep all the sub-image frequencies
c
      if ( nfskip .eq. 0 ) then
c
         nabuff = 2 * nabuff
c
         nakeep = nsa
c
      endif
c
c   If targets are being detected at the sub-image stack then increase
c   the size of the buffer so that we use a larger integration time for
c   this process and put them further into the future.  But, note that
c   for variable integration time modes this is already done.
c
      if ( tgt_si .ne. 0 .and. abs(nfskip) .ge. 0 ) then
         if ( nspt .gt. 0 ) nabuff = nspt 
         nabuff = 3 * nabuff
      endif
c
      write ( lun , '(7x,a16,i8)' )   'Nsa    =        ' , nsa
c
      write ( lun , '(7x,a16,i8)' )   'Nsr    =        ' , nsr
c
      write ( lun , '(7x,a16,i8)' )   'Ncr    =        ' , ncr
c
      write ( lun , '(7x,a16,i8)' )   'Nabuff =        ' , nabuff
c
      write ( lun , '(7x,a16,i8)' )   'Nspt   =        ' , nspt
c
      write ( lun , '(7x,a16,i8)' )   'Nakeep =        ' , nakeep
c
      write ( lun , '(/,7x,a24,/)' )  'Memory Requirements (MW)'
c
      mw = 1.0 / ( 1024.0 ** 2 )
c
      write ( lun , '(7x,a16,f8.3)' ) 'Sub-image Buffer' ,
     .                          mw * float( 2 * nabuff * nsr * nakeep )
c
c   Multiple channels of raw data if nchuse > 1
c
      if ( nchuse .eq. 1 ) then
c
         nchmem = 1
c
      else
c
         nchmem = nchuse + 1
c
      endif
c
      write ( lun , '(7x,a16,f8.3)' ) 'Raw Data Buffer ' ,
     .                          mw * float( 2 * ntr * ntp * nchmem )
c
c   Multiple channels of range compressed data if cmodel > 0
c
      if ( nchuse .gt. 1 .and. cmodel .gt. 0 ) then
c
         nchmem = nchuse + 1
c
      else
c
         nchmem = 1
c
      endif
c
      write ( lun , '(7x,a16,f8.3)' ) 'R-Compr. Buffer ' ,
     .                          mw * float( 2 * ncr * ntp * nchmem )
c
      write ( lun , '(7x,a16,f8.3)' ) 'Fine Res. Image ' ,
     .                  mw * float( ( 2 + 2 * nfocus ) * 2 * nfr * nfa )
c
      write ( lun , * )
c
      return
      end
