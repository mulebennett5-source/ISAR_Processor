C**
C***********************************************************************
C**
      subroutine getacc    ( nf , nr , nt , nspt , dotdot , fstart ,
     .                       rstart , tstart , deltaf , deltar ,
     .                       deltat , athrsh , rcen , fcen , iframe )
C**
C***********************************************************************
C**
c
c                      ISAR Global Motion Model
c
c                            Version 3.32
c
c                            July 27, 2012
c
c                        John R. Bennett, PhD
c                      Kenneth A. Melendez, PhD
c                           David S. Brown
c
c
c  =====================================================================
c
c   Background:
c
c      This subroutine is the start of a chain of routines that
c      implement the Global Motion Model of the RDRTec/Telephonics Inverse
c      SAR processor.
c
c   Purpose:  
c
c      Given acceleration data that has been estimated at values of
c      Doppler frequency, range, and time, the subroutine performs a
c      least squares fit of the data to the model:
c
c         A(F,R,T) = A0 + AT*T + AR*R + AF*F + ART*R*T + AFT*F*T
c
c      and then uses the fit to fill the acceleration array,
c      dotdot(f,r,t), which measures the acceleration in Hz/second.  The
c      acceleration is used in the image formation phase to focus the
c      ISAR image frame.
c
c      The routine also estimates a number of physical parameters which
c      can be used to focus or to interpret the images.
c
c      The algorithm consists of three stages:  Scatterer Screening,
c      Least Squares Fit, and Array Fill.
c        
c  ---------------------------------------------------------------------
c   Passed Variables 
c  ---------------------------------------------------------------------
c
c   nf                              Azimuth cells in an image frame
c
c   nr                              Range cells in a frame
c
c   nt                              Time steps per frame.
c
c   dotdot(nf,nr,nt)                Output array containing acceleration
c                                   in units of Hz/sec
c
c   fstart, rstart, tstart          The coordinates of the first
c                                   element of the dotdot array.
c
c   deltaf, deltar, deltat          Sample spacing between f, r,
c                                   and t cells for the dotdot array
c
c   athrsh                          Maximum magnitude of acceleration
c                                   data to be used in determining the
c                                   acceleration model (Hz/second)
c
c   curve                           Integer flag specifying whether to
c                                   fill acceleration array with zeroes
c                                   (curve=0), or to perform fit and
c                                   fill with acceleration fit estimates
c                                   (curve=1,2,3).
c
c   outlie                          Number of standard deviations from
c                                   the mean for outlier removal
c
c   acoefs                          Number of coefficients to be used
c
c   aghost                          Weight applied to the 'ghost'
c                                   scatterers used to stabilize the
c                                   curve fit in the case of nearly
c                                   co-linear data points
c
c   rcen                            The range to the centroid of the
c                                   target (output)
c
c   fcen                            The frequency to the centroid of the
c                                   target (output)
c
c   editac                          Integer flag for editing points
c
c   lambda                          Radar wavelength (m)
c
c   a0,at,ar,af,art,aft             The fit coefficients
c
c   --------------------------------------------------------------------
c   Other Local Variables 
c   --------------------------------------------------------------------
c
c   a_coeffs (real array)           Vector containing fit coefficients
c                                   of the fit function
c
c   acomax (integer parameter)      Max. Number of fit coefficients
c
c   nlist (integer)                 Number of entries in the target list
c
c   afunc (real array)              Vector containing values of basis
c                                   functions of linear acceleration
c                                   model evaluated at a particular
c                                   f,r,t point.
c
c   n_inside                        Number of target points which have
c                                   f,r,t state vectors within the
c                                   f0,fm,r0,rm,t0,tm boundaries and
c                                   which pass all the editing passes
c
c   sd_acc (real)                   Standard deviation of acceleration
c
c   mean_acc (real)                 Mean of acceleration
c
c   sum_acc (real)                  Sum of acceleration
c
c   sum_sq_acc (real)               Sum of squares of acceleration.
c
c   sd_dw (real)                    Standard deviation of doppler width
c
c   mean_dw (real)                  Mean of doppler width
c
c   sum_dw (real)                   Sum of doppler width
c
c   sum_sq_dw (real)                Sum of squares of doppler width
c
c   minf,minr,mint,minacc (real)    These are the minima of f,r,t,acc
c                                   for those target points which have 
c                                   been selected.
c
c   maxf,maxr,maxt,maxacc (real)    These are the maxima of f,r,t. 
c
c   b,qrrsd (real array)            Arrays needed for sqrlss
c
c   u,qraux,qrwork (real array)     Arrays needed for sqrank,sqrlss
c
c  =====================================================================
c
      implicit none
c
c-----------------------------------------------------------------------
c
c   Use the standard parameters stored in common
c
      include          'sarprm.h'
c
c   Use the global target list stored in common
c
      include          'tglist.h'
c
c   Use Kalman Filter Matrices
c
      include          'kalman.h'
c
      include          'updates.h'
c
c-----------------------------------------------------------------------
c              
      integer           nf , nr , nt , nspt , iframe
c
      real              deltaf , deltat , deltar , dotdot(nf,nr,nt) ,
     .                  fstart , rstart , tstart , sd_acc , mean_acc ,
     .                  sd_dw , mean_dw , athrsh , rcen , fcen ,
     .                  fitmin , fitmax , fitbar , fitsdv , mean_time ,
     .                  domegadt_af , domegadt_art , domegadt_est ,
     .                  omega_lim
c
c   Variables used in sqrank, sqrlss
c
      integer           acomax            !  Maximum no. of acceleration
      parameter       ( acomax = 6 )      !  coefficients
c
      integer           jpvt(acomax) , kr , kk
c
      real              u(nlist+8,acomax) , b(nlist+8) , tmp , dvmax ,
     .                  tol , qraux(acomax) , qrwork(acomax) , rsmall ,
     .                  qrrsd(nlist+8) , qrtol , ri(acomax,acomax) ,
     .                  rit(acomax,acomax) , cov(acomax,acomax)
c
      parameter       ( tol = 1.e-5 )     !  Tolerance for convergence
c
c   Local Variables
c
      real              minf , minr , mint , maxf , maxr , maxt , d_a ,
     .                  minacc , maxacc , tbeg , tend , accdum , lf ,
     .                  tinflu , mean_range , mean_freq , omega_use ,
     .                  vsmall , d_r , d_f , non_physical , lt , lr ,
     .                  min_non_phys , max_non_phys , del , sum_frsq
c
      integer           n_inside , iloop , jloop , kloop , i , j , itgt
c
      real              afunc(acomax) , a_coeffs(acomax)
c
      logical           check , tcheck , pcheck , scheck ! Source check
c
c-----------------------------------------------------------------------
c
c   Variables for diagnostic information of the acceleration field
c   used to focus the image
c
      integer           npass1 , npass2
c
      real              a0mm , armm , afmm , sd_a0mm , sd_at , sd_armm ,
     .                  sd_afmm , sd_art , sd_aft , var_omg , var_armm ,
     .                  xbar , xsdv , ybar , ysdv , ddbar , ddsdv , spin
c
c   AIS local variables
c
      integer           dlength , elength
c
      parameter       ( dlength = 64 , elength = 64 )
c
      integer           didx , eidx
c
      real              sumdwdth2 , curreta , d , dmin , dinc , etamin ,
     .                  d_r0 , d_f0 , eta0 , etainc , pi , degrad ,
     .                  dmax , hough_mat(dlength,elength) , aismajor ,
     .                  aisminor , aismnr , aismnf
     .                  
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Begin Executable Section
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      pi     = atan2( 0.0 , - 1.0 )
c
      degrad = pi / 180.0 
c
c   If no targets pass, then don't update the centroid estimate
c
      rcen   = 0.0
      fcen   = 0.0
c
c   Use an influence time equal to the maximum of the integration time
c   and tlarge.  Also, do not use estimates earlier than half this.
c
      tbar   = tstart + 0.5 * float( nt - 1 ) * deltat
      tinflu = float( nspt ) * deltat
c
      tbeg   = amax1( 0.0 , tbar - 1.0 * tinflu )
      tend   = amax1( 0.0 , tbar + 1.0 * tinflu )
c
c   Zero target selection flag so that we can re-consider all targets
c   for this frame
c
      do iloop = 1 , nlist
c
c   Set the flags to zero unless the target has been killed
c
         if ( iflag(iloop) .eq. 101 ) then
c
            if ( time(iloop) .lt. tbeg ) then
c
               time(iloop)  = - 1.0
c
               iflag(iloop) = 0
c
            endif
c
         else
c
            iflag(iloop) = 0
c
         endif
c
      enddo
c
c   If curve equals zero, fill acceleration array with zeroes or, if in
c   the general strip map mode, fill with the focus function determined
c   from the navigation data.
c
      if ( curve .eq. 0 .or. acoefs .eq. 0 ) then 
c
         if      ( mode .eq. 5 ) then
c
c   Use constant
c
            dotdot(:,:,:) = a_corr
c
         else if ( mode .eq. 6 .or. mode .eq. 7 ) then
c
c   General strip map focus function when the platform and the mo-comp
c   point move at different and time-varying speeds
c
            do kloop = 1 , nt
c
               lt = ( float( kloop - 1 - nt / 2 ) - 0.5 ) * deltat
c
               do jloop = 1 , nr
c
                  lr = ( float( jloop - 1 - nr / 2 ) - 0.5 ) * deltar
c
                  do iloop = 1 , nf 
c
                     lf = ( float( iloop - 1 - nf / 2 ) - 0.5 ) * deltaf
c
                     dotdot(iloop,jloop,kloop) = ( 2.0 / lambda ) * (
     .                  strip(2) * strip(3) + 2.0 * lt * (
     .                  strip(2) * ddt_vspot + strip(3) * ddt_vplat ) )
     .                                         / amax1( 1.0 , strip(1) )
     .            - ( 2.0 / lambda ) * lr * ( strip(2) * sin( strip(5) )
     .                                        / strip(1) ) ** 2
     .            + omegadot_omega * lf
c
                  enddo
c
               enddo
c
            enddo
c
         else if ( mode .eq. 8 ) then
c
c   Special spotlight focus function
c
            do kloop = 1 , nt
c
               lt = ( float( kloop - 1 - nt / 2 ) - 0.5 ) * deltat
c
               do jloop = 1 , nr
c
                  lr = ( float( jloop - 1 - nr / 2 ) - 0.5 ) * deltar
c
                  do iloop = 1 , nf 
c
                     lf = ( float( iloop - 1 - nf / 2 ) - 0.5 ) * deltaf
c
                     dotdot(iloop,jloop,kloop) =
     .                  - ( 2.0 / lambda ) * lr * ( 79.0 / 2995.0 ) ** 2
c
                  enddo
c
               enddo
c
            enddo
c
         else
c
c   Fill with zeroes
c
            dotdot(:,:,:) = 0.0
c
         endif
c
         return
c
      endif
c
      if ( tbar .lt. tinflu .or. nlist .lt. 2 ) then
c
c   Fill with zeroes
c
         dotdot(:,:,:) = 0.0
c
c   Output for MCA
c
         write ( 77 , '(1x,2i8,7f13.3)' ) iframe , 0 , tbar , 0.0 ,
     .                                    - 1.0 , 0.0 , 0.0 , 0.0
c
         return
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Otherwise, fill dotdot acceleration array with model estimates
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                    Stage 1 - Scatterer Screening
c
c   Loop over all elements to find min/max of f,r,t.
c   Only consider elements that fall within the volume and edit the
c   points according to three criteria.
c
c-----------------------------------------------------------------------
c
c                          First Edit Pass
c
c   This stage is always implemented.
c
       do iloop = 1 , nlist
c
         tcheck = ( mod( curve , 2 )     .eq. 1 ) .and.
     .            ( source(iloop) .eq. 't' )
c
         pcheck = ( mod( curve / 2 , 2 ) .eq. 1 ) .and.
     .            ( source(iloop) .eq. 'p' )
c
         scheck = ( curve .ge. 4 )                .and.
     .            ( source(iloop) .eq. 's' )
c
         check  = tcheck .or. pcheck .or. scheck
c
c   Check to make sure this point has an acceleration less than the
c   threshold and has a state vector that lies within the (t,f,r)
c   volume of interest.
c    
         if ( check                               .and.
     .        ( abs( accel(iloop) ) .le. athrsh ) .and.
     .        ( time(iloop)  .ge. tbeg          ) .and.
     .        ( time(iloop)  .le. tend          ) .and.
     .        ( iflag(iloop) .ne. 101           ) ) then
c
            if ( nint( ( freq(iloop)  - fstart ) / deltaf + 1 ) .ge. 1
     .                               .and.
     .           nint( ( range(iloop) - rstart ) / deltar + 1 ) .ge. 1
     .                               .and.
     .           nint( ( freq(iloop)  - fstart ) / deltaf + 1 ) .le. nf
     .                               .and.
     .           nint( ( range(iloop) - rstart ) / deltar + 1 ) .le. nr
     .                               ) then
c
               iflag(iloop) = 1
c
            endif
c
         endif
c
      enddo
c
c   Calculate acceleration statistics for all scatterers which pass this
c   editing stage
c
      call acstat ( 1 , n_inside , nlist , iflag , accel , dwdth ,
     .              range , freq , time , mean_acc , sd_acc , mean_dw ,
     .              sd_dw , mean_range , mean_freq , mean_time )
c
      npass1 = n_inside
c
c-----------------------------------------------------------------------
c
c                          Second Edit Pass
c
c   This stage is always implemented.  It rejects 'outliers' - those
c   scatterers which have accelerations or Doppler widths more than
c   'outlie' standard deviations from their means.  If 'outlie' has been
c   set to a large value then all scatterers will pass.
c
      do iloop = 1 , nlist
c
         if ( iflag(iloop) .eq. 1 ) then
c
            if ( accel(iloop) .le. ( mean_acc + outlie * sd_acc )  .and.
     .           accel(iloop) .ge. ( mean_acc - outlie * sd_acc )  .and.
     .           dwdth(iloop) .le. ( mean_dw  + outlie * sd_dw  ) ) then
c
               iflag(iloop) = 2
c
            endif
c
         endif
c
      enddo
c
c   Calculate acceleration statistics for all scatterers which pass this
c   editing stage
c
      call acstat ( 2 , n_inside , nlist , iflag , accel , dwdth ,
     .              range , freq , time , mean_acc , sd_acc , mean_dw ,
     .              sd_dw , mean_range , mean_freq , mean_time )
c
      npass2 = n_inside
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      if ( acoefs .ge. 3 .and. curve .gt. 3 ) then
c
c   If there is a valid motion model, then use it to remove targets
c   which lie outside the a priori size of the target
c
c   First call to tgtbox - cross section weighted.  This allows targets
c   to be eliminated if they are too far from the main cross section
c   region.
c
         call tgtbox ( n_inside , 0.25 , 1.5 , 2.0 , 5.0 )
c
c   Second call to tgtbox - doppler width weighted
c
         call tgtbox ( n_inside , 0.75 , 1.0 , 2.0 , 5.0 )
c
c   Third call to tgtbox - doppler width weighted; to make sure the box
c   is computed with all targets edited.
c
         call tgtbox ( n_inside , 0.75 , 1.0 , 2.0 , 5.0 )
c
c   Re-calculate the statistics based on the thinned target set
c
         call acstat ( editac , n_inside , nlist , iflag , accel ,
     .                 dwdth , range , freq , time , mean_acc , sd_acc ,
     .                 mean_dw , sd_dw , mean_range , mean_freq ,
     .                 mean_time )
c
      endif
c
c   At this point all the targets have been selected
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                    Stage 2 - Least Squares Fit
c
c   Add the data points to the linear least squares data matrix and
c   obtain minima and maxima of the data within the volume.
c
c   Initialize the mins and maxs and number of points found
c   inside the outer volume.
c
      minf   = 0.0
      minr   = 0.0
      mint   = 0.0
      minacc = 0.0
      maxf   = 0.0
      maxr   = 0.0
      maxt   = 0.0
      maxacc = 0.0
c
      if ( n_inside .gt. 0 ) then
c
         u(:,:) = 0.0
c
         b(:)   = 0.0
c
         itgt   = 0
c
         do iloop = 1 , nlist
c
            if ( iflag(iloop) .eq. editac ) then
c
               itgt = itgt + 1
c
c              Since this point is inside the outer volume,
c              it is eligible to be used for the fit. Feed it 
c              into the fit algorithm by first computing
c              the linear acceleration model basis functions
c              output given the state vector of this point 
c              as input.
c
               call accmodel ( time(iloop) - tbar , range(iloop) ,
     .                         freq(iloop) , afunc , acoefs )
c
c              Add current point to u,b arrays needed for sqrlss.
c
               tmp     = 1.0 / dwdth(iloop) 
c
               do j = 1 , acoefs
c
                  u(itgt,j) = afunc(j) * tmp
c
               enddo
c
               b(itgt) = accel(iloop) * tmp
c
c              Does this point have a state vector for which 
c              freq, range, or time is a max or min?
c
               if ( itgt .eq. 1 ) then
c
                  minf   = freq(iloop)
                  minr   = range(iloop)
                  mint   = time(iloop)
                  minacc = accel(iloop)
c
                  maxf   = freq(iloop)
                  maxr   = range(iloop)
                  maxt   = time(iloop)
                  maxacc = accel(iloop)
c
               else
c
                  if ( freq(iloop)  .lt. minf   ) minf   = freq(iloop)
                  if ( range(iloop) .lt. minr   ) minr   = range(iloop)
                  if ( time(iloop)  .lt. mint   ) mint   = time(iloop)
                  if ( accel(iloop) .lt. minacc ) minacc = accel(iloop)
c
                  if ( freq(iloop)  .gt. maxf   ) maxf   = freq(iloop)
                  if ( range(iloop) .gt. maxr   ) maxr   = range(iloop)
                  if ( time(iloop)  .gt. maxt   ) maxt   = time(iloop)
                  if ( accel(iloop) .gt. maxacc ) maxacc = accel(iloop)
c
               endif
c
            endif
c
         enddo
c
c   Correct to the center of the integration time
c
         mint = mint - tbar
         maxt = maxt - tbar
c
c   The points with state vectors that lie in the specified outer
c   volume have been identified. These will be used for obtaining the
c   parameters which characterize the best fit of the input acceleration
c   data to the linear model for the acceleration.
c
c-----------------------------------------------------------------------
c
c   To stabilize the calculation against co-linear data, add eight
c   points with zero acceleration at the corners of the known data
c   region
c
c   First, make sure the corners are at least one grid point apart
c
         if ( ( maxf - minf ) .le. deltaf ) then
c
            minf = minf - 0.5 * deltaf
            maxf = maxf + 0.5 * deltaf
c
         endif
c
         if ( ( maxr - minr ) .le. deltar ) then
c
            minr = minr - 0.5 * deltar
            maxr = maxr + 0.5 * deltar
c
         endif
c
         if ( ( maxt - mint ) .le. deltat ) then
c
            mint = mint - 0.5 * deltat
            maxt = maxt + 0.5 * deltat
c
         endif
c
c-----------------------------------------------------------------------
c
c   Determine the weight applied to the mean acceleration at 'ghost'
c   points.  Ghost points are the eight dummy points added at the
c   corners of the 3-D area which covers the real observations being
c   used at this time.  If n_inside < 100, then the variance of the
c   error due to the ghost points is a fixed percentage of the total.
c   The parameter, aghost, can be set by the user to reduce the effect
c   of this process or to eliminate it.  Overall (and especially for
c   real data) this process tends to stabilize the calculations since
c   it forces the program to choose the acceleration field with the
c   simplest structure which is consistent with the data.  When the
c   data is nearly co-linear in range and frequency space, the ghost
c   points cause the plane with the smallest slope to be chosen.
c
c   Since there are 8 ghost points, the formula below is set to give
c   them 8 % of the error variance.  However, the fact that they are at
c   the corners of the region increases their effect and the fact that
c   they use the mean doppler width decreases their effect.
c
         tmp    = amax1( 1.0 , 0.1 * sqrt( float( n_inside ) ) ) *
     .            ( amax1( 1.0E-6 , aghost ) / mean_dw )
         accdum = mean_acc * tmp
c
         call accmodel ( mint , minr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+1,j) = afunc(j) * tmp
         enddo
         b(itgt+1) = accdum
c
         call accmodel ( mint , minr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+2,j) = afunc(j) * tmp
         enddo
         b(itgt+2) = accdum
c
         call accmodel ( mint , maxr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+3,j) = afunc(j) * tmp
         enddo
         b(itgt+3) = accdum
c
         call accmodel ( mint , maxr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+4,j) = afunc(j) * tmp
         enddo
         b(itgt+4) = accdum
c
         call accmodel ( maxt , minr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+5,j) = afunc(j) * tmp
         enddo
         b(itgt+5) = accdum
c
         call accmodel ( maxt , minr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+6,j) = afunc(j) * tmp
         enddo
         b(itgt+6) = accdum
c
         call accmodel ( maxt , maxr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+7,j) = afunc(j) * tmp
         enddo
         b(itgt+7) = accdum
c
         call accmodel ( maxt , maxr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u(itgt+8,j) = afunc(j) * tmp
         enddo
         b(itgt+8) = accdum
c
c-----------------------------------------------------------------------
c
c   Solve u*a_coeffs = b for a_coeffs using QR factorization
c
         qrtol = 1.0e-20
c
         do j = 1 , acoefs
c
            do i = 1 , n_inside + 8
c
               if ( abs( u(i,j) ) .gt. qrtol ) qrtol = abs( u(i,j) )
c
            enddo
c
         enddo
c
         do j = 1 , acomax
c
            a_coeffs(j) = 0.0
c
         enddo
c
         qrtol = tol / qrtol
c
c-----------------------------------------------------------------------
c
         if ( dokalm .gt. 0 .and. n_inside .gt. 0 .and.
     .        major_iter_count .gt. 10 ) then
c
            call info_filter ( b , u , acoefs , qrtol )
c 
         else 
c
            if ( major_iter_count .eq. 10 ) then 
c
               state_vec_x(1) = a0
               state_vec_x(2) = at
               state_vec_x(3) = ar
               state_vec_x(4) = af
               state_vec_x(5) = art
               state_vec_x(6) = aft
c
            endif
c
         endif
c
         major_iter_count = major_iter_count + 1
c
c-----------------------------------------------------------------------
c
         call sqrank ( u , nlist + 8 , n_inside + 8 , acoefs , qrtol ,
     .                 kr , jpvt , qraux , qrwork )
c
         call qrcov  ( u , nlist + 8 , acoefs , 1.0 , acomax ,
     .                 cov , ri , rit )
c
         call sqrlss ( u , nlist + 8 , n_inside + 8 , acoefs , kr ,
     .                 b , a_coeffs , qrrsd , jpvt , qraux )
c
         call stats  ( qrrsd , n_inside + 8 , fitmin , fitmax , fitbar ,
     .                 fitsdv )
c
         if ( fitsdv .gt. qrtol ) cov(:,:) = cov(:,:) / fitsdv ** 2
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                       Stage 3 - Array Fill
c
c   Fill variables with the returned coefficients.
c
         if ( dokalm .gt. 0 .and. major_iter_count .gt. 30 ) then
c
            a0  = state_vec_x(1)
            at  = state_vec_x(2)
            ar  = state_vec_x(3)
            af  = state_vec_x(4)
            art = state_vec_x(5)
            aft = state_vec_x(6)
c
         else
c
            a0  = a_coeffs(1)
            at  = a_coeffs(2)
            ar  = a_coeffs(3)
            af  = a_coeffs(4)
            art = a_coeffs(5)
            aft = a_coeffs(6)
c
         endif
c
c-----------------------------------------------------------------------
c
         if ( acoefs .ge. 3 .and. curve .gt. 3 ) then
c
c   Let the PGA algorithm handle all terms which are independent
c   of range and frequency
c
            a0 = 0.0
            at = 0.0
c
         endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   All required calculations completed - fill the acceleration array
c   based on these results
c
      if ( n_inside .gt. 0 ) then
c
c   Fill the dotdot acceleration array based on the fit.
c
         call fill_3d_array ( nf , nr , nt , deltaf , deltar , deltat ,
     .                        minf , maxf , minr , maxr , mint , maxt ,
     .                        fstart , rstart , tstart - tbar , minacc ,
     .                        maxacc , a0 , at , ar , af , art , aft ,
     .                        dotdot )
c
c-----------------------------------------------------------------------
c
c   If no points were found within the limits, then fill with zeros
c
      else
c
         a0            = 0.0
         at            = 0.0
         ar            = 0.0
         af            = 0.0
         art           = 0.0
         aft           = 0.0
c
         dotdot(:,:,:) = 0.0 
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                 END OF BASIC FOCUSING ALGORITHM
c
c   At this point the basic algorithm is complete - every calculation
c   required for the focusing of an image frame is done.  Everything
c   in the following code is for diagnostic purposes.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Apply Dave Brown's formulas for estimating the standard deviations
c   of the coefficients
c
         a0mm     = 0.0
         armm     = 0.0
         afmm     = 0.0
         sd_a0mm  = 0.0
         sd_at    = 0.0
         sd_armm  = 0.0
         var_armm = 0.0
         sd_afmm  = 0.0
         sd_art   = 0.0
         sd_aft   = 0.0
c
         if ( acoefs .gt.  0 ) then
c
            a0mm     = a0
            sd_a0mm  = sqrt( amax1( 0.0 , cov(1,1) ) )
c
         endif
c
         if ( acoefs .gt.  1 ) then
c
            sd_at    = sqrt( abs( cov(2,2 ) ) )
c
         endif
c
         if ( acoefs .gt.  2 ) then
c
            armm     = ar
c
            sd_armm  = sqrt( amax1( 0.0 , cov(3,3) ) )
c
            var_armm = sd_armm ** 2
c
         endif
c
         if ( acoefs .gt.  3 ) then
c
            afmm     = af
            sd_afmm  = sqrt( amax1( 0.0 , cov(4,4) ) )
c
         endif
c
         if ( acoefs .gt. 4 ) then
c
            sd_art   = sqrt( abs( cov(5,5) ) )
c
         endif
c
         if ( acoefs .gt. 5 ) then
c
            sd_aft   = sqrt( abs( cov(6,6) ) )
c
         endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Centering correction - use the minimum of two size corrections so
c   that a target which does not properly belong to the object does not
c   cause an error in centering.  The first estimate is the average of
c   the minimum range/Doppler point and the maximum.  The alternate
c   center is the Doppler-width weighted center.
c
         rcen = 0.5 * ( minr + maxr )
         if ( abs( mean_range ) .lt. abs( rcen ) ) rcen = mean_range
c
         fcen = 0.5 * ( minf + maxf )
         if ( abs( mean_freq )  .lt. abs( fcen ) ) fcen = mean_freq
c
      endif   !   n_inside .gt. 0
c
c-----------------------------------------------------------------------
c
c   Output to GMM file
c
      write ( 85 , '( 2(i5,'',''),(g12.5),9('',''g11.4))' )
     .             npass1 , npass2 , tbar , rcen , fcen ,
     .             a0 , at , ar , af , art , aft , ddsdv
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Automatic Image Selection (AIS) implementation
c
c   First Criterion: Eccentricity calculated from the length and width
c   estimates of tgtbox may redo this in range-frequency coordinates
c
      call ellipse_param ( 0.75 , n_inside , aismajor , aisminor , 
     .                     aismnr , aismnf )
c
      if ( aismajor .ne. 0 ) then
c
         ais(1) = aisminor / aismajor
c
      else
c
         ais(1) = 0
c
      endif
c 
c   Second Criterion: Inverse Harmonic Norm of Doppler Width 
c
c   Third Criterion: Line Detection from Hough Transform 
c
      itgt           = 0
c
      dmin           = 0.0
c
      etamin         = - 180.0 * degrad
c
      etainc         = 360.0 * degrad / ( elength - 1 )
c
      hough_mat(:,:) = 0.0
c
      if ( n_inside .gt. 0 ) then
c
         sumdwdth2 = 0.0
c
         do iloop = 1 , nlist
c
            if ( iflag(iloop) .eq. editac ) then
c
               d_r0 = range(iloop) - mean_range
c
               d_f0 = freq(iloop) - mean_freq
c
               d    =  sqrt( d_r0 ** 2 + d_f0 ** 2 )
c
               if ( d .gt. dmax ) dmax = d
c
            endif
c 
         enddo
c
         dinc = dmax / ( dlength - 1 )
c
         do iloop = 1 , nlist
c
            if ( iflag(iloop) .eq. editac ) then
c
               sumdwdth2 = sumdwdth2 + 1.0 / ( dwdth(iloop) ** 2 )
c
               itgt      = itgt + 1 
c
               d_r0      = range(iloop) - mean_range
c
               d_f0      = freq(iloop)  - mean_freq
c
               eta0      = atan2( d_r0 , d_f0 )
c
               d         = sqrt( d_r0 ** 2 + d_f0 ** 2 )
c
               do jloop = iloop + 1 , nlist
c
                  if ( iflag(jloop) .eq. editac ) then
c
                     d_r     = range(jloop) - range(iloop)
c
                     d_f     = freq(jloop)  - freq(iloop)
c
                     curreta = atan2( d_r , d_f )
c
                     if ( ( d_f0 * d_r - d_r0 * d_f ) .gt. 0.0 ) then
c
                        curreta = curreta - 90.0 * degrad
c
                        if ( curreta .lt. - pi ) 
     .                       curreta = curreta + 2.0 * pi
c
                     else
c
                        curreta = curreta + 90.0 * degrad
c
                        if ( curreta .gt. pi )
     .                       curreta = curreta - 2.0 * pi
c    
                     endif
c                         
                     eidx = nint( ( curreta - etamin ) / etainc ) + 1
c
                     didx = nint( ( d * cos( curreta - eta0 ) - dmin )
     .                            / dinc ) + 1
c
                     if ( eidx .le. 0 ) eidx = eidx + elength
c
                     if ( didx .le. 0 ) didx = didx + dlength
c
                     hough_mat(didx,eidx) = hough_mat(didx,eidx) + 1
c
                  endif
c
               enddo
c
            endif
c
         enddo
c
      endif
c
      if ( itgt .gt. 0 ) then
c
         ais(2) = sqrt( sumdwdth2 / itgt )
c
         ais(3) = 0.0
c
         do didx = 1 , dlength
c
            do eidx = 1 , elength
c
               ais(3) = ais(3) + hough_mat(didx,eidx) ** 2
c
            enddo
c
         enddo
c
         ais(3) =  sqrt( ais(3) )
c
      else
c
         ais(2) = 0.0
c
         ais(3) = 0.0
c
      endif
c
c     Output for MCA
c
      write ( 77 , '(1x,2i8,7f13.3)' ) iframe , n_inside , tbar ,
     .                                 armm , afmm , ais
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      return
      end
C**
C***********************************************************************
C**
      subroutine accmodel ( t , r , f , afunc , acoefs )
C**
C***********************************************************************
C**
c
c       This routine supplies sqrlss with the value of the 
c       basis functions used for the fit of data to :
c   
c         a(r,f,t) = a1 + a2*t + a3*r + a4*f + a5*r*t + a6*f*t
c 
c       so the basis functions are 1,t,r,f
c
c
c   VARIABLES
c
c   t,r,f (real, input)         These scalars contain the point
c                               at which the basis functions are
c                               to be evaluated.
c
c   afunc (real array, output)  This vector will contain on
c                               return the values of the basis
c                               functions evaluated at the point
c                               specified by t,r,f.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      implicit none
c
      integer acoefs
c
      real    t , r , f , afunc(acoefs)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Set all function values to zero
c
      afunc(:) = 0.0
c
c   Depending on the number of coefficients, set the function values
c
      if ( acoefs .gt. 0 ) afunc(1) = 1.0
c
      if ( acoefs .gt. 1 ) afunc(2) = t
c
      if ( acoefs .gt. 2 ) afunc(3) = r
c
      if ( acoefs .gt. 3 ) afunc(4) = f
c
      if ( acoefs .gt. 4 ) afunc(5) = r * t
c
      if ( acoefs .gt. 5 ) afunc(6) = f * t
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      return
      end
C**
C***********************************************************************
C**
      subroutine qrcov ( u , m , n , var , acomax , cov , ri , rit )
C**
C***********************************************************************
C**
c   Computes error covariance matrix from QR factorization as 
c   
c        cov = var * R inverse * R inverse transpose
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      implicit none
c
      integer m                   ! Number rows in matrix u
      integer n                   ! Number columns in matrix u
      real    u(m,n)              ! Upper triangular part = R ,
                                  ! dim(R) = n x n
c
      real    var                 ! Acceleration measurement variance
      integer acomax              ! Max number a_coeffs
      real    cov(acomax,acomax)  ! Error covariance matrix
      real    ri(acomax,acomax)   ! R inverse
      real    rit(acomax,acomax)  ! Transpose of R inverse
c
      integer i , j , k
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Zero R-inverse, R-inverse-transpose, and covariance matrices
c
      ri(:,:)  = 0.0
      rit(:,:) = 0.0
      cov(:,:) = 0.0
c
c   Compute inverse of upper triangular R matrix
c
      do i = 1 , n
c
         ri(i,i) = 1.0 / u(i,i)
c 
      enddo
c
      do j = n , 2 , - 1
c
         do i = j - 1 , 1 , - 1
c
            do k = i + 1 , n
c
               ri(i,j) = ri(i,j) - u(i,k) * ri(k,j)
c
            enddo
c
            ri(i,j) = ri(i,j) / u(i,i)
c
         enddo
c
      enddo
c 
c   Compute transpose of R inverse
c
      do j = 1 , n
c
         do i = 1 , n
c
            rit(i,j) = ri(j,i)
c    
         enddo
c
      enddo
c
c   Compute error covariance
c
      do j = 1 , n
c
         do i = 1 , n
c
            do k = 1 , n
c
               cov(i,j) = cov(i,j) + ri(i,k) * rit(k,j)
c
            enddo
c
            cov(i,j) = cov(i,j) / var
c      
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine acstat ( editlevel , n_inside , nlist , iflag , accel ,
     .                    dwdth , range , freq , time , mean_acc ,
     .                    sd_acc , mean_dw , sd_dw , mean_range ,
     .                    mean_freq , mean_time )
C**
C***********************************************************************
C**
c   Compute Doppler-width weighted acceleration statistics
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      implicit none
c
      integer editlevel , n_inside , nlist , iloop , iflag(nlist) , kk
c
      real    accel(nlist) , dwdth(nlist) , range(nlist) , freq(nlist) ,
     .        time(nlist) , mean_acc , sd_acc , mean_dw , sd_dw ,
     .        mean_range , mean_freq , mean_time , sum_acc ,
     .        sum_sq_acc , sum_dw , sum_sq_dw , sum_range , sum_freq ,
     .        sum_frsq , sum_rsq , sum_time , wt , sumwt
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      n_inside   = 0
c
      sum_acc    = 0.0
      sum_sq_acc = 0.0
c
      sum_dw     = 0.0
      sum_sq_dw  = 0.0
c
      sum_range  = 0.0
      sum_freq   = 0.0
      sum_rsq    = 0.0
      sum_frsq   = 0.0
      sum_time   = 0.0
c
      sumwt      = 0.0
c
      do iloop = 1 , nlist
c
         if ( iflag(iloop) .eq. editlevel ) then
c
            dwdth(iloop) = amax1( 0.1 , dwdth(iloop) )
c
            wt           = 1.0 / dwdth(iloop)
c
            sumwt        = sumwt + wt
c
            sum_acc      = sum_acc    + wt * accel(iloop)
            sum_sq_acc   = sum_sq_acc + wt * accel(iloop) * accel(iloop)
c
            sum_dw       = sum_dw     + dwdth(iloop)
            sum_sq_dw    = sum_sq_dw  + dwdth(iloop) * dwdth(iloop)
c
            sum_range    = sum_range  + wt * range(iloop)
            sum_freq     = sum_freq   + wt * freq(iloop)
            sum_rsq      = sum_rsq    + range(iloop) ** 2
            sum_frsq     = sum_frsq   + freq(iloop) ** 2
            sum_time     = sum_time   + wt * time(iloop)
c
c           Increment Counter Of # Of Points Used For Fit.
c           A point will be used only if its state vector
c           falls within the volume
c
            n_inside     = n_inside + 1
c
         endif
c
      enddo
c
      if ( n_inside .gt. 1 ) then
c
         mean_acc   = sum_acc / sumwt
c
         sd_acc     = sqrt( amax1( 0.0 ,
     .                             ( sum_sq_acc / sumwt )
     .                             - mean_acc ** 2 ) )
c
         mean_dw    = sum_dw / float( n_inside )
         sd_dw      = sqrt( amax1( 0.0 ,
     .                             ( sum_sq_dw / float( n_inside ) )
     .                                   - mean_dw ** 2
     .                      ) )
c
         mean_range = sum_range / sumwt
         mean_freq  = sum_freq  / sumwt
         mean_time  = sum_time  / sumwt
c
      else
c
         if ( n_inside .le. 1 ) sumwt = 1.0
c
         mean_acc   = sum_acc / sumwt
         sd_acc     = 0.0
c
         mean_dw    = sum_dw
         sd_dw      = 0.0
c
         mean_range = sum_range / sumwt
         mean_freq  = sum_freq  / sumwt
         mean_time  = sum_time  / sumwt
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
         subroutine fill_3d_array ( nf , nr , nt , deltaf , deltar ,
     .                              deltat , minf , maxf , minr , maxr ,
     .                              mint , maxt , fstart , rstart ,
     .                              tstart , minacc , maxacc , a0 , at ,
     .                              ar , af , art , aft , dotdot )
c**
c***********************************************************************
c**
c   Fill the 3-D (f,r,t) acceleration array with information from the
c   coefficients computed in the main body of GETACC
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nf , nr , nt , iimin , iimax , jjmin , jjmax , kkmin ,
     .        kkmax , iloop , jloop , kloop , i , j , k , ii , jj , kk
c
      real    a0 , at , ar , af , art , aft , deltaf , deltar , deltat ,
     .        minf , maxf , minr , maxr , mint , maxt , fstart ,
     .        rstart , tstart , minacc , maxacc , dotdot(nf,nr,nt)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      dotdot(:,:,:) = 0.0
c
c   Loop over sub-region (where we have data) to be filled.
c
      kkmin = 1 + nint( ( mint - tstart ) / deltat )
      kkmax = 1 + nint( ( maxt - tstart ) / deltat )
      jjmin = 1 + nint( ( minr - rstart ) / deltar )
      jjmax = 1 + nint( ( maxr - rstart ) / deltar )
      iimin = 1 + nint( ( minf - fstart ) / deltaf )
      iimax = 1 + nint( ( maxf - fstart ) / deltaf )
c
      kkmin = max( 1 , kkmin )
      kkmax = min( nt , kkmax )
      jjmin = max( 1 , jjmin )
      jjmax = min( nr , jjmax )
      iimin = max( 1 , iimin )
      iimax = min( nf , iimax )
c
      kkmin = min( nt , kkmin )
      kkmax = max( kkmax , kkmin )
c
      do kloop = kkmin , kkmax
c
         do jloop = jjmin , jjmax
c
            do iloop = iimin , iimax
c
c   Fill this element with an acceleration value.
c
               dotdot(iloop,jloop,kloop) = 
     .            a0  + 
     .            at  * ( float( kloop - 1 ) * deltat + tstart ) +
     .            ar  * ( float( jloop - 1 ) * deltar + rstart ) +
     .            af  * ( float( iloop - 1 ) * deltaf + fstart ) +
     .            art * ( float( jloop - 1 ) * deltar + rstart ) *
     .                  ( float( kloop - 1 ) * deltat + tstart ) +
     .            aft * ( float( iloop - 1 ) * deltaf + fstart ) *
     .                  ( float( kloop - 1 ) * deltat + tstart )
c
c   Make sure the magnitude of the value is no greater than
c   the max acceleration magnitude or less than the minimum
c   in the target list.
c
                  if ( dotdot(iloop,jloop,kloop) .gt. maxacc ) 
     .                 dotdot(iloop,jloop,kloop) = maxacc 
                  if ( dotdot(iloop,jloop,kloop) .lt. minacc ) 
     .                 dotdot(iloop,jloop,kloop) = minacc 
c
            enddo
c
         enddo
c
      enddo
c
c   Now fill the rest of the dotdot array. 
c
      do k = 1 , nt
c
         do j = 1 , nr
c
            do i = 1 , nf
c
               ii = i
               jj = j
               kk = k
c
               if ( ii .lt. iimin ) ii = iimin
               if ( ii .gt. iimax ) ii = iimax
               if ( jj .lt. jjmin ) jj = jjmin
               if ( jj .gt. jjmax ) jj = jjmax
               if ( kk .lt. kkmin ) kk = kkmin
               if ( kk .gt. kkmax ) kk = kkmax
c
               dotdot(i,j,k) = dotdot(ii,jj,kk)
c
            enddo
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine xyddot ( editac , nlist , iflag , r , rbar , freq ,
     .                    freqbar , accel , accelbar , dwdth , lambda ,
     .                    rsmall , omega_use , xtgt , ytgt , ddot ,
     .                    ddbar , ddsdv , xbar , xsdv , ybar , ysdv )
C**
C***********************************************************************
C**
c   Computes ddot, the time derivative of distance from the scatterers
c   to mean scatterer location and the x and y coordinates of the
c   scatterers
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      implicit none
c
      integer editac , nlist , iflag(nlist) , i , nsum
c
      real    r(nlist) , freq(nlist) , accel(nlist) , dwdth(nlist)
      real    rbar , freqbar , accelbar , sumwt
c
      real    lambda , rsmall , omega_use , xtgt(nlist) , ytgt(nlist) ,
     .        ddot(nlist) , iomegasq
c
      real    ddbar , ddsdv , xbar , xsdv , ybar , ysdv
c
      real    denom , a
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      a        = lambda / 2.0
c
      iomegasq = 1.0 / omega_use ** 2
c
      nsum     = 0
      sumwt    = 0.0
c
      ddbar    = 0.0
      ddsdv    = 0.0
c
      xbar     = 0.0
      xsdv     = 0.0
      ybar     = 0.0
      ysdv     = 0.0
c
      do i = 1 , nlist
c
         if ( iflag(i) .eq. editac ) then
c
            dwdth(i) = amax1( 0.1 , dwdth(i) )
c
            denom    = ( r(i) - rbar ) ** 2 +
     .                 iomegasq * ( a * ( freq(i) - freqbar ) ) ** 2
c
            denom    = sqrt( amax1( denom , rsmall ** 2 ) )
c
            ddot(i)  = ( a / denom ) * ( freq(i) - freqbar ) *
     .                 ( ( r(i) - rbar ) +
     .                   a * iomegasq * ( accel(i) - accelbar ) )
c
            xtgt(i)  = a * ( freq(i) - freqbar ) / omega_use
c
            ytgt(i)  = r(i) - rbar
c
            nsum     = nsum + 1
c
            sumwt    = sumwt + 1.0 / dwdth(i)
c
            ddbar    = ddbar +   ddot(i)        / dwdth(i)
            ddsdv    = ddsdv + ( ddot(i) ** 2 ) / dwdth(i)
c
            xbar     = xbar  +   xtgt(i)        / dwdth(i)
            xsdv     = xsdv  + ( xtgt(i) ** 2 ) / dwdth(i)
c
            ybar     = ybar  +   ytgt(i)        / dwdth(i)
            ysdv     = ysdv  + ( ytgt(i) ** 2 ) / dwdth(i)
c
         else
c
            ddot(i)  = 0.0
c
            xtgt(i)  = 0.0
c
            ytgt(i)  = 0.0
c
         endif
c
      enddo
c
      if ( nsum .gt. 1 ) then
c
         ddbar = ddbar / sumwt
         ddsdv = sqrt( amax1( 0.0 , ( ddsdv / sumwt ) - ddbar ** 2 ) )
c
         xbar  = xbar  / sumwt
         ybar  = ybar  / sumwt
c
         xsdv  = sqrt( amax1( 0.0 , ( xsdv / sumwt ) - xbar ** 2 ) )
         ysdv  = sqrt( amax1( 0.0 , ( ysdv / sumwt ) - ybar ** 2 ) )
c
      else
c
         ddsdv = 0.0
c
         xsdv  = 0.0
         ysdv  = 0.0
c
      endif
c
      return
      end
