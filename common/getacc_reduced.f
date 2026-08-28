C**
C***********************************************************************
C**
      subroutine getacc ( nf , nr , nt , dotdot , fstart , rstart ,
     .                    tstart , deltaf , deltar , deltat , athrsh ,
     .                    rcen , fcen )
C**
C***********************************************************************
C**
c
c                      ISAR Global Motion Model
c
c                             Version 3.32
c
c                             July 16, 2004
c
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
c      implement the Global Motion Model of the SAIC Inverse SAR
c      processor.
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
c-----------------------------------------------------------------------
c              
      integer           nf , nr , nt
c
      real              deltaf , deltat , deltar , dotdot(nf,nr,nt) ,
     .                  fstart , rstart , tstart , sd_acc , mean_acc ,
     .                  sd_dw , mean_dw , athrsh , rcen , fcen ,
     .                  mean_time
c
c   Variables used in sqrank, sqrlss
c
      integer           acomax            !  Maximum no. of acceleration
      parameter       ( acomax = 6 )      !  coefficients
c
      integer           jpvt(acomax) , kr
c
      real              u(nlist+8,acomax) , b(nlist+8) , tmp , dvmax ,
     .                  tol , qraux(acomax) , qrwork(acomax) , rsmall ,
     .                  qrrsd(nlist+8) , qrtol 
c
      parameter       ( tol = 1.e-5 )     !  Tolerance for convergence
c
c   Local Variables
c
      real              minf , minr , mint , maxf , maxr , maxt , 
     .                  minacc , maxacc , tbeg , tend , accdum , lf ,
     .                  tinflu , mean_range , mean_freq , lt , lr
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
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Begin Executable Section
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   If no targets pass, then don't update the centroid estimate
c
      rcen   = 0.0
      fcen   = 0.0
c
c   Zero target selection flag so that we can re-consider all targets
c   for this frame
c
      do iloop = 1 , nlist
c
c   Set the flags to zero
c
         iflag(iloop) = 0
c
      enddo
c
c   Use an influence time equal to the maximum of the integration time
c   and tlarge.
c
      tbar   = tstart + 0.5 * float( nt - 1 ) * deltat
      tinflu = float( nt ) * deltat
c
      tbeg   = amax1( 0.0 , tbar - 1.0 * tinflu )
      tend   = amax1( 0.0 , tbar + 1.0 * tinflu )
c
c   If curve equals zero, fill acceleration array with zeroes and return
c
      if ( ( curve .eq. 0     ) .or. ( acoefs .eq. 0 ) .or.
     .     ( tbar .lt. tinflu ) ) then 
c
c   Fill with zeroes
c
         dotdot(:,:,:) = 0.0
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
c   points according to two criteria.
c
c-----------------------------------------------------------------------
c
c                          First Edit Pass
c
c   This stage is always implemented.
c
      do iloop = 1 , nlist
c
c   Use only data from subroutine TARGET if curve = 1; use only data
c   from subroutine PGA if curve = 2; otherwise use both.  SAIC has
c   found by experience that it is best to ignore the targets from the
c   front-end since they are based on very short integration times.
c   Thus, we recommend that curve should be 2, 5, or 8.  All of these
c   values use only scatterers estimated from the fine resolution image
c   frames.or the sub-image stack
c
         tcheck = ( ( mod( curve , 3 ) .eq. 0 ) .or.
     .              ( mod( curve , 3 ) .eq. 1 ) ) .and.
     .            ( source(iloop) .eq. 't' )
c
         pcheck = ( ( mod( curve , 3 ) .eq. 0 ) .or.
     .              ( mod( curve , 3 ) .eq. 2 ) ) .and.
     .            ( source(iloop) .eq. 'p' )
c
         scheck = ( ( mod( curve , 3 ) .eq. 0 ) .or.
     .              ( mod( curve , 3 ) .eq. 2 ) ) .and.
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
     .        ( time(iloop)  .le. tend          ) ) then
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
            if ( iflag(iloop) .eq. 2 ) then
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
         a_coeffs(:) = 0.0
c
         qrtol = tol / qrtol
c
c-----------------------------------------------------------------------
c
         call sqrank ( u , nlist + 8 , n_inside + 8 , acoefs , qrtol ,
     .                 kr , jpvt , qraux , qrwork )
c
         call sqrlss ( u , nlist + 8 , n_inside + 8 , acoefs , kr ,
     .                 b , a_coeffs , qrrsd , jpvt , qraux )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                       Stage 3 - Array Fill
c
c   Fill variables with the returned coefficients.
c
         a0  = a_coeffs(1)
         at  = a_coeffs(2)
         ar  = a_coeffs(3)
         af  = a_coeffs(4)
         art = a_coeffs(5)
         aft = a_coeffs(6)
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
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      endif   !   n_inside .gt. 0
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Output diagnostic information
c
      write ( 7 , * )
      write ( 7 , '(a22,2f12.3)' ) ' Rcen, Fcen:         ' , rcen , fcen
      write ( 7 , * )
c
c   Output to plot file
c
      write ( 85 , '(a3,'','',2(i5,'',''),(g12.5),8('',''g11.4))' )
     .             '#g#' , npass1 , npass2 , tbar , rcen , fcen ,
     .              a0 , at , ar , af , art , aft
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
      integer editlevel , n_inside , nlist , iloop , iflag(nlist)
c
      real    accel(nlist) , dwdth(nlist) , range(nlist) , freq(nlist) ,
     .        time(nlist) , mean_acc , sd_acc , mean_dw , sd_dw ,
     .        mean_range , mean_freq , mean_time , sum_acc ,
     .        sum_sq_acc , sum_dw , sum_sq_dw , sum_range , sum_freq ,
     .        sum_time , wt , sumwt
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
      sum_time   = 0.0
c
      sumwt      = 0.0
c
      do iloop = 1 , nlist
c
         if ( iflag(iloop) .eq. editlevel ) then
c
            wt         = 1.0 / dwdth(iloop)
c
            sumwt      = sumwt + wt
c
            sum_acc    = sum_acc    + wt * accel(iloop)
            sum_sq_acc = sum_sq_acc + wt * accel(iloop) * accel(iloop)
c
            sum_dw     = sum_dw     + dwdth(iloop)
            sum_sq_dw  = sum_sq_dw  + dwdth(iloop) * dwdth(iloop)
c
            sum_range  = sum_range  + wt * range(iloop)
            sum_freq   = sum_freq   + wt * freq(iloop)
            sum_time   = sum_time   + wt * time(iloop)
c
c           Increment Counter Of # Of Points Used For Fit.
c           A point will be used only if its state vector
c           falls within the volume
c
            n_inside   = n_inside + 1
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
