C**
C***********************************************************************
C**
      subroutine getacc_rt ( nf , nr , nt , dotdot , fstart , rstart ,
     .                       tstart , deltaf , deltar , deltat ,
     .                       athrsh , rcen , fcen , iframe )
C**
C***********************************************************************
C**
c  =====================================================================
c   Purpose:  
c
c              Given acceleration data that has been sampled at
c              irregular values of (f,r,t), perform a singular value
c              decomposition to achieve the best fit of the data to
c              the model:  4 coefficent model for Telephonics ISAR
c
c              A(F,R,T) = A0 + AT*T + AR*R + AF*F 
c
c              and then use the best fit to fill the acceleration
c              array, dotdot(f,r,t).
c
c              This routine also uses a simplified least squares fit to 
c              determine the coefficents defined above.
c              
c  ---------------------------------------------------------------------
c  Passed Variables 
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
c   a0,at,ar,af                     The fit coefficients
c
c   --------------------------------------------------------------------
c   Other Local Variables 
c   --------------------------------------------------------------------
c
c   a_coeffs (real array)           Vector containing fit coefficients
c                                   of the fit function (these are
c                                   returned by the call to svdfitmv).
c                                   These are the A0, AF, AR, AT, ART,
c                                   AFT discussed above
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
c                                   which have magnitudes < athrsh.
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
      include          'realtime.h'    !  Real-time parameters
c
      include          'updates.h'
c
c-----------------------------------------------------------------------
c              
      integer           nf , nr , nt , iframe
c
      real              deltaf , deltat , deltar , dotdot(nf,nr,nt) ,
     .                  fstart , rstart , tstart , sd_acc , mean_acc ,
     .                  mean_time , sd_dw , mean_dw , athrsh , rcen ,
     .                  fcen
c
c   Variables used in sqrank,sqrlss
c
      integer           acomax            !  Maximum no. of acceleration
      parameter       ( acomax = 4 )      !  coefficients
c
      real              b(rt_max_in+8) , tmp , bw(rt_max_in+8) ,
     .                  qrtol , u1((rt_max_in+8)*acomax) ,
     .                  u1w((rt_max_in+8)*acomax) , work1(rt_max_in+8) ,
     .                  work2(rt_max_in+8,rt_max_in+8)
c
c   Local Variables
c
      real              minf , minr , mint , maxf , maxr , maxt ,
     .                  minacc , maxacc , tbeg , tend , accdum ,
     .                  tinflu , mean_range , mean_freq
c
      integer           n_inside , iloop , i , j , itgt
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
      integer           npass1 , npass2 , npass3
c
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c
c   Begin Executable Section
c
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c
c   Don't update the centroid estimate
c
      rcen   = 0.0
      fcen   = 0.0
c
c   Zero target selection flag so that we can re-consider all targets
c   for this frame
c
      do iloop = 1 , nlist
 
         iflag(iloop) = 0
c
      enddo
c
c   If curve equals zero, fill accel array with zeroes and return
c
      if ( curve .eq. 0 .or. acoefs .eq. 0 ) then 
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
c                          Stage 1
c
c   Loop over all elements to find min/max of f,r,t. (slow time)
c   Only consider elements that fall within the volume and edit the
c   points according to three criteria.
c
c   Use an influence time equal to the maximum of the integration time
c   and tlarge.  Also, do not use estimates earlier than half this.
c
      tbar   = tstart + 0.5 * float( nt ) * deltat
      tinflu = amax1( 0.5 , float( nt ) * deltat )
      tbeg   = amax1( 0.0 , tbar - 1.0 * tinflu ) ! 2 * integration time
      tend   = amax1( 0.0 , tbar + 1.0 * tinflu ) ! 2 * integration time
c
c-----------------------------------------------------------------------
c
c                          First Edit Pass
c
      do iloop = 1 , nlist
c
c   Use only data from subroutine TARGET if curve = 1; use only data
c   from subroutine PGA if curve = 2; otherwise use both
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
         if ( check .and. ( abs( accel(iloop) ) .le. athrsh ) .and.
     .      ( time(iloop) .ge. tbeg .and. time(iloop) .le. tend ) ) then
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
      call acstat ( 2 , n_inside , nlist , iflag , accel , dwdth ,
     .              range , freq , time , mean_acc , sd_acc , mean_dw ,
     .              sd_dw , mean_range , mean_freq , mean_time )
c
      npass2 = n_inside
c
      npass3 = npass2        !  No third pass
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                          Stage 2
c
c   Add the data points to the linear least squares data matrix and
c   obtain minima and maxima of the data within the volume.
c
c   Initialize the mins and maxs and number of points found inside the
c   outer volume.
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
      n_inside = min( n_inside , rt_max_in )
c
      if ( n_inside .gt. 0 ) then
c
         itgt = 0
c
         do iloop = 1 , nlist
c
            if ( iflag(iloop) .eq. editac .and. itgt .lt. rt_max_in )
     .                                     then
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
               call accmodel ( time(iloop) , range(iloop) ,
     .                         freq(iloop) , afunc , acoefs )
c
c              Add current point to u,b arrays needed for sqrlss.
c
               tmp     = 1.0 / dwdth(iloop) 
c
               do j = 1 , acoefs
c
                  u1((j-1)*(n_inside+8)+itgt) = afunc(j) * tmp
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
c   The points with state vectors that lie in the specified
c   outer volume have been identified. These will be used for
c   obtaining the fit parameters which characterize the best fit
c   of the input acceleration data to the linear model for the
c   acceleration.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                          Stage 3
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
            u1((j-1)*(n_inside+8)+itgt+1) = afunc(j) * tmp
         enddo
         b(itgt+1) = accdum
c
         call accmodel ( mint , minr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+2) = afunc(j) * tmp
         enddo
         b(itgt+2) = accdum
c
         call accmodel ( mint , maxr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+3) = afunc(j) * tmp
         enddo
         b(itgt+3) = accdum
c
         call accmodel ( mint , maxr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+4) = afunc(j) * tmp
         enddo
         b(itgt+4) = accdum
c
         call accmodel ( maxt , minr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+5) = afunc(j) * tmp
         enddo
         b(itgt+5) = accdum
c
         call accmodel ( maxt , minr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+6) = afunc(j) * tmp
         enddo
         b(itgt+6) = accdum
c
         call accmodel ( maxt , maxr , minf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+7) = afunc(j) * tmp
         enddo
         b(itgt+7) = accdum
c
         call accmodel ( maxt , maxr , maxf , afunc , acoefs )
         do j = 1 , acoefs
            u1((j-1)*(n_inside+8)+itgt+8) = afunc(j) * tmp
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
               if ( u1((j-1)*(n_inside+8) +i) .gt. qrtol ) 
     .              qrtol = u1((j-1)*(n_inside+8)+i)
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
c
c-----------------------------------------------------------------------
c
c   Jeff Richard's fast version of the least squares solver
c
         call fastfit_rt ( u1 , u1w , b , bw , n_inside + 8 , acomax ,
     .                     work1 , work2 )
c
c-----------------------------------------------------------------------
c
c   Fill variables with the returned coefficients.
c
         a0  = b(1)
         at  = b(2)
         if ( acomax .gt. 2 ) ar = b(3)
         if ( acomax .gt. 3 ) af = b(4)
         art = 0.0
         aft = 0.0
c
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
c
c   If no points were found within the limits, then fill with zeros and
c   bail out
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
c   Output to plot file
c
      write ( 85 , '(3(i5,'',''),(g12.5),8('',''g11.4))' )
     .             npass1 , npass2 , npass3 , tbar ,
     .             rcen , fcen , a0 + at * tbar , at ,
     .             ar + art * tbar , af + aft * tbar , art , aft

c
      return
      end
