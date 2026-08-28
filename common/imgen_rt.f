C**
C***********************************************************************
C**
      subroutine imgen_rt ( cframe , nfa , nfr , dff , drf , csbimg ,
     .                      nakeep , nsa , nsr , dfc , drc , fc0 ,
     .                      actual , nabuff , nafill , iisub , cac ,
     .                      mprat , wtaf , work , nwork , lambda ,
     .                      dotdot , mrrat , nspfv , dtv , tframe ,
     .                      nmp , cv , ct , pg , nrcell , nrfcell ,
     .                      ctmp , screen , nafill_a , rc , crc )
C**
C***********************************************************************
C**
c   This routine supports the ISAR movie processor by forming a complex
c   image frame using the Triple-FFT algorithm.  Specifically, the
c   routine implements the third FFT in that method.  The first two
c   stages are range compression and coarse doppler compression
c   (sub-image formation).
c
c   Written in 1994 to support the RDRTec/Telephonics Inverse-SAR (ISAR)
c   processor, IMGEN has grown in uses.  For example, it is a key
c   routine in algorithms for refocusing moving targets and for
c   detecting vibrating targets in standard SAR images.
c
c   There are two stages to the routine.  First, the image is calculated
c   at 'master particles' spaced evenly in the image output frequency
c   coordinate.  This process involves a time integration (sum) of
c   complex sub-image values interpolated to the correct range and
c   frequency appropriate to a test particle following the estimated
c   paths of scatterers in the scene.  Second, the time series of
c   complex values which comprise this sum is Fourier transformed to
c   form estimates of nearby image points.  These estimates are blended
c   into the image by linear interpolation to fill all desired frequency
c   bins.
c
c   To save computer time the above process is done in groups of 'mrrat'
c   range cells which use the same interpolation logic and focus
c   parameters.
c
c   At first glance IMGEN appears to implement a complex algorithm.
c   However, in most practical applications, this algorithm is extremely
c   fast and accurate.  IMGEN achieves its efficiency by carefully
c   separating the terms in the image integration into fast-varying and
c   slowly-varying categories.  Slowly varying terms are approximated as
c   locally constant values.  Since the fast varying terms are largely
c   due to linear phase terms, FFT methods are used to compute them.
c
c-----------------------------------------------------------------------
c
c   MODIFICATION FOR EXACT TIME MODE ( Sept. 26, 1995 )
c
c   In the standard mode, the frame is calculated from a time series of
c   complex values interpolated in range and frequency from the
c   sub-image array but at the same times as the sub-images are defined;
c   the integration time is constrained to be an even integer times
c   the basic time between the sub-images.
c
c   Since a frame may be required at an arbitrary time and for an
c   arbitrary integration time, this routine was modified to form an
c   image frame from a number of 'virtual' sub-images interpolated in
c   time from the sub-image array.
c   
c-----------------------------------------------------------------------
c
c   MODIFICATIONS FOR REAL-TIME ( July 16, 1998 )
c
c
c        1.  Removed the stretch corrections.  Even when the real-time
c            data is collected in stretch mode, we will ignore these
c            terms.
c
c        2.  Removed the small contribution of the focus to the range
c            offset.
c
c        3.  Incorporated a correction to the amplitude of the Master
c            Particle images which depends on the offset from the exact
c            frequency.  The parameter, 'roll_off' is a input to this
c            algorithm declared in 'realtime.h'.
c
c        4.  Removed the 'exact time' option which requires time
c            interpolation of the sub-images
c
c        5.  Broke the algorithm up into two stages - 
c
c               A.  Loading the array of master particle time series
c
c               B.  FFTing this array and adding it into the image
c
c            This prepares the data for an autofocus stage to be
c            inserted between the two stages
c
c        6.  Wrote an autofocus algorithm based on the master particle
c            time series data instead of the complex image.
c
c            This code is done between the two original stages.  The
c            three stages are now:
c
c            Stage 1:  Load the array of master particles (t,f,r)
c
c            Stage 2:  Compute the phase correction by integrating the
c                      Phase Gradient Autofocus sums over range and
c                      coarse frequency
c
c            Stage 3:  Apply the phase correction, FFT the master
c                      particle data with time and add it into the
c                      complex image
c
c        7.  Modified the load and PGA stages to compute the autofocused
c            time series for two frames at once.
c
c            First, the algorithm loads master particle array augmented
c            in time by the time between frames and centered between an
c            odd/even pair of frames.  For example, at 64 sub-images
c            per second and a frame rate of 4 Hz, the normal time series
c            is augmented by 16.
c
c            Second, the augmented time series is autofocused.
c
c            Finally, the first part of the focused time series is used
c            to form the odd frame and the second part of the time
c            series is used to form the even frame.
c
c            The main consequence of this process is that a scatterer in
c            the highest Doppler cells may appear to wiggle back and
c            forth in range between the odd and even frames.  This is
c            because the time used for the range walk calculation is off
c            by one half the time difference between frames.  There is
c            almost no cost to the focus of the two frames since the PGA
c            calculation is valid for the entire augmented time series.
c
c-----------------------------------------------------------------------
c 
      implicit none
c
c***********************************************************************
c
c        Array sizes and raw data length passed as arguments
c
      integer     nsr , nakeep , nsa , nfr , nfa , nwork , nspfv ,
     .            actual , nabuff , nafill , nafill_a
c
c        nsr       :  Range cells in a sub-image
c
c        nakeep    :  Azimuth cells in a sub-image kept
c        nsa       :  Original no. of azimuth cells
c
c        nfr       :  Range cells in a frame
c        nfa       :  Azimuth cells in a frame 
c
c        nwork     :  Real numbers available in work array
c
c-----------------------------------------------------------------------
c
c                  Sub-image time counters and arrray sizes
c
c        nspfv     :  No. of virtual sub-images per frame
c
c        actual    :  Sub-images available for this image (<=nabuff)
c
c        nabuff    :  Size of sub-image buffer
c
c        nafill    :  Sub-images per frame including zero-fill
c
c        nafill_a  :  Sub-images per frame including zero-fill and
c                     augmentation to allow for the computation of two
c                     frames at once
c
c***********************************************************************
c
      integer     ifa , imp , mprat , fminus , fplus , rnear , ts ,
     .            ish , iisub , mrrat , imr , nfrbig , ifrbig , nn(1) ,
     .            ifirst , ifr , nmp , fshift , isbv , jfr , augment
c
c***********************************************************************
c
      real        dff , drf , dfc , drc , rsmin , rsmax , twopi , vf ,
     .            fmp , dfmp , fdelta , delrmp , ltime , fnew , fold ,
     .            wtaf(nspfv) , fc0 , fmp0 , lambda , fmpnom , pprime ,
     .            dtv , tframe , work(nwork) , fprime , rprime , dd ,
     .            dotdot(1+nfa/mprat,nfr/mrrat,nspfv+nafill_a-nafill) ,
     .            pi
c
      complex     cframe(nfa,2*nfr) , csbimg(nsr,nakeep,nabuff) ,
     .            cdfmt , cdfpt , cdelt
c
c   cpfast is a table look-up version of the quantity cexp(j*2*pi*arg)
c
      complex     cpfast
c
c   ac and cac are real and complex names for an equivalenced work
c   vector used in the call to FOURT
c
      complex     cac(nafill_a,nmp,nfr)
c
      real        rc(2,2*nafill_a)             !  Work array for Doppler
      complex     crc(2*nafill_a)              !  compression
c
      include    'realtime.h'                  !  Real-time parameters
c
      integer     r_offset , ifrc , ifrs , t_offset , nrcell(nfr) ,
     .            nrfcell(nfr) , iter , total_cells1 , total_cells2 ,
     .            total_cells3
c
      real        rmpc , cfmax , correl , pwr , pgrms , degrad
c
      complex     covt , cv(nafill_a) , ct(nafill_a,2) ,
     .            ctmp(nafill_a,nfr)
c
      real        pg(nafill_a,4) , screen(nfr,2) , acc , rtgt , ftgt ,
     .            d , dwtgt , atgt , rt_rough
c
      integer     rt_maxit , rt_maxtgt , nacc , middle , np2 , edge ,
     .            nspfv_a , centerdd , idd , na_off
c
      parameter ( rt_maxit = 10 , rt_maxtgt = 20 , middle = 10 ,
     .            rt_rough = 0.175 )
c
      include    'updates.h'
c
c***********************************************************************
c
      pi     = atan2( 0.0 , - 1.0 )
c
      twopi  = 2.0 * pi
c
      degrad = 360.0 / twopi
c
c-----------------------------------------------------------------------
c
c                 Compute Master particle parameters
c
c   Master particles are frequencies in the output complex image at
c   which the image value is evaluated using the full finite difference
c   sum for the contributions of the sub-images.  The first master
c   particle is placed on the most negative frequency in the image.  The
c   separation between them is 'mprat' fine resolution bins.  So that
c   there is no loss of image intensity at the high frequency boundary,
c   the last master particle is placed at beyond the boundary.
c
      fmp0  = - 0.5 * float( nfa ) * dff
      dfmp  = dff * float( mprat )
c
c   Min and max range of sub-images and fine resolution images
c
      rsmin =       - float( nsr / 2 ) * drc
      rsmax = rsmin + float( nsr - 1 ) * drc
c
      nn(1) = nafill                !  FFT length = sub-images per frame
c                                   !  including zero-fill
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Generate the master particle blending function.  The function
c   starts at 0.0 for ifa = 1, increases to 1.0 at ifa = 1 + mprat,
c   and then tapers off symmetrically to 0.0 at ifa = 1 + 2 * mprat.
c   In between 0 and 1 the function generally is larger than the
c   linear trend to account for the roll-off of the sub-image IPR.
c
      if ( mod( rt_img , 2 ) .eq. 0 ) then
c
c   Use Bill Caputi's blending function
c
         do ifa = 2 , mprat
c
c   Linear starting at 1.0 in the middle
c
            work(ifa+mprat)   = 1.001 - float( mprat + 1 - ifa ) /
     .                                  float( mprat )
c
c   Use the linear term to build the full weight function
c
            work(ifa+mprat)   = sin( pi * work(ifa+mprat) ) /
     .                             ( pi * work(ifa+mprat) )
c
            work(mprat+2-ifa) = work(ifa+mprat)
c
c
         enddo
c
         work(1+mprat)   = 1.0
c
         work(1)         = 0.0
c
         work(1+2*mprat) = 0.0
c
      else
c
c   Interpolation weights used for blending frequencies at the master
c   particles into the image frame.  This includes the correction for
c   the roll-off of the sub-image response function.  The raised cosine
c   gives more weight to the nearer point versus the linear formula.
c
         do ifa = 1 , mprat + 1
c
c   Linear starting at 1.0
c
            work(ifa+mprat)   = float( mprat + 1 - ifa ) /
     .                          float( mprat )
c
c   Use the linear term to build the full weight function
c
            work(ifa+mprat)   = 0.5
     .         + 0.5 * cos( 0.5 * twopi * ( 1.0 - work(ifa+mprat) ) )
     .         + roll_off * work(ifa+mprat) * ( 1.0 - work(ifa+mprat) )
c
            work(mprat+2-ifa) = work(ifa+mprat)
c
         enddo
c
      endif
c
c   Offset in the master particle FFT array to the points to be mapped
c   into the image
c
      na_off   = ( nafill / 2 ) - mprat
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Offset in the sub-image range cells for each fine resolution image
c   range cell.  For this algorithm to work well this number should be
c   exactly an integer.
c
      r_offset = nint( drf / drc )
c
c-----------------------------------------------------------------------
c
c   Initialize image to zero
c
      do ifr = 1 , nfr
c
         do ifa = 1 , nfa
c
            cframe(ifa,ifr)     = cmplx( 0.0 , 0.0 )
            cframe(ifa,ifr+nfr) = cmplx( 0.0 , 0.0 )
c
         enddo
c
         do imp = 1 , nmp
c
            do ifa = 1 , nafill_a
c
               cac(ifa,imp,ifr) = cmplx( 0.0 , 0.0 )
c
            enddo
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Assume that nfr is an integer times mrrat
c
c   nfrbig is the number of range regions of mrrat range cells.
c   Over each of these blocks the accelerations and the frequency
c   interpolation parameters are assumed constant.
c
      nfrbig   = nfr / mrrat
c
      vf       = 0.5 * lambda               !  Nominal (non-stretch)
c
      t_offset = max( 0 ,                   !  Number of zeroes at the
     .        ( nafill_a - actual ) / 2 )   !  start and end of time
c
      augment  = nafill_a - nafill          !  Extra sub-images for 2
c                                           !  frames at once
c
      nspfv_a  = min( actual , nspfv + augment )
c
      centerdd = ( nspfv + nafill_a - nafill ) / 2
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                             STAGE 1
c
c-----------------------------------------------------------------------
c
c                Loop over output frame range cells
c
      do 100 ifrbig = 1 , nfrbig
c
c-----------------------------------------------------------------------
c
c   Loop over master particles
c
         do 200 imp = 1 , nmp
c
            fmpnom = fmp0 + dfmp * float( imp )
c
c   Back-track to find the deviations of Doppler frequency and range
c   for this master particle at the start of the coherent integration
c   time from the acceleration field and the values of nominal Doppler
c   and range (the values at the mid-point of the integration time)
c
            fprime = 0.0
            rprime = 0.0
c
            idd = 1 + centerdd
c
            do isbv = 1 , nspfv_a / 2
c
               dd     = 0.5 * ( dotdot(imp+1,ifrbig,idd) +
     .                          dotdot(imp+1,ifrbig,idd+1) )
               fnew   = fprime
               fprime = fprime - dd * dtv
               rprime = rprime - 0.5 * ( fprime + fnew ) * vf * dtv
c
               idd    = idd - 1
c
            enddo
c
c   Form time series of sub-image values for master particles
c
            do 300 isbv = 1 , nspfv_a
c
c   Determine the position in the sub-image array to get the data from
c
               ts     = iisub - nspfv_a + isbv
c
               if ( ts .lt. 1 ) ts = ts + nabuff
c
c   Determine the position in the FFT array to load the data into.
c   Putting the data into the center of the array centers the data on
c   the Nyquist point rather than on the DC.  This makes the complex
c   image with sign flips for alternate points in cross-range.  However,
c   since the autofocus is done already, only the image magnitude is
c   required.  The advantage is that the PGA stage can access the time
c   points continuously rather than having to access the negative times
c   at the end of the array and then swithch to the zero and positive
c   times at the start of the array.
c
               ish    = isbv + t_offset
c
c   Calculate the actual frequency for this time for the master particle
c
               fmp    = fmpnom + fprime
c
c   Calculate local time for use in phase and range walk calculations
c
               ltime  = dtv * float( isbv - 1 - nspfv_a / 2 )

c   Calcculate the phase correction for path curvature
c
               pprime = ( rprime / vf ) - fprime * ltime
c
c-----------------------------------------------------------------------
c
c   Calculate interpolation parameters for frequency
c
               if ( fmp .gt. fc0 ) then
c
                  fminus = 1 + ifix( ( fmp - fc0 ) / dfc )
                  fplus  = fminus + 1
                  fdelta = ( ( fmp - fc0 ) / dfc ) - float( fminus - 1 )
c
               else
c
                  fplus  = 1 + ifix( ( fmp - fc0 ) / dfc )
                  fminus = fplus - 1
                  fdelta = float( 2 - fplus ) - ( ( fc0 - fmp ) / dfc )
c
               endif
c
c   Assume that the sub-image is periodic if all sub-image Doppler cells
c   have been kept.  Otherwise, use the nearest boundary Doppler cell.
c
               if ( fminus .lt. 1 ) then
c
                  if ( nakeep .ne. nsa ) then
c
                     fminus = 1
                     fplus  = 1
c
                  else
c
                     fshift = nsa * ( 1 - fminus / nsa )
                     fminus = 1 + mod( fshift + fminus - 1 , nsa )
                     fplus  = 1 + mod( fshift + fplus - 1 , nsa )
c
                  endif
c
               endif
c
               if ( fplus .gt. nakeep ) then
c
                  if ( nakeep .ne. nsa ) then
c
                     fminus = nakeep
                     fplus  = nakeep
c
                  else   
c
                     fminus = 1 + mod( fminus - 1 , nsa )
                     fplus  = 1 + mod( fplus - 1 , nsa )
c
                  endif
c
               endif
c
c   Calculate the phase correction to the plus and minus frequency
c   bins
c
               cdelt  = cpfast( ltime * fmp + pprime )
c
               cdfmt  = cdelt * ( 1.0 - fdelta )
               cdfpt  = cdelt * fdelta
c
c   Increment the range to each master particle
c
c   In real-time mode, drop the rprime term so that the range
c   interpolation process is independent of range
c
               delrmp = vf * fmpnom * ltime
c
c   Apply these accelerations and frequency interpolation parameters
c   to 'mrrat' range cells at a time
c
c-----------------------------------------------------------------------
c
c         Nearest neighbor interpolation used in real-time mode
c
c              Index to center cell of this range block
c
               ifrc   = mrrat * ( ifrbig - 1 ) + 1 + mrrat / 2
c
c              Range for this master particle at this time
c
               rmpc   = delrmp + drf * float( ifrc - 1 - nfr / 2 )
c
c              Index into sub-image array for this point
c
               ifrs   = 1 + nint( ( rmpc - rsmin ) / drc )
c
               do imr = 1 , mrrat
c
                  ifr   = mrrat * ( ifrbig - 1 ) + imr
c
c   Calculate nearest neighbor for range
c
                  rnear = ifrs + r_offset * ( imr - 1 - mrrat / 2 )
c
c-----------------------------------------------------------------------
c
c   Finally, load the FFT array with the frequency-interpolated value of
c   the sub-image
c
                  cac(ish,imp,ifr) = cdfpt * csbimg(rnear,fplus,ts)
     .                             + cdfmt * csbimg(rnear,fminus,ts)
c
               enddo     !  Loop over mrrat range cells
c
c   Update the deviation of the Doppler and range from their nominal
c   values (except for last time)
c
               if ( isbv .ne. nspfv_a ) then
c
                  dd     = 0.5 * ( dotdot(imp+1,ifrbig,idd) +
     .                             dotdot(imp+1,ifrbig,idd+1) )
                  fold   = fprime
                  fprime = fprime + dd * dtv
                  rprime = rprime + 0.5 * ( fold + fprime ) * vf * dtv
c
                  idd    = idd + 1
c
               endif
c
  300       continue     !  End of loop over time for a master particle
c
  200    continue        !  End of loop over master particles
c
  100 continue           !  End of loop over groups of image range cells
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                             STAGE 2
c
c   End of loops for loading mrrat FFT vectors - now compute autofocus
c   information
c
      if ( rt_pga .ne. 0 ) then    !  Do only if real-time flag set
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                            STAGE 2-A 
c
c   Screen the range and doppler cells to find the most significant ones
c
         total_cells1 = 0
c
         total_cells2 = 0
c
         total_cells3 = 0
c
         do ifr = 1 , nfr
c
            nrfcell(ifr)  = 0
c
            screen(ifr,2) = 0.0
c
         enddo
c
         edge = ( nmp - ( 3 * nmp ) / 4 ) / 2
c
         do ifr = 1 + nfr / 8 , nfr - nfr / 8
c
            cfmax        = 0.0
c
            do imp = 1 + edge , nmp - edge
c
               covt = 0.0
c
               pwr  = 0.0
c
c   Use middle pulse pairs
c
               np2  = middle / 2
c
               do ts = - np2 + nafill_a / 2 , np2 - 1 + nafill_a / 2
c
                  covt = covt + cac(ts,imp,ifr) *
     .                   conjg( cac(ts+1,imp,ifr) )
c
                  pwr  = pwr  + cac(ts,imp,ifr) *
     .                   conjg( cac(ts,imp,ifr) )
c
               enddo
c
               pwr  = pwr  + 0.5 * ( cac(np2+nafill_a/2,imp,ifr) *
     .                        conjg( cac(np2+nafill_a/2,imp,ifr) )
     .                             - cac(-np2+nafill_a/2,imp,ifr) *
     .                        conjg( cac(-np2+nafill_a/2,imp,ifr) ) )
c
               if ( pwr .ne. 0.0 ) then
c
                  correl = cabs( covt )
c
               else
c
                  correl = 0.0
c
               endif
c
               if ( ( ( correl  / pwr ) .gt. 0.67 ) .and.
     .              correl .gt. cfmax ) then
c
                  nrfcell(ifr) = imp
c
                  cfmax        = correl
c
               endif
c
            enddo    !  Master particles
c
            if ( nrfcell(ifr) .gt. 0 ) total_cells1 = total_cells1 + 1
c
            screen(ifr,2) = cfmax
c
         enddo       !  Range cells
c
c-----------------------------------------------------------------------
c
         if ( quiet .gt. 1 )
     .   write ( 6 , * ) total_cells1 , total_cells2 ,
     .                   total_cells3
         write ( 7 , * ) total_cells1 , total_cells2 ,
     .                   total_cells3
c
c   Thin the range cells by eliminating those cells which have scores
c   less than neighboring cells +/- 2 away
c
         if ( total_cells1 .gt. 0 ) then
c
            do ifr = 3 , nfr - 2
c
               if ( screen(ifr,2) .lt. screen(ifr-2,2) .or.
     .              screen(ifr,2) .lt. screen(ifr-1,2) .or.
     .              screen(ifr,2) .lt. screen(ifr+1,2) .or.
     .              screen(ifr,2) .lt. screen(ifr+2,2) ) then
c
                  nrfcell(ifr)  = 0
c
                  screen(ifr,2) = 0.0
c
               endif
c
               if ( nrfcell(ifr) .gt. 0 ) then
c
                  total_cells2 = total_cells2 + 1
c
               endif
c
            enddo
c
         endif
c
         do ifr = 1 , nfr
c
            screen(ifr,1) = screen(ifr,2)
c
         enddo
c
         call sort ( screen , nfr , nrcell )
c
         do ifr = 1 , nfr - 2 * rt_maxtgt
c
            nrfcell(nrcell(ifr))  = 0
c
            screen(nrcell(ifr),2) = 0.0
c
         enddo
c
c   Now the array nrfcell contains a positive number indicating the
c   doppler cell with the best score for that range line.  If the number
c   is zero then no doppler cell passed the first threshold or the cell
c   was thinned by the second phase.  No more than about 5 percent of
c   the range cells should have positive values of nrfcell.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                            STAGE 2-B
c
c              Estimate Accelerations for Future Frames
c
         if ( quiet .gt. 1 )
     .   write ( 6 , * ) total_cells1 , total_cells2 ,
     .                   total_cells3
         write ( 7 , * ) total_cells1 , total_cells2 ,
     .                   total_cells3
c
         if ( total_cells2 .gt. 0 ) then
c
            nacc = nafill_a - 2 * t_offset
c
            do ifr = 3 , nfr - 2
c
               if ( nrfcell(ifr) .gt. 0 ) then   !  Selected range cells
c
                  imp = nrfcell(ifr)
c
                  call accel_rt ( cac(1+t_offset,imp,ifr) , nacc ,
     .                            acc , dtv )
c
c                 Limit the acceleration to two master particles during
c                 the integration time
c
                  if ( abs( acc ) .lt.
     .                 2.0 * dfmp / ( dtv * float( actual ) ) ) then
c
                     rtgt  = drf * float( ifr - 1 - nfr / 2 )
c
                     ftgt  = dfmp * float( imp - 1 - nmp / 2 )
c
                     atgt  = acc +
     .                       dotdot(imp+1,1+(ifr-1)/mrrat,1+nspfv_a/2)
c
                     d     = sqrt( amax1( 0.0 , 1.0 - screen(ifr,2) ) )
c
                     dwtgt = sqrt( 2.0 * d ) / ( twopi * dtv )
c
c                    Use SNR of 1.0 (0 DB) since it is not used
c
                     call addtgt ( tframe , rtgt , ftgt , atgt , dwtgt ,
     .                             1.0 , 'r' )
c
                     total_cells3 = total_cells3 + 1
c
                  else
c
                     nrfcell(ifr)  = 0
c
                     screen(ifr,2) = 0.0
c
                  endif
c
               endif
c
            enddo
c
         endif
c
         do ifr = 1 , nfr
c
            screen(ifr,1) = screen(ifr,2)
c
         enddo
c
         call sort ( screen , nfr , nrcell )
c
         do ifr = 1 , nfr - rt_maxtgt
c
            nrfcell(nrcell(ifr))  = 0
c
            screen(nrcell(ifr),2) = 0.0
c
         enddo
c
         if ( quiet .gt. 1 )
     .   write ( 6 , * ) total_cells1 , total_cells2 ,
     .                   total_cells3
         write ( 7 , * ) total_cells1 , total_cells2 ,
     .                   total_cells3
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                            STAGE 2-C
c
c                       Estimate Phase Error
c
c   Loop over range and use the selected Doppler cell, summing the
c   numerator and denominator
c 
c-----------------------------------------------------------------------
c
c   The real work array pg(nafill_a,4) is used as follows:
c
c      pg(i,1)  :  Phase error delta for this iteration
c      pg(i,2)  :  Numerator for phase gradient calculation
c      pg(i,3)  :  Phase gradient for this iteration
c      pg(i,4)  :  Phase error total for all iterations
c
c   Initialize the total phase error to zero
c
         do ts = 1 + t_offset , nafill_a - t_offset
c
            pg(ts,4) = 0.0   !  Total phase error over all iterations

         enddo
c
c-----------------------------------------------------------------------
c
         if ( total_cells3 .gt. 0 ) then
c
            do iter = 1 , rt_maxit     !  Iterations
c
c-----------------------------------------------------------------------
c
c   The phase gradient is calculated as the ratio of two numbers summed
c   over all relevant cells
c
               do ts = 1 + t_offset , nafill_a - t_offset
c
                  cv(ts)   = 0.0       !  Sum for covariance (numerator)
c
                  pg(ts,2) = 0.0       !  Sum for power (denominator)
c
               enddo
c
c-----------------------------------------------------------------------
c
               do ifr = 3 , nfr - 2
c
                  if ( nrfcell(ifr) .gt. 0 ) then ! Selected range cells
c
c-----------------------------------------------------------------------
c
c   Load temporary array on first iteration
c
                     if ( iter .eq. 1 ) then
c
                        imp = nrfcell(ifr)
c
                        do ts = 2 + t_offset , nafill_a - t_offset - 1
c
c   Load a 'time-roughened' version of the master particle time series.
c   This is done to counteract the fact that the higher frequencies are
c   naturally underestimated because of the filter response due to the
c   sub-image FFT.
c
                           ctmp(ts,ifr) = ( 1.0 + 2.0 * rt_rough ) *
     .                                      cac(ts,imp,ifr)
     .                                    - rt_rough * cac(ts-1,imp,ifr)
     .                                    - rt_rough * cac(ts+1,imp,ifr)
c
                        enddo
c
                     else
c
c   For other than the first iteration correct the phase for the work
c   array
c
                        do ts = 1 + t_offset , nafill_a - t_offset
c
                           ctmp(ts,ifr) = ctmp(ts,ifr) * ct(ts,2)
c
                        enddo
c
                     endif
c
c-----------------------------------------------------------------------
c
c   For this range cell, calculate the contribution to the numerator and
c   the denominator of the phase gradient
c
c   First, calculate the covariance increment and sum it and calculate
c   the contribution to the denominator
c
                     covt = cmplx( 0.0 , 0.0 )
c
                     do ts = 1 + t_offset , nafill_a - t_offset - 1
c
                        ct(ts,1) =        ctmp(ts,ifr) *
     .                             conjg( ctmp(ts+1,ifr) )
c
                        covt     = covt + ct(ts,1)
c
c   Center power estimate on the ts and ts+1 points since this is where
c   the phase gradient estimate applies
c
                        pg(ts,2) = pg(ts,2) + 0.5 * (
     .                          ctmp(ts,ifr) * conjg( ctmp(ts,ifr) )
     .                        + ctmp(ts+1,ifr) * conjg( ctmp(ts+1,ifr) )
     .                                               )
c
                     enddo
c
c   Second, correct the covariance for the mean doppler and add it to
c   the numerator
c
                     if ( covt .ne. 0.0 ) then
c
c   Add the covariance increment for this range cell, corrected for the
c   mean doppler
c
                        covt = conjg( covt ) / cabs( covt )
c
                        do ts = 1 + t_offset , nafill_a - t_offset - 1
c
                           cv(ts) = cv(ts) + ct(ts,1) * covt
c 
                        enddo
c
                     endif       !  covt .ne. 0.0
c  
                  endif          !  Selected range cells
c
               enddo             !  Loop over range cells
c
c   The numerator and denominator versus time are now done
c
c-----------------------------------------------------------------------
c
c   Compute the phase gradient from the ratio
c
               do ts = 1 + t_offset , nafill_a - t_offset - 1
c
                  if ( pg(ts,2) .gt. 0.0 ) then
c
                     pg(ts,3) = - aimag( cv(ts) ) / pg(ts,2)
c
                  else
c
                     pg(ts,3) = 0.0
c
                  endif
c
               enddo
c
c-----------------------------------------------------------------------
c
c   Integrate the phase gradient with time to get the phase error
c
               pg(1+t_offset,1) = 0.0   !  Arbitrary initial phase
c
               do ts = 2 + t_offset , nafill_a - t_offset
c
                  pg(ts,1) = pg(ts-1,1) + 0.9 * pg(ts-1,3)
c
               enddo
c
c-----------------------------------------------------------------------
c
c   Remove phase error exp( -i * ph(t) ) from image by complex
c   multiplication in the time domain.
c
c-----------------------------------------------------------------------
c
c   Compute phase correction for this iteration - it is removed at the
c   top of the iteration loop
c
               pgrms = 0.0
c
               do ts = 1 + t_offset , nafill_a - t_offset
c
                  pg(ts,4) = pg(ts,4) + pg(ts,1)
c
                  pgrms    = pgrms + pg(ts,1) ** 2
c
                  ct(ts,2) = cmplx( cos( pg(ts,1) ) ,
     .                            - sin( pg(ts,1) ) )
c
               enddo
c
               pgrms = degrad *
     .                 sqrt( pgrms / float( nafill_a - 2 * t_offset ) )
c
               if ( quiet .gt. 1 )
     .         write ( 6 , * ) iter , pgrms
               write ( 7 , * ) iter , pgrms
c
c-----------------------------------------------------------------------
c
            enddo             !  Loop over iterations
c
         endif                !  total_cells3 > 0
c
c-----------------------------------------------------------------------
c
c   Compute the phase correction function based on all iterations.
c   Also, apply the weight function.
c
         write ( 7 , * )
         write ( 7 , * ) ' Phase error estimate from real-time PGA'
c
         do ts = 1 + t_offset , nafill_a - t_offset
c
            ct(ts,2) = cmplx( cos( pg(ts,4) ) , - sin( pg(ts,4) ) )
c
            write ( 7 , * ) ts , pg(ts,4)
c
         enddo
c
      endif                   !  Autofocus flag set
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                             STAGE 3
c
c   End of loops for autofocus - now FFT the data and add it into the
c   complex image
c
      do ifr = 1 , nfr
c
         jfr = ifr + nfr       !  Offset to even frame
c
c   Loop over master particles
c
         do imp = 1 , nmp
c
c   Put the odd and even frames into the FFT buffer
c
c   First, zero both frames
c
            do ts = 1 , 2 * nafill
c
               crc(ts) = cmplx( 0.0 , 0.0 )
c
            enddo
c
            if ( rt_pga .ne. 0 ) then
c
c   If the PGA is used, then multiply by the phase error function and
c   load the odd and even arrays with the first and last of the master
c   particle time series.
c
               do ts = 1 + t_offset , nafill - t_offset
c
c                             Odd frame
c
                  crc(ts)        = wtaf(ts-t_offset) *
     .                             cac(ts,imp,ifr) *
     .                             ct(ts,2)
c
c                             Even frame
c
                  crc(ts+nafill) = wtaf(ts-t_offset) *
     .                             cac(ts+augment,imp,ifr) *
     .                             ct(ts+augment,2)
c
               enddo
c
            else
c
c   If no PGA, then just load the odd and even arrays with the first and
c   last of the master particle time series.
c
               do ts = 1 + t_offset , nafill - t_offset
c
c                             Odd frame
c
                  crc(ts)        = wtaf(ts-t_offset) *
     .                             cac(ts,imp,ifr)
c
c                             Even frame
c
                  crc(ts+nafill) = wtaf(ts-t_offset) *
     .                             cac(ts+augment,imp,ifr)
c
               enddo
c
            endif
c
            ifirst = ( imp - 1 ) * mprat
c
c   FFT for the odd frame
c
            call fourt ( rc(1,1) , nn , 1 , + 1 , 1 ,
     .                   work(2*mprat+2) , nwork - 2 * mprat - 1 )
c
c   FFT for the even frame
c
            call fourt ( rc(1,1+nafill) , nn , 1 , + 1 , 1 ,
     .                   work(2*mprat+2) , nwork - 2 * mprat - 1 )
c
c   Use the linear interpolation weights to write these Fourier
c   coefficients into the fine resolution image array
c
            do ifa = 2 , 2 * mprat
c
c                                  Odd frame
c
               cframe(ifa+ifirst,ifr) = cframe(ifa+ifirst,ifr)
     .                    + work(ifa) * crc(ifa+na_off)
c
c                                  Even frame
c
               cframe(ifa+ifirst,jfr) = cframe(ifa+ifirst,jfr)
     .                    + work(ifa) * crc(ifa+na_off+nafill)
c
            enddo
c
         enddo           !  End of loop over master particles
c
      enddo              !  End of loop over image range cells
c
      return
      end
