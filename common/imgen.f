C**
C***********************************************************************
C**
      subroutine imgen ( cframe , nfa , nfr , dff , drf , csbimg ,
     .                   nakeep , nsa , nsr , dfc , drc , fc0 ,
     .                   actual , nabuff , nafill , iisub , dts , ac ,
     .                   cac , mprat , wtaf , work , nwork , br ,
     .                   lambda , clight , dotdot , mrrat , rgwalk ,
     .                   vadd , nspfv , dtv , ti0 , tframe , iwork ,
     .                   strtch , rcen , fcen , reduce )
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
c   MODIFICATION FOR OPTIONAL STRETCH CORRECTIONS ( November 13, 1995 )
c
c   In previous versions it was assumed that the data came from a
c   stretch processing of a linear-swept-frequency range pulse.  This
c   causes two effects which are incorporated into the image formation
c   stage.  First, the range coordinate is 'apparent' range, a value
c   which has a doppler shift error in the range.  Second, the relation
c   between doppler frequency and radial velocity has a correction.
c
c   In this version of the program, these corrections are made optional
c   by adding the variable 'strtch'.  This is done in two code sections
c   labelled 'STRETCH CORRECTION'.  This code is only executed if
c   strtch = 1.
c
c-----------------------------------------------------------------------
c 
      implicit none
c
c***********************************************************************
c
c        Array sizes and raw data length passed as arguments
c
      integer     nfa , nfr , nakeep , nsa , nsr , actual , nabuff ,
     .            nafill , nwork , nspfv
c
c        nsr     :  Range cells in a sub-image
c        nakeep  :  Azimuth cells in a sub-image kept
c        nsa     :  Original no. of azimuth cells
c        actual  :  Sub-images per frame for this image
c        nabuff  :  Size of the sub-images buffer in time
c        nafill  :  Sub-images per frame including zero-fill
c        nfr     :  Range cells in a frame
c        nfa     :  Azimuth cells in a frame 
c        nwork   :  Real numbers available in work array
c        nspfv   :  No. of virtual sub-images per frame
c
c***********************************************************************
c
      integer     ifa , imp , mprat , fminus , fplus , rminus , rplus ,
     .            ish , iisub , mrrat , imr , nfrbig , ifrbig , nn(1) ,
     .            ipos , ineg , ictr , ifr , nmp , fshift , slast ,
     .            iwork(nwork) , isbv , tminus , tplus , sfirst , jmp ,
     .            strtch , rminus2 , rplus2 , rminus3 , rplus3 , ii
c
c***********************************************************************
c
      real        dff , drf , dfc , drc , rsmin , rsmax , rfmin , eps ,
     .            rfmax , vmp , rmp , fmp , dfmp , fdelta , delrmp ,
     .            ltime , rdelta , wtaf(nspfv) , fc0 , dts , fmp0 ,
     .            rbig , br , lambda , clight , fmpnom , vf , fnew ,
     .            work(nwork) , fprime , rprime , drfbig , dd , pprime ,
     .            dotdot(1+nfa/mprat,nfr/mrrat,nspfv) , fold , twopi ,
     .            rgwalk , vadd , dtv , ti0 , tframe , wminus , wplus ,
     .            tsfirst , tslast , tvfirst , tv , tdelta , ltimem ,
     .            ltimep , rcen , fcen , rdelta2 , rdelta3 , r25_128
c
      complex     cframe(nfa,nfr) , csbimg(nakeep,nsr,nabuff) , cdfmtm ,
     .            cdfmtp , cdfptm , cdfptp , cfptm , cfmtm , cfptp ,
     .            cfmtp , cdeltm , cdeltp
c
      complex     c1lag(nspfv)
c
      real        d , dc , dw , pwr , phi
c
      integer     reduce , nsmooth
c
c   cpfast is a table look-up version of the quantity cexp(j*2*pi*arg)
c
      complex     cpfast
c
c   ac and cac are real and complex names for an equivalenced work
c   vector used in the call to FOURT
c
      real        ac(2,nafill,mrrat)
      complex     cac(nafill,mrrat)
c
c   Variables for diagnostic information on the number of points
c   which fell outside the limits of the sub-images
c
      integer     nslow , nfast , nnear , nfar , ntotal
c
c***********************************************************************
c
      twopi = 2.0 * atan2( 0.0 , - 1.0 )
c
c   Set counters for out-of-range particles so that we know whether
c   enough sub-image points have been saved
c
      nfast = 0
      nslow = 0
      nfar  = 0
      nnear = 0
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
      fmp0   = fcen - 0.5 * float( nfa ) * dff
      dfmp   = dff * float( mprat )
c
      nmp    = 1 + ( nfa / mprat )
      if ( mod( nfa , mprat ) .ne. 0 ) nmp = nmp + 1
c
c   Min and max range of sub-images and fine resolution images
c
      rsmin  =       - float( nsr / 2 ) * drc
      rsmax  = rsmin + float( nsr - 1 ) * drc
c
      rfmin  =       - float( nfr / 2 ) * drf
      rfmax  = rfmin + float( nfr - 1 ) * drf
c
      nn(1)  = nafill            !  FFT length = sub-images per frame
c
c   Linear interpolation weights used for blending frequencies at the
c   master particles into the image frame
c
      do ifa = 1 , mprat + 1
c
         work(ifa)          = float( mprat + 1 - ifa ) / float( mprat )
c
         work(nafill+2-ifa) = work(ifa)
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Compute an interpolation table to produce the virtual subimage from
c   the sub-image array
c
c     iisub is the pointer to the last sub-image added to the buffer 
c
      slast   = iisub
      sfirst  = slast - ( actual - 1 )
      if ( sfirst .lt. 1 ) sfirst = sfirst + nabuff
c
      tsfirst = ti0
      tslast  = ti0 + float( actual - 1 ) * dts
      tvfirst = tframe - 0.5 * float ( nspfv - 1 ) * dtv
c
      do isbv = 1 , nspfv
c
         tv = tvfirst + float( isbv - 1 ) * dtv
c
         if ( tv .le. tsfirst ) then
c
            tminus = sfirst
            tplus  = sfirst
            wminus = 1.0
c
         else if ( tv .ge. tslast ) then
c
            tminus = slast
            tplus  = slast
            wminus = 1.0
c
         else
c
            tdelta = ( tv - tsfirst ) / dts
c
            wminus = 1.0 - ( tdelta - float( ifix( tdelta ) ) )
c
            tminus = sfirst + ifix( tdelta )
            if ( tminus .lt. 1      ) tminus = tminus + nabuff
            if ( tminus .gt. nabuff ) tminus = tminus - nabuff
c
            tplus  = tminus + 1
            if ( tplus .lt. 1       ) tplus  = tplus + nabuff
            if ( tplus .gt. nabuff  ) tplus  = tplus - nabuff
c
         endif
c
         iwork(isbv)         = tminus
         iwork(nabuff+isbv)  = tplus
         work(nafill+1+isbv) = wminus
c
      enddo
c
c-----------------------------------------------------------------------
c
c                Loop over output frame range cells
c
c   Assume that nfr is an integer times mrrat
c
c   nfrbig is the number of range regions of mrrat range cells.
c   Over each of these blocks the accelerations and the frequency
c   interpolation parameters are assumed constant.
c
      nfrbig = nfr / mrrat
      drfbig = drf * float( mrrat )
c
      do 100 ifrbig = 1 , nfrbig
c
         rbig = drfbig * float( ifrbig - 1 - nfrbig / 2 )
c
         vf   = 0.5 * lambda  !  Nominal (non-stretch) value
c
c-----------------------------------------------------------------------
c
c                     STRETCH CORRECTION (1 of 2)
c
c   Over very large swath and for high chirp rate there is an 'extra'
c   term which can be accounted for by making the ratio of velocity
c   to Doppler frequency a slow function of range.  This is rarely
c   important but is easy to include.
c
         if ( strtch .eq. 1 ) then
c
            eps  = ( 2.0 * br * 1.0E+12 * lambda * rbig ) /
     .                       clight ** 2
            vf   = 0.5 * lambda / ( 1.0 - eps )
c
         endif
c
c-----------------------------------------------------------------------
c
c
c   Loop over master particles
c
         do 200 imp = 1 , nmp
c
            jmp    = min( imp , 1 + nfa / mprat )
c
            fmpnom = fmp0 + dfmp * float( imp - 1 )
c
c   Back-track to find the deviations of Doppler frequency and range
c   for this master particle at the start of the coherent integration
c   time from the acceleration field and the values of nominal Doppler
c   and range (the values at the mid-point of the integration time)
c
            fprime = 0.0
            rprime = 0.0
c
            do isbv = 1 , nspfv / 2
c
               dd     = 0.5 * ( dotdot(jmp,ifrbig,1+nspfv/2-isbv) +
     .                          dotdot(jmp,ifrbig,2+nspfv/2-isbv) )
               fnew   = fprime
               fprime = fprime - dd * dtv
               rprime = rprime - 0.5 * ( fprime + fnew ) * vf * dtv
c
            enddo
c
c   Form time series of sub-image values for master particles
c
            do imr  = 1 , mrrat
c
               do isbv = 1 , nafill
c
                  cac(isbv,imr) = cmplx( 0.0 , 0.0 )
c
               enddo
c
            enddo
c
            do 300 isbv = 1 , nspfv
c
c   Interpolation coefficients to compute virtual sub-image from the
c   sub-image array
c
               tminus = iwork(isbv)
               tplus  = iwork(nabuff+isbv)
c
               wminus = work(nafill+1+isbv)
               wplus  = 1.0 - wminus
c
c   Calculate the actual frequency for this time for the master particle
c
               fmp    = fmpnom + fprime
c
c   Calculate local time for use in phase and range walk calculations
c
               ltime  = dtv * float( isbv - 1 - nspfv / 2 )
c
c   Calculate the phase correction for path curvature
c
               pprime = ( rprime / vf ) - fprime * ltime
c
c-----------------------------------------------------------------------
c
c                     STRETCH CORRECTION (2 of 2)
c
               if ( strtch .eq. 1 ) pprime = pprime -
     .                ( clight / ( lambda * 1.0E+12 * br ) ) * fprime
c
c-----------------------------------------------------------------------
c
c   Determine the correct position in the FFT array to load the point
c   into
c
               ish = isbv - nspfv / 2
               if ( ish .lt. 1 ) ish = ish + nafill
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
                  nslow = nslow + 1
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
                  nfast = nfast + 1
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
               ltimem = ltime - wplus  * dts
               ltimep = ltime + wminus * dts
c
               cdeltm = cpfast( ltimem * fmp + pprime )
               cdeltp = cpfast( ltimep * fmp + pprime )
c
               cdfmtm = cdeltm * ( 1.0 - fdelta )
               cdfptm = cdeltm * fdelta
c
               cdfmtp = cdeltp * ( 1.0 - fdelta )
               cdfptp = cdeltp * fdelta
c
c   Increment the range to each master particle
c
c   The parameters rgwalk and vadd are normally equal to 1.0 and 0.0;
c   however, they can be used to disable the nominal value of range walk
c   or to add a constant to it.  This is only for research purposes and
c   for non-standard applications such as vibration echo studies.
c
               vmp    = rgwalk * vf * fmpnom + vadd
c
               delrmp = rcen + vmp * ltime + rprime
c
c   Apply these accelerations and frequency interpolation parameters
c   to 'mrrat' range cells at a time
c
               do imr = 1 , mrrat
c
c   Index of range cell in fine resolution image
c
                  ifr = mrrat * ( ifrbig - 1 ) + imr
c
c   Calculate range for this master particle at this time
c
                  rmp = delrmp + drf * float( ifr - 1 - nfr / 2 )
c
c   Fill with zeros if local range following this particle lies
c   outside the sub-images
c
                  if ( rmp .le. rsmin + 4.0 * drc
     .            .or. rmp .ge. rsmax - 4.0 * drc ) then
c
                     cac(ish,imr) = cmplx( 0.0 , 0.0 )
c
                     if ( rmp .le. rsmin ) then
                        nnear = nnear + 1
                     else
                        nfar  = nfar  + 1
                     endif
c
                  else
c
c   Calculate linear interpolation parameters for range
c
                     rminus = 1 + ifix( ( rmp - rsmin ) / drc )
                     rplus  = rminus + 1
                     rdelta = ( ( rmp - rsmin ) / drc )
     .                        - float( rminus - 1 )
c
c-----------------------------------------------------------------------
c
c   Over-ride with nearest neighbor interpolation to test real-time
c   approximations
c
c                    if ( rdelta .gt. 0.5 ) then
c
c                       rdelta = 1.0
c
c                    else
c
c                       rdelta = 0.0
c
c                    endif
c
c-----------------------------------------------------------------------
c
c   csbfp and csbfm are the range interpolated sub-image values at
c   the plus and minus frequency points
c
                     cfptm  = rdelta * csbimg(fplus,rplus,tminus) +
     .                  ( 1.0 - rdelta ) * csbimg(fplus,rminus,tminus)
                     cfmtm  = rdelta * csbimg(fminus,rplus,tminus) +
     .                  ( 1.0 - rdelta ) * csbimg(fminus,rminus,tminus)
c
                     cfptp  = rdelta * csbimg(fplus,rplus,tplus) +
     .                  ( 1.0 - rdelta ) * csbimg(fplus,rminus,tplus)
                     cfmtp  = rdelta * csbimg(fminus,rplus,tplus) +
     .                  ( 1.0 - rdelta ) * csbimg(fminus,rminus,tplus)
c
c-----------------------------------------------------------------------
c
c               Higher-order corrections to range interpolation
c
                     rminus2 = rminus - 1
                     rplus2  = rplus  + 1
c
                     rdelta2 = ( 1.0 + rdelta ) / 3.0
c
                     r25_128 = 25.0 / 128.0
c
                     cfptm   = ( 75.0 / 64.0 ) * cfptm - r25_128 * (
     .                       rdelta2   * csbimg(fplus,rplus2,tminus) +
     .               ( 1.0 - rdelta2 ) * csbimg(fplus,rminus2,tminus) )
c
                     cfmtm   = ( 75.0 / 64.0 ) * cfmtm - r25_128 * (
     .                       rdelta2   * csbimg(fminus,rplus2,tminus) +
     .               ( 1.0 - rdelta2 ) * csbimg(fminus,rminus2,tminus) )
c
                     cfptp   = ( 75.0 / 64.0 ) * cfptp - r25_128 * (
     .                       rdelta2   * csbimg(fplus,rplus2,tplus) +
     .               ( 1.0 - rdelta2 ) * csbimg(fplus,rminus2,tplus) )
c
                     cfmtp   = ( 75.0 / 64.0 ) * cfmtp - r25_128 * (
     .                       rdelta2   * csbimg(fminus,rplus2,tplus) +
     .               ( 1.0 - rdelta2 ) * csbimg(fminus,rminus2,tplus) )
c
                     rminus3 = rminus - 2
                     rplus3  = rplus  + 2
c
                     rdelta3 = ( 2.0 + rdelta ) / 5.0
c
                     cfptm   = cfptm + ( 3.0 / 128.0 ) * (
     .                       rdelta3   * csbimg(fplus,rplus3,tminus) +
     .               ( 1.0 - rdelta3 ) * csbimg(fplus,rminus3,tminus) )
c
                     cfmtm   = cfmtm + ( 3.0 / 128.0 ) * (
     .                       rdelta3   * csbimg(fminus,rplus3,tminus) +
     .               ( 1.0 - rdelta3 ) * csbimg(fminus,rminus3,tminus) )
c
                     cfptp   = cfptp + ( 3.0 / 128.0 ) * (
     .                       rdelta3   * csbimg(fplus,rplus3,tplus) +
     .               ( 1.0 - rdelta3 ) * csbimg(fplus,rminus3,tplus) )
c
                     cfmtp   = cfmtp + ( 3.0 / 128.0 ) * (
     .                       rdelta3   * csbimg(fminus,rplus3,tplus) +
     .               ( 1.0 - rdelta3 ) * csbimg(fminus,rminus3,tplus) )
c
c-----------------------------------------------------------------------
c
c   Finally, load the FFT array with the time-weighted and
c   frequency-interpolated value of the sub-image
c
                     cac(ish,imr) = wtaf(isbv) * (
     .                    wminus * ( cdfptm * cfptm + cdfmtm * cfmtm )
     .                  + wplus  * ( cdfptp * cfptp + cdfmtp * cfmtp ) )
c
                  endif
c
               enddo     !  Loop over mrrat range cells
c
c   Update the deviation of the Doppler and range from their nominal
c   values (except for last time)
c
               if ( isbv .ne. nspfv ) then
c
                  dd     = 0.5 * ( dotdot(jmp,ifrbig,isbv) +
     .                             dotdot(jmp,ifrbig,isbv+1) )
                  fold   = fprime
                  fprime = fprime + dd * dtv
                  rprime = rprime + 0.5 * ( fold + fprime ) * vf * dtv
c
               endif
c
  300       continue  !  End of loop over time for a master particle
c
c   End of loop for loading mrrat FFT vectors - now FFT them and add
c   them into the complex image
c
            do imr = 1 , mrrat
c
               ifr = mrrat * ( ifrbig - 1 ) + imr
c
c   ictr = Fine resolution index for this master particle
c
               ictr = 1 + ( imp - 1 ) * mprat
c
               if ( mprat .eq. 1 ) then
c
                  if ( ictr .le. nfa ) then
c
                     cframe(ictr,ifr) = cmplx( 0.0 , 0.0 )
c
                     do isbv = 1 , nafill
c
                        cframe(ictr,ifr) = cframe(ictr,ifr) +
     .                                     cac(isbv,imr)
c
                     enddo
c
                  endif
c
               else
c
c-----------------------------------------------------------------------
c
                  if ( mod(reduce,2) .eq. 1 ) then
c
c                         Implement filtering by Doppler width
c
c                    Load temporary array with data FFT_Shifted back to
c                    the original time order
c
                     do ii = 1 , nspfv / 2
c
                        c1lag(ii+nspfv/2) = cac(ii,imr)
                        c1lag(ii)         = cac(ii+(nafill-nspfv)/2,imr)
c
                     enddo
c
c                    Remove Taylor weight
c
                     do ii = 1 , nspfv
c
                        c1lag(ii) = c1lag(ii) / wtaf(ii)
c
                     enddo
c
c                    Smooth to remove Doppler beyond the local pair of
c                    master particles
c
                     nsmooth = nint( float( nafill / ( 2 * mprat ) ) )
                     
                     if ( nspfv .gt. 4 * nsmooth ) then
c
                        call c1smooth ( c1lag , nspfv , nsmooth )
c
c                       Estimate the Doppler time series parameters
c
                        call onelag ( c1lag , nspfv , dtv , nsmooth ,
     .                                dc , dw , d , pwr , phi )
c
c                       Downweight the complex data based on the local Doppler
c                       width (non-dimensional version)
c
                        do ii = 1 , nafill
                           cac(ii,imr) = cac(ii,imr) / sqrt( 0.01 + d )
                        enddo
c                    
                     endif
c
                  endif
c
c-----------------------------------------------------------------------
c
                  call fourt ( ac(1,1,imr) , nn , 1 , + 1 , 1 ,
     .                         work(2*nafill+2) , nwork - nafill - 1 )
c
c   Use the linear interpolation weights to write these Fourier
c   coefficients into the fine resolution image array
c
                  if ( imp .ne. 1 ) then
c
c   First, for all but the most negative frequency add the negative
c   Fourier component into the image frame
c
                     do ifa = 2 , mprat
c
                        ineg = ictr + 1 - ifa
                        if ( ineg .le. nfa )
     .                       cframe(ineg,ifr) = cframe(ineg,ifr)
     .                     + work(nafill+2-ifa) * cac(nafill+2-ifa,imr)
c
                     enddo
c
                  endif
c
c   Next, for all but the most positive frequency add the positive
c   Fourier component into the image frame
c
                  if ( imp .ne. nmp ) then
c
                     do ifa = 1 , mprat
c
                        ipos = ictr - 1 + ifa
                        if ( ipos .le. nfa )
     .                  cframe(ipos,ifr) = work(ifa) * cac(ifa,imr)
c
                     enddo
c
                  endif
c
               endif
c
            enddo  !  End of loop over mrrat range cells
c
  200    continue  !  End of loop over master particles
c
  100 continue     !  End of loop over groups of image range cells
c
c   If points outside of sub-image were used, write out diagnostic
c
      ntotal = nnear + nfar + nslow + nfast
c
      if ( ntotal .ne. 0 ) then
c
c        write ( 6 , '(/,i12,a48)' )
c    .      ntotal , ' Points outside of the sub-image range were used'
c         write ( 6 , '(a32,4i8)' ) '    nnear, nfar, nslow, nfast =' ,
c    .                                   nnear , nfar , nslow , nfast
c        write ( 7 , '(/,i12,a48)' )
c    .      ntotal , ' Points outside of the sub-image range were used'
c        write ( 7 , '(a32,4i8)' ) '     nnear, nfar, nslow, nfast =' ,
c    .                                   nnear , nfar , nslow , nfast
c
      endif
c
      return
      end
