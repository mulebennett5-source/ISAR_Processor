C**
C***********************************************************************
C**
      subroutine imgenu ( cframe , nfa , nfr , dff , drf , csbimg ,
     .                    nakeep , nsr , dfc , drs , fc0 , actual ,
     .                    nabuff , nafill , iisub , dts , ac , cac ,
     .                    mprat , wtaf , work , nwork , br , lambda ,
     .                    clight , dotdot , mrrat , rgwalk , vadd ,
     .                    nspfv , dtv , ti0 , tframe , iwork , strtch ,
     .                    vplat , slant0 )
C**
C***********************************************************************
C**
c
c      Image Formation Processor for Ultra-Wide Band Strip-map SAR
c
c
c   This routine supports the ISAR movie processor by forming a complex
c   image frame using the Triple-FFT algorithm.  Specifically, the
c   routine implements the third FFT in that method.
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
c   labelled 'STRETCH CORRECTION'.
c
c-----------------------------------------------------------------------
c
c   MODIFICATION FOR UWB STRIP-MAP PROCESSING ( November 22, 1995 )
c
c   In the Ultra Wide Band strip map option, the sub-image array was
c   computed in 'broadside range' instead of in actual range.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c***********************************************************************
c
c        Array sizes and raw data length passed as arguments
c
      integer     nfa , nfr , nakeep , nsr , actual , nabuff ,
     .            nafill , nwork , nspfv
c
c        nsr     :  Range cells in a sub-image
c        nakeep  :  Azimuth cells in a sub-image kept
c        actual  :  Sub-images per frame for this image
c        nabuff  :  Number of sub-images
c        nafill  :  Sub-images per frame including zero-fill
c        nfr     :  Range cells in a frame
c        nfa     :  Azimuth cells in a frame 
c        nwork   :  Real numbers available in work array
c        nspfv   :  No. of virtual sub-images per frame
c
c***********************************************************************
c
      integer     ifa , imp , mprat , fminus , fplus , rminus , rplus ,
     .            ish , iisub , mrrat , imr , nfrbig , ifrbig , nn(2) ,
     .            ipos , ineg , ictr , ifr , nmp , isbv , iwork(nwork) ,
     .            tminus , tplus , sfirst , slast , rplusp , rminusp ,
     .            strtch
c
c***********************************************************************
c
      real        dff , drf , dfc , drs , rsmin , rsmax , rfmin , eps ,
     .            rfmax , vmp , fmp , dfmp , fdelta , delrmp , ltime ,
     .            wtaf(nspfv) , fc0 , dts , fmp0 , rbig , br , lambda ,
     .            clight , fmpnom , vf , fnew , work(nwork) , fprime ,
     .            rprime , drfbig , dd , pprime , rfm , rfp , rdelta ,
     .            rdeltp , dotdot(1+nfa/mprat,nfr/mrrat,nspfv) , fold ,
     .            rgwalk , vadd , dtv , ti0 , tframe , wminus , wplus ,
     .            twopi , tsfirst , tslast , tvfirst , tv , tdelta ,
     .            ltimem , ltimep , vplat , slant0 , fstrip , rslant ,
     .            xdist , rsbig , rstrip , rleft , rright , delrmpl ,
     .            delrmpr , eps1 , foo , fgridl , fgridr , dtr , tcent ,
     .            dtdf , rbb , fcmax
c
      complex     cframe(nfa,nfr) , csbimg(nakeep,nsr,nabuff) , cdfmtm ,
     .            cdfmtp , cdfptm , cdfptp , cfptm , cfmtm , cfptp ,
     .            cfmtp , cdeltm , cdeltp
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
      logical     toofar , first
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
      fmp0   = - 0.5 * float( nfa ) * dff
      dfmp   = dff * float( mprat )
      nmp    = 1 + ( nfa / mprat )
c
      fcmax  = fc0 + float( nakeep - 1 ) * dfc
c
c   Min and max range of sub-images and fine resolution images
c
      rsmin  = - float( nsr / 2 ) * drs
      rsmax  = rsmin + float( nsr - 1 ) * drs
      rfmin  = - float( nfr / 2 ) * drf
      rfmax  = rfmin + float( nfr - 1 ) * drf
c
      nn(1)  = nafill            !  FFT length = sub-images per frame
c
c   Linear interpolation weights used for blending frequencies at the
c   master particles into the image frame
c
      do ifa = 1 , mprat + 1
c
         work(ifa) = float( mprat + 1 - ifa ) / float( mprat )
         work(nafill+2-ifa) = work(ifa)
c
      enddo
c
c   Calculate the times for the sub-images
c
c   iisub is the pointer to the last sub-image added to the buffer 
c
      slast   = iisub
      sfirst  = slast - ( actual - 1 )
      if ( sfirst .lt. 1 ) sfirst = sfirst + nabuff
c
      tsfirst = ti0
      tslast  = ti0 + float( actual - 1 ) * dts
c
      first   = .true.
      tcent   = tframe
c
      dtdf    = - lambda * slant0 / ( 2.0 * vplat ** 2 )
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
         rbig   = drfbig * float( ifrbig - 1 - nfrbig / 2 )
         dtr    = dtv
c
         rsbig  = slant0 + rbig + 0.5 * drfbig
c
         fstrip = 0.0
         rstrip = 0.0
c
         dtr    = ( rsbig / slant0 ) * dtv
c
c-----------------------------------------------------------------------
c
         vf     = 0.5 * lambda
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
            rbb = sqrt( rsbig ** 2 +
     .                  ( 0.5 * float( nspfv ) * dtr * vplat ) ** 2 )
            rbb = 0.5 * ( rsbig + rbb ) - slant0
            eps = ( 2.0 * br * 1.0E+12 * lambda * rbb ) / clight ** 2
            vf  = vf / ( 1.0 - eps )
c
         endif
c
c-----------------------------------------------------------------------
c
c   Loop over master particles
c
         do 200 imp = 1 , nmp
c
            fmpnom  = fmp0 + dfmp * float( imp - 1 )
            fstrip  = fmpnom
c
c-----------------------------------------------------------------------
c
c   Compute a time interpolation table to produce the virtual subimage
c   from the sub-image array
c
            tcent   = tframe + fmpnom * dtdf
c
            tvfirst = tcent - 0.5 * float ( nspfv - 1 ) * dtr
c
            do isbv = 1 , nspfv
c
               tv = tvfirst + float( isbv - 1 ) * dtr
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
                  tminus = sfirst + ifix( tdelta )
                  wminus = 1.0 - ( tdelta - float( ifix( tdelta ) ) )
                  if ( tminus .lt. 1 ) tminus = tminus + nabuff
                  if ( tminus .gt. nabuff ) tminus = tminus - nabuff
c
                  tplus  = tminus + 1
                  if ( tplus .lt. 1 ) tplus = tplus + nabuff
                  if ( tplus .gt. nabuff ) tplus = tplus - nabuff
c
               endif
c
               iwork(isbv)         = tminus
               iwork(nabuff+isbv)  = tplus
               work(nafill+1+isbv) = wminus
c
            enddo
c
c   Back-track to find the deviations of Doppler frequency and range
c   for this master particle at the start of the coherent integration
c   time from the acceleration field and the values of nominal Doppler
c   and range (the values at the mid-point of the integration time)
c
c-----------------------------------------------------------------------
c
            fprime = 0.0
            rprime = 0.0
c
            do isbv = 1 , nspfv / 2
c
               dd      = 0.5 * ( dotdot(imp,ifrbig,1+nspfv/2-isbv) +
     .                           dotdot(imp,ifrbig,2+nspfv/2-isbv) )
c
               fnew   = fprime
               fprime = fprime - dd * dtr
               rprime = rprime - 0.5 * ( fprime + fnew ) * vf * dtr
c
            enddo
c
            xdist  = - vplat * dtr * ( float( nspfv / 2 ) - 0.5 )
            rslant = sqrt( rsbig ** 2 + xdist ** 2 )
            fstrip = vplat * xdist / ( rslant * vf )
            rstrip = rslant - slant0
c
c-----------------------------------------------------------------------
c
c   Form time series of sub-image values for master particles
c
            do imr = 1 , mrrat
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
               wminus = work(nafill+1+isbv)
               wplus  = 1.0 - wminus
c
c   Calculate the actual frequency for this time for the master particle
c
               fmp    = fprime + fstrip
c
c   Calculate local time for use in phase and range walk calculations
c
               ltime  = dtr * float( isbv - 1 - nspfv / 2 )
c
c   Calculate the phase correction for path curvature
c
               pprime = ( ( rprime + rstrip ) / vf ) - fmp * ltime
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
               ish    = isbv - nspfv / 2
               if ( ish .lt. 1 ) ish = ish + nafill
c
c   Calculate interpolation parameters for frequency
c
               if ( fmp .gt. fc0 .and. fmp .lt. fcmax ) then
c
                  toofar = .false.
                  fminus = 1 + ifix( ( fmp - fc0 ) / dfc )
                  fplus  = min( fminus + 1 , nakeep )
                  fdelta = ( ( fmp - fc0 ) / dfc ) - float( fminus - 1 )
c
               else
c
                  toofar = .true.
c
               endif
c
c   Calculate the phase correction to the plus and minus frequency
c   bins
c
c              cdelta  = cexp( cmplx( 0.0 ,
c    .                         twopi * ( ltime * fmp + pprime ) ) )
c
               ltimem  = ltime - wplus * dts
               ltimep  = ltime + wminus * dts
c
c   Calculate the range offset for each frequency bin separately
c
               eps1    = ( vf / vplat ) ** 2
               foo     = 1.0 - eps1 * fmp ** 2
               fgridl  = fmp - fdelta * dfc
               fgridr  = fgridl + dfc
c
               rleft   = rsbig * ( sqrt( ( 1.0 - eps1 * fgridl ** 2 )
     .                                    / foo ) - 1.0 )
               rright  = rsbig * ( sqrt( ( 1.0 - eps1 * fgridr ** 2 )
     .                                    / foo ) - 1.0 )
c
c-----------------------------------------------------------------------
c
               cdeltm  = cpfast( ltimem * fmp + pprime )
               cdeltp  = cpfast( ltimep * fmp + pprime )
c
               cdfmtm  = cdeltm * ( 1.0 - fdelta )
               cdfptm  = cdeltm * fdelta
               cdfmtp  = cdeltp * ( 1.0 - fdelta )
               cdfptp  = cdeltp * fdelta
c
c   Increment the range to each master particle
c
c   The parameters rgwalk and vadd are normally equal to 1.0 and 0.0;
c   however, they can be used to disable the nominal value of range walk
c   or to add a constant to it.  This is only for research purposes and
c   for non-standard applications such as vibration echo studies.
c
               vmp     = rgwalk * vf * fmpnom + vadd
               delrmp  = vmp * ltime + rprime
c
               delrmpl = delrmp
               delrmpr = delrmp
c
               delrmpl = delrmp + rleft
               delrmpr = delrmp + rright
c
c-----------------------------------------------------------------------
c
c   Apply these accelerations and frequency interpolation parameters
c   to 'mrrat' range cells at a time
c
               do 402 imr = 1 , mrrat
c
c   Index of range cell in fine resolution image
c
                  ifr = mrrat * ( ifrbig - 1 ) + imr
c
c   Calculate range for this master particle at this time
c
                  rfm = delrmpl + drf * float( ifr - 1 - nfr / 2 )
                  rfp = delrmpr + drf * float( ifr - 1 - nfr / 2 )
c
c   Fill with zeros if local range following this particle lies
c   outside the sub-images
c
                  if ( rfm .le. rsmin .or. rfm .ge. rsmax .or.
     .                 rfp .le. rsmin .or. rfp .ge. rsmax .or.
     .                 ( toofar ) ) then
c
                     cac(ish,imr) = cmplx( 0.0 , 0.0 )
c
                     if ( rfm .le. rsmin .or. rfp .lt. rsmin ) then
                        nnear = nnear + 1
                     else
                        nfar  = nfar  + 1
                     endif
c
                  else
c
c   Calculate linear interpolation parameters for range
c
                     rminus  = 1 + ifix( ( rfm - rsmin ) / drs )
                     rplus   = rminus + 1
                     rdelta  = ( ( rfm - rsmin ) / drs )
     .                         - float( rminus - 1 )
c
                     rminusp = 1 + ifix( ( rfp - rsmin ) / drs )
                     rplusp  = rminusp + 1
                     rdeltp  = ( ( rfp - rsmin ) / drs )
     .                         - float( rminusp - 1 )
c
c   csbfp and csbfm are the range interpolated sub-image values at
c   the plus and minus frequency points
c
                     cfptm   = rdeltp * csbimg(fplus,rplusp,tminus) +
     .               ( 1.0 - rdeltp ) * csbimg(fplus,rminusp,tminus)
                     cfmtm   = rdelta * csbimg(fminus,rplus,tminus) +
     .               ( 1.0 - rdelta ) * csbimg(fminus,rminus,tminus)
c
                     cfptp   = rdeltp * csbimg(fplus,rplusp,tplus) +
     .               ( 1.0 - rdeltp ) * csbimg(fplus,rminusp,tplus)
                     cfmtp   = rdelta * csbimg(fminus,rplus,tplus) +
     .               ( 1.0 - rdelta ) * csbimg(fminus,rminus,tplus)
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
  402          continue
c
c   Update the deviation of the Doppler and range from their nominal
c   values (except for last time)
c
               if ( isbv .ne. nspfv ) then
c
                  xdist  = xdist + vplat * dtr
                  rslant = sqrt( rsbig ** 2 + xdist ** 2 )
                  fstrip = vplat * xdist / ( rslant * vf )
                  rstrip = rslant - slant0
c
                  dd     = 0.5 * ( dotdot(imp,ifrbig,isbv) +
     .                             dotdot(imp,ifrbig,isbv+1) )
                  fold   = fprime
                  fprime = fprime + dd * dtr
                  rprime = rprime + 0.5 * ( fold + fprime ) * vf * dtr
c
               endif
c
  300       continue  !  End of loop over time for a master particle
c
c   End of loop for loading mrrat FFT vectors - now FFT them and add
c   them into the complex image
c
            do 302 imr = 1 , mrrat
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
                  call fourt ( ac(1,1,imr) , nn , 1 , + 1 , 1 ,
     .                         work(nafill+2) , nwork - nafill - 1 )
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
  302       continue
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
         write ( 6 , '(/,i12,a48)' )
     .      ntotal , ' Points outside of the sub-image range were used'
         write ( 6 , '(a32,4i8)' ) '     nnear, nfar, nslow, nfast =' ,
     .                                   nnear , nfar , nslow , nfast
         write ( 7 , '(/,i12,a48)' )
     .      ntotal , ' Points outside of the sub-image range were used'
         write ( 7 , '(a32,4i8)' ) '     nnear, nfar, nslow, nfast =' ,
     .                                   nnear , nfar , nslow , nfast
c
      endif
c
      return
      end
