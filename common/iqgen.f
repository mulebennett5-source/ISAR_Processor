C**
C***********************************************************************
C**
      subroutine iqgen ( crc , ntr , ipulse , spulse , nerror )
C**
C***********************************************************************
C**
c
c     This routine generates a synthetic pulse of uncompressed signal
c     history as a function of fast time according to a test pattern
c     of targets
c
c
c     INPUTS
c
c         ntr:     Number of un-compressed samples
c         ipulse : Present pulse number
c         spulse : Start pulse number
c
c
c     OUTPUTS
c
c         crc    : Uncompressed signal history
c         nerror : Error type, 0 means no error
c
c
c   This routine is called for simulated cases ( nfmt < 0 ).  Positive
c   values of nfmt are reserved for real data.
c
c   There are 6 values of nfmt (defined in 'sarprm.h'):
c
c          - 1  :  A simple linear pattern of point targets
c
c          - 2  :  A user-specified pattern of point targets which
c                  is allowed only some geometrically simple motion
c                  patterns, some of which may not be realizable, but
c                  which serve to test the motion model used in the
c                  processor
c
c          - 3  :  A 3-D rigid body (ship) which undergoes a variety
c                  of periodic and general motions
c
c          - 4  :  A strip-map scene using the exact expression for
c                  the phase and using a simple test pattern of targets
c
c          - 5  :  A strip-map scene using the quadratic expansion
c                  of the phase with time
c
c          - 6  :  A strip-map scene using the exact expression for
c                  the phase and using a general pattern of targets
c
c***********************************************************************
c
c Revision History:
c
c     Bart Schade   Sep. 12, 1995   Modified target file format for
c                                   multiple types of targets and added
c                                   periodic ship target.
c
c     Bart Schade   Sep. 19, 1995   Added sea clutter.
c                       
c     John Bennett  Oct. 31, 1995   Added stripmap modes (nfmt=-4,-5)
c
c     John Bennett  Jan. 31, 1997   Added mo-comp errors to strip-map
c
c     John Bennett  Jan. 30, 2008   Added step-chirp pulses
c
c***********************************************************************
c
      implicit none
c
c   Fundamental parameters held in common for use in multiple routines
c
      include     'sarprm.h'
c
      include     'updates.h'
c
      integer      ntr , i , ipulse , iseed , k , middle , spulse ,
     .             tgt_type(ntestm) , ndot , lastdot , nclut , nerror ,
     .             ntrband , iband , ii
c
      real         rnor , uni , r1 , r2 , tfast , tslow , ph , bc , pi ,
     .             twopi , amp , vn , epsvn , ptime , rmin , rmax ,
     .             vprime , slant , cross , tprime , toff , rsprime,
     .             beamwidth , tslow_sc , tslow_last , tfoff
c
c   ***** Start Ship Simulation Variables ******
c
      real         inyaw , inpitch , inroll , avelyaw , avelpitch ,
     .             avelroll , aacyaw , aacpitch , aacroll , pyawa ,
     .             pyawp , pyawt , ppitcha , ppitchp , ppitcht ,
     .             prolla , prollp , prollt
c
      real         x_scat(ntestm) , y_scat(ntestm) , z_scat(ntestm) ,
     .             newxsc(ntestm) , newysc(ntestm) , newzsc(ntestm) ,
     .             rcs(ntestm) , oldrange(ntestm) , svnois(ntestm) ,
     .             svtime(ntestm) , svfreq(ntestm) , svamp(ntestm) ,
     .             spfreq(ntestm) , spdur(ntestm) , spamp(ntestm) ,
     .             spstart(ntestm) , vpert(ntestm) , rpert(ntestm) ,
     .             svactl(ntestm) , acluts(ntestm) , dcluts(ntestm) ,
     .             refrange(ntestm)
c
c   ***** End Ship Simulation Variables ******
c
c   Strip Map cross-range offset (meters)
c
      real         r0(ntestm) , x0(ntestm)
c
      complex      crc(ntr) , cpfast
c
      character    iqfile*80
c
      save         !  Make sure that temporary variables last over
c                  !  multiple calls
c
c***********************************************************************
c
      middle    = 1 + ntest / 2
c
      pi        = atan2( 0.0 , - 1.0 )
c
      twopi     = 2.0 * pi
c
c-----------------------------------------------------------------------
c
c   2-degree beamwidth for illumination function in strip-map mode
c
      beamwidth = 2.0 * pi / 180.0
c
c-----------------------------------------------------------------------
c
c   Consistency checks - ntr divisible by nbands and ntr/nbands even
c
      if ( mod( ntr , nbands ) .ne. 0 ) stop 'ntr, nbands inconsistent'
c
      ntrband = ntr / nbands
c
      if ( mod( ntrband , 2 ) .ne. 0 ) stop 'ntr/nbands not even'
c
c-----------------------------------------------------------------------
c
c   On the first call, calculate the fundamental parameters of the
c   test targets and initialize the random number generator for
c   the Gaussian noise
c
      if ( ipulse .eq. spulse ) then
c
         iseed      = 1234567
c
         tslow_sc   = 0.0
c
         tslow_last = tslow_sc 
c
         r1         = rnor( iseed )
         vn         = 0.0
c
         rmin       = - clight / ( 4.0 * 1.0E+12 * abs( br ) * dtr )
         rmax       = - rmin
c
         if ( nfmt .eq. - 1 ) then
c
c-----------------------------------------------------------------------
c
c                    ISAR Linear Test Pattern
c
            do k = 1 , ntest
c
               r(k)        = rbar + float( k - middle ) * drtest
               v0(k)       = vbar + slipv - float( k - middle ) * dvtest
               v(k)        = v0(k)
               tgdvdt(k)   = dvdt + 0.5 * lambda * (
     .                       ( daccdf * v0(k) * 2.0 / lambda ) +
     .                       ( daccdr * r(k) ) )
c
               oldrange(k) = r(k)
               refrange(k) = 0.0
c
               if ( k .eq. middle ) then
c
                  tgtamp(k) = 4.0
c
               else
c
                  tgtamp(k) = 1.0
c
               endif
c
            enddo
c
         else if ( nfmt .eq. - 2 ) then
c
c-----------------------------------------------------------------------
c
c                    Read in Test Pattern from File
c
            ndot   = lastdot( ifile ) - 1
            iqfile = ifile(1:ndot) // '.iqg'
c
            open ( 9 , file = iqfile , form = 'FORMATTED' ,
     .                 status = 'OLD' , err = 2000 )
c
            read ( 9 , * , end = 1000 , err = 2000 ) ntest
            if ( ntest .gt. ntestm ) ntest = ntestm
c
            do k = 1 , ntest
c
               read ( 9 , * , end = 1000 , err = 2000 )
     .                r(k) , v0(k) , tgtamp(k)
c
               r(k)        = r(k) + rbar
               v(k)        = v0(k) + vbar + slipv
               tgdvdt(k)   = dvdt + 0.5 * lambda * (
     .                       ( daccdf * v0(k) * 2.0 / lambda ) +
     .                       ( daccdr * r(k) ) )
c
               oldrange(k) = r(k)
               refrange(k) = 0.0
c
            enddo
c
            close ( 9 )
c
         else if ( nfmt .eq. - 3 ) then
c
c-----------------------------------------------------------------------
c
c                    SHIP ROTATION SIMULATION
c
            ndot   = lastdot( ifile ) - 1
            iqfile = ifile(1:ndot) // '.sim'
c
c           GET USER PROVIDED RUN PARAMETERS FROM INPUT TEXT FILE
c
            call getuparms ( inyaw , inpitch , inroll , avelyaw ,
     .                       avelpitch , avelroll , aacyaw ,
     .                       aacpitch , aacroll , pyawa , pyawp ,
     .                       pyawt , ppitcha , ppitchp , ppitcht ,
     .                       prolla , prollp , prollt , iqfile , 9 ,
     .                       nerror )
c
            if ( nerror .ne. 0 ) return
c
            iqfile = ifile(1:ndot) // '.shp'
c
c           GET FIXED (X,Y,Z,RADAR XSECTION) COORDS OF SCATTERERS ON SHIP
c
            do k = 1 , ntestm
c
               svnois(k)  = 0.0
               svtime(k)  = 0.0
               rpert(k)   = 0.0
               vpert(k)   = 0.0
               svactl(k)  = 0.0
               svamp(k)   = 0.0
               svfreq(k)  = 0.0
               spdur(k)   = 0.0
               spamp(k)   = 0.0
               spstart(k) = 0.0
c
            enddo
c
            call getshpcoord ( tgt_type , x_scat , y_scat , z_scat ,
     .                         rcs , svnois , svtime , svfreq , svamp ,
     .                         spfreq , spdur , spamp , spstart ,
     .                         ntestm , ntest , iqfile , 9 , nerror )
c
            if ( nerror .ne. 0 ) return
c
c           ROTATE THE SHIP TO POSITION AT ONE NEGATIVE TIME STEP
c           SO THAT WE WILL BE ABLE TO GET A VELOCITY ESTIMATE (NEED
c           TWO TIME SAMPLES OF THE RANGE TO GET A VELOCITY)
c
            call rot_ship ( ntest , inyaw , inpitch , inroll ,
     .                      avelyaw , avelpitch , avelroll ,
     .                      aacyaw , aacpitch , aacroll , pyawa ,
     .                      pyawp , pyawt , ppitcha , ppitchp ,
     .                      ppitcht , prolla , prollp , prollt ,
     .                      x_scat , y_scat , z_scat ,
     .                      newxsc , newysc , newzsc , - dtp )
c
c           SAVE OLD SCATTERING COORDS FOR VELOCITY COMPUTATION
c           (ONLY X IS IMPORTANT BECAUSE WE ARE MAKING THE APPROX
c           THAT RANGE = X)
c
            do k = 1 , ntest
c
               rpert(k)   = 0.0
               vpert(k)   = 0.0
               svactl(k)  = 0.0
               svamp(k)   = 0.0
c
               newxsc(k)   = newxsc(k) + rbar
c
               oldrange(k) = newxsc(k)
c
            enddo
c
c           ROTATE THE SHIP TO POSITION AT TIME = 0.0
c
            call rot_ship ( ntest , inyaw , inpitch , inroll , avelyaw ,
     .                      avelpitch , avelroll , aacyaw , aacpitch ,
     .                      aacroll , pyawa , pyawp , pyawt , ppitcha ,
     .                      ppitchp , ppitcht , prolla , prollp ,
     .                      prollt , x_scat , y_scat , z_scat , newxsc ,
     .                      newysc , newzsc , 0.0 )
c
            do k = 1 , ntest
               write( 11 , '(3f8.4)' ) newxsc(k) , newysc(k) , newzsc(k)
            enddo
c
c           CALCULATE RANGE, LINE OF SIGHT VEL, RCS FOR EACH SCATTERER
c
            do k = 1 , ntest
c
               r(k)      = newxsc(k)
c
               v0(k)     = ( newxsc(k) - oldrange(k) ) / dtp
c
               tgtamp(k) = sqrt( abs( rcs(k) ) )
               v(k)      = v0(k)
               tgdvdt(k) = dvdt + 0.5 * lambda * (
     .                     ( daccdf * v0(k) * 2.0 / lambda ) +
     .                     ( daccdr * r(k) ) )
c
c              INITIALIZE PERTURBATION VELOCITY (VELOCITY NOT ASSOCIATED
c              WITH ROTATION) AND CUMULATIVE POSITION PERTURBATION DUE
c              TO VELOCITY PERTURBATION
c
               vpert(k)  = 0.0
               rpert(k)  = 0.0
               svactl(k) = 0.0
c
            enddo
c
c           SAVE OLD SCATTERING COORDS FOR VELOCITY COMPUTATION
c           NEXT TIME
c
            do k = 1 , ntest
c
               oldrange(k) = newxsc(k)
               refrange(k) = 0.0
c
            enddo
c
         else if ( nfmt .eq. - 4 .or. nfmt .eq. - 5 ) then
c
c-----------------------------------------------------------------------
c
c                    Strip-Map Test Patterns
c
            do k = 1 , ntest
c
               if ( k .eq. middle ) then
c
                  tgtamp0(k) = .4
c
               else
c
                  tgtamp0(k) = .1
c
               endif
c
            enddo
c
            rsprime = 0.0
c
         else if ( nfmt .eq. - 6 ) then
c
c-----------------------------------------------------------------------
c
c                    Read in Strip Map Test Pattern from File
c
            ndot   = lastdot( ifile ) - 1
            iqfile = ifile(1:ndot) // '.iqs'
c
            open ( 9 , file = iqfile , form = 'FORMATTED' ,
     .                 status = 'OLD' , err = 2000 )
c
            read ( 9 , * , end = 1000 , err = 2000 ) ntest
            if ( ntest .gt. ntestm ) ntest = ntestm
c
            do k = 1 , ntest
c
               read ( 9 , * , end = 1000 , err = 2000 )
     .                r0(k) , x0(k) , tgtamp0(k)
c
            enddo
c
            close ( 9 )
c
            rsprime = 0.0
c
         else
c
            do k = 1 , ntest
c
               r(k)        = rbar
               v0(k)       = vbar
               v(k)        = vbar
               tgtamp(k)   = 0.0
c
               oldrange(k) = rbar
               refrange(k) = rbar
c
            enddo
c
         endif
c             
         if ( pclut .ne. 0.0 ) then
c
            nclut = nint( ntr * pclut )
c
            do k = 1 , nclut
c
               cr(k)     = uni( 0 )
               cr(k)     = rmin + cr(k) * ( rmax - rmin )
               cv(k)     = vclut + rnor( 0 ) * vdclut
               dcluts(k) = uni( 0 )
               acluts(k) = sin( dcluts(k) * pi ) * aclut
               dcluts(k) = dcluts(k) * dclut
c
            enddo
c
         else
c
            nclut = 0
c
         endif
c
      endif      
c
c***********************************************************************
c
c                    End of initialization
c
c***********************************************************************
c
      tslow_sc = tslow_sc + dtp  !  Increment slow-time
c
      do iband = 1 , nbands
c
c   Offset slow-time for each of the step-chirps
c
         tslow = tslow_sc + ( dtp / float( nbands ) ) *
     .                    ( float( iband ) - 0.5 * float( nbands + 1 ) )
c
c   If there are targets, pre-compute numbers used in generating
c   the synthetic pulse
c
         if ( ( ntest .gt. 0 ) .or. ( nclut .gt. 0 ) ) then
c
            bc       = ( br * 1.0E+12 ) / ( clight ** 2 )
c
c   Compute velocity noise realization using noise standard deviation
c   and decay time. [vnoise and vntime]
c
            if ( vntime .lt. dtp ) then
c
               epsvn = 1.0
c
            else
c
               epsvn = dtp / vntime
c
            endif
c
            vn      = ( 1.0 - epsvn ) * vn + epsvn * vnoise * rnor( 0 )
c
c   Global amplitude modulation function - linear and periodic in time.
c
            amp     = 1.0 + dadt * tslow +
     .                aamp * sin( twopi * afreq * tslow )
c
c   Linear change of acceleration with time and periodic velocity
c   change with time
c
            vprime  = daccdt * 0.5 * lambda * 0.5 * tslow * tslow +
     .                  vamp  * sin( twopi * vfreq * tslow ) +
     .                  vamp2 * sin( twopi * vfreq2 * tslow )
c
c   Add random and deterministic time-dependent velocity to strip-map
c
            if ( nfmt .eq. - 4 .or. nfmt .eq. - 5 ) rsprime = rsprime +
     .                          ( vprime + vn ) * ( tslow - tslow_last )
c
c   Compute clutter parameters
c     
            if ( nclut .ne. 0 ) then
c
               do k = 1 , nclut
c
                  if ( dcluts(k) - tslow .lt. 0.0 ) then
c
                     dcluts(k) = tslow + dclut
                     cr(k)     = uni( 0 )
                     cr(k)     = rmin + cr(k) * ( rmax - rmin )
c
                  endif   
c
                  acluts(k) = aclut * sin( pi *
     .                               ( ( dcluts(k) - tslow ) / dclut ) )
                  cr(k)     = cr(k) + cv(k) * ( tslow - tslow_last )
c
               enddo
c
            endif         
c
c   Check to see whether we are using ship rotation simulation.
c
            if ( nfmt .eq. - 1 .or. nfmt .eq. - 2 ) then
c
c        We are NOT using ship rotation simulation.
c
               do k = 1 , ntest
c
                  r(k) = r(k) + v(k) * ( tslow - tslow_last )
                  v(k) = v0(k) + vn + tgdvdt(k) * tslow + vprime
c
               enddo
c
            else if ( nfmt .eq. - 3 ) then
c
c        We are using ship rotation simulation.
c
c           Rotate the ship to its position at the current time
c
               call rot_ship2 ( ntest , inyaw , inpitch , inroll ,
     .                          avelyaw , avelpitch , avelroll ,
     .                          aacyaw , aacpitch , aacroll , pyawa ,
     .                          pyawp , pyawt , ppitcha , ppitchp ,
     .                          ppitcht , prolla , prollp , prollt ,
     .                          x_scat , y_scat , z_scat ,
     .                          newxsc , newysc , newzsc , tslow )
c
               do k = 1 , ntest
c
c              The range change is the superposition of the change due to
c              pure rotation and that due to other velocity perturbations
c              like translation.
c
                  r(k)        = newxsc(k)
c
                  v(k)        = ( newxsc(k) - oldrange(k) ) /
     .                          ( tslow - tslow_last ) 
c
                  oldrange(k) = newxsc(k)
c
c              Calculate the perturbation velocity (the velocity not
c              associated with rotation) and the range perturbation.
c
                  vpert(k)    = vbar + slipv + tgdvdt(k) * tslow +
     .                          vn + vprime
c
c   Target type 2 is used to mark a target which does not
c   have vibration.  This is used to simulate a mast which can
c   vibrate independently of the other scatterers.
c
                  if ( tgt_type(k) .eq. 2 ) then
c
                     if ( svtime(k) .lt. dtp ) then
c
                        epsvn = 1.0
c
                     else
c
                        epsvn = dtp / svtime(k)
c
                     endif
c
                     svactl(k) = ( 1.0 - epsvn ) * svactl(k) +
     .                           epsvn * svnois(k) * rnor( 0 )
c
                     vpert(k)  = vpert(k) + svactl(k) + svamp(k) *
     .                           sin( twopi * svfreq(k) * tslow )

                  else if ( tgt_type(k) .eq. 3 ) then
c
c                    tgt type 3 is a periodic scatterer with cosine
c                    amplitude function over the duration.
c
                     ptime = mod( tslow + spstart(k) , 1.0 / spfreq(k) )
c
                     if ( tslow .gt. spstart(k) .and.
     .                    ptime .lt. spdur(k) ) then
c
c   May 13, 1998:  Modified this formula to make the amplitude positive
c   definite - added 1 to the cosine.
c
                        tgtamp(k) = sqrt( abs( spamp(k) ) ) * ( 1.0 +
     .                              cos( pi * ( ptime * 2.0 - spdur(k) )
     .                                         / spdur(k) ) )
c
                     else
c
                        tgtamp(k) = 0.0
c
                     endif
c
                  endif
c
                  rpert(k) = rpert(k) + vpert(k) *
     .                                  ( tslow - tslow_last )
c
                  v(k)     = v(k) + vpert(k)
c
               enddo
c
            else if ( nfmt .eq. - 4 ) then
c
c   Strip-map Simple Test Pattern
c
               do k = 1 , ntest
c
                  slant       = slant0 + rsprime + float( k - middle )
     .                                           * drtest
c
                  if ( abs( k - middle ) .le. 3 ) then
                     toff = 0.0
                  else
                     toff = ( slant - slant0 ) / vplat
                  endif
c
c   Added extra second so the target does not repeat so fast.KAM 7/21/99
c
                  tprime      = amod( tinteg + 1.0 + tslow + toff ,
     .                                tinteg + 1.0 )
     .                                - 0.5 * ( tinteg + 1.0 )
c
c   Weight by an effective illumination function so that the transition
c   at high doppler doesn't cause sidelobe type effects 
c
c   The illumination function assumes 2 degree beamwidth (0.035 rad)
c   with a footprint of slant0 * 0.035 converted to an integration time
c   with vplat.  Hence the sigma = 0.85 * beamwidth * slant0 / vplat.
c   This is sigma where 0.5 = exp( - 0.5 * x^2 / sigma^2 ) and where
c   x = 0.035 * slant0.
c                           KAM 7/21/99
c 
                  tgtamp(k)   = tgtamp0(k) * exp( - 0.5 * (
     .                          ( tprime * vplat ) /
     .                          ( 0.85 * beamwidth * slant0 ) ) ** 2 )
c
                  cross       = vplat * tprime
c
                  r(k)        = sqrt( slant ** 2 + cross ** 2 ) - slant0
c
                  v(k)        = vprime + vplat * vplat * tprime /
     .                                   ( r(k) + slant0 )
c
                  refrange(k) = 0.0  !  Slant Range
c
               enddo
c
            else if ( nfmt .eq. - 5 ) then
c
c   Strip-map Idealized Test Pattern (Quadratic Approximation to Range)
c
               do k = 1 , ntest
c
                  slant       = slant0 + rsprime + float( k - middle )
     .                                           * drtest
c
                  toff        = ( slant - slant0 ) / vplat
c
                  tprime      = amod( tinteg + tslow + toff , tinteg )
     .                          - 0.5 * tinteg
c
c   Weight by an effective illumination function so that the transition
c   at high doppler doesn't cause sidelobe type effects
c
c              Old illumination function left as a legacy.
c
                  tgtamp(k)   = tgtamp0(k) *
     .                          exp( - 8.0 * ( tprime / tinteg ) ** 2 )
c
                  r(k)        = slant - slant0 + 0.5 *
     .                          vplat * vplat * tprime * tprime / slant0
c
                  v(k)        = vprime + vplat * vplat * tprime / slant0
c
                  refrange(k) = 0.0  !  slant
c
               enddo
c
            else if ( nfmt .eq. - 6 ) then
c
c   Strip-map Simple Test Pattern
c
               do k = 1 , ntest
c
                  slant       = slant0 + rsprime + r0(k)
c
                  toff        = x0(k) / vplat
c
c   Added extra second so the target does not repeat so fast.KAM 7/21/99
c
                  tprime      = amod( tinteg + 1.0 + tslow + toff ,
     .                                tinteg + 1.0 )
     .                                - 0.5 * ( tinteg + 1.0 )
c
c   Weight by an effective illumination function so that the transition
c   at high doppler doesn't cause sidelobe type effects 
c
c   The illumination function assumes 2 degree beamwidth (0.035 rad)
c   with a footprint of slant0 * 0.035 converted to an integration time
c   with vplat.  Hence the sigma = 0.85 * beamwidth * slant0 / vplat.
c   This is sigma where 0.5 = exp( - 0.5 * x^2 / sigma^2 ) and where
c   x = 0.035 * slant0.
c                           KAM 7/21/99
c 
                  tgtamp(k)   = tgtamp0(k) * exp( - 0.5 * (
     .                          ( tprime * vplat ) /
     .                          ( 0.85 * beamwidth * slant0 ) ) ** 2 )
c
                  cross       = vplat * tprime
c
                  r(k)        = sqrt( slant ** 2 + cross ** 2 ) - slant0
c
                  v(k)        = vprime + vplat * vplat * tprime /
     .                                   ( r(k) + slant0 )
c
                  refrange(k) = 0.0  !  Slant Range
c
               enddo
c
            endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Pre-compute the constants used in calculating the phase history
c   for the scatterers & clutter
c
            do k = 1 , ntest
c
               p(k) = - 2.0 * ( ( ( r(k) - refrange(k) - slipv * tslow )
     .                          / lambda ) - bc * r(k) ** 2 )
               c(k) = - 2.0 * ( ( r(k) * br * 1.0E+12 / clight ) +
     .                          ( v(k) / lambda ) )
c
            enddo
c
            do k = 1 , nclut
c
               cp(k) = - 2.0 * ( ( ( cr(k) - slipv * tslow ) / lambda )
     .                           - bc * cr(k) ** 2 )
               cc(k) = - 2.0 * ( ( cr(k) * br * 1.0E+12 / clight ) +
     .                           ( cv(k) / lambda ) )
c
            enddo
c
         endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Finally, all preparation is done.  Form the uncompressed range pulse
c   from the information on the noise power, the test targets and the
c   clutter.
c
c   Fast-time offset for each band
c
         tfoff   = ( ( dtr * float( ntrband ) )
     .             - ( rfdelt / ( br * 1.0E+12 ) ) ) *
     .             ( float( iband ) - 0.5 * float( nbands + 1 ) )
c
         do i = 1 , ntrband
c
            ii      = i + ntrband * ( iband - 1 )
c
            tfast   = float( ii - 1 - ntr / 2 ) * dtr - tfoff
c
c   First, compute a carpet of Gaussian noise
c
            r1      = rnor( 0 )
            r2      = rnor( 0 )
            crc(ii) = noise * cmplx( r1 , r2 )
c
c   Next, if there are test targets add them to the noise
c
            if ( ntest .gt. 0 ) then
c
               do k = 1 , ntest
c
                  ph      = p(k) + c(k) * tfast
                  crc(ii) = crc(ii) + tgtamp(k) * amp * cpfast( ph )
c
               enddo
c
            endif
c
c   If there are clutter noise features, add them in
c         
            if ( nclut .gt. 0 ) then
c
               do k = 1 , nclut
c
                  ph      = cp(k) + cc(k) * tfast
                  crc(ii) = crc(ii) + acluts(k) * cpfast( ph )
c
               enddo
c
            endif
c
c   To test the adaptive phase correction, shift the middle band by 90
c   degrees
c
            if ( iband .eq. 1 + nbands / 2 )
     .         crc(ii) = crc(ii) * cmplx( 0.0 , 1.0 )
c
         enddo  !  i = 1 , ntrband
c
         tslow_last = tslow
c
      enddo     !  iband = 1 , nbands
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      nerror = 0  !  Normal return
c
      return
c
 1000 nerror = 1  !  Error return - end of file encountered
c
      write ( 6 , * ) ' End of file encountered in IQGEN'
      write ( 7 , * ) ' End of file encountered in IQGEN'
c
      return
c
 2000 nerror = 2  !  Error return - file I/O error on open or read
c
      write ( 6 , * ) ' File I/O error encountered in IQGEN'
      write ( 7 , * ) ' File I/O error encountered in IQGEN'
c
      return
      end
