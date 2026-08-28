C**
C***********************************************************************
C**
      subroutine gettgt_s_b ( csbimg , nakeep , nsr , nabuff , pnoise ,
     .                        pgatgt , nrcell , dts , iisub , isub ,
     .                        work , nptgt , rscore , dopcen , curtime )
C**
C***********************************************************************
C**
      implicit none
c
      integer nakeep , nsr , nabuff , i , j , jindex , nptgt , iisub ,
     .        isub , sfirst , slast , isb , actual , idop , nbest
c
      complex csbimg(nakeep,nsr,nabuff) , work(nabuff,nakeep) , cdeltf ,
     .        covt , cpfast
c
      real    pnoise , pgatgt(nsr,5) , vel , d0 , dd , aavg , npower ,
     .        ddw , rscore(nsr) , peak , fcon , tgtime , tgfreq , dts ,
     .        pgatmp(5) , fdop , ltime , time0 , db , rtgt , athrsh ,
     .        cfmax , correl , pwr , dopcen , dopsum , dwsum , curtime
c
      integer nrcell(nsr) , itgt , itg , nrfcell , nrlose , jlocal ,
     .        nka , nkd , nks , nkf , nktot , kill
c
c-----------------------------------------------------------------------
c
c   Include files
c
      include     'sarprm.h'      !  Standard ISAR-T parameters
c
      include     'updates.h'     !  Updates to parameters from the
c                                    first major release of the code
c
      include     'realtime.h'    !  Real-time parameters
c
      include     'tglist.h'
c
c-----------------------------------------------------------------------
c
c   Get Targets: Selects targets and computes track file;
c   pgatgt: centroid, width, pwr, accel, snr vs range line
c
c   Loop over range lines and compute track file
c
      actual = min( isub , nabuff / 3 )
      time0  = curtime - 0.5 * float(actual-1) * dts
c
      ddw    = 0.2 * dff
      athrsh = 5.0 / ( ( 0.25 * float( actual ) * dts ) ** 2 )
c
c-----------------------------------------------------------------------
c
c   Don't do the algorithm unless there are at leat 16 sub-images in the
c   half-buffer used for target detection
c
      if ( actual .lt. 16 ) return
c
c-----------------------------------------------------------------------
c
      slast  = iisub
c
      sfirst = slast - actual + 1
c
      if ( sfirst .lt. 1 ) sfirst = sfirst + nabuff
c
c   Allow buffer room at top and bottom for range-walk.  Use double the
c   nominal maximum number of cells.
c
      nrlose = 2 * nint( ( dfc * float( nakeep / 2 ) * 0.5 * lambda *
     .                           float( actual / 2 ) * dts ) / drs )
c
      do j = 1 + nrlose , nsr - nrlose
c
         do idop = 1 , nakeep
c
c   Doppler frequency for this sub-image cell
c
            fdop = dfc * float( idop - 1 - nakeep / 2 )
c
c   Loop over time
c
            isb  = sfirst
c
            do i = 1 , actual
c
c   Correct phase to motion-compensate to the center of this cell
c
               ltime        = dts * ( float( i - 1 - actual / 2 ) )
c
               cdeltf       = cpfast( ltime * fdop )
c
c   Use nearest neighbor range cell, accounting for range-walk
c
               jlocal       = j + nint( ltime * fdop * 0.5 * lambda
     .                                  / drs )
c
c   Fetch time series value for this cell
c
               work(i,idop) = csbimg(idop,jlocal,isb) * cdeltf
c
c   Update time index
c
               isb          = isb + 1
c
               if ( isb .gt. nabuff ) isb = isb - nabuff
c
            enddo
c
         enddo
c
         cfmax = 0.0
c
         do idop = 1 , nakeep
c
            covt = 0.0
c
            pwr  = 0.0
c
c   Use middle pulse pairs
c
            do i = 1 , actual - 1
c
               covt = covt + work(i,idop) *
     .                conjg( work(i+1,idop) )
c
               pwr  = pwr  + work(i,idop) *
     .                conjg( work(i,idop) )
c
            enddo
c
            pwr = pwr  + 0.5 * ( work(actual,idop) *
     .                    conjg( work(actual,idop) )
     .                         - work(1,idop) *
     .                    conjg( work(1,idop) ) )
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
     .           correl .gt. cfmax ) then
c
               nrfcell = idop
c
               cfmax   = correl
c
            endif
c
         enddo    !  Doppler cells
c
c   Estimate the spectral parameters
c
         if ( nrfcell .gt. 0 ) then
c
            if ( rt_pga .eq. 0 ) then
c
c   Normal version of dechirp
c
               call dchirp    ( work(1,nrfcell) , actual , pgatmp(4) ,
     .                          vel , pgatmp(3) , d0 , dd , dts , - 1 ,
     .                          aavg , lambda , pgatmp(1) , pgatmp(2) )
c
            else
c
c   Real-time version of dechirp
c
               call dchirp_rt ( work(1,nrfcell) , actual , pgatmp(4) ,
     .                          vel , pgatmp(3) , d0 , dd , dts , - 1 ,
     .                          aavg , lambda , pgatmp(1) , pgatmp(2) )
c
            endif
c
c   Add in the center frequency for this Doppler cell
c
            pgatmp(1)   = pgatmp(1) +
     .                    float( nrfcell - 1 - nakeep / 2 ) * dfc
c
c   Choose the highest score for each range line
c
            pgatmp(5)   = abs( pgatmp(3) / ( pgatmp(2) + ddw ) )
c
            pgatgt(j,:) = pgatmp(:)
c
         else
c
            pgatgt(j,:) = 0.0
c
         endif
c
      enddo  !  Range cells
c
c   Fill in the ends by extrapolation
c
      do j = 1 , nrlose
c
         pgatgt(j,:)            = pgatgt(1+nrlose,:)
         pgatgt(nsr-nrlose+j,:) = pgatgt(nsr-nrlose,:)
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Sort scores to determine the noise floor
c
      rscore(:) = pgatgt(:,3)
c
      call sort ( rscore , nsr , nrcell )
      jindex    = nint( 0.01 * pnoise * float( nsr ) )
      if ( jindex .eq. 0 ) jindex = 1
c
c   Bullet-proof against very small or zero noise values by limiting it
c   to a small fraction of the peak cross section or, if the peak is
c   zero, set it to one.
c 
      npower    = amax1( rscore(jindex) , 1.0E-5 * rscore(nsr) )
      if ( npower .eq. 0.0 ) npower = 1.0
c
c-----------------------------------------------------------------------
c
c   Scale SNR by noise level and remove large acceleration, large
c   Doppler width points and points more than 40 dB down from the
c   peak value
c
      peak      = rscore(nsr) / ( npower * 2.0 * dts )
c
      fcon      = 0.15
c
      nka       = 0
      nkd       = 0
      nks       = 0
      nkf       = 0
c
      do j = 1 , nsr
c
         kill = 0
         pgatgt(j,5) = pgatgt(j,5) / ( npower * 2.0 * dts )
c
         if ( abs( pgatgt(j,4) ) .ge. athrsh ) then
            kill = 1
            nka  = nka + 1
         endif
c            
         if ( pgatgt(j,2) .gt. fcon * float( npass ) * dff ) then
            kill = 1
            nkd  = nkd + 1
         endif
c            
         if ( pgatgt(j,5) .lt. 1.0E-4 * peak ) then
            kill = 1
            nks  = nks + 1
         endif
c
         if ( abs( pgatgt(j,1) ) .gt. 0.4 * float( nakeep ) * dfc ) then
            kill = 1
            nkf  = nkf + 1
         endif
c            
         if ( kill .gt. 0 ) pgatgt(j,5) = - abs( pgatgt(j,5) )
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Mark targets with the same Dopplers as the strongest targets for 
c   deletion
c
c-----------------------------------------------------------------------
c
      nbest = nint( overrg * drf / drc )
c
      call best ( pgatgt(1,5) , rscore , snrmin , nsr , nrcell , nptgt ,
     .            nbest )
c
      if ( nptgt .gt. 0 ) then
c
         if ( rt_pga .eq. 0 ) then
c
            if ( quiet .gt. 1 )
     .      write ( 6 ,'(/,10x,a)')
     .      'Sub-Image Targets (i,snr,f,da,a0,dw,r)'
c
            write ( 7 ,'(/,10x,a)')
     .      'Sub-Image Targets (i,snr,f,da,a0,dw,r)'
c
         endif
c
         dwsum  = 0.0
c
         dopsum = 0.0
c
         do itgt = 1 , nptgt
c
c   PGA sorted the SNR values from smallest to largest; nrcell contains
c   the true range index of the targets
c
            itg  = nrcell(nsr+1-itgt)
c
            rtgt = rsmin + drc * float( itg - 1 )
c
c   Target parameters returned from PGA:
c
c        k   -   pgatgt(itg,k)
c
c        1   -   Doppler centroid (Hz)
c        2   -   Doppler width (Hz)
c        3   -   Cross section
c        4   -   Acceleration (Hz/sec)
c        5   -   SNR
c
c   This information must be reported to routine ADDTGT so that future
c   frames can benefit from this focus information.
c
            if ( rt_pga .eq. 0 ) then
c
               if ( quiet .gt. 1 )
     .         write ( 6 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                   pgatgt(itg,1) , pgatgt(itg,4) , 0.0 ,
     .                   pgatgt(itg,2) , rtgt
c
               write ( 7 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                   pgatgt(itg,1) , pgatgt(itg,4) , 0.0 ,
     .                   pgatgt(itg,2) , rtgt
c
            endif
c
            tgtime = time0
            tgfreq = pgatgt(itg,1)
c
            call addtgt ( tgtime , rtgt , tgfreq , pgatgt(itg,4) ,
     .                    pgatgt(itg,2) , pgatgt(itg,5) , 's' )
c
            dopsum = dopsum + tgfreq / pgatgt(itg,2)
c
            dwsum  = dwsum  + 1.0 / pgatgt(itg,2)
c
         enddo
c
         dopcen = dopsum / dwsum
c
      else
c
         dopcen = 0.0
c
      endif
c
c-----------------------------------------------------------------------
c
      return
      end
