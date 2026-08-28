C**
C***********************************************************************
C**
      subroutine gettgt_s ( csbimg , nakeep , nsr , nabuff , pnoise ,
     .                      pgatgt , nrcell , dts , iisub , isub ,
     .                      dcen , work , nptgt , rscore , curtime ,
     .                      arbar , afbar )
C**
C***********************************************************************
C**
      implicit none
c
      integer nakeep , nsr , nabuff , i , j , jindex , nptgt , iisub ,
     .        isub , sfirst , slast , isb , actual , idop , nbest
c
      complex csbimg(nakeep,nsr,nabuff) , work(nabuff,nakeep) , cdeltf ,
     .        covt , cpfast , dcov
c
      real    pnoise , pgatgt(nsr,5) , vel , d0 , dd , aavg , npower ,
     .        ddw , rscore(nsr) , peak , tgtime , tgfreq , dts , dft ,
     .        pgatmp(5) , fdop , ltime , time0 , db , rtgt , athrsh ,
     .        covmax , correl , curtime , pwr , cova , covamax ,
     .        covabs , cormax , d , pi , cycles , dcen , arbar ,
     .        afbar , rlocal , dotdot
c
      integer nrcell(nsr) , itgt , itg , nrfcell , nrlose , jlocal ,
     .        k , nka , nkd , nks , nkf , nktot , kill
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
c-----------------------------------------------------------------------
c
      actual = min( isub , nabuff / 3 )
      actual = actual-mod(actual,4)
c
c   Don't do the algorithm unless there are at leat 16 sub-images in the
c   half-buffer used for target detection
c
      if ( actual .lt. 16 ) return
c
c-----------------------------------------------------------------------
c
      time0  = curtime - 0.5 * float(actual-1) * dts
c
      dft    = 1.0 / ( float( actual ) * dts )
      ddw    = 0.1 * dft
      athrsh = 5.0 / ( ( 0.25 * float( actual ) * dts ) ** 2 )
c
      pi     = atan2( 0.0 , - 1.0 )
c
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
      nrlose = 2 * nint( ( ( dcen + dfc * float( nakeep / 2 ) ) * 0.5 *
     .                   lambda * float( actual / 2 ) * dts ) / drs )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     Loop over range
c
      do j = 1 + nrlose , nsr - nrlose
c
         rlocal  = ( j - 1 - nsr / 2 ) * drc
c
c        Phase 1: Load complex array for time and Doppler
c
         nrfcell = 0
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c        Loop over Doppler
c
         do idop = 1 , nakeep
c
c           Doppler frequency for this sub-image cell
c
            fdop   = dcen + dfc * float( idop - 1 - nakeep / 2 )
c
            dotdot = ( 2.0 / lambda ) * arbar * rlocal + afbar * fdop
c
c           Loop over time
c
            isb  = sfirst
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c           Loop over time
c
            do i = 1 , actual
c
c              Correct phase to motion-compensate to the center of this cell
c
               ltime        = dts * ( float( i - 1 - actual / 2 ) )
c
               cycles       = ltime * ( fdop + 0.5 * dotdot * ltime )
c
               cdeltf       = cpfast( cycles )
c
c              Use nearest neighbor range cell, accounting for range-walk
c
               jlocal       = j + nint( cycles * 0.5 * lambda
     .                                  / drs )
c
c              Fetch time series value for this cell
c
               work(i,idop) = csbimg(idop,jlocal,isb) * cdeltf
c
c              Update time index, considering the circular buffer
c
               isb          = isb + 1
c
               if ( isb .gt. nabuff ) isb = isb - nabuff
c
            enddo  !  Time
c
         enddo     !  Doppler
c
c        Phase 2: Compute covariance and statistics
c
         covmax  = 0.0
c
         cormax  = 0.0
c
         covamax = 0.0
c
c        Second loop over Doppler - find the primary Doppler cell
c
         do idop = 1 , nakeep
c
            covt   = 0.0
c
            cova   = 0.0
c
            pwr    = 0.0
c
c           Second loop over time
c
            do i = 1 , actual - 1
c
               dcov = work(i,idop) * conjg( work(i+1,idop) )
c
               covt = covt + dcov
c
               cova = cova + cabs(dcov)
c
               pwr  = pwr  + work(i,idop) * conjg( work(i,idop) )
c
            enddo  !  Time
c
            pwr    = pwr  + 0.5 * ( work(actual,idop) *
     .                       conjg( work(actual,idop) )
     .                            - work(1,idop) *
     .                       conjg( work(1,idop) ) )
c
            if ( pwr .ne. 0.0 ) then
c
               covabs = cabs( covt )
               correl = covabs / pwr
c
            else
c
               covabs = 0.0
               correl = 0.0
c
            endif
c
            if ( covabs .gt. covmax ) then
c
               covmax  = covabs
               covamax = cova
c
               if ( (correl .gt. 0.67) .and. (correl .gt. cormax) ) then
c
                  nrfcell = idop
                  cormax  = correl
c
               endif
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
     .                          vel , pgatmp(3) , d0 , dd , dts , - 2 ,
     .                          aavg , lambda , pgatmp(1) , pgatmp(2) )
c
            else
c
c   Real-time version of dechirp
c
               call dchirp_rt ( work(1,nrfcell) , actual , pgatmp(4) ,
     .                          vel , pgatmp(3) , d0 , dd , dts , - 2 ,
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
            pgatmp(5)   = abs( pgatmp(3) / ( pgatmp(2) + ddw ) ) /
     .                    ( 2.0 * dts )
c
            pgatgt(j,:) = pgatmp(:)
c
         else
c
            pgatgt(j,1) = 0.0
            d           = amax1( 0.0 , 1.0 - covmax ) / covamax
            pgatgt(j,2) = sqrt( 2.0 * d ) / ( 2.0 * pi * dts )
            pgatgt(j,3) = covmax
            pgatgt(j,4) = 2.0 * athrsh
            pgatgt(j,5) = covmax / ( pgatgt(j,2) + ddw ) / ( 2.0 * dts )
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
      rscore(:) = pgatgt(:,5)
c
      call sort ( rscore , nsr , nrcell )
      jindex    = nint( 0.01 * pnoise * float( nsr ) )
      if ( jindex .eq. 0 ) jindex = 1
c
c   Bullet-proof against very small or zero noise values by limiting it
c   to a small fraction of the peak cross section or, if the peak is
c   zero, set it to one.
c 
      npower    = amax1( rscore(jindex) , 1.0E-6 * rscore(nsr) )
      if ( npower .eq. 0.0 ) npower = 1.0
c
      write( 97 , * ) 'Max to noise power  ' ,
     .                jindex , time0 , rscore(nsr) / npower
c
c-----------------------------------------------------------------------
c
c   Scale SNR by noise level and remove large acceleration, large
c   Doppler width points and points more than 40 dB down from the
c   peak value
c
      peak      = rscore(nsr) / npower
c
      nka       = 0
      nkd       = 0
      nks       = 0
      nkf       = 0
c
      do j = 1 , nsr
c
         kill        = 0
         pgatgt(j,5) = pgatgt(j,5) / npower
c
         if ( abs( pgatgt(j,4) ) .ge. athrsh ) then
            kill = 1
            nka  = nka + 1
         endif
c            
         if ( pgatgt(j,2) .gt. 0.075 * float( npass ) * dft ) then
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
      nbest = nint( overrg * drf / drs )
c
      call best_f ( pgatgt(1,5) , rscore , snrmin , nsr , nrcell ,
     .              nptgt , nbest , nrlose )
c
      if ( nptgt .gt. 0 ) then
c
         if ( quiet .gt. 1 )
     .   write ( 6 ,'(/,10x,a)')
     .   'Sub-Image Targets (i,snr,f,da,a0,dw,r)'
c
         write ( 7 ,'(/,10x,a)')
     .   'Sub-Image Targets (i,snr,f,da,a0,dw,r)'
c
         do itgt = 1 , nptgt
c
c   PGA sorted the SNR values from smallest to largest; nrcell contains
c   the true range index of the targets
c
            itg  = nrcell(nsr+1-itgt)
c
            rtgt = rsmin + drs * float( itg - 1 )
c
c   Target parameters returned from dchirp:
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
c
            if ( quiet .gt. 1 )
     .      write ( 6 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                pgatgt(itg,1) , pgatgt(itg,4) , 0.0 ,
     .                pgatgt(itg,2) , rtgt
c
            write ( 7 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                pgatgt(itg,1) , pgatgt(itg,4) , 0.0 ,
     .                pgatgt(itg,2) , rtgt
c
            tgtime = time0
            tgfreq = pgatgt(itg,1)
c
            call addtgt ( tgtime , rtgt , tgfreq , pgatgt(itg,4) ,
     .                    pgatgt(itg,2) , pgatgt(itg,5) , 's' )
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c
      return
      end
C**
C***********************************************************************
C**
      subroutine best_f ( score , rscore , scmin , nr , nrcell , nptgt ,
     .                    nbest , nrlose )
C**
C***********************************************************************
C**
c   Select points which are isolated maxima and above a threshold - then
c   sort them in ascending order
c
      implicit none
c
      integer nr , nrcell(nr) , j , k , nptgt , nbest , nrlose , nrlim
c
      real    score(nr) , rscore(nr) , scmin
c
c     Ignore the end regions
c
      nrlim = nbest + nrlose
      do j = 1 , nrlim
c
         score(j)      = - abs( score(j) )
         score(1+nr-j) = - abs( score(1+nr-j) )
c
      enddo
c
c     Use only points which are above the threshold and which are local
c     maxima
c
      do j = 1 + nrlim , nr - nrlim
c
         if ( abs( score(j) ) .le. scmin ) then
c
            score(j)  = - abs( score(j) )
c
         else
c
            do k = 1 , nbest
c
               if ( ( abs( score(j) ) .lt. abs( score(j+k) ) ) .or.
     .              ( abs( score(j) ) .lt. abs( score(j-k) ) ) )
     .                       score(j)  = - abs( score(j) )
c
            enddo
c
         endif
c
      enddo
c
c   For SNR quality compute number of targets
c
      nptgt = 0
c
      do j = 1 , nr
c
         rscore(j) = score(j)
         if ( score(j) .gt. scmin ) nptgt = nptgt + 1
c
      enddo
c
c   Now sort them by score
c
      call sort ( rscore , nr , nrcell )
c
      return
      end

