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
     .        isub , sfirst , slast , isb , actual , idop , id , nbest
c
      complex csbimg(nakeep,nsr,nabuff) , work(nabuff,nakeep) , cdeltf ,
     .        covt , cpfast , dcov
c
      real    pnoise , pgatgt(nsr,2,5) , vel , d0 , dd , aavg , npower ,
     .        ddw , rscore(2*nsr) , peak , tgtime , tgfreq , dts , dft ,
     .        pgatmp(5) , fdop , ltime , time0 , db , rtgt , athrsh ,
     .        covmax  , curtime , pwr , cova , covamax , covabs ,
     .        cormax , correl(nakeep), d , pi , dcen , arbar , afbar ,
     .        rlocal , dotdot , cycles , dflim
c
      integer nrcell(2*nsr) , itgt , itg , nrfcell1 , nrfcell2 , k ,
     .        nrlose , jlocal , nka , nkd , nks , nkf , nktot , kill ,
     .        ir , i2 , nrfcellm , nrfcellp , tgt_si_use
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
      if ( tgt_si .eq. 4 ) then     ! MCA compatibility mode
         dflim      = nakeep * dfc
         tgt_si_use = 1
      else                          ! Normal mode
         dflim      = 2 * dfc
         tgt_si_use = tgt_si
      endif 
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
         rlocal  = ( j - 1 - nsr / 2 ) * drs
c
c        Phase 1: Load complex array for time and Doppler
c
         nrfcell1 = 0
         nrfcell2 = 0
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
            isb    = sfirst
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
         do id = 1 , nakeep
c
            if ( mod(j,2) .eq. 0 .or. tgt_si .eq. 4 ) then
                idop   = id
            else
                idop   = nakeep + 1 - id
            endif
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
               covabs       = cabs( covt )
               correl(idop) = covabs / pwr
c
            else
c
               covabs       = 0.0
               correl(idop) = 0.0
c
            endif
c
            if ( covabs .gt. covmax ) then
c
               covmax  = covabs
               covamax = cova
c
               if ( ( correl(idop) .gt. 0.67 )     .and.
     .              ( correl(idop) .gt. cormax ) ) then
c
                  nrfcell1 = idop
                  cormax   = correl(idop)
                  
               endif
c
            endif
c           
         enddo    !  Doppler cells
c
         if ( nrfcell1 .gt. 0 .and. tgt_si_use .gt. 1 ) then
c
c           Find secondary Doppler cell
c
            cormax   = 0.0
c
            nrfcellm = max(1,nrfcell1-1)
            nrfcellp = min(nakeep,nrfcell1+1)
c
            do idop = 1 , nakeep
                
               if ( ( correl(idop) .gt. 0.67 )     .and.
     .              ( correl(idop) .gt. cormax ) ) then
                  
                  if ( idop .ne. nrfcell1 .and. idop .ne. nrfcellm .and. 
     .                 idop .ne. nrfcellp ) then
                     
                     nrfcell2 = idop
                     cormax   = correl(idop)

                  endif
               
               endif
            
            enddo
c
         endif
c
c        Estimate the spectral parameters
c
         if ( nrfcell1 .gt. 0 ) then
c
            call dchirp    ( work(1,nrfcell1) , actual , pgatmp(4) ,
     .                       vel , pgatmp(3) , d0 , dd , dts , - 2 ,
     .                       aavg , lambda , pgatmp(1) , pgatmp(2) )
c
c           Add in the center frequency for this Doppler cell
c
            pgatmp(1)     = pgatmp(1) + dcen +
     .                      float( nrfcell1 - 1 - nakeep / 2 ) * dfc
c
c           Calculate the score
c
            pgatmp(5)     = abs( pgatmp(3) / ( pgatmp(2) + ddw ) ) /
     .                      ( 2.0 * dts )
c
            pgatgt(j,1,:) = pgatmp(:)
c
         else
c
c           Fill with dummy info that will allow the noise level to be
c           determined - this point will be edited later
c
            pgatgt(j,1,1) = dcen
            d             = amax1( 0.0 , 1.0 - covmax ) / covamax
            pgatgt(j,1,2) = sqrt( 2.0 * d ) / ( 2.0 * pi * dts )
            pgatgt(j,1,3) = covmax
            pgatgt(j,1,4) = 2.0 * athrsh  !  Guaranteed to be edited
            pgatgt(j,1,5) = covmax / ( pgatgt(j,1,2) + ddw )
     .                             / ( 2.0 * dts )
c
         endif
c
         if ( nrfcell2 .gt. 0 ) then
c
            call dchirp    ( work(1,nrfcell2) , actual , pgatmp(4) ,
     .                       vel , pgatmp(3) , d0 , dd , dts , - 2 ,
     .                       aavg , lambda , pgatmp(1) , pgatmp(2) )
c
c           Add in the center frequency for this Doppler cell
c
            pgatmp(1)     = dcen + pgatmp(1) +
     .                      float( nrfcell2 - 1 - nakeep / 2 ) * dfc
c
c           Calculate the score
c
            pgatmp(5)     = abs( pgatmp(3) / ( pgatmp(2) + ddw ) ) /
     .                      ( 2.0 * dts )
c
            pgatgt(j,2,:) = pgatmp(:)
c
         else
c
c           Fill with dummy info that will allow the noise level to be
c           determined - this point will be edited later
c
            pgatgt(j,2,:) = pgatgt(j,1,:)
c
         endif
c
      enddo  !  Range cells
c
c     Fill in the ends by extrapolation
c
      do j = 1 , nrlose
c
         pgatgt(j,:,:)            = pgatgt(1+nrlose,:,:)
         pgatgt(nsr-nrlose+j,:,:) = pgatgt(nsr-nrlose,:,:)
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Sort scores to determine the noise floor
c
      rscore(1:nsr)       = pgatgt(:,1,5)
      rscore(1+nsr:2*nsr) = pgatgt(:,2,5)
c
c     Use only the primary set of targets (the first nsr)
c
      call sort ( rscore , nsr , nrcell )
      jindex              = nint( 0.01 * pnoise * float( nsr ) )
      if ( jindex .eq. 0 ) jindex = 1
c
c   Bullet-proof against very small or zero noise values by limiting it
c   to a small fraction of the peak cross section or, if the peak is
c   zero, set it to one.
c 
      npower              = amax1( rscore(jindex) ,
     .                             1.0E-6 * rscore(nsr) )
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
         do k = 1 , 2
c
            kill          = 0
            pgatgt(j,k,5) = pgatgt(j,k,5) / npower
c
            if ( abs( pgatgt(j,k,4) ) .ge. athrsh ) then
               kill = 1
               nka  = nka + 1
            endif
c            
            if ( pgatgt(j,k,2) .gt. 0.075 * float( npass ) * dft ) then
               kill = 1
               nkd  = nkd + 1
            endif
c            
            if ( pgatgt(j,k,5) .lt. 1.0E-4 * peak ) then
               kill = 1
               nks  = nks + 1
            endif
c
            if ( abs( pgatgt(j,k,1) ) .gt.
     .           0.4 * float( nakeep ) * dfc ) then
               kill = 1
               nkf  = nkf + 1
            endif
c            
            if ( kill .gt. 0 ) pgatgt(j,k,5) = - abs( pgatgt(j,k,5) )
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     Thin the targets in range based on the oversampling factor and
c     proximity in Doppler, then sort the full set of targets
c
      nbest = nint( overrg * drf / drs )
c
      call best_f ( pgatgt(1,1,5) , rscore , snrmin , nsr , nrcell ,
     .              pgatgt(1,1,1) , dflim , nptgt , nbest , nrlose ,
     .              tgt_si_use )
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
            itg  = nrcell(2*nsr+1-itgt)
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
            if ( itg .le. nsr ) then
               ir = itg
               i2 = 1
            else
               ir = itg-nsr
               i2 = 2
            endif
            
            rtgt = rsmin + drs * float( ir - 1 )
c
            if ( quiet .gt. 1 )
     .      write ( 6 ,'(i4,6f9.3)' ) itgt , db( pgatgt(ir,i2,5) ) ,
     .                pgatgt(ir,i2,1) , pgatgt(ir,i2,4) , 0.0 ,
     .                pgatgt(ir,i2,2) , rtgt
c
            write ( 7 ,'(i4,6f9.3)' ) itgt , db( pgatgt(ir,i2,5) ) ,
     .                pgatgt(ir,i2,1) , pgatgt(ir,i2,4) , 0.0 ,
     .                pgatgt(ir,i2,2) , rtgt
c
            tgtime = time0
            tgfreq = pgatgt(ir,i2,1)
c
            call addtgt ( tgtime , rtgt , tgfreq , pgatgt(ir,i2,4) ,
     .                    pgatgt(ir,i2,2) , pgatgt(ir,i2,5) , 's' )
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
      subroutine best_f ( score , rscore , scmin , nr , nrcell ,
     .                    dopplr , df , nptgt , nbest , nrlose ,
     .                    tgt_si )
C**
C***********************************************************************
C**
c   Select points which are isolated maxima and above a threshold - then
c   sort them in ascending order
c
      implicit none
c
      integer nr , nrcell(2*nr) , j , k , nptgt , nbest , nrlose ,
     .        nrlim , tgt_si
c
      real    score(nr,2) , rscore(2*nr) , dopplr(nr,2) , df , scmin
c
      logical tplus , tminus
c
c     Ignore the end regions
c
      nrlim = nbest + nrlose
      do j = 1 , nrlim
c
         score(j,:)      = - abs( score(j,:) )
         score(1+nr-j,:) = - abs( score(1+nr-j,:) )
c
      enddo
c
c     Use only points which are above the threshold and which are local
c     maxima
c
      do j = 1 + nrlim , nr - nrlim
c
         do k = 1 , 2
            
             if ( abs( score(j,k) ) .le. scmin ) then
c
               score(j,k) = - abs( score(j,k) )

c
            endif
         
         enddo
c
      enddo
c
      do j = 1 + nrlim , nr - nrlim
c
         do k = 1 , nbest
c
c-----------------------------------------------------------------------
c
c           Primary targets
c
            tplus  = ( abs( score(j,1) ) .lt. abs( score(j+k,1) ) )
     .            .and. ( abs( dopplr(j,1) - dopplr(j+k,1) ) .lt. df )
c
            tminus = ( abs( score(j,1) ) .lt. abs( score(j-k,1) ) )
     .            .and. ( abs( dopplr(j,1) - dopplr(j-k,1) ) .lt. df )
c
            if ( tplus .or. tminus ) score(j,1) = - abs( score(j,1) )
c         
c-----------------------------------------------------------------------
c
c           Secondary targets
c
            if ( tgt_si .gt. 1 ) then
                
               tplus  = ( abs( score(j,2) ) .lt. abs( score(j+k,2) ) )
     .            .and. ( abs( dopplr(j,2) - dopplr(j+k,2) ) .lt. df )
c
               tminus = ( abs( score(j,2) ) .lt. abs( score(j-k,2) ) )
     .            .and. ( abs( dopplr(j,2) - dopplr(j-k,2) ) .lt. df )
c
               if ( tplus .or. tminus ) score(j,2) = - abs( score(j,2) )
         
               tplus  = ( abs( score(j,2) ) .lt. abs( score(j+k,1) ) )
     .            .and. ( abs( dopplr(j,2) - dopplr(j+k,1) ) .lt. df )
c
               tminus = ( abs( score(j,2) ) .lt. abs( score(j-k,1) ) )
     .            .and. ( abs( dopplr(j,2) - dopplr(j-k,1) ) .lt. df )
c
               if ( tplus .or. tminus ) score(j,2) = - abs( score(j,2) )
c         
            endif
c
c-----------------------------------------------------------------------
c
         enddo
c
      enddo
c
c     Choose which targets are uses - primary, secondary or both
c
      if ( tgt_si .eq. 1 ) score(:,2) = - abs( score(:,2) )
      if ( tgt_si .eq. 2 ) score(:,1) = - abs( score(:,1) )
c
c     For SNR quality compute number of targets
c
      nptgt  = 0
c
      do j = 1 , nr
c
         rscore(j)    = score(j,1)
         if ( score(j,1) .gt. scmin ) nptgt = nptgt + 1
         
         rscore(j+nr) = score(j,2)
         if ( score(j,2) .gt. scmin ) nptgt = nptgt + 1
c
      enddo
c
c     Now sort them by score
c
      call sort ( rscore , 2*nr , nrcell )
c
      return
      end

