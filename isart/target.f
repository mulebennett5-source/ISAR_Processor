C**
C***********************************************************************
C**
      subroutine target ( c , nr , np , nr1 , nre , np1 , npe , dr ,
     .                    dt , cw , lambda , nalias , ntg , sg0 , d0 ,
     .                    dd , a , b , dc , dw , tgsnr , tgdw , tgvel ,
     .                    tgrvel , tgacc , tgr , tgdadt , vrbar ,
     .                    tgrmin , tgrmax , tgrbar , curtime , snoise ,
     .                    vr0 , dr0 , snrmin , athrsh , notch , vnotch ,
     .                    dvntch , rnotch , drntch , rwindo , vwindo ,
     .                    eslip )
C**
C***********************************************************************
C**
C     This routine detects point targets in range-compressed signal
C     history and estimates their paths.
C
C
C     Definitions:
C
C         **********  INPUT VARIABLES  **********
C
C       c(np,nr)       :  Range-compressed signal history (complex)
C
C       nr             :  Number of range cells
C
C       np             :  Number of pulses
C
C       dr             :  Range cell separation
C
C       dt             :  Pulse repetition interval
C
C       cw(np)         :  Complex work array
C
C       athrsh         :  Maximum acceleration
C
C       lambda         :  Wavelength (meters)
C
C       nalias         :  Number of alias paths to search
C
C       curtime        :  Time of current sub-image (seconds)
C
C       snoise         :  Noise floor
C
C       snrmin         :  SNR threshold for target detection
C
C       athrsh         :  Acceleration limit
C
C         **********  OUTPUT VARIABLES  **********
C
C       ntg            :  Number of targets identified
C
C       sg0(nr)        :  Radar cross section
C
C       d0(nr)         :  [1-covar] for raw data
C
C       a(nr)          :  Trend of cross section with time
C
C       dc(nr)         :  Doppler centroid (Hz)
C
C       tgsnr(nr)      :  Signal-to-noise-ratio of targets ( sg0*(1-d) )
C
C       tgdw(nr)       :  Doppler width for targets 
C
C       tgvel(n)       :  Target absolute velocities (m/s)
C
C       tgrvel(nr)     :  Target relative velocities (m/s)
C
C       tgacc(nr)      :  Target accelerations (Hz/s)
C
C       tgr(nr)        :  Average range of targets (meters)
C
C       tgdadt(nr)     :  Rate of change of amplitude with time
C
C       vrbar          :  Estimated average radial velocity
C
C       tgrbar         :  Estimated average range
C
C***********************************************************************
C
      implicit none
c
      integer nr , np , nalias , ntg , ir , ip , irp , irp1 , it , k ,
     .        ialias , np1 , npe , nr1 , nre , ntot , nguard , notch ,
     .        nend , nfirst , ntmid , nmid
c
      complex c(np,nr) , cw(np)
c
      real    tgsnr(nr) , tgdw(nr) , tgvel(nr) , tgrvel(nr) ,
     .        tgacc(nr) , tgr(nr) , tgdadt(nr)
c
      real    b(nr) , rvel , sg0(nr) , d0(nr) , a(nr) , dc(nr) ,
     .        dw(nr) , dd(nr) , dt , lambda , snrtmp , athrsh , abar ,
     .        dr , db , valias , del , rbar , r , delr , vel , pi ,
     .        rswath , snr , oldrv , dvel , vrbar ,  tgrbar , tgrmin ,
     .        tgrmax , ddw , n0 ,  snoise , snrmin , rwloss , snrmax ,
     .        vnotch , dvntch , rnotch , drntch , vr0 , dr0 , rwindo ,
     .        vwindo , n0mid , pmid , eslip , curtime , dplus , dminus ,
     .        snrlim
c
      complex skill(2,3)
c
      logical boss
c
c***********************************************************************
c
      include 'updates.h'
c
c***********************************************************************
c
      pi     = atan2( 0.0 , - 1.0 )
      valias = 0.5 * ( lambda / dt ) / prfrat
      rswath = dr * float( nr )
      ntot   = npe - np1 + 1
      ddw    = 0.1 / ( dt * float( ntot ) )
c
c   To keep alias targets from being detected, a guard band of cells
c   is required.  The value is the minimum number of grid cells allowed
c   between targets.  It is proportional to the maximum alias velocity
c   and the integration time.
c
      nguard = 8 + ifix( 0.5 *
     .         float( nalias / 2 ) * valias * float( ntot ) * dt / dr )
c
c***********************************************************************
c
c   First pass - select all range cells for which the SNR meets the
c   SNR threshold
c
c   Loop over range lines
c
c   In the first pass use only the middle fraction of the pulses to
c   save computer time
c
      pmid   = 0.25
c
c   Noise estimates for full record and middle
c
      n0     = snoise * float( ntot )
      n0mid  = n0 * pmid
c
c   Calculate loss of SNR due to range walk of incorrect aliases
c
      if ( nalias .eq. 1 ) then
         rwloss = 1.0
      else
         rwloss = dr /
     .            ( valias * float( nalias / 2 ) * float( ntot ) * dt )
      endif
c
      rwloss = amin1( 1.0 , rwloss / pmid )
c
      nfirst = np1 + ifix( float( ntot ) * ( 1.0 - pmid ) * 0.5 )
      nend   = nfirst - 1 + ifix( float( ntot ) * pmid )
      ntmid  = nend + 1 - nfirst
c
c***********************************************************************
c***********************************************************************
c-----------------  Calculate basic SNR values  ------------------------
c
c   This section requires the range-compressed signal history array but
c   no global information - it can therefore be split into overlapping
c   sections for computation be separate processors
c
      do ir = nr1 , nre
c
c   Compute first order model of time series - power, centroid,
c   covariance amplitude, acceleration, amplitude trend
c
         do ip = nfirst , nend
c
            cw(ip) = c(ip,ir)
c
         enddo
c
         call dchirp ( cw(nfirst) , ntmid , b(ir) , rvel , sg0(ir) ,
     .                 d0(ir) , dd(ir) , dt , 0 , a(ir) , lambda ,
     .                 dc(ir) , dw(ir) )
c
         tgsnr(ir) = snr( sg0(ir) , dw(ir) + ddw / pmid , dt , n0mid )
c
      enddo
c
c   Indentify those range lines which meet the SNR threshold and which
c   have doppler widths smaller than their neighbors
c
      ntg = 0
      do ir = nr1 + nguard , nre - nguard
c
         if ( tgsnr(ir) .gt. rwloss * snrmin ) then
c
            boss = .true.
            do irp = 1 , nguard
c
               if ( tgsnr(ir) .le. tgsnr(ir+irp) .or.
     .              tgsnr(ir) .le. tgsnr(ir-irp) ) boss = .false.
c
            enddo
c
            if ( boss ) then
c
               ntg        = ntg + 1
               tgdw(ir)   = dw(ir) + ddw
               tgvel(ir)  = 0.5 * lambda * dc(ir)
               tgrvel(ir) = tgvel(ir)
               tgacc(ir)  = b(ir)
               tgr(ir)    = dr * float( ir - 1 - nr / 2 )
               tgdadt(ir) = a(ir)
c
            else
c
               tgsnr(ir)  = 0.0
c
            endif
c
         endif
c
      enddo
c
      if ( quiet .gt. 1 ) write ( 7  , * ) 'First pass targets:  ' , ntg
c
      if ( ntg .eq. 0 ) then
c
         vrbar  = 0.0
         tgrbar = 0.0
c
         return
c
      endif
c
c***********************************************************************
c
c   Second pass - form time series along alias paths corresponding to
c   estimated centroids and repeat first order time series model
c
      do ir = nr1 + nguard , nre - nguard
c
         if ( tgsnr(ir) .gt. snrmin ) then
c
            tgsnr(ir) = 0.0
c
            vel   = tgvel(ir) - float( nalias / 2 ) * valias
            oldrv = tgrvel(ir)
c
c   r and rbar are range from the first range cell, not MCP
c
            rbar  = 0.5 * rswath + tgr(ir)
c
            do ialias = 1 , nalias
c
c   Form time series of range-compressed signal history along all
c   alias paths having the estimated Doppler frequency
c
               delr = vel * dt
               r    = rbar - 0.5 * delr * float( ntot )
c
c   Assume data is periodic in range
c
               if ( r .lt. 0.0 ) r = r + rswath
               if ( r .ge. rswath ) r = r - rswath
c
               do ip = np1 , npe
c
                  irp  = 1 + ifix( r / dr )
                  irp1 = irp + 1
                  if ( irp1 .gt. nr ) irp1 = irp1 - nr
                  del = ( r / dr ) - float( irp - 1 )
c
                  cw(ip) = del * c(ip,irp1) + ( 1.0 - del ) * c(ip,irp)
c
                  r      = r + delr
c
                  if ( r .lt. 0.0    ) r = r + rswath
                  if ( r .ge. rswath ) r = r - rswath
c
               enddo
c
               call dchirp ( cw(np1) , ntot , b(ir) , rvel , sg0(ir) ,
     .                       d0(ir) , dd(ir) , dt , - 1 , a(ir) ,
     .                       lambda , dc(ir) , dw(ir) )
c
               snrtmp  = snr( sg0(ir) , dw(ir) + ddw , dt , n0 )
c
               if ( snrtmp .gt. tgsnr(ir) ) then
c
                  tgsnr(ir)  = snrtmp
                  tgdw(ir)   = dw(ir) + ddw
c
c   Correct the new velocity for the new centroid estimate
c
                  dvel       = rvel - oldrv
                  if ( dvel .gt. + 0.5 * valias ) dvel = dvel - valias
                  if ( dvel .lt. - 0.5 * valias ) dvel = dvel + valias
c
                  tgvel(ir)  = vel + dvel
                  tgrvel(ir) = rvel
                  tgacc(ir)  = b(ir)
                  tgr(ir)    = dr * float( ir - 1 - nr / 2 )
                  tgdadt(ir) = a(ir)
c
               endif
c
               vel = vel + valias
c
            enddo
c
         endif
c
      enddo
c
c-----------------  Calculate basic SNR values  ------------------------
c***********************************************************************
c***********************************************************************
c-------------------  Sort target detections  --------------------------
c
c   This section does not require the range-compressed signal history
c   but requires global information over all range cells
c
      if ( notch .ne. 0 ) then
c
         if ( notch .eq. - 1 .or. notch .le. - 3 ) then
c
c   Exclude velocity notch
c
            do ir = nr1 + nguard , nre - nguard
c
               if ( abs( tgvel(ir) + vr0 - vnotch ) .lt. 0.5 * dvntch )
     .            tgsnr(ir) = - abs( tgsnr(ir) )
c
            enddo
c
         endif
c
c   Include velocity notch
c
         if ( notch .eq. + 1 .or. notch .ge. + 3 ) then
c
            do ir = nr1 + nguard , nre - nguard
c
               if ( abs( tgvel(ir) + vr0 - vnotch ) .gt. 0.5 * dvntch )
     .            tgsnr(ir) = - abs( tgsnr(ir) )
c
            enddo
c
         endif
c
c   Exclude range notch
c
         if ( notch .eq. - 2 .or. notch .le. - 3 ) then
c
            do ir = nr1 + nguard , nre - nguard
c
               if ( abs( tgr(ir) + dr0 - rnotch ) .lt. 0.5 * drntch )
     .            tgsnr(ir) = - abs( tgsnr(ir) )
c
            enddo
c
         endif
c
c   Include range notch
c
         if ( notch .eq. + 2 .or. notch .ge. + 3 ) then
c
            do ir = nr1 + nguard , nre - nguard
c
               if ( abs( tgr(ir) + dr0 - rnotch ) .gt. 0.5 * drntch )
     .            tgsnr(ir) = - abs( tgsnr(ir) )
c
            enddo
c
         endif
c
      endif
c
c   Zero out aliases
c
      call alias0 ( nr , tgsnr , tgvel , tgacc , tgr , abar ,
     .              rwindo , vwindo , nr1 + nguard , nre - nguard ,
     .              snrmin , snrmax )
c
c             Thin targets by SNR and acceleration thresholds and
c             estimate the average velocity and range
c
c   All detected targets are treated equally.  The average range is
c   estimated as the mean of the max and min ranges of detected targets.
c
      vrbar  = 0.0
      tgrmax = - 0.5 * rswath
      tgrmin = + 0.5 * rswath
c
      if ( quiet .gt. 1 ) then
c
         write ( 6  ,'(/,10x,a39)')
     .                        'Candidate Targets (i,snr,v,a,dw,r,dadt)'
         write ( 6  ,'(10x,a18,f16.3)' ) 'Current time =    ' , curtime
c
      endif
c
      write ( 7  ,'(/,10x,a39)')
     .                        'Candidate Targets (i,snr,v,a,dw,r,dadt)'
      write ( 7  ,'(10x,a18,f16.3)' ) 'Current time =    ' , curtime
c
      ntg    = 0
c
c   Limit the snr values to 0.0001 times the maximum to eliminate
c   sidelobes of big targets
c
      snrlim = amax1( snrmin , 1.0E-4 * snrmax )
c
      do ir = nr1 + nguard , nre - nguard
c
c   Apply minimum SNR constraint and both absolute and relative
c   acceleration limits
c
         if ( tgsnr(ir) .gt. snrlim .and.
     .        abs( tgacc(ir) ) .le. 4.0 * athrsh .and.
     .        abs( tgacc(ir) - abar ) .le. athrsh ) then
c
            ntg    = ntg + 1
            vrbar  = vrbar + tgvel(ir)
            rbar   = tgr(ir)
            tgrmax = amax1( tgrmax , rbar )
            tgrmin = amin1( tgrmin , rbar )
c
            if ( quiet .gt. 1 )
     .         write ( 6 ,'(i10,6f10.3)' ) ntg , db( tgsnr(ir) ) ,
     .                      tgvel(ir) , tgacc(ir) , tgdw(ir) , rbar ,
     .                      tgdadt(ir)
c
            write ( 7 ,'(i10,6f10.3)' ) ntg , db( tgsnr(ir) ) ,
     .                   tgvel(ir) , tgacc(ir) , tgdw(ir) , rbar ,
     .                   tgdadt(ir)
c
         else
c
            tgsnr(ir) = - abs( tgsnr(ir) )
c
         endif
c
      enddo
c
      if ( ntg .gt. 0 ) then
c
         vrbar  = vrbar / float( ntg )
         tgrbar = 0.5 * ( tgrmin + tgrmax )
c
      endif
c
      if ( ntg .eq. 0 ) then
c
         vrbar  = 0.0
         tgrbar = 0.0
c
      endif
c
      do ir = nr1 + nguard , nre - nguard
c
         if ( tgsnr(ir) .gt. snrmin ) then
c
            call addtgt ( curtime , tgr(ir) ,
     .                    2.0 * tgvel(ir) / lambda , tgacc(ir) ,
     .                    tgdw(ir) ,  tgsnr(ir) , 't' )
c
         endif 
c
      enddo
c
      if ( quiet .gt. 1 ) write ( 7  , * ) 'Second pass targets: ' , ntg
c
c-------------------  Sort target detections  --------------------------
c***********************************************************************
c***********************************************************************
c
c---------------------- SLIP VELOCITY ESTIMATION -----------------------
c
c   To correct for inconsistencies between the doppler velocity and the
c   actual velocity for some systems such as the NRL Advanced Profile
c   radar, it is necessary to estimate the 'slip' velocity.
c
c   Since this is not expected to be necessary on most systems, this
c   code is isolated here.  Also, it is assumed that the doppler target
c   detector above is the primary mode - thus, the slip velocity is
c   estimated only after the final set of targets has been selected.
c
      eslip = 0.0
c
      if ( ntg .gt. 0 ) then
c
         do 1000 ir = nr1 + nguard , nre - nguard
c
         if ( tgsnr(ir) .le. snrmin ) go to 1000
c
c   Six skill scores are computed.  These represent the scores for the
c   first and second halves of the time series for the main path, the
c   path at minus one range cell and the path at plus one range cell.
c
c                  (1,3)    |    (2,3)     ^
c                                          |
c                  (1,2)    |    (2,2)     r
c                                          |
c                  (1,1)    |    (2,1)     |
c
c                       --- t --->
c
            skill(:,:) = cmplx( 0.0 , 0.0 )
c
            do k = 1 , 3
c
               vel   = tgvel(ir)
c
c   r and rbar are range from the first range cell, not MCP
c
               rbar  = 0.5 * rswath + tgr(ir) + dr * float( k - 2 )
c
c   Form time series of range-compressed signal history along all
c   alias paths having the estimated Doppler frequency
c
               delr  = vel * dt
               r     = rbar - 0.5 * delr * float( ntot )
c
c   Assume data is periodic in range
c
               if ( r .lt. 0.0 )    r = r + rswath
               if ( r .ge. rswath ) r = r - rswath
c
               do ip = np1 , npe
c
                  irp  = 1 + ifix( r / dr )
                  irp1 = irp + 1
                  if ( irp1 .gt. nr ) irp1 = irp1 - nr
c
                  del    = ( r / dr ) - float( irp - 1 )
c
                  cw(ip) = del * c(ip,irp1) + ( 1.0 - del ) * c(ip,irp)
c
                  r      = r + delr
c
                  if ( r .lt. 0.0 )    r = r + rswath
                  if ( r .ge. rswath ) r = r - rswath
c
               enddo
c
c   Compute skill scores for first and second half of time series for
c   the k-th range cell
c
               nmid = ( npe + np1 - 1 ) / 2
               do ip = np1 , npe - 1
c
                  if ( ip .lt. nmid ) then
c
                     skill(1,k) = skill(1,k) +
     .                            cw(ip) * conjg( cw(ip+1) )
c
                  else if ( ip .gt. nmid ) then
c
                     skill(2,k) = skill(2,k) +
     .                            cw(ip) * conjg( cw(ip+1) )
c
                  endif
c
               enddo
c
            enddo
c
            if ( cabs( skill(1,2) ) .gt. cabs( skill(1,1) ) .and.
     .           cabs( skill(1,2) ) .gt. cabs( skill(1,3) ) ) then
               dminus = 0.5 * dr *
     .                 ( cabs( skill(1,3) ) - cabs( skill(1,1) ) ) /
     .                 ( 2.0 * cabs( skill(1,2) ) - cabs( skill(1,1) )
     .                                            - cabs( skill(1,3) ) )
            else if ( cabs( skill(1,1) ) .gt. cabs( skill(1,3) ) ) then
               dminus = - 0.5 * dr
            else if ( cabs( skill(1,1) ) .lt. cabs( skill(1,3) ) ) then
               dminus = + 0.5 * dr
            else
               dminus = 0.0
            endif
c
            if ( cabs( skill(2,2) ) .gt. cabs( skill(2,1) ) .and.
     .           cabs( skill(2,2) ) .gt. cabs( skill(2,3) ) ) then
               dplus = 0.5 * dr *
     .                 ( cabs( skill(2,3) ) - cabs( skill(2,1) ) ) /
     .                 ( 2.0 * cabs( skill(2,2) ) - cabs( skill(2,1) )
     .                                            - cabs( skill(2,3) ) )
            else if ( cabs( skill(2,1) ) .gt. cabs( skill(2,3) ) ) then
               dplus = - 0.5 * dr
            else if ( cabs( skill(2,1) ) .lt. cabs( skill(2,3) ) ) then
               dplus = + 0.5 * dr
            else
               dplus = 0.0
            endif
c
            eslip = eslip + ( dplus - dminus ) /
     .                      ( float( ntot / 2 ) * dt )
c
 1000    continue
c
         eslip = eslip / float( ntg )
c
      endif
c
c---------------------- SLIP VELOCITY ESTIMATION -----------------------
c
c***********************************************************************
c
      return
      end
C**
C***********************************************************************
C**
      subroutine alias0 ( nr , tgsnr , tgvel , tgacc , tgr , abar ,
     .                    rwindo , vwindo , nfirst , nlast , snrmin ,
     .                    snrmax )
C**
C***********************************************************************
C**
c  This routine eliminates alias detections by recursively removing the
c  target with velocity furthest from the snr-weighted mean velocity
c
      implicit none
c
      integer nr , ir , ivmin , ivmax , irmin , irmax , nfirst , nlast
c
      real    tgsnr(nr) , tgvel(nr) , tgacc(nr) , tgr(nr) , abar ,
     .        rwindo , vwindo , vbar , qbar , vmin , vmax , rbar ,
     .        rmin , rmax , snrmin , snrmax
c
      logical done
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c             Stage 1 - Thin to eliminate Doppler outliers
c
      done = .false.
c
      do while ( .not. done )
c
         vbar  = 0.0
         qbar  = 0.0
         vmin  = + 1000.0 * vwindo
         vmax  = - 1000.0 * vwindo
         ivmin = 0
         ivmax = 0
c
         do ir = nfirst , nlast
c
            if ( tgsnr(ir) .gt. snrmin ) then
c
               if ( tgvel(ir) .le. vmin ) then
                  vmin  = tgvel(ir)
                  ivmin = ir
               endif
c
               if ( tgvel(ir) .ge. vmax ) then
                  vmax  = tgvel(ir)
                  ivmax = ir
               endif
c
               vbar   = vbar + tgvel(ir) * tgsnr(ir)
               qbar   = qbar + tgsnr(ir)
c
            endif
c
         enddo
c
         if ( qbar .gt. snrmin ) then
c
            vbar = vbar / qbar
c
            if ( ( vmax - vmin ) .ge. vwindo ) then
c
c   Eliminate the point furthest from the mean
c
               if ( ( vmax - vbar ) .gt. ( vbar - vmin ) ) then
                  tgsnr(ivmax) = - abs( tgsnr(ivmax) )
               else
                  tgsnr(ivmin) = - abs( tgsnr(ivmin) )
               endif
c
            else
c
               done = .true.
c
            endif
c
         else
c
            done = .true.
c
         endif
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c             Stage 2 - Thin to eliminate range outliers
c
      done = .false.
c
      do while ( .not. done )
c
         vbar   = 0.0
         rbar   = 0.0
         abar   = 0.0
         qbar   = 0.0
         snrmax = 0.0
         rmin   = + 1000.0 * rwindo
         rmax   = - 1000.0 * rwindo
         irmin  = 0
         irmax  = 0
c
         do ir = nfirst , nlast
c
            if ( tgsnr(ir) .gt. snrmin ) then
c
               if ( tgr(ir) .le. rmin ) then
                  rmin  = tgr(ir)
                  irmin = ir
               endif
c
               if ( tgr(ir) .ge. rmax ) then
                  rmax  = tgr(ir)
                  irmax = ir
               endif
c
               vbar   = vbar + tgvel(ir) * tgsnr(ir)
               rbar   = rbar + tgr(ir) * tgsnr(ir)
               abar   = abar + tgacc(ir) * tgsnr(ir)
               qbar   = qbar + tgsnr(ir)
c
               snrmax = amax1( snrmax , tgsnr(ir) )
c
            endif
c
         enddo
c
         if ( qbar .gt. snrmin ) then
c
            vbar = vbar / qbar
            rbar = rbar / qbar
            abar = abar / qbar
c
            if ( ( rmax - rmin ) .gt. rwindo ) then
c
c   Eliminate the point furthest from the mean
c
               if ( ( rmax - rbar ) .gt. ( rbar - rmin ) ) then
                  tgsnr(irmax) = - abs( tgsnr(irmax) )
               else
                  tgsnr(irmin) = - abs( tgsnr(irmin) )
               endif
c
            else
c
               done = .true.
c
            endif
c
         else
c
            done = .true.
c
         endif
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      real function noise0 ( pwr , n_sort , nr , nperct )
C**
C***********************************************************************
C**
c
c     Extract the "nperct" percentile value, return as noise power
c     level.
c
      implicit none
c
      integer  index , nr , n_sort(nr)
c
      real     pwr(nr) , nperct
c
      index = nint( nperct * float( nr ) / 100.0 )
c
      if ( index .eq. 0 ) index = 1
c
      call sort ( pwr , nr , n_sort )
c
      noise0 = pwr(index)
      noise0 = amax1( noise0 , 1.0e-8 * pwr(nr) )
c
      return
      end
