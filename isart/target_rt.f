C**
C***********************************************************************
C**
      subroutine target_rt ( c , nr , np , nr1 , nre , np1 , npe , dr ,
     .                       dt , cw , lambda , nalias , ntg , sg0 ,
     .                       d0 , dd , a , b , dc , dw , tgsnr , tgdw ,
     .                       tgvel , tgrvel , tgacc , tgr , tgdadt ,
     .                       vrbar , tgrmin , tgrmax , tgrbar ,
     .                       curtime , snoise , vr0 , dr0 , snrmin ,
     .                       athrsh , notch , vnotch , dvntch , rnotch ,
     .                       drntch , rwindo , vwindo , eslip )
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
      integer nr , np , nalias , ntg , ir , ip , irp , np1 , npe , nr1 ,
     .        nre , ntot , nguard , notch , nend , nfirst , ntmid
c
      complex c(np,nr) , cw(np)
c
      real    tgsnr(nr) , tgdw(nr) , tgvel(nr) , tgrvel(nr) ,
     .        tgacc(nr) , tgr(nr) , tgdadt(nr)
c
      real    b(nr) , rvel , sg0(nr) , d0(nr) , a(nr) , dc(nr) ,
     .        dw(nr) , dd(nr) , dt , lambda , athrsh , abar , dr ,
     .        valias , rbar , pi , rswath , snr , vrbar ,  tgrbar ,
     .        tgrmin , tgrmax , ddw , n0 ,  snoise , snrmin , rwloss ,
     .        snrmax , vnotch , dvntch , rnotch , drntch , vr0 , dr0 ,
     .        rwindo , vwindo , n0mid , pmid , eslip , curtime , snrlim
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
         call dchirp_rt ( cw(nfirst) , ntmid , b(ir) , rvel , sg0(ir) ,
     .                    d0(ir) , dd(ir) , dt , 0 , a(ir) , lambda ,
     .                    dc(ir) , dw(ir) )
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
      if ( ntg .eq. 0 ) then
c
         vrbar  = 0.0
         tgrbar = 0.0
c
         return
c
      endif
c
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
c     write ( 6  ,'(/,10x,a39)')
c    .                        'Candidate Targets (i,snr,v,a,dw,r,dadt)'
c     write ( 6  ,'(10x,a18,f16.3)' ) 'Current time =    ' , curtime
c     write ( 7  ,'(/,10x,a39)')
c    .                        'Candidate Targets (i,snr,v,a,dw,r,dadt)'
c     write ( 7  ,'(10x,a18,f16.3)' ) 'Current time =    ' , curtime
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
c           write ( 6 ,'(i10,6f10.3)' ) ntg , db( tgsnr(ir) ) ,
c    .                   tgvel(ir) , tgacc(ir) , tgdw(ir) , rbar ,
c    .                   tgdadt(ir)
c           write ( 7 ,'(i10,6f10.3)' ) ntg , db( tgsnr(ir) ) ,
c    .                   tgvel(ir) , tgacc(ir) , tgdw(ir) , rbar ,
c    .                   tgdadt(ir)
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
c***********************************************************************
c
      eslip = 0.0
c
      return
      end
