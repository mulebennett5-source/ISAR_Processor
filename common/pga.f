C**
C***********************************************************************
C**
      subroutine pga ( cim , nf , nr , nrcell , nrfcell , worknf ,
     .                 wt , cx , x , cpp , workft , pg , denom , ph ,
     .                 icount , dtf , corg , cnew , aveamp , nptgt ,
     .                 pnoise , pgatgt , rscore , dotdot , mprat ,
     .                 mrrat , nafill , actual , dfmp , time0 , rcen ,
     .                 fcen , nfocus )
C**
C***********************************************************************
C**
c   Phase Gradient Autofocus: Removes phase and amplitude error from 
c   complex image cim which is common to all range lines.  Estimates
c   phase gradient pg(t) from complex image cim using nrkeep of the
c   range lines with the brightest targets.  The phase error is removed
c   from cim through complex multiplication by
c
c                  exp{ -i * integral[ pg(t) * dt ] }.
c
c   The total phase error removed from cim, possibly over multiple
c   iterations, is returned in ph(t).  Also computes track file
c   (centroid, Doppler width, total pwr, acceleration, SNR) for each
c   range line.
c
c  Inputs...
c
c    integer  nf           Number of freq bins in image
c    integer  nr           Number of range lines in image
c    complex  cim(nf,nr,1) Phase-corrupted complex image
c    complex  cim(nf,nr,2) Work array for effective signal history
c    integer  npass        Freq passband (odd # bins) about DC
c    integer  ntaper       Freq taper length (bins); npass + 2 * ntaper
c                                                   <= nf
c    integer  maxit        Max iterations for convergence
c    real     dtf          Sample rate in time domain
c    real     pnoise       Noise level (percentile of range cells)
c    real     snrmin       SNR detection threshold
c
c  Work Arrays...
c
c    real     worknf(nf)   Real work array of length nf
c    real     wt(nf)       Freq domain window
c    complex  cx(nf)       data array; equivalenced with x for FOURT
c    real     x(2,nf)
c    complex  cpp(nf)      Phase array
c    real     workft(2,nf) used by FOURT
c    real     pg(nf)       Phase gradient
c    real     denom(nf)    Denominator
c    real     aveamp(nf)   Average amplitude in time domain
c
c  Outputs...
c
c    complex  cim(nf,nr,1) Phase-corrected image
c    complex  cim(nf,nr,2) Phase-corrected image in real space
c    real     ph(nf)       Total phase (rad) removed from cim
c                          (FOURT storage)
c    integer  icount       Number of iterations completed
c    real     cnew         Contrast of phase-corrected image
c    real     pgatgt(nr,5) Track file: centroid, Doppler width,
c                          total pwr, accel, snr vs range line
c    integer  nrkeep       Max. no. range lines to be used
c    integer  nptgt        Actual no. of lines used in pg estimate
c    integer  nrcell(nr)   Range line with smallest to largest SNR
c    integer  nrfcell(nr)  Freq. bin with largest SNR for range line
C**
C***********************************************************************
C**
      implicit none
c
      integer nn(1) , nf , nr , nrcell(nr) , nrfcell(nr) , j , k , m ,
     .        n , icount , nptgt , nfocus
c
      complex cim(nf,nr,2+2*nfocus) , cx(nf) , cpp(nf)
c
      real    wt(nf) , x(2,nf) , pg(nf) , denom(nf) , ph(nf) , db ,
     .        worknf(nf) , dtf , ctrast , corg , cold , cnew , ctc ,
     .        aveamp(nf) , pnoise , pgatgt(nr,5) , rcen , fcen ,
     .        acc_sd , workft(2,nf) , rscore(nr) , ct(1+2*nfocus) ,
     .        ac_use , ar_use
c
      integer mprat , mrrat , nafill , actual
c
      real    dotdot(1+nf/mprat,nr/mrrat,actual) , dfmp , time0
c
      logical all
c
c-----------------------------------------------------------------------
c
c   Include files
c
      include     'sarprm.h'      !  Standard SAR parameters
c
      include     'updates.h'     !  Updated SAR parameters
c
      nn(1) = nf
c
      ph(:) = 0.0
c
c   Calculate initial contrast
c
      corg = ctrast( cim , nf , nr )
      cnew = corg
      cold = corg
c
      if ( quiet .gt. 1 ) write ( 6 , '(/,1x,a40,f10.2)' )
     .   ' Contrast into PGA:                     ' , corg
c
      write ( 7 , '(/,1x,a40,f10.2)' )
     .   ' Contrast into PGA:                     ' , corg
c
c   Convert back to effective signal history
c
      do j = 1 , nr
c
         call cshift ( cim(1,j,1) , nf , nf / 2 , cx )
c
         call fourt ( x , nn , 1 , - 1 , 1 , workft , 2 * nf )
c
         call cshift ( cx , nf , nf / 2 , cim(1,j,2) )
c
         cim(:,j,2) = cim(:,j,2) / float( nf )
c
      enddo
c
c   Select targets
c
      call gettgt ( cim , nf , nr , pnoise , pgatgt , nrcell , dtf ,
     .              cx , nptgt , rscore , dotdot , mprat , mrrat ,
     .              nafill , actual , dfmp , time0 , rcen , fcen ,
     .              acc_sd )
c
      write ( 6 , '(1x,a,i8,f10.3)' ) 'nptgt, acc_sd: ' , nptgt , acc_sd
c
c   Second part - remove range independent phase errors
c
      if ( nptgt .eq. 0 ) then
c
         if ( quiet .gt. 1 ) then
c
            write ( 6 , * ) 'No lines passed SNR threshold in PGA'
            write ( 6 , * )
c
         endif
c
         write ( 7 , * ) 'No lines passed SNR threshold in PGA'
         write ( 7 , * )
c
      else
c
c   Estimate phase gradient
c
         if ( quiet .gt. 1 ) then
c
            write ( 6 , '(1x,i6,a34)' ) nptgt ,
     .                             ' Lines passed SNR threshold in PGA'
            write ( 6 , * )
c
         endif
c
         write ( 7 , '(1x,i6,a34)' ) nptgt ,
     .                             ' Lines passed SNR threshold in PGA'
         write ( 7 , * )
c
c   Do at least one interation of PGA - then quit if minimum
c   contrast change is attained
c
         do icount = 1 , maxit
c
            call estpg ( cim , nf , nr , nrcell , nrfcell ,
     .                   nptgt , wt , npass , ntaper , cx ,
     .                   x , cpp , workft , pg , denom ,
     .                   worknf , dtf , aveamp )
c
c   Integrate phase gradient
c
            call integ ( pg , nf , dtf  , worknf )
c
c   Update total phase array
c
            pg(:) = worknf(:)
c
            ph(:) = ph(:) + pg(:)
c
c   Remove phase error
c
c   Subroutine remove will only re-compute all image lines at the last
c   of the iterations.  Otherwise, it re-computes the phase for all
c   lines and re-computes the image for only those range lines that
c   contain targets.
c
            all = ( icount .eq. maxit )
c 
            call remove ( cim , nf , nr , pg  , cx , x , cpp ,
     .                    workft , nrcell , nptgt , all , 2 , 1 )
c
            if ( all .and. ( nfocus .gt. 0 ) ) then
c
               do j = 1 , nr
c
                  cim(:,j,2+2*nfocus) = cim(:,j,2)  !  Nominal S.H.
c
                  cim(:,j,1+nfocus)   = cim(:,j,1)  !  Nominal image
c
               enddo
c
c   Fill all the other focus images with nominal signal history
c
               do k = 1 , nfocus
c
                  do j = 1 , nr
c
                     cim(:,j,k)          = cim(:,j,2+2*nfocus)
c
                     cim(:,j,k+1+nfocus) = cim(:,j,2+2*nfocus)
c
                  enddo
c
               enddo
c
c   Refocus the images
c
               do k = 1 , nfocus
c
                  ac_use = float( k ) * accorr
c
                  ar_use = float( k ) * arcorr
c
                  n      = 1 + nfocus + k
c
                  m      = 1 + nfocus - k
c
                  call remove_r ( cim , nf , nr , pg  , cx , x , cpp ,
     .                            workft , m , m , drf , dtf ,
     .                            - ac_use , - ar_use )
c
                  call remove_r ( cim , nf , nr , pg  , cx , x , cpp ,
     .                            workft , n , n , drf , dtf ,
     .                            ac_use , ar_use )
c
               enddo
c
            endif
c
         enddo  !  Iterations
c
         do k = 1 , 1 + 2 * nfocus
c
            ct(k) = ctrast( cim(1,1,k) , nf , nr )
c
         enddo
c
         ctc = ct(1+nfocus)
c
c-----------------------------------------------------------------------
c
      endif     !  Lines passed selection criteria
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Compute final contrast
c
      cnew = ctrast( cim(1,1,1+nfocus) , nf , nr )
c
      write ( 6 , '(/,1x,a40,f10.2)' )
     .   ' Contrast after PGA:                    ' , cnew 
c
      write ( 7 , '(/,1x,a40,f10.2)' )
     .   ' Contrast after PGA:                    ' , cnew
c
      return
      end
C**
C***********************************************************************
C**
      subroutine pquad ( pg , quad , nf , dtf )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , i
c
      real    pg(nf) , quad , dtf , t , pi
c
      pi = atan2( 0.0 , - 1.0 )
c
      do i = 1 , nf
c
         t     = float( i - 1 - nf / 2 ) * dtf
c
         pg(i) = pi * quad * t ** 2
c
      enddo
c
      return
      end

C**
C***********************************************************************
C**
      subroutine estpg ( cim , nf , nr , nrcell , nrfcell , nptgt, wt ,
     .                   npass , ntaper , cx , x , cxdot , work , pg ,
     .                   denom , worknf , dtf , aveamp )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , nr , nrcell(nr) , nrfcell(nr) , nptgt ,
     .        npass , ntaper , i
c
      complex cim(nf,nr) , cx(nf) , cxdot(nf)
c
      real    wt(nf) , x(2,nf) , work(2,nf) , pg(nf) , denom(nf) ,
     .        worknf(nf) , dtf , aveamp(nf)
c
c   Estimates phase gradient pg(t) from complex image cim using nptgt
c   of the range lines with the brightest targets.
c
c   For each target range line, select freq cell with largest intensity.
c
      call select ( cim , nf , nr , nrcell , nrfcell , nptgt )
c
c   Compute window centered at bin 1.
c
      call hwindo ( wt , npass + 2 * ntaper , ntaper , nf , worknf )
c
c   From nptgt of the range lines with the brightest targets, estimate
c   the phase gradient pg(t) and average amplitude aveamp(t).
c
      aveamp(:) = 0.0
      pg(:)     = 0.0
      denom(:)  = 0.0
c
      do i = 1 , nptgt
c
         call rshift ( wt , nf , nrfcell(nrcell(nr-i+1))-1 , worknf )
c
         call tfns ( cim , nf , nr , nrcell(nr-i+1) ,
     .               nrfcell(nrcell(nr-i+1)) , cx , x , cxdot ,
     .               worknf , work , dtf )
c
         aveamp(:) = aveamp(:) + abs( cx(:) )
c
         pg(:)     = pg(:) + imag( conjg( cx(:) ) * cxdot(:) )
c
         denom(:)  = denom(:) + cx(:) * conjg( cx(:) )
c
      enddo
c
      aveamp(:) = aveamp(:) / float( nptgt )
      pg(:)     = pg(:) / denom(:)
c
      return
      end
C**
C***********************************************************************
C**
      subroutine select ( cim , nf , nr , nrcell , nrfcell , nptgt )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , nr , nrcell(nr) , nrfcell(nr) , i , j , nptgt ,
     .        index
c
      complex cim(nf,nr)
c
      real    old , new
c
c   For each target range line, select freq cell with largest intensity.
c   Link from range line to largest freq cell is nrfcell.
c
      do j = 1 , nptgt
c
         old            = 0.0
         index          = nrcell( nr - j + 1 )
         nrfcell(index) = 1
c
         do i = 1 , nf
c
            new = cim(i,index) * conjg( cim(i,index) )
c
            if ( new .gt. old ) then
c
               old            = new
               nrfcell(index) = i
c
            endif
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine tfns ( cim , nf , nr , nrcell , nfcell , cx , x ,
     .                  cwork , wt , work , dtf )
C**
C***********************************************************************
C**
      implicit none
c
      integer nn(1) , nf , nr , nrcell , nfcell , i , idist , iplus ,
     .        iminus , nwrap
c
      complex cim(nf,nr) , cx(nf) , cwork(nf)
c
      real    x(2,nf) , work(2,nf) , wt(nf) , dtf , cpeak , aplus ,
     .        aminus , delta , lpeak
c
      nn(1) = nf
c
c   Computes time functions x(t) and xdot(t) for shifted, windowed target
c   for selected range line.
c
      cx(:) = wt(:) * cim(:,nrcell)
c
c   Shift windowed target to DC
c
      call cshift ( cx , nf , - ( nfcell - 1 ) , cwork )
c
c   Zero out all points which wrapped around
c
      nwrap = min( nfcell , nf + 1 - nfcell )
c
      do i = 1 + nwrap , 1 + nf / 2
c
         cwork(i)      = cmplx( 0.0 , 0.0 )
         cwork(nf+2-i) = cmplx( 0.0 , 0.0 )
c
      enddo
c
c   Linearly interpolate to peak frequency
c
c   Compute delta to interpolate to peak
c
      if ( cwork(1) .eq. cwork(nf) .and.
     .     cwork(1) .eq. cwork(2) ) then
c
         delta = 0.0
c
      else
c
         delta = 0.5 * ( cabs( cwork(2) ) - cabs( cwork(nf) ) )
     .         / ( 2.0 * cabs( cwork(1) ) - cabs( cwork(nf) )
     .                                    - cabs( cwork(2) ) )
c
      endif
c
      do i = 1 , nf
c
         if ( delta .ge. 0.0 ) then
c
            iplus  = i + 1
            if ( iplus .gt. nf ) iplus  = iplus - nf
            cx(i) = delta * cwork(iplus) + ( 1.0 - delta ) * cwork(i)
c
         else
c
            iminus = i - 1
            if ( iminus .lt. 1 ) iminus = iminus + nf
            cx(i) = - delta * cwork(iminus) + ( 1.0 + delta ) * cwork(i)
c
         endif
c
      enddo
c
c   Apply target-centered window to x, clipping peaks which are greater
c   than three cells away and greater than 0.1 times the peak
c
      cpeak = cabs( cx(1) )
c
      do i = 1 , nf
c
         idist = min( abs( i - 1 ) , abs( i - 1 - nf ) )
c
         if ( idist .gt. 1 ) then
c
c   The mask function is 1.0, 1.0, 0.9, 0.45, 0.3, 0.225, 0.18, 0.15, ..
c
            lpeak = amax1( 0.1 , ( 0.9 / float( idist - 1 ) ) ) * cpeak
c
            if ( cabs( cx(i) ) .gt. lpeak ) then
c
               cx(i) = cx(i) * ( lpeak / cabs( cx(i) ) )
c
c   Under-relaxation: take 90% of true value
c
            else
c
               cx(i) = 0.9 * cx(i)
c
            endif
c
         endif
c
      enddo
c
c   Use the minimum of the left or right side of the spectrum to
c   eliminate effects of real scatterers which appear on only one side
c
      do i = 2 , ( nf / 2 ) - 1
c
         aplus  = cabs( cx(i) )
         aminus = cabs( cx(nf+2-i) )
c
         if ( aplus .gt. aminus ) then
c
            cx(i)      = amax1( 0.25 , aminus / aplus ) * cx(i)
c
         else if ( aminus .gt. aplus ) then
c
            cx(nf+2-i) = amax1( 0.25 , aplus / aminus ) * cx(nf+2-i)
c
         endif
c
      enddo
c
c   IFFT x and normalize
c
      call fourt ( x , nn , 1 , - 1 , 1 , work , 2 * nf )
c
      cwork(:) = cx(:) / float( nf )
c
c   Reorder x to increasing time
c
      call cshift ( cwork , nf , nf / 2 , cx )
c
c   Compute time derivative of x
c
      do i = 2 , nf - 1
c
         cwork(i) = ( cx(i+1) - cx(i-1) ) / ( 2.0 * dtf )
c
      enddo
c
      cwork(1)  = ( cx(2) - cx(1) ) / dtf
      cwork(nf) = ( cx(nf) - cx(nf-1) ) / dtf
c
      return
      end
C**
C***********************************************************************
C**
      subroutine integ ( ydot , n , dt , y )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i
c
      real    ydot(n) , dt , y(n)
c
c   Set y(t) = integral( ydot(t) * dt ) using Trapezoidal Rule
c
      y(1) = 0.0
c
      do i = 2 , n
c
         y(i) = y(i-1) + 0.5 * ( ydot(i-1) + ydot(i) ) * dt
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine remove ( cim , nf , nr , ph , cx , x , cp , work ,
     .                    nrcell , nptgt , all , sh , im )
C**
C***********************************************************************
C**
      implicit none
c
      integer nn(1) , nf , nr , j , nrcell(nr) , nptgt , jtgt , sh , im
c
      complex cim(nf,nr,sh) , cx(nf) , cp(nf)
c
      real    ph(nf) , x(2,nf) , work(2,nf)
c
      logical all
c
      nn(1) = nf
c
c   Remove phase error exp( -i * ph(t) ) from image by complex
c   multipication in time domain.
c
c   Compute phase correction
c
      cp(:) = cmplx( cos( ph(:) ) , - sin( ph(:) ) )
c
c   Remove phase error over all range lines
c 
      do j = 1 , nr
c
         cim(:,j,sh) = cim(:,j,sh) * cp(:)
c
      enddo
c
c   Restore image
c
      if ( all ) then
c
c   Do all range lines
c
         do j = 1 , nr
c
            call cshift ( cim(1,j,sh) , nf , nf / 2 , cx )
c
            call fourt ( x , nn , 1 , + 1 , 1 , work , 2 * nf )
c
            call cshift ( cx , nf , nf / 2 , cim(1,j,im) )
c
         enddo
c
      else
c
c   Do only those range lines containing targets
c
         do jtgt = 1 , nptgt
c
            j = nrcell(nr+1-jtgt)
c
            call cshift ( cim(1,j,sh) , nf , nf / 2 , cx )
c
            call fourt ( x , nn , 1 , + 1 , 1 , work , 2 * nf )
c
            call cshift ( cx , nf , nf / 2 , cim(1,j,im) )
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine remove_r ( cim , nf , nr , ph , cx , x , cp , work ,
     .                      sh , im , drf , dtf , accorr , arcorr )
C**
C***********************************************************************
C**
      implicit none
c
      integer nn(1) , nf , nr , j , sh , im
c
      complex cim(nf,nr,sh) , cx(nf) , cp(nf)
c
      real    ph(nf) , x(2,nf) , work(2,nf) , drf , accorr , arcorr ,
     .        dtf , quad
c
      nn(1) = nf
c
c   Remove phase error exp( -i * ph(t) ) from image by complex
c   multipication in time domain.
c
c   Remove phase error over all range lines
c 
      do j = 1 , nr
c
c   Compute phase correction
c
         quad        = accorr + arcorr * drf * float( j - 1 - nr / 2 )
c
         call pquad ( ph , quad , nf , dtf )
c
         cp(:)       = cmplx( cos( ph(:) ) , - sin( ph(:) ) )
c
         cim(:,j,sh) = cim(:,j,sh) * cp(:)
c
      enddo
c
c   Restore image
c
      do j = 1 , nr
c
         call cshift ( cim(1,j,sh) , nf , nf / 2 , cx )
c
         call fourt ( x , nn , 1 , + 1 , 1 , work , 2 * nf )
c
         call cshift ( cx , nf , nf / 2 , cim(1,j,im) )
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine gettgt ( cim , nf , nr , pnoise , pgatgt , nrcell ,
     .                    dtf , work , nptgt , rscore , dotdot , mprat ,
     .                    mrrat , nafill , actual , dfmp , time0 ,
     .                    rcen , fcen , acc_sd )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , nr , i , j , nrcell(nr) , index , nptgt , offset
c
      complex cim(nf,nr,2) , work(nf)
c
      real    pnoise , pgatgt(nr,5) , dtf , vel , d0 , dd , aavg ,
     .        npower , ddw , rscore(nr) , peak , fcon , dtdf , tgtime ,
     .        tgfreq , rcen , fcen , acc_bar , acc_sd
c
      integer mprat , mrrat , nafill , actual , irm , itgt , itg ,
     .        ileft , ileftp , ioff , ncent , ntot , nbest
c
      real    dotdot(1+nf/mprat,nr/mrrat,actual) , dfmp , accorg ,
     .        time0 , db , rtgt , f , delta , athrsh , dts
c
c-----------------------------------------------------------------------
c
c   Variables for Sidelobe Cancellation
c
c   Added by KAM 6-29-98,Mod. 7-30-98
c
      real    dopscore(nr) , rdscore(nr) , strong_targs_dop
c
      integer target_index , k , dopcell(nr) , rdcell(nr) ,
     .        slbkilled
c
c-----------------------------------------------------------------------
c
c   Include files
c
      include     'sarprm.h'      !  Standard SAR parameters
c
      include     'updates.h'     !  Standard SAR parameters
c
c   Get Targets: Selects targets and computes track file;
c   pgatgt: centroid, width, pwr, accel, snr vs range line
c
c   Loop over range lines and compute track file
c
c   Form a time series over the middle of the integration period.  The
c   length of the time series is half the integration time.  It is also
c   chosen to be a multiple of 4 so that dchirp will break it cleanly
c   into 4 disjoint sections.
c
      ntot   = int( 3 * ( nf * actual ) / nafill )/ 4 ! 3/4 actual inte-
c                                                     ! gration time
c
      ntot   = 4 * ( max( 1 , ntot / 4 ) )            ! Multiple of 4,
c                                                     ! at least 4
c
      ioff   = ( nf - ntot ) / 2                      ! Offset to start
c                                                     ! of middle ntot
c
      ddw    = 0.2 * dff * ( float( nafill ) / float( actual ) )
      dts    = float( ntot / 4 ) * dtf
      athrsh = 5.0 / ( dts ** 2 )
c
      do j = 1 , nr
c
         do i = 1 , ntot
c
            work(i) = cim(i+ioff,j,2)
c
         enddo
c
         call dchirp ( work , ntot , pgatgt(j,4) , vel ,
     .                 pgatgt(j,3) , d0 , dd , dtf , - 1 , aavg ,
     .                 lambda , pgatgt(j,1) , pgatgt(j,2) )
c
c   Downweight range cells with large Doppler width
c
c        if ( reduce .eq. 2 .or. reduce .eq. 3 ) then
c
c           do i = 1 , nf
c              cim(i,j,2) = cim(i,j,2) / ( 0.01 + dd )
c           enddo
c
c        endif
c
      enddo
c
      do  j = 1 , nr
c
         rscore(j)   = pgatgt(j,3)
         pgatgt(j,5) = pgatgt(j,3) / ( pgatgt(j,2) + ddw )
c
      enddo
c
      call sort ( rscore , nr , nrcell )
      index = nint( 0.01 * pnoise * float( nr ) )
      if ( index .eq. 0 ) index = 1
c
c   Bullet-proof against very small or zero noise values by limiting it
c   to a small fraction of the peak cross section or, if the peak is
c   zero, set it to one.
c 
      npower = amax1( rscore(index) , 1.0E-5 * rscore(nr) )
      if ( npower .eq. 0.0 ) npower = 1.0
c
c   Scale SNR by noise level and remove large acceleration, large
c   Doppler width points and points more than 40 dB down from the
c   peak value
c
      peak = rscore(nr) / ( npower * 2.0 * dtf )
c
      if ( mode .eq. 1 ) then
c
         fcon = 0.15
c
      else
c
         fcon = 4.0
c
      endif
c
      do j = 1 , nr
c
         pgatgt(j,5) = pgatgt(j,5) / ( npower * 2.0 * dtf )
c
         if ( abs( pgatgt(j,4) ) .ge. athrsh .or.
     .        pgatgt(j,2) .gt. fcon * float( npass ) * dff .or.
     .        pgatgt(j,5) .lt. 1.0E-4 * peak .or.
     .        abs( pgatgt(j,1) ) .gt. 0.4 * float( nf ) * dff ) then
c
            pgatgt(j,5) = - abs( pgatgt(j,5) )
c
         endif
c
      enddo
c
c-----------------------------------------------------------------------
c
c    Sidelobe Cancellation
c    Sort the targets based upon rscore to find the strongest targets
c
c    - Added by KAM 6-29-98
c

c    Sort the targets based upon Doppler bin
c
      if ( slbkil .gt. 0 ) then
c
         slbkilled = 0
c
         do j = 1 , nr
c
            dopscore(j) = pgatgt(j,1)
c
            dopcell(j)  = j
c
         enddo
c      
         call sort ( dopscore , nr , dopcell )

         do  j = 1 , nr
c
            rdscore(j) = pgatgt(dopcell(j),5) 
            rdcell(j)  = j
c
         enddo
c      
         call sort ( rdscore , nr , rdcell )
c 
         k = 0
c
         do  j = 1 , 5
c
            target_index = rdcell(nr-j-k+1)
c    
            do while ( ( ( nr - j - k + 1 ) .gt. 0 ) .and. 
     .                 ( target_index .gt. 0 )       .and.
     .                 ( pgatgt(dopcell(target_index),5) .le. 0.0 ) )
c
               k            = k + 1 
c
               target_index = rdcell(nr-j-k+1)
c
            enddo
c
            if ( ( nr - j - k + 1 ) .lt. 1 ) go to 101 
c
            if ( target_index .gt. 0 ) then
c
               i = 1
c
               strong_targs_dop = dopscore(target_index) + 1.5 * dff
c
               do while ( ( target_index + i .le. nr ) .and. 
     .                    ( dopscore(target_index+i) .le. 
     .                      strong_targs_dop ) )
c
                  pgatgt(dopcell(target_index+i),5) 
     .                      = - abs( pgatgt(dopcell(target_index+i),5) )
c
                  slbkilled = slbkilled + 1
c
                  i         = i + 1
c
               enddo
c     
               strong_targs_dop = dopscore(target_index) - 1.5 * dff
c
               i = 1
c
               do while ( ( target_index - i .gt. 0 ) .and.
     .                    ( dopscore(target_index-i) .ge. 
     .                      strong_targs_dop ) )
c
                  pgatgt(dopcell(target_index-i),5) 
     .                  = - abs( pgatgt(dopcell(target_index-i),5) )
c
                  slbkilled = slbkilled + 1
c
                  i         = i + 1
c
               enddo
c          
            endif     !  target_index .gt. 0
c
         enddo        !  do  j = 1 , 5
c
         write ( 6 , * ) ' @ Sidelobes killed: ' , slbkilled
c
         write ( 7 , * ) ' @ Sidelobes killed: ' , slbkilled
c
      endif
c
 101  continue
c
c-----------------------------------------------------------------------
c
c   Mark targets with the same Dopplers as the strongest targets for 
c   deletion
c
c-----------------------------------------------------------------------
c
      nbest = 1
c
      call best ( pgatgt(1,5) , rscore , snrmin , nr , nrcell , nptgt ,
     .            nbest )
c
      acc_bar = 0.0
c
      acc_sd  = 0.0
c
      if ( nptgt .gt. 0 ) then
c
         if ( quiet .gt. 1 )
     .   write ( 6 ,'(/,10x,a32)') 'PGA Targets (i,snr,f,da,a0,dw,r)'
c
         write ( 7 ,'(/,10x,a32)') 'PGA Targets (i,snr,f,da,a0,dw,r)'
c
         if ( mode .eq. 2  .or. mode .eq. 3 )
     .        dtdf = - 0.5 * lambda * slant0 / ( vfocus ** 2 )
c
         do itgt = 1 , nptgt
c
c   PGA sorted the SNR values from smallest to largest; nrcell contains
c   the true range index of the targets
c
            itg = nrcell(nr+1-itgt)
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
c   frames can benefit from this focus information.  However, the
c   acceleration value must correct for the acceleration which has
c   already been applied.  This is estimated from the dotdot array
c   evaluated at the range and frequency at which the target was
c   observed.
c
            f      = ( pgatgt(itg,1) + 0.5 * dff * float( nf ) ) / dfmp
            ileft  = 1 + ifix( f )
            delta  = f  - float( ileft - 1 )
c
            ileft  = min( ileft , 1 + nf / mprat )
            ileft  = max( ileft , 1 )
            ileftp = ileft + 1
            ileftp = min( ileftp , 1 + nf / mprat )
            ileftp = max( ileftp , 1 )
            ncent  = 1 + actual / 2
c
            irm    = 1 + ( itg - 1 ) / mrrat
            rtgt   = rfmin + drf * float( itg - 1 )
            accorg = 0.5 * ( ( 1.0 - delta ) *
     .                             ( dotdot(ileft,irm,ncent) +
     .                               dotdot(ileft,irm,ncent-1) ) +
     .                       delta *
     .                             ( dotdot(ileftp,irm,ncent) +
     .                               dotdot(ileftp,irm,ncent-1) ) )
c
            if ( quiet .gt. 1 )
     .      write ( 6 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                pgatgt(itg,1) , pgatgt(itg,4) , accorg ,
     .                pgatgt(itg,2) , rtgt
c
            write ( 7 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                pgatgt(itg,1) , pgatgt(itg,4) , accorg ,
     .                pgatgt(itg,2) , rtgt
c
            if ( mode .eq. 1 ) then
c
               tgtime = time0
               tgfreq = pgatgt(itg,1)
c
            else
c
               tgtime = time0 + pgatgt(itg,1) * dtdf
               tgfreq = 0.0
c
            endif
c
            call addtgt ( tgtime , rtgt + rcen ,
     .                    tgfreq + fcen , accorg + pgatgt(itg,4) ,
     .                    pgatgt(itg,2) , pgatgt(itg,5) , 'p' )
c
            acc_bar = acc_bar + pgatgt(itg,4)
c
            acc_sd  = acc_sd  + pgatgt(itg,4) ** 2
c
         enddo
c
         acc_bar = acc_bar / float( nptgt )
c
         acc_sd  = sqrt( amax1( 0.0 ,
     .                   acc_sd / float( nptgt ) - acc_bar ** 2 ) )
c
      endif
c
c   If there are plenty of targets, use only the middle scores - this
c   reduces the impact of the atypical huge scatterers in estimating
c   the vibration errors.  Of the scatterers which need to be ignored
c   eliminate 1/8 at the top and 7/8 at the bottom.
c
      if ( nptgt .gt. nrkeep ) then
c
         offset = ( nptgt + 1 - nrkeep ) / 8
c
         do j = 1 , nrkeep
c
            nrcell(nr+1-j) = nrcell(nr+1-j-offset)
c
         enddo
c
         nptgt = nrkeep
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine best ( score , rscore , scmin , nr , nrcell , nptgt ,
     .                  nbest )
C**
C***********************************************************************
C**
c   Select points which are isolated maxima and above a threshold - then
c   sort them in ascending order
c
      implicit none
c
      integer nr , nrcell(nr) , j , k , nptgt , nbest
c
      real    score(nr) , rscore(nr) , scmin
c
c   Ignore the end-points
c
      do j = 1 , nbest
c
         score(j)      = - abs( score(j) )
         score(1+nr-j) = - abs( score(1+nr-j) )
c
      enddo
c
c   Use only points which are above the threshold and which are local
c   maxima
c
      do j = 1 + nbest , nr - nbest
c
         if ( abs( score(j) ) .le. scmin ) then
c
            score(j)  = - abs( score(j) )
c
         else
c
            do k = 1 , nbest
c
               if ( ( abs( score(j) ) .le. abs( score(j+k) ) ) .or.
     .              ( abs( score(j) ) .le. abs( score(j-k) ) ) )
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
