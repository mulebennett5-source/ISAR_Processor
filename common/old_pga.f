C**
C***********************************************************************
C**
      subroutine pga ( cim , nf , nr , nrcell , nrfcell , worknf ,
     .                 wt , cx , x , cpp , cpaccum , workft , pg ,
     .                 denom , maxitu , ph , icount , dtf , corg ,
     .                 cnew , taywt , aveamp , nptgt , pnoise ,
     .                 pgatgt , rscore , dotdot , mprat , mrrat ,
     .                 nafill , actual , dfmp , time0 , tint ,
     .                 rcen , fcen )
C**
C***********************************************************************
C**
c   Phase Gradient Autofocus: Removes phase and amplitude error from 
c   complex image cim which is common to all range liness.  Estimates
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
c    integer  nf           Number of freq bins in image (power of 2)
c    integer  nr           Number of range lines in image
c    complex  cim(nf,nr,1) Phase-corrupted complex image
c    integer  npass        Freq passband (odd # bins) about DC
c    integer  ntaper       Freq taper length (bins); npass + 2 * ntaper
c                                                   <= nf
c    real     percon       Contrast diff (percent) between interations
c                          for convergence
c    integer  ampcor       Control parameter for amplitude correction
c                          0 --> None
c                          1 --> Adaptive correction
c
c    integer  maxit        Max iterations for convergence
c    real     dtf          Sample rate in time domain
c    real     taywt        Sibelobe level in dB for Taylor weight
c    real     pnoise       Noise level (percentile of range cells)
c    real     snrmin       SNR detection threshold
c    real     tint         Image integration time (sec)
c
c  Work Arrays...
c
c    real     worknf(nf)   Real work array of length nf
c    real     wt(nf)       Freq domain window
c    complex  cx(nf)       data array; equivalenced with x for FOURT
c    real     x(2,nf)
c    complex  cpp(nf)      Phase array
c    complex  cpaccum(nf)  Accumulation Phase Array
c    real     workft(2,nf) used by FOURT
c    real     pg(nf)       Phase gradient
c    real     denom(nf)    Denominator
c    real     aveamp(nf)   Average amplitude in time domain
c
c  Outputs...
c
c    complex  cim(nf,nr,1) Amp/phase-corrected image
c    complex  cim(nf,nr,2) Amp/phase-corrected image in real space
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
      integer nn(1) , nf , nr , nrcell(nr) , nrfcell(nr) , i , j ,
     .        icount , nptgt , maxitu
c
      complex cim(nf,nr,3) , cx(nf) , cpp(nf), cpaccum(nf)
c
      real    wt(nf) , x(2,nf) , pg(nf) , denom(nf) , ph(nf) ,
     .        worknf(nf) , dtf , ctrast , corg , cold , cnew , taywt ,
     .        aveamp(nf) , pnoise , pgatgt(nr,5) , rcen ,
     .        fcen , workft(2,nf) , rscore(nr) , tint
c
      integer mprat , mrrat , nafill , actual , inner
c
      real    dotdot(1+nf/mprat,nr/mrrat,actual) , dfmp , time0
c
      logical all
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
c   Number of inner iterations over target range lines before the whole
c   image is re-computed and the convergence criterion is tested.
c
      integer     ninner
      parameter ( ninner = 3 )
c
      integer     ninner_use
c
      if ( rt_pga .eq. 0 ) then
c
         ninner_use = ninner
c
      else
c
         ninner_use = 10
c
      endif
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
         call fourt ( x , nn , 1 , - 1 , 1 , workft , nf )
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
     .              nafill , actual , dfmp , time0 , rcen , fcen )
c
      if ( maxitu .eq. 0 ) return
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
         do icount = 1 , maxitu
c
            if ( rt_pga .eq. 0 ) then
c
c   Store old image in case contrast decreases
c
               cim(:,:,3) = cim(:,:,1)
c
            endif
c
c   To save computer time and to allow the algorithm to overcome local
c   minima in contrast, only compute the full image and check the
c   convergence criteria every 'ninner' iterations.
c
            do inner = 1 , ninner_use
c
c   ** KAM 3/3/99 - Multiple types of PGA kernels per Jakowatz paper **
c
               if ( pgtype .lt. 2 ) then
c
                  call estpg ( cim , nf , nr , nrcell , nrfcell ,
     .                         nptgt , wt , npass , ntaper , cx ,
     .                         x , cpp , workft , pg , denom ,
     .                         worknf , dtf , aveamp )
c
c   Integrate phase gradient
c
                  call integ ( pg , nf , dtf  , worknf )
c
               else
c
                  call estpg_jako ( cim , nf , nr , nrcell , nrfcell ,
     .                              nptgt , wt , npass , ntaper , cx ,
     .                              x , cpp , cpaccum , workft , pg ,
     .                              denom , worknf , dtf , aveamp ,
     .                              pgtype )
c
                  worknf(1) = pg(1)
c
                  do i = 2, nf
c
                     worknf(i) = worknf(i-1) + pg(i) 
c
                  enddo
c
               endif
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
c   of the 'ninner' iterations.  Otherwise, it re-computes the phase
c   for all lines and re-computes the image for only those range lines
c   which contain targets.
c
               all = ( inner .eq. ninner_use ) 
               call remove ( cim , nf , nr , pg  , cx , x , cpp ,
     .                       workft , nrcell , nptgt , all )
c
            enddo
c
c   Compute new image contrast
c
            cnew = ctrast( cim , nf , nr )
c
            if ( ( rt_pga .eq. 0 ) .and. ( cold .gt. cnew ) .and.
     .           ( ampcor .eq. 0 ) ) then
c
c   Replace new image with last iteration (except for real-time mode)
c
               cim(:,:,1) = cim(:,:,3)
c
               cnew = cold
c
            else
c
               if ( quiet .gt. 1 ) write ( 6 , '(1x,a30,i10,f10.2)' )
     .            ' Contrast after iteration:    ' , icount , cnew
c
               write ( 7 , '(1x,a30,i10,f10.2)' )
     .            ' Contrast after iteration:    ' , icount , cnew
c
            endif
c
c   Iterate if necessary; remove amp error on last iteration
c
            if ( ( rt_pga .eq. 0 ) .and.
     .           ( cnew - cold ) * 100.0 / cold .le. percon ) then
c
               if ( ampcor .ne. 0 ) then
c
c   Store old image in case contrast decreases
c
                  cold = cnew
c
                  if ( rt_pga .eq. 0 ) then
c
                     cim(:,:,3) = cim(:,:,1)
c
                  endif
c
                  call reamp ( cim , nf , nr , aveamp , wt , taywt ,
     .                         cx , x , workft , worknf , nptgt ,
     .                         nrcell , nrfcell , dtf , cpp , tint )
c
                  cnew = ctrast( cim , nf , nr )
c
                  if ( ( rt_pga .eq. 0 ) .and. ( cold .gt. cnew ) ) then
c
c   Replace new image with last iteration
c
                     cim(:,:,1) = cim(:,:,3)
c
                     cnew = cold
c
                  endif
c
                  if ( quiet .gt. 1 ) write ( 6 , '(a40,f10.2)' )
     .               ' Contrast after amplitude correction:   ' , cnew
c
                  write ( 7 , '(a40,f10.2)' )
     .               ' Contrast after amplitude correction:   ' , cnew
c
               endif
c
               return
c
            else
c
               cold = cnew
c
            endif
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
      subroutine estpg ( cim , nf , nr , nrcell , nrfcell , nptgt, wt ,
     .                   npass , ntaper , cx , x , cxdot , work , pg ,
     .                   denom , worknf , dtf , aveamp )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , nr , nrcell(nr) , nrfcell(nr) , nptgt ,
     .        npass , ntaper , i , j
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
      subroutine estpg_jako ( cim , nf , nr , nrcell , nrfcell , nptgt ,
     .                        wt , npass , ntaper , cx , x , cxdot ,
     .                        cpaccum , work , pg , denom , worknf ,
     .                        dtf , aveamp , pgtype )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , nr , nrcell(nr) , nrfcell(nr) , nptgt ,
     .        npass , ntaper , i , j , pgtype
c
      complex cim(nf,nr) , cx(nf) , cxdot(nf), cpaccum(nf)
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
      aveamp(:)  = 0.0
c
      pg(:)      = 0.0
c
      denom(:)   = 0.0
c
      cpaccum(:) = cmplx( 0.0 , 0.0 )
c
      do i = 1 , nptgt
c
         call rshift ( wt , nf , nrfcell(nrcell(nr-i+1))-1 , worknf )
c
         call tfns ( cim , nf , nr , nrcell(nr-i+1) ,
     .               nrfcell(nrcell(nr-i+1)) , cx , x , cxdot ,
     .               worknf , work , dtf )
c
         if ( pgtype .eq. 2 ) then
c
            cpaccum(1) = cpaccum(1) + conjg( cx(1) ) * cx(2)
c
            aveamp(1)  = aveamp(1)  + abs( cx(1) )
c
            denom(1)   = denom(1)   + cx(1) * conjg( cx(1) )
c
            do j = 2 , nf - 1 
c
               aveamp(j)  = aveamp(j) + abs( cx(j) )
c
               cpaccum(j) = cpaccum(j) + conjg( cx(j) ) * cx(j+1) 
c
               denom(j)   = denom(j) + cx(j) * conjg( cx(j) )
c
            enddo
c
            cpaccum(nf) = cpaccum(nf) - conjg( cx(nf) ) * cx(nf-1)
c
            aveamp(nf)  = aveamp(nf)  + abs( cx(nf))
c
            denom(nf)   = denom(nf)   + cx(nf) * conjg( cx(nf) )
c
         else
c
c   Place the more general eigenvector calculation cases here if
c   pgtype = 2 works
c
         endif
c
      enddo
c
      pg(:)     = imag( log( cpaccum(:) ) ) 
c
      aveamp(:) = aveamp(:) / float( nptgt )
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
      call fourt ( x , nn , 1 , - 1 , 1 , work , nf )
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
     .                    nrcell , nptgt , all )
C**
C***********************************************************************
C**
      implicit none
c
      integer nn(1) , nf , nr , i , j , nrcell(nr) , nptgt , jtgt
c
      complex cim(nf,nr,3) , cx(nf) , cp(nf)
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
         cim(:,j,2) = cim(:,j,2) * cp(:)
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
            call cshift ( cim(1,j,2) , nf , nf / 2 , cx )
c
            call fourt ( x , nn , 1 , + 1 , 1 , work , nf )
c
            call cshift ( cx , nf , nf / 2 , cim(1,j,1) )
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
            call cshift ( cim(1,j,2) , nf , nf / 2 , cx )
c
            call fourt ( x , nn , 1 , + 1 , 1 , work , nf )
c
            call cshift ( cx , nf , nf / 2 , cim(1,j,1) )
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
      subroutine reamp ( cim , nf , nr , aveamp , wt , taywt , cx ,
     .                   x , work , worknf , nptgt , nrcell , nrfcell ,
     .                   dtf , cwork , tint )
C**
C***********************************************************************
C**
      implicit none
c
c-----------------------------------------------------------------------
c
c   Frequency passband (Hz) of amplitude modulation window
c
      real        fwidth
      parameter ( fwidth = 15.0 )
c
c-----------------------------------------------------------------------
c
      integer     nn(1) , nf , nr , i , j , npass_a , ntaper_a ,
     .            nrcell(nr) , nrfcell(nr) , nptgt , nstart , nend
c
      complex     cim(nf,nr,3) , cx(nf) , cwork(nf)
c
      real        aveamp(nf) , wt(nf) , taywt , x(2,nf), work(2,nf) ,
     .            ave , worknf(nf) , dtf , tint
c
c-----------------------------------------------------------------------
c
c   Inside passband of amplitude modulation window, remove amplitude
c   modulation error which is common to range cells from
c   real-space image.  The following procedure is used to remove the
c   amplitude modulation.
c
c     1. Compute Taylor Window Wt(t) and time average <Wt(t)>.  Compute
c        ratio r(t) = Wt(t) / <Wt(t)>.
c
c     2. For each target in (range,freq) image, window (fwidth Hz) about
c        target, shift to dc, and back transform to (range,time) image.
c
c     3. Model i-th target cell time series as constant signal + noise
c
c          yi(t) = ai*[1+e(t)]*W(t) = ai*W(t)+e(t)*W(t)
c
c        where ai = amplitude, e(t) = modulation, W(t) = weight
c
c     4. For each target cell, estimate "amplitude" time average
c
c          ai*< W(t) > = sum_over_time( yi(t) ) / nf
c
c     5. For each target cell, remove "amplitude" from time series
c
c          yi'(t) = yi(t) / ( ai*< W(t) > ) = [1+e(t)]*W(t)/<W(t)>
c
c     6. For each time t, estimate modulation correction across targets
c
c          c(t) = sum_over_targets( yi'(t) ) / ntargets
c              
c               = <[1+e(t)]*W(t)/<W(t)>
c
c     7. Apply modulation correction over all range cells
c
c          [yi(t)/c(t)] * r(t) = [ai*<W(t)>] * Wt(t)/<Wt(t)>
c
c        which is approximately equal to ai*W(t) as desired.
c
c-----------------------------------------------------------------------
c
c     cim(nf,nr,3)  :  Complex image; 1,3 freq domain, 2 time (in/out)
c     nf            :  Number of freq bins in image (in)
c     nr            :  Number of range bins in image (in)
c     aveamp(nf)    :  Estimated amplitude correction (out)
c     wt(nf)        :  Weight function work array
c     taywt         :  Sidelobe level in dB for Taylor weight (in)
c     cx(nf)        :  Work array
c     x(2,nf)       :  Work array
c     work(nf)      :  Work array
c     worknf(nf)    :  Work array
c     nptgt         :  Number of point targets (in)
c     nrcell(nr)    :  Range line with smallest to largest SNR (in)
c     nrfcell(nr)   :  Freq bin with largest SNR for range line (in)
c     dtf           :  Time separation between samples (in)
c     cwork(nf)     :  Work array
c     tint          :  Non-zero image integration time (sec) (in)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                      Frequency Domain Window
c
c   Compute amplitude modulation window with 10% taper; but use at
c   least 3 points.  Should be an odd number.
c
      npass_a   = nint( 2.0 * fwidth * float( nf ) * dtf )
      ntaper_a  = max( 3 , nint( 0.1 * float( npass_a ) ) )
      npass_a   = npass_a + 2 * ntaper_a
c
      if ( mod( npass_a , 2 ) .eq. 0 ) npass_a = npass_a + 1
c
      call hwindo ( wt , npass_a , ntaper_a , nf , worknf )
c
c-----------------------------------------------------------------------
c
c                        Time Domain Window
c
c   Compute index of start and end of non-zero time series
c
      nstart    = max(  1 , nint( 0.5 * ( float( nf ) - tint / dtf ) ) )
c
      nend      = min( nf , nint( 0.5 * ( float( nf ) + tint / dtf ) ) )
c
c-----------------------------------------------------------------------
c
      nn(1)     = nf    !  Fourier transform size
c
c-----------------------------------------------------------------------
c
c   Zero amplitude array for accumulation over targets
c
      aveamp(:) = 0.0
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Loop over range lines containing targets, estimating correction
c
      do j = 1 , nptgt
c
c   Center amplitude modulation window on target and apply
c
         call rshift ( wt , nf , nrfcell(nrcell(nr-j+1)) - 1 , worknf )
c
         cwork(:) = worknf(:) * cim(:,nrcell(nr-j+1),1)
c
c   Shift windowed target to DC
c
         call cshift ( cwork , nf , - ( nrfcell(nrcell(nr-j+1) ) - 1 ) ,
     .                 cx )
c
c   IFFT x and normalize
c
         call fourt ( x , nn , 1 , - 1 , 1 , work, nf )
c
         cx(:) = cx(:) / float( nf )
c
c   Reorder x to increasing time
c
         call cshift ( cx , nf , nf / 2 , cwork )
c
c   Compute the time average of < amp * w(t) >
c
         ave = 0.0
c
         do i = nstart , nend
c
            ave = ave + cabs( cwork(i) )
c
         enddo
c
         ave = ave / float( nend - nstart + 1 )
c
c   Accumulate amplitude time correction
c
         do i = nstart , nend
c
            aveamp(i) = aveamp(i) + abs( cwork(i) ) / ave
c
         enddo
c
      enddo  !  Loop over target range cells
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Compute average amplitude time correction over all targets
c
      do i = nstart , nend
c
         aveamp(i) = aveamp(i) / float( nptgt )
c
      enddo
c
c   Compute normalized Taylor weight function
c
      call taylor ( nend - nstart + 1 , taywt , wt )
c
      ave = 0.0
c
      do i = 1 , nend - nstart + 1
c
         ave = ave + wt(i)
c
      enddo
c
      ave = ave / float( nend - nstart + 1 )
c
      do i = 1 , nend - nstart + 1
c
         wt(i) = wt(i) / ave
c
      enddo
c
c   Compute Taylor-weighted amplitude correction over time.  This step
c   consists of dividing out the weight estimated from the data and
c   multiplying by the desired Taylor weight function.
c
      do i = nstart , nend
c
         wt(i-nstart+1) = wt(i-nstart+1) / aveamp(i)
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                 Implement amplitude correction
c
c   Loop over all range lines
c
      do j = 1 , nr
c
c   Remove amplitude error in the time domain
c
         do i = nstart , nend
c
            cim(i,j,2) = wt(i-nstart+1) * cim(i,j,2)
c
         enddo
c
c   Restore image lines by FFT and shift operations
c
         call cshift ( cim(1,j,2), nf , nf / 2 , cx )
c
         call fourt ( x , nn , 1 , + 1 , 1 , work , nf )
c
         call cshift ( cx , nf , nf / 2 , cim(1,j,1) )
c
      enddo
c
c-----------------------------------------------------------------------
c
      return
      end
C**
C***********************************************************************
C**
      subroutine gettgt ( cim , nf , nr , pnoise , pgatgt , nrcell ,
     .                    dtf , work , nptgt , rscore , dotdot , mprat ,
     .                    mrrat , nafill , actual , dfmp , time0 ,
     .                    rcen , fcen )
C**
C***********************************************************************
C**
      implicit none
c
      integer nf , nr , i , j , nrcell(nr) , index , nptgt , offset
c
      complex cim(nf,nr,3) , work(nf)
c
      real    pnoise , pgatgt(nr,5) , dtf , vel , d0 , dd , aavg ,
     .        npower , ddw , rscore(nr) , peak , fcon , dtdf , tgtime ,
     .        tgfreq , rcen , fcen
c
      integer mprat , mrrat , nafill , actual , irm , itgt , itg ,
     .        ileft , ileftp , ioff , ncent , ntot , nbest
c
      real    dotdot(1+nf/mprat,nr/mrrat,actual) , dfmp , accorg ,
     .        time0 , db , rtgt , f , delta , athrsh , dts , rbig ,
     .        rsmall
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
      include     'sarprm.h'      !  Standard ISAR-T parameters
c
      include     'updates.h'     !  Updates to parameters from the
c                                    first major release of the code
c
      include     'realtime.h'    !  Real-time parameters
c
      include     'tglist.h'
c
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
      ntot = int( 3 * ( nf * actual ) / nafill )/ 4 ! 3/4 actual inte-
c                                                   ! gration time
c
      ntot = 4 * ( max( 1 , ntot / 4 ) )            ! Multiple of 4, at
c                                                   ! least 4
c
      ioff = ( nf - ntot ) / 2                      ! Offset to start of
c                                                   ! middle ntot points
c
      do j = 1 , nr
c
         do i = 1 , ntot
c
            work(i) = cim(i+ioff,j,2)
c
         enddo
c
         if ( rt_pga .eq. 0 ) then
c
            call dchirp    ( work , ntot , pgatgt(j,4) , vel ,
     .                       pgatgt(j,3) , d0 , dd , dtf , - 1 , aavg ,
     .                       lambda , pgatgt(j,1) , pgatgt(j,2) )
c
         else
c
            call dchirp_rt ( work , ntot , pgatgt(j,4) , vel ,
     .                       pgatgt(j,3) , d0 , dd , dtf , - 1 , aavg ,
     .                       lambda , pgatgt(j,1) , pgatgt(j,2) )
c
         endif
c
      enddo
c
      ddw    = 0.2 * dff * ( float( nafill ) / float( actual ) )
      dts    = float( ntot / 4 ) * dtf
      athrsh = 5.0 / ( dts ** 2 )
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
c    Side Lobe Cancellation
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
      nbest = nint( 2.0 * overrg )
c
      call best ( pgatgt(1,5) , rscore , snrmin , nr , nrcell , nptgt ,
     .            nbest )
c
      if ( nptgt .gt. 0 ) then
c
         if ( rt_pga .eq. 0 ) then
c
            if ( quiet .gt. 1 )
     .      write ( 6 ,'(/,10x,a32)') 'PGA Targets (i,snr,f,da,a0,dw,r)'
c
            write ( 7 ,'(/,10x,a32)') 'PGA Targets (i,snr,f,da,a0,dw,r)'
c
         endif
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
            if ( rt_pga .eq. 0 ) then
c
               if ( quiet .gt. 1 )
     .         write ( 6 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                   pgatgt(itg,1) , pgatgt(itg,4) , accorg ,
     .                   pgatgt(itg,2) , rtgt
c
               write ( 7 ,'(i4,6f9.3)' ) itgt , db( pgatgt(itg,5) ) ,
     .                   pgatgt(itg,1) , pgatgt(itg,4) , accorg ,
     .                   pgatgt(itg,2) , rtgt
c
            endif
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
c   If we are not using sub-image targets then add these to the target
c   list.  Even if we are using sub-image targets make the call to
c   addtgt if we have read the targets from disk - in this case addtgt
c   will pass the target to the output diagnostic file without using it
c   for focusing.
c
c   Only add targets to list for Doppler autofocus (maxit>0)
c
            if ( ( maxit .gt. 0 ) .and.
     .           ( ( tgt_si .eq. 0 ) .or. ( rd_tgt .ne. 0 ) ) )
     .      call addtgt ( tgtime , rtgt + rcen ,
     .                    tgfreq + fcen , accorg + pgatgt(itg,4) ,
     .                    pgatgt(itg,2) , pgatgt(itg,5) , 'p' )
c
         enddo
c
      endif
c
c   If there is a valid rotation rate then eliminate all targets which
c   fall outside the ranges which comprise the target
c
      if ( omega_valid .ne. 0 ) then
c
         rsmall = corner(1,1,1)
         rbig   = corner(1,1,1)
c
         do k = 2 , 4
c
            rsmall = amin1( corner(k,1,1) , rsmall )
            rbig   = amax1( corner(k,1,1) , rbig   )
c
         enddo
c
         do itgt = 1 , nr
c
            rtgt = rfmin + drf * float( itg - 1 )
c
            if ( ( rtgt .gt. rbig ) .or.
     .           ( rtgt .lt. rsmall ) ) then
c
               pgatgt(itg,5) = - abs( pgatgt(itg,5) )
c
            endif
c
         enddo
c
         nbest = nint( 2.0 / overrg )
c
         call best ( pgatgt(1,5) , rscore , snrmin , nr , nrcell ,
     .               nptgt , nbest )
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
