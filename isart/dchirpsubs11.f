C**
C***********************************************************************
C**
      subroutine dchirp ( z , n , bavg , vel , sg0 , d0 , dd , dt ,
     .                    nterms , aavg , lambda , dc , dw )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , nterms , i , nlag , m , m2 , m22 , m222 , mlim
c
      complex z(n) , cj
c
      real    b1 , b2 , b3 , bavg , dt , a1 , a2 , a3 , aavg , lambda ,
     .        vel , d0 , dd , sg0 , pi , arg , ph , dc , dw , dcnew ,
     .        dwnew , sg0new , d0new , argnew , snr , snrold , snrnew ,
     .        bavg0 , dtm2
c
      logical done
c
c   Computes estimates of target acceleration, signal + noise power,
c   centroid, and Doppler width from a time series for a single range
c   cell.  Recursively improves estimates using filtered versions of
c   time series centered on the target.  All processing is done in the
c   time domain.
c
c    z(n)   :  Complex time series
c    n      :  Number of samples
c    bavg   :  Estimated target acceleration ( Hz/s )
c    vel    :  Estimated target velocity ( m/s )
c    sg0    :  Estimated target + noise power
c    d0     :  Estimated square of Doppler width ( non-dimensional )
c    dd     :  Initial square of Doppler width ( non-dimensional )
c    dt     :  Time separation between samples
c
c    nterms :  abs(nterms) =0,1,2 means despin time series with nothing,
c              current velocity estimate, current velocity and
c              acceleration estimates before making new
c              parameter estimates
c
c              nterms < 0 means perform recursive estimate improvement
c
c    aavg   :  Estimated target amplitude time rate of change
c    lambda :  Wavelength ( m )
c    dc     :  Estimated target Doppler centroid ( Hz )
c    dw     :  Estimated target Doppler width ( Hz )
c  
      pi   = atan2( 0.0 , - 1.0 )
      cj   = cmplx( 0.0 , 1.0 )
c
      nlag = 1
c
c   Estimate spectral parameters dc, dw, d0, sg0, and arg from signal
c   covariance using entire time series
c
      call onelag ( z , n , dt , nlag , dc , dw , d0 , sg0 , arg )
c
      vel  = dc * lambda / 2.0
      dd   = d0
c
c   If desired, remove amplitude variation, velocity, and possibly 
c   acceleration from time series to improve estimates of spectral
c   parameters.  Removing velocity by itself does not improve spectral
c   estimates, but having the target at dc will be necessary for the
c   next step
c
      if ( nterms .eq. 0 ) then ! don't estimate/remove amp, vel, accel
c
         aavg = 0.0
         bavg = 0.0
c
      else
c
c        Estimate accel and amp variation using three subsequences
c
         call accel ( z(1)     , n / 2 , b1 , dt , a1 )
         call accel ( z(1+n/4) , n / 2 , b2 , dt , a2 )
         call accel ( z(1+n/2) , n / 2 , b3 , dt , a3 )
c
         aavg = ( a1 + a2 + a3 ) / 3.0 ! ave amp variation
         bavg = ( b1 + b2 + b3 ) / 3.0 ! ave accel over 3 subsequences
c
c        Remove vel, amp variation, and possibly acceleration
c
         do i = 1 , n
c
            if ( iabs( nterms ) .eq. 1 ) ph = arg * float( i ) ! vel
c
            if ( iabs( nterms ) .eq. 2 ) ph = arg * float( i )
     .             + bavg * pi * ( dt * float( i - n / 2 ) ) ** 2
c
            z(i) = z(i) * cexp(
     .             cj * ph - 0.5 * aavg * float( i - n / 2 ) * dt )
c
         enddo
c
c        Compute improved spectral estimates
c
         call onelag ( z , n , dt , nlag , dcnew , dw , dd , sg0 , arg )
c
         dc  = dc + dcnew
         vel = dc * lambda / 2.0
c
      endif
c
c   Compute how many times the filtered time series centered on the
c   target can be cut in half and still end up with 16 points in the
c   smoothed time series.  16 points represent the minimum desired
c   length time series that can be cleanly broken into four disjoint
c   sections allowing three independent acceleration estimates
c
      m    = n
      mlim = 0
c
      do while ( m .ge. 32 )
c
        m    = m / 2
        mlim = mlim + 1
c
      enddo
c
c   Recursively improve spectral estimates using filtered versions
c   of time series centered on the target.  All processing is done in
c   the time domain.  Nyquist is effectively halved on each iteration.
c
      if ( nterms .lt. 0 .and. mlim .gt. 0 ) then
c
c        Make more passes using filtered versions of the time series
c
         m      = 0
         snrold = snr( sg0 , dw , dt , 1.0 )
         bavg0  = bavg
         done   = .false.
c
         do while ( m .lt. mlim .and. ( .not. done ) )
c
            m    = m + 1
            m2   = 2 ** m
            m22  = 2 * m2
            m222 = 2 * m22
c
            dtm2 = float( m2 ) * dt
c
c           Decimate 1/2 and smooth time series
c
            do i = 1 , n / m2
c
               z(i) = 0.5 * ( z(2*i-1) + z(2*i) )
c
            enddo
c
c           Compute improved spectral estimates
c
            call onelag ( z , n / m2 , dtm2 , nlag , dcnew ,
     .                    dwnew , d0new , sg0new , argnew )
c
            snrnew = snr( float( m2 ) * sg0new , dwnew , dt , 1.0 )
c          
c           If snr increases and Doppler width decreases, update
c           spectral estimates, despin the target to dc, and iterate
c
            if ( snrnew .gt. snrold .and. dwnew .lt. dw ) then
c
               snrold = snrnew
               sg0    = float( m2 ) * sg0new
               dw     = dwnew
               d0     = d0new
               arg    = argnew
c
               call accel ( z(1)        , n / m22 , b1 , dtm2 , a1 )
               call accel ( z(1+n/m222) , n / m22 , b2 , dtm2 , a2 )
               call accel ( z(1+n/m22)  , n / m22 , b3 , dtm2 , a3 )
c
               if ( nterms .eq. - 2 ) then
c
                  bavg = bavg0 + ( b1 + b2 + b3 ) / 3.0
c
               else
c
                  bavg = ( b1 + b2 + b3 ) / 3.0
c
               endif
c
               dc  = dc + dcnew
               vel = dc * lambda / 2.0
c
               do i = 1 , n / m2
c
                  z(i) = z(i) * cexp( cj * arg * float( i ) )
c
               enddo
c
            else
c
               done = .true.
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
      subroutine accel ( z , n , b , dt , a )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i
c
      complex z(n) , zsum1 , zsum2
c
      real    b , dt , pi , arg1 , arg2 , argdif , a
c
      pi = atan2( 0.0 , - 1.0 )
c
c  
c  Estimates signal acceleration and amplitude rate for time series
c  from velocities and amplitudes from the first and second half
c  of time series.
c
c     z(n)  :  Complex time series
c     n     :  Number of samples
c     b     :  Acceleration of signal ( Hz/s )
c     dt    :  Time separation between samples
c     a     :  Amplitude rate of signal
c
c  Compute signal covariance for first and second half of time series
c
      zsum1 = cmplx( 0.0 , 0.0 )
      zsum2 = cmplx( 0.0 , 0.0 )
c
      do i = 1 , n / 2 - 1
c
         zsum1 = zsum1 + z(i) * conjg( z(i+1) )
c
      enddo
c
      do i = 1 + n / 2 , n - 1
c
         zsum2 = zsum2 + z(i) * conjg( z(i+1) )
c
      enddo
c
c  Compute phase of signal covariance (proportional to velocity)
c  for first and second half of time series 
c
      if ( zsum1 .eq. cmplx( 0.0 , 0.0 ) ) then
c
         arg1 = 0.0
c
      else
c
         arg1 = atan2( aimag( zsum1 ) , real( zsum1 ) )
c
      endif
c
      if ( zsum2 .eq. cmplx( 0.0 , 0.0 ) ) then
c
         arg2 = 0.0
c
      else
c
         arg2 = atan2( aimag( zsum2 ) , real( zsum2 ) )
c
      endif
c
c  Compute phase difference (proportional to velocity difference) from
c  first to second half of time series; if necessary, dealias phase diff
c
      argdif = arg2 - arg1
      if ( argdif .gt.   pi ) argdif = argdif - 2.0 * pi
      if ( argdif .lt. - pi ) argdif = argdif + 2.0 * pi
c
c  Estimate fractional amplitude rate (1/amp) * d[amp]/dt from the
c  first to second half of the time series. Use the equivalent
c  expression d[ln(amp)]/dt and the ratio amp2/amp1 to compute the rate.
c  The 0.1 terms were added to save the logarithm if the amp in the
c  first or second half (but not both) of the time series was zero.
c
      if ( cabs( zsum1 ) .eq. 0.0 .or. cabs( zsum2 ) .eq. 0.0 ) then
c
         a = 0.0
c
      else
c
         a = 0.5 * alog( ( 0.1 * cabs( zsum1 ) + cabs( zsum2 ) ) /
     .                   ( 0.1 * cabs( zsum2 ) + cabs( zsum1 ) ) ) /
     .                       ( float( n ) * dt / 2.0 )
c
      endif
c
c  Estimate acceleration in Hz/s as d[phi]/dt / (2*pi*dt)
c
      b = argdif / ( 2.0 * pi * dt * dt * float( n / 2 ) )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine onelag ( c , n , dt , nlag , dc , dw , d , pwr , phi )
C**
C***********************************************************************
C**
c
c   Estimate the spectral parameters from the signal covariance
c   calculated at a single lag
c
c      c(n)  :  Complex time series
c      n     :  Number of samples
c      dt    :  Time separation between samples
c      nlag  :  Lag at which covariance is to be estimated
c      dc    :  Estimated Doppler centroid ( Hz )
c      dw    :  Estimated Doppler width ( Hz )
c      d     :  Square of Doppler width ( non-dimensional )
c      pwr   :  Total signal + noise power
c      phi   :  Phase of covariance ( radians )
c
      implicit none
c
      integer n , nlag , i
c
      complex c(n) , cov , dcov
c
      real    dt , dc , dw , d , cova , pwr , phi , pi
c
      pi   = atan2( 0.0 , - 1.0 )
      cov  = cmplx( 0.0 , 0.0 )
      cova = 0.0
      pwr  = 0.0
c
      do i = 1 , n - nlag
c
         dcov = c(i) * conjg( c(i+nlag) )           ! lagged product
         cov  = cov  + dcov                         ! signal covariance
         cova = cova + cabs( dcov )                 ! signal power
         pwr  = pwr  + cabs( c(i) * conjg( c(i) ) ) ! signal + noise
c
      enddo
c
c   Compute total signal + noise power
c
      pwr  = pwr + 0.5 * ( cabs( c(n) * conjg( c(n) ) ) -
     .                     cabs( c(1) * conjg( c(1) ) ) )
      pwr  = amax1( 0.0 , pwr ) ! possibly trap above minus sign
c
c   Compute phase of signal covariance
c
      if ( cov .eq. cmplx( 0.0 , 0.0 ) ) then
c
         phi = 0.0
c
      else
c
         phi = atan2( aimag( cov ) , real( cov ) )
c
      endif
c
c   Estimate Doppler centroid from covariance phase
c
      dc   = phi / ( 2.0 * pi * float( nlag ) * dt )
c
c   Estimate square of Doppler width from the signal covariance
c   magnitude and signal power
c
      if ( cova .eq. 0.0 ) then
c
         d = 1.0
c
      else
c
         d = amax1( 0.0 , 1.0 - cabs( cov ) / cova )
c
      endif
c
c   Compute Doppler width in Hz
c
      dw   = sqrt( 2.0 * d ) / ( 2.0 * pi * float( nlag ) * dt )
c
      return
      end
