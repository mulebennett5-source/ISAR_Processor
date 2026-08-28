C**
C***********************************************************************
C**
      subroutine mocomp ( vr0 , ctgt , ncr , ip , isub , isub0 ,
     .                    naskip , nfa , nr1 , nre , nr1fm , nrefm ,
     .                    dr0 , vrsave , ntp , sigsq , ipulse , spulse ,
     .                    cstr , vwt )
C**
C***********************************************************************
C**
c   Gently adjust the motion compensation point to try to follow the
c   centroid.  This is intended to reduce both low frequency (drift)
c   errors and high frequency (vibration) errors.
c 
c   There are 4 options for finemc:
c
c      1   :  Cross section weighting, un-smoothed time series
c
c      2   :  'Intelligent' weighting, un-smoothed time series
c
c      3   :  Cross section weighting, smoothed time series
c
c      4   :  'Intelligent' weighting, smoothed time series
c
c***********************************************************************
c
      implicit none
c
      include    'sarprm.h'
c
      integer     kp , kpwr , k0 , kplus , kminus , klast , know , ir ,
     .            nr1 , nre , ncr , ndwts , ip , naskip , isub , isub0 ,
     .            nfa , nr1fm , nrefm , nvc , ipulse , spulse , ntp 
c
      real        vwtsum , phbar , vwt(ncr) , u0 , uminus , uplus ,
     .            udel , vminus , vplus , vdel , delsum , dwt , dwtsum ,
     .            alpha , xmid , b2 , vr0 , vcent , power , sigsq ,
     .            sigsqp , dd , dr0 , pi , vrsave , vcbar , vcsq
c
      complex     c0 , tzero , tminus , tplus , cvsum , cvbar , ctemp ,
     .            ctgt(ncr,13) , cstr(ntp,ncr) , czero
c
      logical     no_vnotch
c
      integer     q
      parameter ( q = 2 )
c
      save
c
c   Initialize variables on first pulse
c
      if ( ip .eq. 1 ) then
c
c   The motion compensation algorithm is a simple fading estimate of
c   the average velocity estimated over the target region.
c
c   Fade constant for Mo-comp velocity estimates
c
         alpha  = alpha0 / float( naskip )
c
c   Fade constant for Mo-comp raw statistics
c
         b2     = amax1( beta ** 2 , 2.0 * alpha )
c
         pi     = atan2( 0.0 , - 1.0 )
c
         czero  = cmplx( 0.0 , 0.0 )
c
         nvc    = 0
         vcbar  = 0.0
         vcsq   = 0.0
c
c   Define pointers to the vectors for power and the three covariance
c   values
c
         kpwr   = 8   !  Power - no lag in time or range
         k0     = 9   !  Lag only in time
         kminus = 10  !  Lag in time and -q in range
         kplus  = 11  !  Lag in time and +q in range
         klast  = 12  !  Lag in time, possibly smoothed
         know   = 13  !  Present time, possibly smoothed
c
c   Apply 'intelligent' average over the middle half of the doppler
c   range of the fine resolution image
c
         xmid   = 6.28 * 0.25 * float( nlag ) * float( nfa ) * dff * dtp
c
c   Zero the covariance values
c
         do kp = 8 , 13
c
            ctgt(:,kp) = czero
c
         enddo
c
         if ( isub .eq. isub0 ) then
c
c   Initialize the Mo-comp variance estimate (single measurement) to
c   the most pessimistic value - half the alias velocity squared.
c
            sigsqp = ( lambda / ( 4.0 * dtp ) ) ** 2
            sigsq  = sigsqp
c
         endif
c
      else if ( ip .gt. nlag ) then
c
c-----------------------------------------------------------------------
c
c   Update the running estimates of the time series values
c
         if ( isub   .eq. isub0 .or.
     .        abs(finemc) .eq. 1     .or. abs(finemc) .eq. 2 ) then
c
c   Un-smoothed time series
c
            do ir = nr1 , nre
c
               ctgt(ir,klast) = cstr(ip-nlag,ir)
               ctgt(ir,know)  = cstr(ip,ir)
c
            enddo
c
         else
c
c   Smoothed time series
c
            do ir = nr1 , nre
c
               ctgt(ir,klast) = ( 1.0 - beta ) * ctgt(ir,klast)
     .                                  + beta * cstr(ip-nlag,ir)
c
               ctgt(ir,know)  = ( 1.0 - beta ) * ctgt(ir,know)
     .                                  + beta * cstr(ip,ir)
c
            enddo
c
         endif
c
c-----------------------------------------------------------------------
c
c   Update the running estimates of the covariance values
c
         do ir = nr1 + q , nre - q
c
            ctgt(ir,kpwr)   = ( 1.0 - b2 ) * ctgt(ir,kpwr) +
     .                  b2 * ctgt(ir,know) * conjg( ctgt(ir,know) )
c
            ctgt(ir,k0)     = ( 1.0 - b2 ) * ctgt(ir,k0) +
     .                  b2 * ctgt(ir,klast) * conjg( ctgt(ir,know) )
c
         enddo
c
      endif
c
c   After the covariance estimates are based on a long enough time
c   average, start estimating the local mo-comp velocity
c
      if ( ip .ge. nlag + ( naskip / 2 ) ) then
c
         tzero  = czero
         tminus = czero
         tplus  = czero
c
         vwtsum = 0.0
         ndwts  = 0
         dwtsum = 0.0
         cvsum  = czero
         delsum = 0.0
c
         vwt(:) = 0.0
c
         do ir = nr1fm + q , nrefm - q
c
c   Compute the drop in the covariance over one pulse and then, for the
c   intelligent weighting scheme compute the weight for the velocity
c   estimates.
c
            power  = cabs( ctgt(ir,kpwr) )
c
            if ( power .ne. 0.0 ) then
c
               dd = 1.0 - cabs( ctgt(ir,k0) ) / power
c
               if ( mod( abs(finemc) , 2 ) .eq. 0 ) then
c
                  vwt(ir) = 1.0 / sqrt( power * amax1( 0.0001 , dd ) )
c
               else
c
                  vwt(ir) = 1.0
c
               endif
c
            else
c
               vwt(ir) = 0.0
c
            endif
c
            vwtsum = vwtsum + vwt(ir)
            cvsum  = cvsum + vwt(ir) * ctgt(ir,k0)
c
         enddo
c
         if ( cvsum .ne. czero .and. vwtsum .gt. 0.0 ) then
c
            cvbar  = cvsum / vwtsum
            vwt(:) = vwt(:) / vwtsum
            phbar  = atan2( aimag( cvbar ) , real( cvbar ) )
c
         else
c
            phbar  = 0.0
c
         endif
c
c   Estimate of Doppler centroid for this pulse
c
         vcent = phbar * ( lambda / ( 4.0 * pi * float( nlag ) * dtp ) )
c
      endif   !  Done with all statistic updates
c
c-----------------------------------------------------------------------
c
c   For the first pass, set the mo-comp velocity at the naskip pulse
c
      no_vnotch = .not. ( notch .eq. 1 .or. notch .eq. 3 )
c
      if ( no_vnotch .and. ( isub .le. 0 .and. ip .eq. naskip ) ) then
c
         vr0    = vr0 + vcent
         vrsave = vr0
         vcent  = 0.0
c
c   Correct covariance estimates for this jolt in mo-comp velocity
c
         ctemp = cexp( cmplx( 0.0 , - phbar ) )
c
         do  ir = nr1 , nre
c
            ctgt(ir,k0)     = ctemp * ctgt(ir,k0)
c
            ctgt(ir,kminus) = ctemp * ctgt(ir,kminus)
c
            ctgt(ir,kplus)  = ctemp * ctgt(ir,kplus)
c
         enddo
c
c   Update the Mo-Comp velocity and its variance estimate
c
      else if ( ipulse .gt. spulse + naskip ) then
c
         vr0   = vr0 + alpha * vcent 
         sigsq = sigsqp * ( alpha ** 2 ) +
     .           sigsq * ( ( 1.0 - alpha ) ** 2 )
c
      endif
c
      vcbar = vcbar + ( vr0 + vcent )
      vcsq  = vcsq  + ( vr0 + vcent ) ** 2
      nvc   = nvc + 1
c
c   Update the single pulse variance estimate of the centroid
c
      if ( ip .eq. ntp ) then
c
         vcbar  = vcbar / float( nvc )
         sigsqp = ( vcsq / float( nvc ) ) - vcbar ** 2
         nvc    = 0
         vcbar  = 0.0
         vcsq   = 0.0
c
      endif
c
      return
      end
