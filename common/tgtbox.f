C**
C***********************************************************************
C**
      subroutine tgtbox ( n_inside , dw_wt , wiggle , kill , sdfact )
C**
C***********************************************************************
C**
c   Computes (f,r) corners of box around target. If omega is not known,
c   corners represent target extent.  If omega is known, use given
c   length x width box to edit outliers in (x,r) and return box (f,r)
c   corners.
c
c     n_inside  :  Total targets used
c
c     dw_wt     :  Relative weight applied to Doppler width
c
c     wiggle    :  Fraction of half-length or half-width allowed before
c                  a target is considered outside and therefore not used
c                  for the motion model for this frame
c
c     kill      :  Fraction outside box before target is permanently
c                  eliminated
c
c-----------------------------------------------------------------------
c
      implicit none
c
      include 'tglist.h'
c
      include 'sarprm.h'
c
      include 'updates.h'
c
      integer  n_inside , iloop, i , j , k
c
      real     dw_wt , wiggle , kill , sdfact
c
      integer  nsleep , nkill
c
      real     mean_range , mean_xrange , sum_range , sum_xrange , pi ,
     .         wt , sumwt , pc(2,2) , a , snrwt , sumsnrwt , raddeg ,
     .         snr_range , freq_t0 , snr_xrange , sumboth , sn_wt ,
     .         grace_x , toofar_x , grace_y , toofar_y , lw2 , xpmin ,
     .         xpmax , ypmin , ypmax , xpbar , ypbar , range_t0 ,
     .         eigx , eigy , temp , a_apriori , xpc_ap , ypc_ap , sd_x ,
     .         sd_y , sumwt2
c
      real     l2  , w2   !  Half-length and half-width
c
      real     xp  , yp   !  Coords. of current center in rotated frame
c
      real     xpc , ypc  !  Coords. of current corner in rotated frame
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Stewart Strait's F-statistic calculation for future diagnostic
c
      real     fstat , fmean , fvar  !  F-statistic for length = width
c                                    !  and mean and variance of that
c                                    !  statistic
c
      real     fz                    !  Z-score of the F-statistic
c
      real     dof                   !  Approximate degrees of freedom
c                                    !  for F-statistic
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      pi          = atan2( 0.0 , - 1.0 )
      raddeg      = pi / 180.0
c
c-----------------------------------------------------------------------
c
c   Initialize counters for targets allowed to 'sleep' for this frame
c   and those which are killed (zapped for this and all future frames).
c
      nsleep      = 0
c
      nkill       = 0
c
c   Compute half-length and half-width
c
      l2          = length / 2.0
      w2          = width  / 2.0
c
      lw2         = sqrt( abs( l2 * w2 ) )
c
c   Corners are set to negative 10 ** 6 if there is no basis for
c   computing the box due to insufficient targets or no valid rotation
c   rate.
c
c-----------------------------------------------------------------------
c
      do i = 1 , 4
c
         do j = 1 , 2
c
            corner(i,j,1) = - 1000000.0
            corner(i,j,2) = - 1000000.0
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   If insufficient targets reported, don't do algorithm
c
      if ( n_inside .lt. 3 ) return
c
c-----------------------------------------------------------------------
c   The sum of the Doppler width and the snr weights should be 1.
c
      sn_wt = 1.0 - dw_wt
c
c-----------------------------------------------------------------------
c
c   Compute mean range and frequency
c
c-----------------------------------------------------------------------
c
c   Doppler width weights
c
      sum_range   = 0.0    ! Initialize for statistics
      sum_xrange  = 0.0
      sumwt       = 0.0
c
c   SNR weights
c
      snr_range   = 0.0    ! Initialize for statistics
      snr_xrange  = 0.0
      sumsnrwt    = 0.0
c
      do iloop = 1 , nlist
c
         if ( iflag(iloop) .eq. editac ) then
c
c   Compute weight functions
c
            wt         = 1.0 / dwdth(iloop)
c
            snrwt      = snr(iloop)
c
            sumwt      = sumwt + wt
c
            sumsnrwt   = sumsnrwt + snrwt
c
c   Correct range and frequency to center of aperture time
c
            range_t0   = range(iloop) + 0.5 * lambda * freq(iloop) *
     .                                        ( tbar - time(iloop) )
c
            freq_t0    = freq(iloop)  + accel(iloop) *
     .                                  ( tbar - time(iloop) )
c
            sum_range  = sum_range    
     .                 + wt * range_t0
c
            sum_xrange = sum_xrange   
     .                 + wt * a * freq_t0
c
            snr_range  = snr_range    
     .                 + snrwt * range_t0
c
            snr_xrange = snr_xrange   
     .                 + snrwt * a * freq_t0
c
         endif
c
      enddo
c
      mean_range  = dw_wt * ( sum_range  / sumwt ) + 
     .              sn_wt * ( snr_range  / sumsnrwt )
c
      mean_xrange = dw_wt * ( sum_xrange / sumwt ) +
     .              sn_wt * ( snr_xrange / sumsnrwt )
c
      return
      end
