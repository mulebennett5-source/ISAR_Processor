C**
C***********************************************************************
C**
      subroutine ellipse_param ( dw_wt , n_inside , major_axis , 
     .                           minor_axis , mean_r , mean_f )
C**
C***********************************************************************
C**
c   Computes the major and minor axes of the 1-sigma concentration
c   ellipses formed by assuming a correlated gaussian scatter model for
c   the current target list.
c   
c
c     n_inside    :  Total targets used
c
c     dw_wt       :  Relative weight applied to Doppler width
c
c     major_axis  :  The major and minor axis of the ellipse
c     minor_axis     in pixel coordinates
c   
c     mean_r      :  Mean position vector of the center of the ellipse
c     mean_f         in pixels
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
      integer  n_inside , iloop
c
      real     dw_wt 
c
      real     sum_range , sum_xrange , pi , wt , sumwt , pc(2,2) , a ,
     .         snrwt , sumsnrwt , raddeg , snr_range , freq_t0 ,
     .         snr_xrange , sumboth , sn_wt , range_t0 , eigx , eigy , 
     .         major_axis , minor_axis , sumwt2 , mean_r , mean_f 
c
c-----------------------------------------------------------------------
c
      if ( n_inside .lt. 3 ) then
c
         major_axis = - 1.0
c
         minor_axis =   0.0
c
         return
c
      endif
c
c-----------------------------------------------------------------------
c
      pi       = atan2( 0.0 , - 1.0 )
      raddeg   = pi / 180.0
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      a        = drf / dff                ! Calculate in pixels
c
c   The sum of the Doppler width and the snr weights should be 1.
c
      sn_wt    = 1.0 - dw_wt
c
c-----------------------------------------------------------------------
c
c   Compute mean range and frequency
c
c-----------------------------------------------------------------------
c
c   Doppler width weights
c
      sum_range  = 0.0    ! Initialize for statistics
      sum_xrange = 0.0
      sumwt      = 0.0
c
c   SNR weights
c
      snr_range  = 0.0    ! Initialize for statistics
      snr_xrange = 0.0
      sumsnrwt   = 0.0
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
            sum_xrange = sum_xrange   
     .                 + wt * a * freq_t0
c
            snr_range  = snr_range    
     .                 + snrwt * range_t0
            snr_xrange = snr_xrange   
     .                 + snrwt * a * freq_t0
c
         endif
c
      enddo
c
      mean_r  = dw_wt * ( sum_range  / sumwt ) + 
     .          sn_wt * ( snr_range  / sumsnrwt )
c
      mean_f  = dw_wt * ( sum_xrange / sumwt ) +
     .          sn_wt * ( snr_xrange / sumsnrwt )
c
      return
      end
