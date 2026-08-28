C**
C***********************************************************************
C**
      subroutine uwbsub ( csbimg , nsa , nsr , drs , drc , lambda ,
     .                    slant0 , cstr , ncr , dfc , naskip , wtac ,
     .                    dtp , vplat , rphase , cphase , br , strtch ,
     .                    clight )
C**
C***********************************************************************
C**
c
c   For Ultra-Wide-Band strip-map mode the sub-image formation must
c   account for range walk and quadratic focus.  This routine does this
c   by directly integrating the range-compressed signal history on a
c   linear path for each value of frequency and broadside range.
c
c   INPUT:
c
c      nsa      :  Number of sub-image frequency bins
c
c      nsr      :  Number of sub-image range bins
c
c      drs      :  Broadside range separation between sub-image samples
c
c      drc      :  Slant range separation for range-compressed data
c
c      lambda   :  Radar wavelength
c
c      slant0   :  Broadside range to sub-image center
c
c      cstr     :  Range-compressed signal history data
c
c      ncr      :  Number of range lines in range-compressed data
c
c      naskip   :  Number of time samples between calls to this routine
c
c      wtac     :  Array of Taylor weights
c
c      dtp      :  Pulse repetition interval (1/PRF)
c
c      vplat    :  Platform velocity
c
c   OUTPUT:
c
c      csbimg   :  Sub-image array
c
C**
C***********************************************************************
C**
      implicit none
c
      integer     nsa , nsr , naskip , ncr , i , j , k , nacent , jj ,
     .            nlast
c
      real        slant0 , rbroad , drs , drc , lambda , dtp , vplat ,
     .            fc , wtac(3*naskip) , dfc , eps , rdot , rcmin ,
     .            rloc , rcen , rcenr0 , twolam , dotdot , t , rcmax ,
     .            br , epsbr , clight , drcent , drhalf
c
      complex     cstr(4*naskip,ncr) , csbimg(nsa,nsr) , cpfast
c
      complex     cphase(3*naskip)
c
      real        rphase(3*naskip)
c
      integer     strtch
c
      integer     jbig
      parameter ( jbig = 32 )
c
c-----------------------------------------------------------------------
c
c   Min and max slant range of the range-compressed data in units of
c   the range cell size
c
      rcmin  = ( slant0 / drc ) - float( ncr / 2 )
      rcmax  = rcmin + float( ncr - 2 )
c
      drhalf = float( jbig / 2 ) * drc
c
c   Note: rcmin, rcmax, rphase, rloc and rcen are non-dimensional (in
c         units of the slant range cell size.  However, the variable
c         rbroad is in units of meters and measures the slant range at
c         broadside.
c
      twolam = 2.0 * drc / lambda
c
      nlast  = 3 * naskip
      nacent = nlast / 2
c
c   Loop over coarse doppler frequency
c
      do k = 1 , nsa
c
         fc     = dfc * ( float( k ) - 1 - nsa / 2 )
         rdot   = 0.5 * lambda * fc
         eps    = ( rdot / vplat ) ** 2
         rcenr0 = 1.0 / sqrt( 1.0 - eps )
c
c   Loop over sub-image broadside range
c
         do j = 1 , nsr
c
            rbroad = slant0 + float( j - 1 - nsr / 2 ) * drs
c 
c   To decrease computation, pre-compute the phase function every
c   16 range cells using the range at the center of each region in
c   the calculation of the quadratic term.
c
            if ( mod( j , jbig ) .eq. 1 ) then
c
               if ( strtch .eq. 1 ) then
c
                  drcent = rcenr0 * ( rbroad + drhalf ) - slant0
                  epsbr  = ( 2.0 * br * 1.0E+12 * lambda * drcent ) /
     .                               clight ** 2
                  rdot   = 0.5 * lambda * fc / ( 1.0 - epsbr )
                  eps    = ( rdot / vplat ) ** 2
                  rcenr0 = 1.0 / sqrt( 1.0 - eps )
c
                  twolam = ( 2.0 * drc / lambda ) * ( 1.0 - epsbr )
c
               endif
c
               dotdot = ( vplat ** 2 - rdot ** 2 ) /
     .                  ( rcenr0 * ( rbroad + drhalf ) )
c
               do i = 1 , nlast
c
                  t         = ( float( i - nacent ) - 0.5 ) * dtp
                  rphase(i) = ( rdot * t + 0.5 * t * t * dotdot ) / drc
                  cphase(i) = wtac(i) * cpfast( twolam * rphase(i) )
c
               enddo
c
            endif
c
c   Compute a complex time series along the path for this doppler
c   frequency; then sum the complex weighted values of range-compressed
c   signal history.
c
            csbimg(k,j) = cmplx( 0.0 , 0.0 )
c
c   Only use line segments which fall completely in the range bins of
c   the range-compressed data.
c
            rcen   = rcenr0 * rbroad / drc
c
            if ( rcen + rphase(1) .gt. rcmin .and.
     .           rcen + rphase(nlast) .lt. rcmax ) then
c
c   Loop over time
c
               do i = 1 , nlast
c
                  rloc   = rcen + rphase(i)
c
                  jj     = 1 + nint( rloc - rcmin )
c
                  csbimg(k,j) = csbimg(k,j) +
     .                          cphase(i) * cstr(i+naskip/2,jj)
c
               enddo
c
            endif
c
         enddo    ! Done with range cells for this coarse frequency
c
      enddo       ! Done with coarse frequencies
c
      return
      end
