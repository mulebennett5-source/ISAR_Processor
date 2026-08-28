      subroutine choose_tinteg( tinteg_1 , ti_con , tframe , deltat ,
     .                          nlim , rmsdop , rmsacc , ntot )
c
      implicit none
c
c   Compute a recommended image integration time from Doppler statistics
c   of the sub-image target reports
c
c   Use the global target list stored in common
c
      include 'tglist.h'
c
c   Local and output variables
c
      integer nt , nlim , iloop , ntot
c
      real    ti_con , tframe , deltat , tinflu , tbeg , tend ,
     .        dopsq , accsq , tinteg_1 , rmsdop , rmsacc
c
c   Set the time window for this frame
c
      tinflu = float( nlim ) * deltat
c
      tbeg   = amax1( 0.0 , tframe - 1.0 * tinflu )
      tend   = amax1( 0.0 , tframe + 1.0 * tinflu )
c
c   Compute the Doppler statistics for the selected targets
c
      ntot   = 0
      dopsq  = 0.0
      accsq  = 0.0
      
      do iloop = 1 , nlist
c
         if ( ( time(iloop)  .ge. tbeg ) .and.
     .        ( time(iloop)  .le. tend ) .and.
     .        ( iflag(iloop) .ne. 101 ) ) then
c		   
            ntot   = ntot+1
            dopsq  = dopsq + freq(iloop)**2
            accsq  = accsq + accel(iloop)**2
c
         endif
c
      enddo
c
c   If there are enough target reports then compute the statistics and
c   the recommended integration time. Otherwise just set the output
c   values to -1
c
      if ( ntot .ge. 100 ) then
c
         rmsdop   = sqrt( dopsq / ntot )
         rmsacc   = sqrt( accsq / ntot )
c
         tinteg_1 = abs( ti_con ) * sqrt( rmsdop / rmsacc )
c
      else
c
         rmsdop   = -1.0
         rmsacc   = -1.0
c
         tinteg_1 = -1.0
c
      endif
c
      return
      end
