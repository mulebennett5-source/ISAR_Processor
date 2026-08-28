C**
C***********************************************************************
C**
      subroutine quant ( qfirst , iqlsb , rc , nr )
C**
C***********************************************************************
C**
      implicit none
c
      logical     qfirst
c
      integer     nr , ir , ncall , navg
c
      real        iqlsb , rc(nr) , sqtot , iqsq , iqquant , delta
c
      integer     twindo         ! Time window for averaging stats
      parameter ( twindo = 100 )
c
      save        ncall , sqtot
c
c-----------------------------------------------------------------------
c
c     Initialization
c
      if ( qfirst ) then
c
         sqtot  = 0.0
c
         ncall  = 0
c
         qfirst = .false.
c
      endif
c
      ncall = ncall + 1
c
c-----------------------------------------------------------------------
c
c     Compute RMS value for this time
c
      iqsq  = 0.0
c
      do ir = 1 , nr
c
         iqsq = iqsq + rc(ir) ** 2
c
      enddo
c
c-----------------------------------------------------------------------
c
c     Apply time window to get smoothed RMS
c
      if ( ncall .le. twindo ) then
c
         sqtot = sqtot + iqsq
c
         navg  = ncall
c
      else
c
         delta = 1.0 / float( twindo )
c
         sqtot = ( 1.0 - delta ) * sqtot + iqsq
c
         navg  = twindo
c
      endif
c
c-----------------------------------------------------------------------
c
c     Apply quantization
c
      iqquant = iqlsb * sqrt( sqtot / float( nr * navg ) )
c
      do ir = 1 , nr
c
         rc(ir) = iqquant * float( nint( rc(ir) / iqquant ) )
c
      enddo
c
      return
	end