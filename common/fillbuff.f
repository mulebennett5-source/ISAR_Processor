C**
C***********************************************************************
C**
      subroutine fillbuff ( nt1 , ntp , spulse , npulse , isub , ntr ,
     .                      nchuse , nerror , naskip , nwork , rawfil ,
     .                      crc_raw , rc_raw , first , auxil , work )
C**
C***********************************************************************
C**
c   This routine fills a buffer with raw radar data from simulation or
c   from disk files
c
      implicit none
c
      include     'sarprm.h'
c
      include     'updates.h'
c
      include     'realtime.h'
c
      integer      nt1 , ntr , ntp , isub , nchuse , nerror , jpulse ,
     .             ir , ip , ich , loc , naskip , ndot , lastdot ,
     .             spulse , npulse , nwork , jpulse_prime , j , jj
c
c   Data buffer
c
      real         rc_raw(2,ntr,ntp,nchuse+1)  !  Raw data array -
      complex      crc_raw(ntr,ntp,nchuse+1)   !  possibly multi-channel
c
      real         auxil(10) , work(nwork) , inew , qnew
c
      character    rawfil*80
c
      logical      first , cfirst
c
      complex      ctmp
c
c-----------------------------------------------------------------------
c
c   Fill the buffer only for pulses nt1 --> nte
c
      do ip = nt1 , ntp
c
c   Calculate the pulse number required
c
         if ( npulse .gt. 0 ) then
c
c   Normal case, read pulses forwards in file
c
            jpulse = naskip * max( isub - 1 , 0 ) +
     .               ip + ( spulse - 1 )
c
         else
c
c   If npulse < 0, read pulses backwards in file (Ed Stockburger's
c   suggestion to focus the first section of the data by using target
c   motion derived from later data).
c
            jpulse = - naskip * max( isub - 1 , 0 )
     .               - ip + ( npulse + 1 )
c
         endif
c
c   Get a pulse either by generating it with one of the internal test
c   cases ( nfmt < 0 ) or by reading it from a data file ( nfmt >= 0 )
c
         if ( nfmt .lt. 0 ) then
c
            call iqgen ( crc_raw(1,ip,1) , ntr , jpulse , spulse ,
     .                   nerror )
c
            if ( nerror .ne. 0 ) return
c
         else
c
c   Read data
c
            do ich = 1 , nchuse           !  Possibly multi-channel data
c
               if ( nchuse .eq. 1 ) then
c
c   If single channel, it is stored in layer 1.
c
                  loc = 1
c
               else
c
c   If multiple channels, shift layers by one to make room for the
c   cancelled channel in layer 1.
c
                  loc = ich + 1
c
               endif
c
c   Up to three channels, with different file names
c
               if ( ich .eq. 1 ) then
c
                  ndot   = lastdot( ifile ) - 1
                  rawfil = ifile(1:ndot)  // '.raw'
c
               else if ( ich .eq. 2 ) then
c
                  ndot   = lastdot( ifile2 ) - 1
                  rawfil = ifile2(1:ndot) // '.raw'
c
               else if ( ich .eq. 3 ) then
c
                  ndot   = lastdot( ifile3 ) - 1
                  rawfil = ifile3(1:ndot) // '.raw'
c
               endif
c
c   Offset each pulse to line up the channels
c
               jpulse_prime = jpulse + pulse_offset(2,ich)
c
c   Global initialization - use temporary variable cfirst
c
               cfirst       = first
c
c   Finally, get a pulse of radar data for this channel
c
               call rawdat ( cfirst , rawfil , nfmt , jpulse_prime , 
     .                       ntr , hist , rc_raw(1,1,ip,loc) , auxil ,
     .                       nerror , work , nwork , nrcent )
c
c-----------------------------------------------------------------------
c
c   Bail out with error message if there is a problem with reading data
c
               if ( nerror .ne. 0 ) then
c
                  write ( 6 , * ) ' I/O Error: ' , nerror
                  write ( 7 , * ) ' I/O Error: ' , nerror
c
                  if ( nerror .eq. 1 ) then
c
                     write ( 6 , '(/,a25,/)' )
     .                     ' End of data file reached'
                     write ( 7 , '(/,a25,/)' )
     .                     ' End of data file reached'
c
                  else if ( nerror .eq. 2 ) then
c
                     write ( 6 , '(/,a27,a)' )
     .                     ' Error in reading data file'
                     write ( 7 , '(/,a27,a)' )
     .                     ' Error in reading data file'
c
                  endif
c
                  return
c
               endif
c
c-----------------------------------------------------------------------
c
c   Shift pulses in fast-time if required
c
               if ( pulse_offset(1,ich) .gt. 0 ) then
c
                  do j = 1 + pulse_offset(1,ich) , ntr
c
                     jj                 = j - pulse_offset(1,ich)
c
                     crc_raw(jj,ip,loc) = crc_raw(j,ip,loc)
c
                     crc_raw(j,ip,loc)  = cmplx( 0.0 , 0.0 )
c
                  enddo
c
               else if ( pulse_offset(1,ich) .lt. 0 ) then
c
                  do j = ntr + pulse_offset(1,ich) , 1 , - 1
c
                     jj                 = j - pulse_offset(1,ich)
c
                     crc_raw(jj,ip,loc) = crc_raw(j,ip,loc)
c
                     crc_raw(j,ip,loc)  = cmplx( 0.0 , 0.0 )
c
                  enddo
c
               endif
c
c-----------------------------------------------------------------------
c
c   Modify data to fit the ISAR-T format
c
c   Three transformations:
c
c         1.  Bit 1:  Conjugate data
c
c         2.  Bit 2:  Switch right and left sides to rotate the middle
c                     to DC
c
c         3.  Bit 3:  Flip signs of odd points to put zero range in
c                     middle
c
               if ( mod( moddat , 2 ) .eq. 1 ) then
c
c   Conjugate raw data
c
                  do ir = 1 , ntr
c
                     crc_raw(ir,ip,loc) = conjg( crc_raw(ir,ip,loc) )
c
                  enddo
c
               endif
c
               if ( ( mod( moddat , 4 ) .eq. 2 ) .or.
     .              ( mod( moddat , 4 ) .eq. 3 ) ) then
c
c   Shift to put DC in the middle
c
                  do ir = 1 , ntr / 2
c
                     ctmp                     = crc_raw(ir,ip,loc)
c
                     crc_raw(ir,ip,loc)       = crc_raw(ir+ntr/2,ip,loc)
c
                     crc_raw(ir+ntr/2,ip,loc) = ctmp
c
                  enddo
c
               endif
c
               if ( moddat .gt. 4 ) then
c
c   Flip sign of odd points
c
                  do ir = 1 , ntr , 2
c
                     crc_raw(ir,ip,loc) = - crc_raw(ir,ip,loc)
c
                  enddo
c
               endif
c
c-----------------------------------------------------------------------
c
               if ( rt_img .ne. 0 ) then
c
c   If in Real-time mode, add I/Q imbalance terms to simulate A/D errors
c
                  do ir = 1 , ntr
c
                     inew               = rc_raw(1,ir,ip,loc) +
     .                                    iqmix * rc_raw(2,ir,ip,loc)
c
                     qnew               = rc_raw(2,ir,ip,loc) +
     .                                    iqmix * rc_raw(1,ir,ip,loc)
c
                     qnew               = iqratio * qnew
c
                     crc_raw(ir,ip,loc) = cmplx( inew + ibar ,
     .                                           qnew + qbar )
c
                  enddo
c
               endif
c
            enddo     !  Loop over channels
c
         endif        !  Choice of data or simulation
c
         if ( addvib .ne. 0 .or. addamp .ne. 0 .or. fixamp .ne. 0 ) then
c
            call addvibamp( crc_raw(1,ip,loc) , ntr , addvib , addamp ,
     .                      fixamp , jpulse , dtp )
c
         endif
c
      enddo           !  Loop over pulses
c
      return
      end
C**
C***********************************************************************
C**
      subroutine addvibamp ( crc_raw , ntr , addvib , addamp , fixamp ,
     .                       jpulse , dtp )
C**
C***********************************************************************
C**
c   This routine adds vibration and amplitude modulation to data
c
c
      implicit none
c
      integer ntr , addvib , addamp , fixamp , jpulse , j
c
      real    dtp , t , pi , amp , vibAmp , vibFreq , ampAmp , ampFreq ,
     .        ptv , ampfix , totpow4 , totamp
c
      complex crc_raw(ntr) , cvib
c
      parameter ( vibAmp = 0.1 , vibFreq = 15.0 ,
     .            ampAmp = 0.1 , ampFreq = 20.0 )
c
      pi         = atan2( 0.0 , - 1.0 )
      t          = float(jpulse-1) * dtp
      ptv        = vibAmp*sin(2*pi*vibFreq*t)
c
      if ( addvib .ne. 0 ) then
         cvib   = cexp( cmplx( 0.0 , ptv ) )
      else
         cvib   = 1.0
      endif
c
      if ( addamp .ne. 0 ) then
         amp    = 1.0+ampAmp*sin(2.0*pi*ampFreq*t)
      else
         amp    = 1.0
      endif
c
      crc_raw(:) = crc_raw(:) * amp * cvib
c
      if ( fixamp .ne. 0 ) then
         totpow4 = 0.0
         totamp  = 0.0
         do j = 1 , ntr
            totpow4 = totpow4+cabs(crc_raw(j))**4
            totamp  = totamp +cabs(crc_raw(j))
         enddo
         totpow4 = totpow4/float(ntr)
         totamp  = totamp /float(ntr)
         ampfix  = 1.0/sqrt(sqrt(totpow4-totamp**4))
      else
         ampfix = 1.0
      endif
c
      crc_raw(:) = crc_raw(:) * ampfix
c
      return
      end
