c**
c***********************************************************************
c**
      subroutine rfcomb0 ( cac , crc , nsa , nakeep , ntr , nbands ,
     .                     rfdelt , drfreq , nrfsmth , isub , ndop )
c**
c***********************************************************************
c**
c   Purpose:  To combine RF bands from a step-chirp system
c
c   cac     :  On input this is the data in Doppler and fast-time with
c              the data as bands that were collected sequentially.  This
c              is implicitly dimensioned as [nsa, ntr/nbands, nbands].
c
c              On output this data is merged by time-correcting the
c              Doppler bins and blending the fast-time data together.
c
c   crc     :  Complex work array
c
c   nsa     :  Number of Doppler cells - assumed to be the entire FFT of
c              the slow-time data with the DC at 1+nsa/2.
c
c   nakeep  :  The Number of Doppler cells to be used - cells outside of
c              this center are ignored (i.e., not corrected).
c
c   ntr     :  Total number of fast-time samples for all bands - must be
c              an even number times nbands.
c
c   nbands  :  Number of step-chirp bands
c
c   rfdelt  :  Difference between center frequencies of bands [Hz]
c
c   drfreq  :  Frequency step [Hz] for 'fast-time'
c
c   nrfsmth :  Number of smoothing passes done at the first stage -
c              usually zero.
c
c   isub    :  Sequence number for input/output data, used to implement
c              a low-pass filter over the adaptive phase estimates.
c              Setting isub=1 resets the routine for a new data set.
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer     nsa , ntr , nbands , iband , n2 , ovrlap2 , joff ,
     .            iovrlap , ntrband , nctr , lctr , rctr , isub ,
     .            ileft , iright , i , j , nftotal , nakeep , nalose ,
     .            jj , iter , ntotal , middle , nbdiv , nrfsmth ,
     .            dopwidth , ndop
c
      real        rfc(nbands) , rfi(nbands) , pi , rfdelt , drfreq ,
     .            ridelt , delta , phasec , tfade , wtall(ntr) ,
     .            wt(ntr/nbands,nbands)
c
      complex     cac(ndop,ntr) , crc(ntr) , czero , covrf(nbands-1) ,
     .            c , power(nbands-1) , covcor , cleft , cright ,
     .            covrfd(nbands-1)
c
c-----------------------------------------------------------------------
c
c   Total number of bands allowable and number of calls used for time
c   averaging of adaptive weights
c
      integer     nbandmax , ntfilter
c
      parameter ( nbandmax = 32 , ntfilter = 10 )
c
      complex     covrfsave(nbandmax)
c
      save        covrfsave  !  Accumulate adaptive weights over time
c
c   Rarely used option to use only the average value over the blended
c   region for computing the adaptive weights.  This would weight only
c   those targets near the range center.
c
      integer     dcweight
c
      parameter ( dcweight = 0 )
c
c-----------------------------------------------------------------------
c
      if ( nbands .lt. 2 ) return  !  Nothing to do
c
c   Consistency checks - ntr divisible by nbands and ntr/nbands even
c
      if ( mod( ntr , nbands ) .ne. 0 ) stop 'ntr, nbands inconsistent'
c
      if ( mod( ntr / nbands , 2 ) .ne. 0 ) stop 'ntr/nbands not even'
c
      if ( nbands .gt. nbandmax ) stop 'Too many bands in rfcomb0'
c
c-----------------------------------------------------------------------
c
c   Constants
c
      czero      = cmplx( 0.0 , 0.0 )
c
      pi         = atan2( 0.0 , - 1.0 )
c
c-----------------------------------------------------------------------
c
c   Initialize the running total for the adaptive phase estimator
c
      if ( isub .eq. 1 ) then
c
         covrfsave(:) = czero
c
      endif
c
c   Offset to first Doppler cell to be corrected
c
      nalose     = max( ( nsa - nakeep ) / 2 , 0 )
c
c   Middle (DC) for Doppler
c
      if ( nakeep .eq. nsa * nbands ) then
c
         middle = 1 + ( nbands / 2 ) * nsa + nsa / 2
c
      else
c
         middle = 1 + nsa / 2
c
      endif
c
c   Assume that adaptive algorithm for band phase difference can use
c   only the middle quarter (+1) of the cells being kept
c
      dopwidth   = min( nsa , nakeep ) / 8
c
c   Number of fast-time samples per RF band
c
      ntrband    = ntr / nbands
c
c   Center cell in each band
c
      nctr       = 1 + ntrband / 2
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 1
c
c   Smooth the RF data to isolate the target at DC
c
      if ( nrfsmth .gt. 0 ) then
c
         do i = 1 + nalose , nakeep + nalose
c
            do iband = 1 , nbands
c
               do iter = 1 , nrfsmth
c
c                 Even points
c
                  do j = 2 , ntrband - 2 , 2
c
                     jj        = ( iband - 1 ) * ntrband + j
c
                     cac(i,jj) = 0.50 * cac(i,jj) +
     .                           0.25 * ( cac(i,jj-1) + cac(i,jj+1) )
c
                  enddo
c
c                 Odd points
c
                  do j = 3 , ntrband - 1 , 2
c
                     jj        = ( iband - 1 ) * ntrband + j
c
                     cac(i,jj) = 0.50 * cac(i,jj) +
     .                           0.25 * ( cac(i,jj-1) + cac(i,jj+1) )
c
                  enddo
c
               enddo
c
            enddo
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 2
c
c   Shift each band to make them all line up on the same pickets and to
c   be consistent with the new interpolated frequencies.
c
c
c   nftotal = Nearest integer number of points between the centers of
c             the first and last band
c
c   nbdiv   = Number of convenient half-bands to split this up into
c
      call ntrcomb ( nftotal , nbdiv , nbands , rfdelt , drfreq )
c
c   Array of center frequencies for the bands in the original data
c
      do iband = 1 , nbands
c
         rfc(iband) = rfdelt * ( float( iband - 1 - nbands / 2 )
     .                - 0.5 * float( 1 - mod( nbands , 2 ) ) )
c
      enddo
c
c   Compute the interpolated center frequencies
c
      ridelt  = float( nftotal ) * drfreq / float( nbands - 1 )
c
      do iband = 1 , nbands
c
         rfi(iband) = ridelt * ( float( iband - 1 - nbands / 2 )
     .                - 0.5 * float( 1 - mod( nbands , 2 ) ) )
c
      enddo
c
c  Finally, interpolate all bands to the nicely spaced pickets
c
      do iband = 1 , nbands
c
         delta = ( rfi(iband) - rfc(iband) ) / drfreq
c
         call dshift ( cac , ndop , nakeep , nalose , ntr / nbands ,
     .                 nbands , iband , delta )
c
      enddo
c
c   The original data has now been adjusted so that there are exactly
c   2*nftotal/nbdiv points between the centers of the bands
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 3
c
c   Form blending weights.  Here we use just 0 or 1 but we could use
c   an interpolation in the overlap regions.
c
c   Frequency cells from midpoint of a band to the midpoint between
c   adjacent bands
c
      n2      = nftotal / nbdiv
c
c   Number of overlap cells at each end of a band
c
      ovrlap2 = ( ntrband / 2 ) - n2
c
c   Use a weight of 1 for the non-overlap regions
c
      wt(:,:) = 1.0
c
c   Zero the end (overlap) regions
c
      do iband = 2 , nbands        !  All except first band
c
         do j = 1 , ovrlap2
c
            wt(j,iband)           = 0.0
c
         enddo
c
      enddo
c
      do iband = 1 , nbands - 1    !  All except last band
c
         do j = 1 , ovrlap2
c
            wt(ntrband+1-j,iband) = 0.0
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 4
c
c   Correct the phase gradient in Doppler due to firing sequence of
c   bands using the center band as the reference
c
      do iband = 1 , nbands
c
         delta = ( float( iband - 1 - nbands / 2 ) -
     .             0.5 * float( 1 - mod( nbands , 2 ) ) )
     .           / float( nbands )
c
         do i = 1 + nalose , nakeep + nalose
c
            c = cexp( cmplx( 0.0 ,
     .                       pi * delta * float( i - 1 - nsa / 2 ) /
     .                                    float( nsa / 2 ) ) )
c
            do j = 1 , ntrband
c
               jj        = j + ( iband - 1 ) * ntrband
c
               cac(i,jj) = c * cac(i,jj)
c
            enddo
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 5
c
c   Compute covariance between bands for adaptive phase correction
c
      covrf(:)  = czero
c
      covrfd(:) = czero
c
      power(:)  = czero
c
      do iband = 1 , nbands - 1
c
c   Indexes of the center frequencies for the left and right bands
c
         lctr   = nctr + ( iband - 1 ) * ntrband
c
         rctr   = nctr +         iband * ntrband
c
c   Use only the middle range of Doppler cells
c
         do i = middle - dopwidth , middle + dopwidth
c
            cleft  = czero
c
            cright = czero
c
            do iovrlap = n2 - ovrlap2 / 2 , n2 + ovrlap2 / 2
c
c   Indexes of the pixels in each band.  By construction (Stage 2) the
c   bands line up at +n2 from the center of the left band and -n2 from
c   the right band.
c
               ileft        = lctr + iovrlap
c
               iright       = rctr - 2 * n2 + iovrlap
c
               cleft        = cleft  + cac(i,ileft)
c
               cright       = cright + cac(i,iright)
c
               covrf(iband) = covrf(iband) + cac(i,ileft) *
     .                                conjg( cac(i,iright) )
c
               power(iband) = power(iband) +
     .                        cac(i,ileft)  * conjg( cac(i,ileft) ) *
     .                        cac(i,iright) * conjg( cac(i,iright) )
c
            enddo
c
            covrfd(iband) = covrfd(iband) + cleft * conjg( cright )
c
         enddo
c
      enddo
c
c   Output the correlation coefficients and phases as diagnostic info
c
      ntotal = ( 1 + 2 * ( ovrlap2 / 2 ) ) * ( 1 + 2 * ( dopwidth ) )
c
c     write ( 49 , * ) ( cabs( covrf(j) ) /
c    .                   sqrt( float( ntotal ) * cabs( power(j) ) ) ,
c    .                   phasec( covrf(j) ) , j = 1 , nbands - 1 )
c
c-----------------------------------------------------------------------
c
c  Optionally, use the covariance computed from only the average values
c  over fast-time.  This emphasizes the target near the middle range at
c  the expense of the noise and other scatterers further out.
c
      if ( dcweight .ne. 0 ) then
c
         covrf(:) = covrfd(:)
c
      endif
c
c-----------------------------------------------------------------------
c
c   Implement a fading filter in time
c
      if ( isub .le. ntfilter ) then
c
         do j = 1 , nbands - 1
c
            covrfsave(j) = covrfsave(j) + covrf(j)
c
         enddo
c
      else
c
         tfade = 1.0 / float( ntfilter )
c
         do j = 1 , nbands - 1
c
            covrfsave(j) = ( 1.0 - tfade ) * covrfsave(j) +
     .                               tfade * covrf(j)
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c
c   Make unit vectors
c
      do j = 1 , nbands - 1
c
         if ( cabs( covrfsave(j) ) .gt. 0.0 ) then
c
            covrf(j) = covrfsave(j) / cabs( covrfsave(j) )
c
         else
c
            covrf(j) = czero
c
         endif
c
      enddo
c
c     write ( 49 , * ) ( cabs( covrf(j) ) ,
c    .                   phasec( covrf(j) ) , j = 1 , nbands - 1 )
c
c     write ( 49 , * )
c
      write ( 96 , '(1x,i6,2f10.3)' ) isub ,
     .              ( phasec( covrf(j) ) , j = 1 , nbands - 1 )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 6
c
c   Make the middle band the pivot, correcting F2F phase and applying
c   blending weights
c
      wtall(:)   = 0.0                !  Diagnostic - sum of all weights
c
c-----------------------------------------------------------------------
c
c   Pivot band - no correction required
c
      iband = 1 + nbands / 2
c
      do i = 1 + nalose , nakeep + nalose
c
         do j = 1 , ntrband
c
            jj             = j + ( iband - 1 ) * ntrband
c
            cac(i,jj)      = wt(j,iband) * cac(i,jj)
c
            joff           = ( iband - 1 ) * 2 * ovrlap2
c
            wtall(jj-joff) = wtall(jj-joff) + wt(j,iband)
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
c   From the pivot band to higher
c
      if ( nbands .gt. 2 ) then
c
         covcor = cmplx( 1.0 , 0.0 )
c
         do iband = ( nbands / 2 ) + 2 , nbands
c
            covcor = covcor * covrf(iband-1)
c
            do i = 1 + nalose , nakeep + nalose
c
               do j = 1 , ntrband
c
                  jj             = j + ( iband - 1 ) * ntrband
c
                  cac(i,jj)      = wt(j,iband) * cac(i,jj) * covcor
c
                  joff           = ( iband - 1 ) * 2 * ovrlap2
c
                  wtall(jj-joff) = wtall(jj-joff) + wt(j,iband)
c
               enddo
c
            enddo
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c
c   From the pivot band to lower
c
      covcor = cmplx( 1.0 , 0.0 )
c
      do iband = ( nbands / 2 ) , 1 , - 1
c
         covcor = covcor * conjg( covrf(iband) )
c
         do i = 1 + nalose , nakeep + nalose
c
            do j = 1 , ntrband
c
               jj             = j + ( iband - 1 ) * ntrband
c
               cac(i,jj)      = wt(j,iband) * cac(i,jj) * covcor
c
               joff           = ( iband - 1 ) * 2 * ovrlap2
c
               wtall(jj-joff) = wtall(jj-joff) + wt(j,iband)
c
            enddo
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c                           STAGE 7
c
c   Load output array with phase-corrected data
c
      do i = 1 + nalose , nakeep + nalose
c
         crc(:) = czero
c
         do iband = 1 , nbands
c
            joff = ( iband - 1 ) * 2 * ovrlap2
c
            do j = 1 , ntrband
c
               jj           = j + ( iband - 1 ) * ntrband
c
               crc(jj-joff) = crc(jj-joff) + cac(i,jj)
c
            enddo
c
         enddo
c
         cac(i,:) = crc(:)
c
      enddo
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      wtall(:) = wtall(:) / float( nakeep )
c
      if ( isub .eq. 1 ) then
c
c        write ( 66 , '(8f10.3)' ) wtall !  Diagnostic - sum of weights
c
         write ( 66 , * ) 'ovrlap2 , n2 , nbdiv , nbdiv * n2:'
     .                   , ovrlap2 , n2 , nbdiv , nbdiv * n2
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      subroutine ntrcomb ( nftotal , nbdiv , nbands , rfdelt , drfreq )
c**
c***********************************************************************
c**
c   Purpose:  To compute the number of adjusted frequency bins between
c             the first and last band for a step-chirp system
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nftotal , nbdiv , nbands , nfsmidg
c
      real    rfdelt , drfreq
c
c-----------------------------------------------------------------------
c
c   Outputs:

c      nftotal = Number of points between the centers of the first and
c                last band after adjustment
c
c      nbdiv   = Number of convenient half-bands to split this up into
c
c   Inputs:
c
c      nbands  = Number of frequency bands
c
c      rfdelta = Original frequency separation between bands
c
c      drfreq  = Frequency step per pixel
c
c-----------------------------------------------------------------------
c
c   Initial value of nftotal
c
      nftotal = nint( float( nbands - 1 ) * rfdelt / drfreq )
c
      nbdiv   = 2 * ( nbands - 1 )
c
c   nfsmidg = Number of original extra points - correct this by adding
c             or subtracting to get a number divisible by nbdiv
c
      nfsmidg = mod( nftotal , nbdiv )
c
c   To make the smallest correction, either add a bit or subtract a bit
c   to make nftotal divisible by nbdiv
c
      if ( nfsmidg .gt. ( nbands - 1 ) ) then
c
         nftotal = nftotal + nbdiv - nfsmidg
c
      else
c
         nftotal = nftotal - nfsmidg
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      subroutine dshift ( cac , ndop , nakeep , nalose , nt , nb , ib ,
     .                    delta )
c**
c***********************************************************************
c**
c   Purpose:  To shift an array a small amount in one coordinate by
c             interpolation
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer ndop , nakeep , nalose , nt , nb , ib , i , j , nzero ,
     .        jshift , jj
c
      real    delta , deltap
c
      complex cac(ndop,nt,nb)
c
c-----------------------------------------------------------------------
c
      if ( abs( delta ) .gt. float( nt - 2 ) ) then
c
         write ( 6 , * ) delta ,  'Delta too large in DSHIFT'
         read  ( 5 , * )
c
	endif
c
c-----------------------------------------------------------------------
c
      nzero  = 1 + ifix( abs( delta ) )
c
      jshift = nzero - 1
c
c   Shift array by replacing data.  The sign of the shift determines the
c   direction of scan.
c
      if ( delta .gt. 0.0 ) then
c
c   Shift right
c
         deltap = delta - float( jshift )
c
         do i = 1 + nalose , nakeep + nalose
c
            do j = 1 , nt - nzero
c
               jj          = j + jshift
c
               cac(i,j,ib) = ( 1.0 - deltap ) * cac(i,jj,ib) +
     .                                 deltap * cac(i,jj+1,ib)
c
            enddo
c
            do j = nt - nzero + 1 , nt
c
               cac(i,j,ib) = cmplx( 0.0 , 0.0 )
c
            enddo
c
         enddo
c
      else
c
c   Shift left
c
         deltap = delta + float( jshift )
c
         do i = 1 + nalose , nakeep + nalose
c
            do j = nt , 1 + nzero , - 1
c
               jj          = j - jshift
c
               cac(i,j,ib) = ( 1.0 + deltap ) * cac(i,jj,ib) -
     .                                 deltap * cac(i,jj-1,ib)
c
            enddo
c
            do j = 1 , nzero
c
               cac(i,j,ib) = cmplx( 0.0 , 0.0 )
c
            enddo
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c
      return
      end
