c**
c***********************************************************************
c**
      subroutine subimg0 ( mpass , isub , csbimg , nsa , nsr , ntp ,
     .                     ncr , naskip , nakeep , str , cstr , wtac ,
     .                     work , ac , cac , nwork , nerror , wtr ,
     .                     ntr , rc , crc , curtime , dopcen , ntrcom ,
     .                     nrfsmth , ndop )
c**
c***********************************************************************
c**
c   Purpose:  To return a coarse resolution complex image either by
c             computing it from the range-compressed data or by reading
c             it from disk.
c
c             The routine also handles disk buffering of the sub-image
c             array for multi-pass operation.  This is controlled by the
c             parameter 'pass', which has the following interpretation.
c
c             mpass = 0 :  Single pass method - don't buffer data
c
c             mpass = 1 :  First pass of multi-pass method.  Calculate
c                          the sub-image and write it to disk.
c
c             mpass = 2 :  Second pass of multi-pass method.  Do not
c                          calculate sub-image but read it from disk.
c
c   This version implements the Doppler FFT before the range FFT to
c   support step-chirp systems - systems that need to make corrections
c   that are functions of Doppler before the range-compression is done.
c
c-----------------------------------------------------------------------
c
      implicit none
c
      include  'sarprm.h'      !  Standard ISAR-T parameters
c
      include  'updates.h'     !  Updates to parameters
c
      integer   mpass , isub , nsa , nsr , ntp , ncr , naskip , nakeep ,
     .          ir , ia , isr , nacent , nalose , nwork , nerror , ntr ,
     .          active , iactive , ntrcom , nrfsmth , ib , ndop , ncen
c
      real      wtac(3*naskip) , str(2,ntp,ncr) , ac(2,ndop,ntr) ,
     .          wtr(ntr) , rc(2,ncr) , work(nwork) , dopcen , avgsg0 ,
     .          curtime , drfreq , trigr , ctrast
c
      complex   csbimg(nakeep,nsr) , cstr(ntp,ncr) , crc(ncr) , czero ,
     .          cac(ndop,ntr)
c
      integer   reclda , newntr , ndot
c
      character nrfile*80
c
      save      newntr
c
c-----------------------------------------------------------------------
c
      czero  = cmplx( 0.0 , 0.0 )
c
c   Sampling in the frequency domain
c
      drfreq = br * 1.0E+12 * dtr
c
c-----------------------------------------------------------------------
c
c   Compress in time to form coarse Doppler (sub-image)
c
      if ( mpass .lt. 2 ) then
c
c   For ISAR mode the sub-image is simply the Fourier transform with
c   respect to time at fixed range
c
c-----------------------------------------------------------------------
c
c   Load Doppler work array with fast-time samples versus slow-time
c
         cac(:,:) = czero
c
         do ir = 1 , ntr
c
            do ia = 1 , 3 * naskip
c
               cac(ia,ir) = cstr(ia+naskip/2,ir)
c
            enddo
c
         enddo
c
c-----------------------------------------------------------------------
c
c   Keep only the middle Doppler cells to save memory
c
         nalose = max( 0 , ( nsa - nakeep ) / 2 )
c
c   FFT in slow-time to convert to coarse Doppler
c
         nacent = ( 3 * naskip ) / 2
c
         do ir = 1 , ntr          !  All fast-time samples
c
            call wtfft ( ac(1,1,ir) , wtac , ac(1,1,ir) , work , + 1 ,
     .                   - nacent , 3 * naskip , nsa )
c
            if ( nalose .gt. 0 ) then
c
               do ia = 1 , nalose
c
                  cac(ia,ir)               = czero
c
                  cac(ia+nalose+nakeep,ir) = czero
c
               enddo
c
            endif
c
            if ( nbands .gt. 1 .and. nakeep .eq. ( nsa * nbands ) ) then
c
               do ib = 2 , nbands
c
                  do ia = 1 , nsa
c
                     cac((ib-1)*nsa+ia,ir) = cac(ia,ir)
c
                  enddo
c
               enddo
c
            endif
c
         enddo
c
c-----------------------------------------------------------------------
c
c   Compress in range.  Allow for up-chirp (br>0) or down-chirp (br<0)
c
         if ( nbands .gt. 1 ) call rfcomb0 ( cac , crc , nsa , nakeep ,
     .                             ntr , nbands , rfdelt , drfreq ,
     .                             nrfsmth , isub , ndop )
c
         do ia = 1 , nakeep
c
            crc(:) = czero
c
            do ir = 1 , ntr
c
               crc(ir) = cac(ia+nalose,ir)
c
            enddo
c
            if ( br .gt. 0.0 ) then
c
               call wtfft ( rc , wtr , rc , work , + 1 , - ntrcom / 2 ,
     .                      ntrcom , ncr )
c
            else
c
               call wtfft ( rc , wtr , rc , work , - 1 , - ntrcom / 2 ,
     .                      ntrcom , ncr )
c
            endif
c
            do isr = 1 , nsr
c
               ir             = isr + ( ncr - nsr ) / 2
c
               csbimg(ia,isr) = crc(ir)
c
            enddo
c
         enddo
c
c-----------------------------------------------------------------------
c
c   Estimate the centroid for use in adaptive motion compensation
c
         avgsg0  = 0.0
c
         ncen    = min( nakeep , nsa ) / 4
c
         do isr = 1 , nsr
c
            do ia = 1 + ( ndop / 2 ) - ncen , 1 + ( ndop / 2 ) + ncen
c
               avgsg0 = avgsg0 + csbimg(ia,isr) *
     .                    conjg( csbimg(ia,isr) )
c
            enddo
c
         enddo
c
         avgsg0  = sqrt( avgsg0 / float( nakeep * nsr ) )
c
         trigr   = 4.0 * avgsg0
c
         active  = 0
c
         iactive = 0
c
         do isr = 2 + nsr / 8 , nsr - nsr / 8
c
            do ia = 1 + ( ndop / 2 ) - ncen , 1 + ( ndop / 2 ) + ncen
c
               if ( cabs( csbimg(ia,isr) ) .gt. trigr ) then
c
                  active  = active  + 1
c
                  iactive = iactive + ia
c
               endif
c
            enddo
c
         enddo
c
         if ( active .eq. 0 ) then
c
            dopcen = 0.0
c
         else
c
            dopcen = dfc * ( ( float( iactive ) / float( active ) ) -
     .                       float( nakeep / 2 + 1 ) )
c
         endif
c
      endif  !  mpass .lt. 2
c
c-----------------------------------------------------------------------
c
c   Write sub-image to disk file or read it if required
c
      if ( mpass .eq. 1 ) then
c
         write ( 86 , rec = isub , err = 10 ) csbimg
c
      else if ( mpass .eq. 2 ) then
c
         read  ( 86 , rec = isub , err = 10 ) csbimg
c
      endif
c
      nerror = 0
c
      return
c
   10 nerror = 1
c
      return
      end
