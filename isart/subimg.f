c**
c***********************************************************************
c**
      subroutine subimg  ( mpass , isub , csbimg , nsa , nsr , ntp ,
     .                     ncr , naskip , nakeep , str , cstr , wtac ,
     .                     work , ac , cac , nwork , nerror )
c**
c***********************************************************************
c**
c   Purpose:  To return a coarse resolution complex image either by
c             computing it from the range-compressed data or by reading
c             it from disk.
c
c             There are two algorithms for producing the sub-image:
c
c             mode = 1  :  Normal ISAR mode - simply FFT the data at
c                          fixed range
c
c             mode = 2  :  Ultra-Wide Band mode - do an accurate
c                          direct integration on a linear path in
c                          range versus time
c
c             mode = 3  :  Telephonics Stripmap
c
c             mode = 4  :  Telephonics Stripmap - realtime mode
c
c             The routine also handles disk buffering of the sub-image
c             array for multi-pass operation.  This is controlled by the
c             parameter 'mpass', which has the following interpretation.
c
c             mpass = 0 :  Single pass method - don't buffer data
c
c             mpass = 1 :  First pass of multi-pass method.  Calculate
c                          the sub-image and write it to disk.
c
c             mpass = 2 :  Second pass of multi-pass method.  Do not
c                          calculate sub-image but read it from disk.
c
c-----------------------------------------------------------------------
c
      implicit none
c
      include 'sarprm.h'      !  Standard ISAR-T parameters
c
      include 'updates.h'     !  Updates to parameters
c
      integer  mpass , isub , nsa , nsr , ntp , ncr , naskip , nakeep ,
     .         ir , ia , isr , nacent , nalose , nwork , nerror
c
      real     wtac(3*naskip) , str(2,ntp,ncr) , ac(2,nsa) , work(nwork)
c
      complex  csbimg(nakeep,nsr) , cstr(ntp,ncr) , cac(nsa)
c
c-----------------------------------------------------------------------
c
c   Compress in time to form coarse Doppler (sub-image)
c
      if ( mpass .lt. 2 ) then
c
         nacent = ( 3 * naskip ) / 2
c
         if ( uwb .eq. 0 ) then
c
c   For ISAR mode the sub-image is simply the Fourier transform with
c   respect to time at fixed range
c
            do isr = 1 , nsr          !  All range cells in sub-image
c
               ir = isr + ( ncr - nsr ) / 2
c
               call wtfft ( str(1,1+naskip/2,ir) , wtac , ac , work ,
     .                      + 1 , - nacent , 3 * naskip , nsa )
c
c   Keep only the middle Doppler cells to save memory
c
               nalose = ( nsa - nakeep ) / 2
c
               do ia = 1 , nakeep
c
                  csbimg(ia,isr) = cac(ia+nalose)
c
               enddo
c
            enddo
c
         else
c
c   For Ultra-Wide-Band strip-map mode the sub-image formation must
c   account for range walk and quadratic focus
c
            call uwbsub ( csbimg , nakeep , nsr , drs , drc , lambda ,
     .                    slant0 , cstr , ncr , dfc , naskip , wtac ,
     .                    dtp , vfocus , work , cac , br , strtch ,
     .                    clight )
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c
c   Write sub-image to disk file or read it if required
c
      if      ( mpass .eq. 1 .and. isub .gt. 0 ) then
c
         write ( 86 , rec = isub )            csbimg
c
      else if ( mpass .eq. 2 .and. isub .gt. 0 ) then
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
