c
c-----------------------------------------------------------------------
c
c
c             IMGENLIB - Complex Image Formation and Modification
c
c
c      subroutine imgeni ( ciorg , csbimg , nfa , nfr , mpr , cw , rw ,
c      subroutine fliplr ( ci , nfa , nfr )
c      subroutine fliplrb ( b , nfa , nfr )
c      subroutine flipdu ( ci , nfa , nfr )
c      subroutine flipdub ( ci , nfa , nfr )
c      subroutine fliptr ( ci , co , nfa , nfr )
c      subroutine cutout ( ifile , ofile , nxin , nyin , nxout , nyout ,
c      subroutine addtwo ( fil1 , fil2 , ofile , nx1 , ny1 , nx2 , ny2 ,
c      subroutine sumdiff ( fil1 , fil2 , ofile , nx1 , ny1 , nx2 , ny2
c      subroutine upsamp ( ifile , ofile , nxin , nyin , nxout , nyout ,
c      subroutine fftup ( r , nin , nout , work , nwork )
c      subroutine flipit ( ifile , ofile , nxin , nyin , cin , option ,
c      subroutine fftimg ( ifile , ofile , nx , ny , nz , cin , rin ,
c      subroutine fft1d ( cd , nn , ndim , isign , iform , work , nwork
c      subroutine fft2d ( cd , nn , ndim , isign , iform , work , nwork
c      subroutine qswap ( cin , cout , nx , ny )
c      subroutine ambig ( ifile , ofile , nx , ny , cin , cout , rout ,
c      subroutine ciflat ( c , nx , ny , w , eps )
c      subroutine removewt ( ifile , ofile , c , nx , ny , w , eps ,
c      subroutine diff_images ( c1 , c2 , csum , cdif , nx , ny , phlag
c      subroutine csumdif ( c1 , c2 , csum , cdif , cmin , nx , ny ,
c      subroutine addimg ( ci , cibig , nxci , nyci , nxcbig , nycbig ,
c
c-----------------------------------------------------------------------
c
c
c                        ISARTLIB - General SAR
c
c
c                        Subroutines
c
c      subroutine accel ( z , n , b , dt , a )
c      subroutine cshift ( cin , n , nshift , cout )
c      subroutine dchirp ( z , n , bavg , vel , sg0 , d0 , dd , dt ,
c      subroutine hwindo ( w , n , nt , nf , work )
c      subroutine onelag ( c , n , dt , nlag , dc , dw , d , pwr , phi )
c      subroutine rshift ( rin , n , nshift , rout )
c      subroutine sort ( ra , n , index )
c      subroutine taylor ( nx , sll_db , weights )
c      subroutine rcomp ( br , clight , dr0 , rc , crc , pt , ntr , ncr ,
c      subroutine chirp_compress ( rawdat , compressed , chirp , ntr ,
c      subroutine iqstats ( rdata , ntr , imean , qmean , iqcor ,
c
c
c                        Functions
c
c      real function ctrast ( ci , nx , ny )
c      real function snr ( sg0 , dw , dt , n0 )
c      complex function cpfast ( p )
c      integer function twopwr ( k )
c      real function mflops ( work , nft )
c      real function alog2 ( x )
c      real function rtrast ( ri , nx , ny )
c
c-----------------------------------------------------------------------
c
c
c                        IOLIB - Input/Output Routines
c
c
c                        Subroutines
c
c      subroutine readci ( ci , nfa , nfrb , ifile , ftype , xoffset ,
c      subroutine iocimg ( lunit , irec0 , nfa , nfr , option , ci , ne
c      subroutine getfrm ( iframe , pframe , ci , nfa , nfr , mode ,
c      subroutine oframe ( ci , mi , bi , iframe , ifile , nx , ny ,
c      subroutine outdbs ( ofile , image , imageb , mxpixx , mxpixy ,
c      subroutine dabyte ( option , array , ni , nj , ofile , nerror )
c      subroutine apack6 ( rdwrt , r , npts , iofile , nerror ,
c      subroutine putiq4 ( c , nx , ny , i2_b , ofile , kframe )
c      subroutine put16bit ( r , nx , i2_b , scale , outunit , recnum )
c      subroutine sun_flt ( nx , ny , work , infile , outfile )
c      subroutine numfile ( infile , ifile , k )
c      subroutine decode_double ( raw , float8 )
c      subroutine sun2flt ( ibyt , obyt , n )
c      subroutine hex2int ( hex , nhex , a )
c
c
c                        Functions
c
c      integer function ibmi2 ()
c      integer function lastdot ( ch80var )
c
c-----------------------------------------------------------------------
c
C**
C***********************************************************************
C**
      subroutine imgeni ( ciorg , csbimg , nfa , nfr , mpr , cw , rw ,
     .                    rwt , nsr , double )
C**
C***********************************************************************
C**
c
c   This routine converts a fine resolution SAR image (doppler,range)
c   into a three-dimensional complex field (doppler,range,time).  This
c   process can be thought of as the inverse of the image focusing
c   routine, IMGEN, part of the RDRTec/Telephonics Inverse-SAR processor.
c
c   Arguments:
c
c      ciorg(nfa,nfr)             :  Input fine resolution image
c
c      csbimg(1+nfa/mpr,nfr,2*mpr):  Output freq-range-time field
c
c      nfa                        :  No. of input doppler bins
c
c      nfr                        :  No. of input and output range cells
c
c      mpr                        :  Master particle ratio, the ratio
c                                    of input fine resolution cells to
c                                    the output coarse resolution cells
c
c      cw,rw                      :  Work space - equivalenced to real
c                                    and complex
c
      implicit none
c
      integer nfa , nfr , mpr , nn(1) , imp , nmp , ictr , ipos , ineg ,
     .        ifa , i , j , it , js , nsr
c
      complex ciorg(nfa,nfr) , csbimg(1+nfa/mpr,nsr,2*mpr) , pcor
c
      real    rwt(1+4*mpr) , fmp , ltime , twopi
c
      logical double
c
c-----------------------------------------------------------------------
c
c   cw and rw are complex and real names for the same work storage area
c
c   The usage of this space is:
c
c        cw(1-->2*mpr)          :  Complex array loaded with a small
c                                  region of the complex image so that
c                                  it can be FFT'd
c
c        rw(remainder)          :  Work space passed to subroutine FOURT
c
c
      complex cw(4*mpr)
c
      real    rw(2,4*mpr)
c
c-----------------------------------------------------------------------
c
      twopi = 2.0 * atan2( 0.0 , - 1.0 )
c
      nmp   = 1 + nfa / mpr
c
      nn(1) = 2 * mpr
c
c   Load linear weight array
c
      do i = 1 , mpr + 1
c
c   Commented out option to weight the image points
c
c        rwt(i)         = float( mpr - i + 1 ) / float( mpr )
c
         rwt(i)         = 1.0
         rwt(2+2*mpr-i) = rwt(i)
c
      enddo
c
c   Loop over range lines
c
      do j = 1 , nfr
c
c   If doubling is on then load the sub-image array into alternate rows
c   so that interpolation can be used to fill the others
c
         if ( .not. double ) then
            js = j
         else
            js = 2 * j - 1
         endif
c
c   Loop over master particles
c
         do imp = 1 , nmp
c
            ictr = 1 + ( imp - 1 ) * mpr
c
            do it = 1 , 2 * mpr
c
               cw(it) = cmplx( 0.0 , 0.0 )
c
            enddo
c
c   Use the linear interpolation weights to write these Fourier
c   coefficients into the fine resolution image array
c
            if ( imp .ne. 1 ) then
c
c   First, for all but the most negative frequency extract the negative
c   Fourier component from the image frame
c
               do ifa = 2 , mpr
c
                  ineg            = ictr + 1 - ifa
                  cw(2*mpr+2-ifa) = ciorg(ineg,j) * rwt(2*mpr+2-ifa)
c
               enddo
c
            endif
c
c   Next, for all but the most positive frequency extract the positive
c   Fourier component from the image frame
c
            if ( imp .ne. nmp ) then
c
               do ifa = 1 , mpr
c
                  ipos    = ictr - 1 + ifa
                  cw(ifa) = ciorg(ipos,j) * rwt(ifa)
c
               enddo
c
            endif
c
            call fourt ( rw , nn , 1 , - 1 , 1 , rwt(2+2*mpr) ,
     .                   2 * mpr )
c
            do it = 1 , mpr
c
               csbimg(imp,js,it+mpr) = cw(it)
               csbimg(imp,js,it)     = cw(it+mpr)
c
            enddo
c
c   Zero Nyquist frequency
c
            csbimg(imp,js,1) = cmplx( 0.0 , 0.0 )
c
c   Correct phase so that it is consistent with IMGEN definition
c
            fmp = float( imp - nmp / 2 - 1 )
c
            do it = 2 , 2 * mpr
c
               ltime             = 0.5 * float( it - mpr - 1 )
               pcor              = ( 1.0 / float( 2 * mpr ) ) * cexp(
     .                           cmplx( 0.0 , - twopi * fmp * ltime ) )
               csbimg(imp,js,it) = pcor * csbimg(imp,js,it)
c
            enddo
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Option to double the size of the output array by interpolating the
c   even rows from the odd ones
c
      if ( double ) then
c
c   Fill the even rows by interpolating the odd rows
c
         do it  = 1 , 2 * mpr
c
            do imp = 1 , nmp
c
c   Fill the last row with the last odd row
c
               csbimg(imp,nsr,it)   = csbimg(imp,nsr-1,it)
c
c  Two-point interpolation for the next-to-last even row
c
               csbimg(imp,nsr-2,it) = 0.5 * ( csbimg(imp,nsr-1,it)
     .                                      + csbimg(imp,nsr-3,it) )
c
c  Two-point interpolation for the second row
c
               csbimg(imp,2,it)     = 0.5 * ( csbimg(imp,1,it)
     .                                      + csbimg(imp,3,it) )
c
c   Four-point interpolation for all other even rows
c
               do js = 4 , nsr - 4 , 2
c
                  csbimg(imp,js,it) = (
     .                               9.0 * ( csbimg(imp,js+1,it) +
     .                                       csbimg(imp,js-1,it) )
     .                                   - ( csbimg(imp,js+3,it) +
     .                                       csbimg(imp,js-3,it) )
     .                             ) / 16.0
c
               enddo
c
            enddo
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine addbyt ( r , b , nx , ny , iflag )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nx , ny , iflag , i , j
c
      real      r(nx,ny)
c
      character b(nx,ny)*1
c
      if ( iflag .eq. 0 ) then
c
         do j = 1 , ny
c
            do i = 1 , nx
c
               r(i,j) = 0.0
c
            enddo
c
         enddo
c
      else if ( iflag .eq. 1 ) then
c
         do j = 1 , ny
c
            do i = 1 , nx
c
               r(i,j) = r(i,j) + float( ichar( b(i,j) ) )
c
            enddo
c
         enddo
c
c
      else if ( iflag .eq. 2 ) then
c
         do j = 1 , ny
c
            do i = 1 , nx
c
               r(i,j) = r(i,j) +
     .                10.0 ** ( 0.025 * float( ichar( b(i,j) ) - 64 ) )
c
            enddo
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine fliplr ( ci , nfa , nfr )
C**
C***********************************************************************
C**
      implicit none
c
      integer      nfa , nfr , i , j
c
      complex      ci(nfa,nfr) , ctmp
c
      do j = 1 , nfr
c
         do i = 1 , nfa / 2
c
            ctmp          = ci(nfa+1-i,j)
            ci(nfa+1-i,j) = ci(i,j)
            ci(i,j)       = ctmp
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine fliplrb ( b , nfa , nfr )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nfa , nfr , i , j
c
      character b(nfa,nfr)*1 , btmp*1
c
      do j = 1 , nfr
c
         do i = 1 , nfa / 2
c
            btmp          = b(nfa+1-i,j)
            b(nfa+1-i,j)  = b(i,j)
            b(i,j)        = btmp
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine fliptr ( ci , co , nfa , nfr )
C**
C***********************************************************************
C**
c     Simple transpose of a complex array
c
      implicit none
c
      integer      nfa , nfr , i , j , ij , ijtr
c
      complex      ci(nfa*nfr) , co(nfa*nfr)
c
      do j = 1 , nfr
c
         do i = 1 , nfa
c
            ij       = ( j - 1 ) * nfa + i
            ijtr     = ( i - 1 ) * nfr + j
c
            co(ijtr) = ci(ij)
c
         enddo
c
      enddo
c
c   Transfer back to original array
c
      do j = 1 , nfr
c
         do i = 1 , nfa
c
            ij     = ( j - 1 ) * nfa + i
c
            ci(ij) = co(ij)
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine flipdu ( ci , nfa , nfr )
C**
C***********************************************************************
C**
      implicit none
c
      integer      nfa , nfr , i , j
c
      complex      ci(nfa,nfr) , ctmp
c
      do i = 1 , nfa
c
         do j = 1 , nfr / 2
c
            ctmp          = ci(i,nfr+1-j)
            ci(i,nfr+1-j) = ci(i,j)
            ci(i,j)       = ctmp
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine flipdub ( b , nfa , nfr )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nfa , nfr , i , j
c
      character b(nfa,nfr)*1 , btmp*1
c
      do i = 1 , nfa
c
         do j = 1 , nfr / 2
c
            btmp         = b(i,nfr+1-j)
            b(i,nfr+1-j) = b(i,j)
            b(i,j)       = btmp
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine flipit ( ifile , ofile , nxin , nyin , cin , option ,
     .                    cout , mode )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nxin , nyin , i , j , option , mode , reclda
c
      character ifile*80 , ofile*80
c
      complex   cin(nxin,nyin) , cout(nxin,nyin)
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nxin ) )
c
         do j = 1 , nyin
c
            read ( 8 , rec = j ) ( cin(i,j) , i = 1 , nxin )
c
         enddo
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c
      if ( mode .eq. 1 ) then
c
         call fliplr ( cin , nxin , nyin )         !  Left/Right Flip
c
      else if ( mode .eq. 2 ) then
c
         call flipdu ( cin , nxin , nyin )         !  Down/Up Flip
c
      else
c
         call fliptr ( cin , cout , nxin , nyin )  !  Transpose
c
      endif
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nxin ) )
c
         do j = 1 , nyin
c
            write ( 8 , rec = j ) ( cin(i,j) , i = 1 , nxin )
c
         enddo
c
         close ( 8 )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine cutout ( ifile , ofile , nxin , nyin , nxout , nyout ,
     .                    imin , jmin , itot , jtot , cin , cout , opt )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nxin , nyin , nxout , nyout , imin , jmin , itot ,
     .          jtot , nxlim , nylim , i , j , reclda , jrec , opt
c
      character ifile*80 , ofile*80
c
      complex   cin(nxin,nyin) , cout(nxout,nyout) , cw(nxin)
c
c-----------------------------------------------------------------------
c
      if ( opt .eq. 1 ) then
c
         close ( 8 )
c
         open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nxin ) ,
     .              status = 'OLD' )
c
      endif
c
c-----------------------------------------------------------------------
c
c   Initialize with zeroes
c
      do j = 1 , nyout
c
         do i = 1 , nxout
c
            cout(i,j) = cmplx( 0.0 , 0.0 )
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
      if ( imin .gt. nxin .or. jmin .gt. nyin ) then
c
         write ( 6 , * ) ' Imin or Jmin outside array --> zero-filling'
c
      else if ( nxout .eq. itot .and. nyout .eq. jtot ) then
c
         nxlim = min( nxout , nxin + 1 - imin )
         nylim = min( nyout , nyin + 1 - jmin )
c
         do j = 1 , nylim
c
            jrec = j + jmin - 1
c
            if ( opt .eq. 1 ) then
c
               read ( 8 , rec = jrec ) cw
c
            else
c
               cw(:) = cin(:,jrec)
c
            endif
c
            do i = 1 , nxlim
c
               cout(i,j) = cw(i+imin-1)
c
            enddo
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c
      if ( opt .eq. 1 ) then
c
         close ( 8 )
c
         open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nxout ) )
c
         do j = 1 , nyout
c
            write ( 8 , rec = j ) ( cout(i,j) , i = 1 , nxout )
c
         enddo
c
         close ( 8 )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine addtwo ( fil1 , fil2 , ofile , nx1 , ny1 , nx2 , ny2 ,
     .                    xoffset , yoffset , c1 , c2 , wgt , option )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nx1 , ny1 , nx2 , ny2 , i , j , ii , jj , option ,
     .          reclda , xoffset , yoffset
c
      character fil1*80 , fil2*80 , ofile*80
c
      complex   c1(nx1,ny1) , c2(nx2,ny2)
c
      real      wgt , wgt1 , wgt2 , awt
c
c-----------------------------------------------------------------------
c
c   If OPTION = 0, use memory; = 1, read and write to disk
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = fil1 , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx1 ) )
c
         do j = 1 , ny1
c
            read ( 8 , rec = j ) ( c1(i,j) , i = 1 , nx1 )
c
         enddo
c
         close ( 8 )
c
         open ( 8 , file = fil2 , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx2 ) )
c
         do j = 1 , ny2
c
            read ( 8 , rec = j ) ( c2(i,j) , i = 1 , nx2 )
c
         enddo
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c
      wgt1 = 0.0
      wgt2 = 0.0
c
      do j = 1 , ny2
c
         do i = 1 , nx2
c
            ii = i + xoffset
            jj = j + yoffset
c
            if ( ii .gt. 0 .and. ii .le. nx1 .and.
     .           jj .gt. 0 .and. jj .le. ny1 ) then
c
               wgt1 = wgt1 + cabs( c1(ii,jj) ) ** 2
               wgt2 = wgt2 + cabs( c2(i,j) )   ** 2
c
            endif
c
         enddo
c
      enddo
c
      if ( wgt2 .ne. 0.0 ) then
c
         awt = wgt * sqrt( wgt1 / wgt2 )
c
      else
c
         awt = 0.0
c
      endif
c
      do j = 1 , ny2
c
         do i = 1 , nx2
c
            ii = i + xoffset
            jj = j + yoffset
c
            if ( ii .gt. 0 .and. ii .le. nx1 .and.
     .           jj .gt. 0 .and. jj .le. ny1 ) c1(ii,jj) =
     .                                         c1(ii,jj) + awt * c2(i,j)
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx1 ) )
c
         do j = 1 , ny1
c
            write ( 8 , rec = j ) ( c1(i,j) , i = 1 , nx1 )
c
         enddo
c
         close ( 8 )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine sumdiff ( fil1 , fil2 , ofile , nx , ny , c1 , c2 ,
     .                     csum , cdif , iq8bit , option )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nx , ny , i , j , option , reclda , ndot , lastdot
c
      character fil1*80 , fil2*80 , ofile*80 , sumfile*80 , diffile*80 ,
     .          logfil*80
c
      character iq8bit(nx,ny)*1
c
      complex   c1(nx,ny) , c2(nx,ny) , csum(nx,ny) , cdif(nx,ny)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Maximum lag for correlating the two images
c
      integer      maxlag
c
      parameter  ( maxlag = 1 )
c
c   Arrays for storing sum and difference image statistics
c
      real         rcovs(1+2*maxlag,1+2*maxlag) ,
     .             rcovd(1+2*maxlag,1+2*maxlag)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Cancel option parameter
c
      integer      cancel_option
c
      parameter  ( cancel_option = - 2 )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Lag used for calculating the phase gradient
c
      integer      phlag
c
      parameter  ( phlag = 128 )
c
c-----------------------------------------------------------------------
c
c   If OPTION = 0, use memory; = 1, read and write to disk
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = fil1 , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
         do j = 1 , ny
c
            read ( 8 , rec = j ) ( c1(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
         open ( 8 , file = fil2 , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
         do j = 1 , ny
c
            read ( 8 , rec = j ) ( c2(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c
      ndot   = lastdot( ofile ) - 1
c
      logfil  = ofile(1:ndot) // '.cancel.txt'
      open ( 7  , file = logfil , status = 'unknown' ,
     .            form = 'formatted' )
c
      call diff_images ( c1 , c2 , csum , cdif , nx , ny , phlag ,
     .                   maxlag , cancel_option , rcovs , rcovd ,
     .                   iq8bit )
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         ndot    = lastdot( ofile ) - 1
c
         sumfile = ofile(1:ndot) // '.sum.flt'
c
         open ( 8 , file = sumfile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
         do j = 1 , ny
c
            write ( 8 , rec = j ) ( csum(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
         diffile = ofile(1:ndot) // '.dif.flt'
c
         open ( 8 , file = diffile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
         do j = 1 , ny
c
            write ( 8 , rec = j ) ( cdif(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine upsamp ( ifile , ofile , nxin , nyin , nxout , nyout ,
     .                    cin , cout , cac , ac , work , nwork ,
     .                    option )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nxin , nyin , nxout , nyout , nwork , i , j ,
     .          option , reclda
c
      character ifile*80 , ofile*80
c
      complex   cin(nxin,nyin) , cout(nxout,nyout)
c
      complex   cac(nwork)
      real      ac(2,nwork) , work(nwork)
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nxin ) )
c
         do j = 1 , nyin
c
            read ( 8 , rec = j ) ( cin(i,j) , i = 1 , nxin )
c
         enddo
c
         call remove_phase_gradient_2 ( cin , nxin , nyin )
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c
c   Initialize with zeroes
c
      do j = 1 , nyout
c
         do i = 1 , nxout
c
            cout(i,j) = cmplx( 0.0 , 0.0 )
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
      if ( nxout .lt. nxin .or. nyout .lt. nyin ) then
c
         write ( 6 , * ) ' Output smaller than input --> zero-filling'
c
      else
c
c   First, up-sample in y
c
         do i = 1 , nxin
c
            do j = 1 , nyin
c
               cac(j) = cin(i,j)
c
            enddo
c
            if ( nyout .gt. nyin ) then
c
               call fftup ( ac , nyin , nyout , work , nwork )
c
            endif
c
            do j = 1 , nyout
c
               cout(i,j) = cac(j)
c
            enddo
c
         enddo
c
         if ( nxout .gt. nxin ) then
c
            do j = 1 , nyout
c
               do i = 1 , nxin
c
                  cac(i) = cout(i,j)
c
               enddo
c
               call fftup ( ac , nxin , nxout , work , nwork )
c
               do i = 1 , nxout
c
                  cout(i,j) = cac(i)
c
               enddo
c
            enddo
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nxout ) )
c
         do j = 1 , nyout
c
            write ( 8 , rec = j ) ( cout(i,j) , i = 1 , nxout )
c
         enddo
c
         close ( 8 )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine fftup ( r , nin , nout , work , nwork )
C**
C***********************************************************************
C**
      implicit none
c
      integer nin , nout , i , nwork , nn(1)
c
      real    r(2,nout) , scale
c
      real    work(nwork)
c
      nn(1) = nin
      call fourt ( r , nn , 1 , - 1 , 1 , work , nwork )
c
      scale = 1.0 / float( nin )
c
      do i = 1 , nin / 2
c
         r(1,i)        = scale * r(1,i)
         r(2,i)        = scale * r(2,i)
c
         r(1,nout+1-i) = scale * r(1,nin+1-i)
         r(2,nout+1-i) = scale * r(2,nin+1-i)
c
      enddo
c
      if ( mod( nin , 2 ) .eq. 1 ) then
c
         r(1,1+nin/2)  = scale * r(1,1+nin/2)
         r(2,1+nin/2)  = scale * r(2,1+nin/2)
c
      endif
c
      do i = 1 + nin / 2 + mod( nin , 2 ) , nout - nin / 2
c
         r(1,i) = 0.0
         r(2,i) = 0.0
c
      enddo
c
c   Divide the Nyquist between the + and - frequencies in the new series
c
      if ( mod( nin , 2 ) .eq. 0 ) then
c
         r(1,nout+1-nin/2) = 0.5 * r(1,nout+1-nin/2)
         r(2,nout+1-nin/2) = 0.5 * r(2,nout+1-nin/2)
c
         r(1,1+nin/2)      = r(1,nout+1-nin/2)
         r(2,1+nin/2)      = r(2,nout+1-nin/2)
c
      endif
c
      nn(1) = nout
      call fourt ( r , nn , 1 , + 1 , 1 , work , nwork )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine qswap ( cin , cout , nx , ny )
C**
C***********************************************************************
C**
c
c   Note:  Unless the compiler is told to assume arguments can be
c   aliased then cin and cout must be separate arrays.
c
      implicit none
c
      integer i , j , nx , ny
c
      complex cin(nx,ny) , cout(nx,ny) , ctemp
c
      do j = 1 , ny / 2
c
         do i = 1 , nx / 2
c
            ctemp               = cin(nx/2+i,ny/2+j)
            cout(nx/2+i,ny/2+j) = cin(i,j)
            cout(i,j)           = ctemp
c
            ctemp               = cin(i,ny/2+j)
            cout(i,ny/2+j)      = cin(nx/2+i,j)
            cout(nx/2+i,j)      = ctemp
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine fftimg ( ifile , ofile , nx , ny , nz , cin , rin ,
     .                    work , ctemp , rtemp , option , whichxy )
C**
C***********************************************************************
C**
      implicit none
c
      character ifile*80 , ofile*80
c
      character whichxy*1
c
      integer   nx , ny , nz , option , nn(2) , reclda , i , j , k , jj
c
      complex   cin(nx,ny) , ctemp(ny)
c
      real      rin(2,nx,ny) , work(2*(nx+ny)) , rtemp(2,ny)
c
c-----------------------------------------------------------------------
c
      do k = 1 , nz
c
         if ( option .eq. 1 .and.
     .        ( whichxy .eq. 'b' .or. whichxy .eq. 'B' .or.
     .          whichxy .eq. 'y' .or. whichxy .eq. 'Y' ) ) then
c
c   Read the 2-D layer from disk
c
            close ( 8 )
c
            open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .                 access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
            do j = 1 , ny
c
               jj = j + ( k - 1 ) * ny
c
               read ( 8 , rec = jj ) ( cin(i,j) , i = 1 , nx )
c
            enddo
c
c           call remove_phase_gradient_2 ( cin , nx , ny )
c
            close ( 8 )
c
         endif
c
c-----------------------------------------------------------------------
c
         nn(1) = nx
         nn(2) = ny
c
         if ( whichxy .eq. 'b' .or. whichxy .eq. 'B' ) then
c
            call fft2d ( cin , nn , 2 , - 1 , 1 , work ,
     .                   2 * ( nx + ny ) )
c
         else if ( whichxy .eq. 'y' .or. whichxy .eq. 'Y' ) then
c
            do i = 1 , nx
c
               do j = 1 , ny , 2
c
                  ctemp(2*j-1) = cin(i,2*j-1)
c
                  ctemp(j)     = - cin(i,j)
c
               enddo
c
               nn(1) = ny
c
               call fourt ( rtemp , nn , 1 , - 1 , 1 , work ,
     .                      2 * ny )
c
               do j = 1 , ny , 2
c
                  cin(i,2*j-1) = ctemp(2*j-1)
c
                  cin(i,j)     = - ctemp(j)
c
               enddo
c
            enddo
c
         else if ( whichxy .eq. 'x' .or. whichxy .eq. 'X' ) then
c
            close ( 8 )
c
            open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .                 access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
            open ( 9 , file = ofile , form = 'UNFORMATTED' ,
     .                 access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
            do j = 1 , ny
c
               jj = j + ( k - 1 ) * ny
c
               read ( 8 , rec = jj ) ( cin(i,1) , i = 1 , nx )
c
               if ( whichxy .eq. 'x' ) then
c
                  do i = 1 , nx , 2
c
                     rin(1,i,1) = - rin(1,i,1)
c
                     rin(2,i,1) = - rin(2,i,1)
c
                  enddo
c
               endif
c
               call fourt ( rin(1,1,1) , nn , 1 , - 1 , 1 , work ,
     .                      2 * nx )
c
               do i = 1 , nx , 2
c
                  rin(1,i,1) = - rin(1,i,1)
c
                  rin(2,i,1) = - rin(2,i,1)
c
               enddo
c
               jj = j + ( k - 1 ) * ny
c
               write ( 9 , rec = jj ) ( cin(i,1) , i = 1 , nx )
c
            enddo
c
            close ( 8 )
c
            close ( 9 )
c
         endif
c
c-----------------------------------------------------------------------
c
         if ( option .eq. 1 .and.
     .        ( whichxy .eq. 'b' .or. whichxy .eq. 'B' .or.
     .          whichxy .eq. 'y' .or. whichxy .eq. 'Y' ) ) then
c
c   Put the 2-D layer back on disk
c
            open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .                 access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
            do j = 1 , ny
c
               jj = j + ( k - 1 ) * ny
c
               write ( 8 , rec = jj ) ( cin(i,j) , i = 1 , nx )
c
            enddo
c
            close ( 8 )
c
         endif
c
      enddo
c
      return
      end
c**
c***********************************************************************
c**
      subroutine fft1d ( cd , nn , ndim , isign , iform , work , nwork )
c**
c***********************************************************************
c**
      implicit none
c
      integer ndim , nn(ndim) , isign , iform , nwork , i
c
      complex cd(nn(1))
c
      real    work(nwork)
c
      if ( ndim .eq. 1 ) then
c
         do i = 2 , nn(1) , 2
c
            cd(i) = - cd(i)
c
         enddo
c
      endif
c
      call fourt ( cd , nn , ndim , isign , iform , work , nwork )
c
      if ( ndim .eq. 1 ) then
c
         do i = 2 , nn(1) , 2
c
            cd(i) = - cd(i)
c
         enddo
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      subroutine fft2d ( cd , nn , ndim , isign , iform , work , nwork )
c**
c***********************************************************************
c**
      implicit none
c
      integer ndim , nn(ndim) , isign , iform , nwork , i , j
c
      complex cd(nn(1),nn(2))
c
      real    work(nwork)
c
      if ( ndim .eq. 2 ) then
c
         do j = 1 , nn(2)
c
            do i = 1 , nn(1)
c
               if ( mod( i + j , 2 ) .eq. 0 ) cd(i,j) = - cd(i,j)
c
            enddo
c
         enddo
c
      endif
c
      call fourt ( cd , nn , ndim , isign , iform , work , nwork )
c
      if ( ndim .eq. 2 ) then
c
         do j = 1 , nn(2)
c
            do i = 1 , nn(1)
c
               if ( mod( i + j , 2 ) .eq. 0 ) cd(i,j) = - cd(i,j)
c
            enddo
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine ambig ( ifile , ofile , nx , ny , cin , cout , cwork ,
     .                   rout , work , bytes , lag1 , lag2 , option )
C**
C***********************************************************************
C**
      implicit none
c
      character ifile*80 , ofile*80 , phasefil*80
c
      integer   nx , ny , option , nn(2) , reclda , i , j , k , k2 ,
     .          lag1 , lag2 , iph , kk
c
      complex   cin(nx,ny) , cout(nx,ny) , cwork(nx,ny)
c
      character bytes(nx,ny)*1
c
      real      rout(2,nx,ny) , work(nx,ny) , ph , pi
c
c-----------------------------------------------------------------------
c
      if ( option .eq. 1 ) then
c
         open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
         do j = 1 , ny
c
            read ( 8 , rec = j ) ( cin(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c
      pi = atan2( 0.0 , - 1.0 )
c
      do k = lag1 , lag2
c
         kk = k + 1 - lag1
c
         k2 = k / 2
c
         do j = 1 , ny
c
            do i = 1 , nx
c
               if ( i .gt. k2 .and. i .le. nx - k + k2 ) then
c
                  cout(i,j) = cin(i-k2,j) * conjg( cin(i+k-k2,j) )
c
               else
c
                  cout(i,j) = cmplx( 0.0 , 0.0 )
c
               endif
c
            enddo
c
         enddo
c
         ofile = 'covar'
         call oframe ( cout , work , bytes , - kk , ofile , nx , ny ,
     .                 0 , 0 , 0.25 , 0 , 0 , 0 , 0 , 0.0 , 0 , 0 )
c
         do j = 1 , ny
c
            do i = 1 , nx
c
               if ( cout(i,j) .eq. cmplx( 0.0 , 0.0 ) ) then
c
                  bytes(i,j) = char( 0 )
c
               else
c
                  ph         = atan2( aimag( cout(i,j) ) ,
     .                                 real( cout(i,j) ) )
                  iph        = min( 255 , max( 0 ,
     .                           128 + int( ph * 128 / pi ) ) )
                  bytes(i,j) = char( iph )
c
               endif
c
            enddo
c
         enddo
c
         phasefil = 'phase' // '.rmd'
c
         open ( unit = 8 , file = phasefil , form = 'unformatted' ,
     .          status = 'unknown' , access = 'direct' ,
     .          recl = reclda( nx , 0 ) )
c
         do j = 1 , ny
c
            write ( 8 , rec = j + ny * ( kk - 1 ) )
     .            ( bytes(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
         nn(1) = nx
         nn(2) = ny
c
         call fourt ( rout , nn , 2 , - 1 , 1 , work , nx * ny )
c
         cwork(:,:) = cout(:,:)
c
         call qswap ( cwork , cout , nx , ny )
c
         ofile = 'ambig'
         call oframe ( cout , work , bytes , - kk , ofile , nx , ny ,
     .                 0 , 0 , 0.25 , 0 , 0 , 0 , 0 , 0.0 , 0 , 0 )
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine ciflat ( c , nx , ny , w , eps )
C**
C***********************************************************************
C**
      implicit none
c
      integer nx , ny , i , j
c
      complex c(nx,ny)
c
      real    eps , w(nx*ny) , a , wcen
c
c   Zero weights
c
      do j = 1 , ny + nx
c
         w(j) = 0.0
c
      enddo
c
c   Compute weights
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            a       = sqrt( cabs( c(i,j) ) )
c
            w(j)    = w(j) + a
c
            w(i+ny) = w(i+ny) + a
c
         enddo
c
      enddo
c
c   Normalize weights
c
      wcen = w(1+ny/2)
c
      do j = 1 , ny
c
         w(j) = ( w(j) / wcen ) ** 2
c
         if ( w(j) .lt. eps ) w(j) = 1.0
c
      enddo
c
      wcen = w(ny+1+nx/2)
c
      do i = 1 , nx
c
         w(i+ny) = ( w(i+ny) / wcen ) ** 2
c
         if ( w(i+ny) .lt. eps ) w(i+ny) = 1.0
c
      enddo
c
c   Divide by weights
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            c(i,j) = c(i,j) / ( w(j) * w(i+ny) )
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine removewt ( ifile , ofile , c , nx , ny , w , eps ,
     .                      xyboth )
C**
C***********************************************************************
C**
      implicit none
c
      character ifile*80 , ofile*80
c
      character xyboth*1
c
      integer   nx , ny , i , j , reclda
c
      complex   c(nx)
c
      real      eps , w(nx+ny) , a , wcen
c
c   Open input and output files
c
      open ( 1 , file = ifile       , form = 'unformatted' ,
     .           status = 'old'     , access = 'direct'    ,
     .           recl = reclda( 0 , 2 * nx ) )
c
      open ( 2 , file = ofile       , form = 'unformatted' ,
     .           status = 'unknown' , access = 'direct'    ,
     .           recl = reclda( 0 , 2 * nx ) )
c
c   Zero weights
c
      do j = 1 , ny + nx
c
         w(j) = 0.0
c
      enddo
c
c   Compute weights
c
      do j = 1 , ny
c
         read ( 1 , rec = j ) c
c
         do i = 1 , nx
c
            a       = sqrt( cabs( c(i) ) )
c
            w(j)    = w(j) + a
c
            w(i+ny) = w(i+ny) + a
c
         enddo
c
      enddo
c
c   Normalize weights
c
      wcen = w(1+ny/2)
c
      do j = 1 , ny
c
         w(j) = ( w(j) / wcen ) ** 2
c
         if ( w(j) .lt. eps ) w(j) = 1.0
c
      enddo
c
      wcen = w(ny+1+nx/2)
c
      do i = 1 , nx
c
         w(i+ny) = ( w(i+ny) / wcen ) ** 2
c
         if ( w(i+ny) .lt. eps ) w(i+ny) = 1.0
c
      enddo
c
c   Divide by weights
c
      do j = 1 , ny
c
         read ( 1 , rec = j ) c
c
         if      ( xyboth .eq. 'x' .or. xyboth .eq. 'X' ) then
c
            do i = 1 , nx
c
               c(i) = c(i) / w(i+ny)
c
            enddo
c
         else if ( xyboth .eq. 'y' .or. xyboth .eq. 'Y' ) then
c
            do i = 1 , nx
c
               c(i) = c(i) / w(j)
c
            enddo
c
         else if ( xyboth .eq. 'b' .or. xyboth .eq. 'B' ) then
c
            do i = 1 , nx
c
               c(i) = c(i) / ( w(j) * w(i+ny) )
c
            enddo
c
         endif
c
         write ( 2 , rec = j ) c
c
      enddo
c
      return
      end
c**
c***********************************************************************
c**
      subroutine diff_images ( c1 , c2 , csum , cdif , nx , ny , phlag ,
     .                         maxlag , cancel_option , rcovs , rcovd ,
     .                         iq8bit )
c**
c***********************************************************************
c**
      implicit none
c
      integer   nx , ny , i , j , phlag , maxlag , cancel_option ,
     .          ilag , jlag , total , ilagm , jlagm , ilagp , jlagp ,
     .          im , jm , ip , jp , ir , jr , border , ilagbest ,
     .          jlagbest , nn(2)
c
      complex   c1(nx,ny) , c2(nx,ny) , csum(nx,ny) , cdif(nx,ny) ,
     .          cov12 , covsd , cov , cmin , czero , cminbest
c
      real      rcovs(1+2*maxlag,1+2*maxlag) , rcov1 , rcov2 ,
     .          rcovd(1+2*maxlag,1+2*maxlag) , raddeg , covbest , scale
c
      character iq8bit(nx,ny)*1
c
      character file*80
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      raddeg = 180.0 / atan2( 0.0 , - 1.0 )
c
      czero  = cmplx( 0.0 , 0.0 )
c
c   Buffer region around edges
c
      border = ( maxlag + 1 ) / 2
c
c   Size of non-zero region in sum and difference images
c
      total  = ( nx - 2 * border ) * ( ny - 2 * border ) 
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      if ( cancel_option .lt. 0 ) then
c
c   Image Domain: Estimate and remove the phase gradient between the
c                 two images
c
         write ( 7 , * ) ' Removing image-domain phase gradients'
c
         if ( cancel_option .lt. - 1 )
     .        call equal2imgs ( c1 , c2 , nx , ny , 128 )
c
         call remove_phase_gradient_2_2 ( c1 , c2 , nx , ny , phlag )
c
c   Time Domain:  Estimate and remove the phase gradient between the
c                 two images
c
         if ( cancel_option .lt. - 2 ) then
c
            nn(1)  = nx
c
            nn(2)  = ny
c
            call fourt2c ( c1 , nn , 2 , + 1 , 0 , csum , 2 * nx * ny )
c
            call fourt2c ( c2 , nn , 2 , + 1 , 0 , csum , 2 * nx * ny )
c
            write ( 7 , * ) ' Removing time-domain phase gradients'
c
            file = 'file1fft'
c
            call oframe ( c1 , csum , iq8bit , 0 , file , nx , ny , 0 ,
     .                    0 , 0.25 , 0 , 0 , 0 , 0 , 0.0 , 0 , 0 )
c
            file = 'file2fft'
c
            call oframe ( c2 , csum , iq8bit , 0 , file , nx , ny , 0 ,
     .                    0 , 0.25 , 0 , 0 , 0 , 0 , 0.0 , 0 , 0 )
c
            if ( cancel_option .lt. - 3 )
     .           call equal2imgs ( c1 , c2 , nx , ny , 128 )
c
            call remove_phase_gradient_2_2 ( c1 , c2 , nx , ny , phlag )
c
            call fourt2c ( c1 , nn , 2 , - 1 , 0 , csum , 2 * nx * ny )
c
            call fourt2c ( c2 , nn , 2 , - 1 , 0 , csum , 2 * nx * ny )
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Loop over all positive or negative lags up to maxlag, splitting plus
c   and minus lags to use the center of the image and make the smallest
c   border region
c
      do jlag = - maxlag , + maxlag
c
         jlagp = jlag / 2
         jlagm = jlagp - jlag
c
         jr    = jlag + 1 + maxlag        !  Array subscript
c
         do ilag = - maxlag , + maxlag
c
            ilagp = ilag / 2
            ilagm = ilagp - ilag
c
            ir    = ilag + 1 + maxlag     !  Array subscript
c
c-----------------------------------------------------------------------
c-------------------  Calculation for each lag  ------------------------
c
c   Calculate the covariance of the two images at this lag
c
            rcov1 = 0.0
c
            rcov2 = 0.0
c
            cov12 = czero
c
            do j = 1 + border , ny - border
c
               jp = j + jlagp
               jm = j + jlagm
c
               do i = 1 + border , nx - border
c
                  ip        = i + ilagp
                  im        = i + ilagm
c
                  rcov1 = rcov1  + c1(ip,jp) * conjg( c1(ip,jp) )
c
                  rcov2 = rcov2  + c2(im,jm) * conjg( c2(im,jm) )
c
                  cov12 = cov12  + c2(im,jm) * conjg( c1(ip,jp) )
c
               enddo
c
            enddo
c
            rcov1 = rcov1 / float( total )
c
            rcov2 = rcov2 / float( total )
c
            cov12 = cov12 / float( total )
c
c   Compute optimal cancellation coefficient
c
            if      ( cancel_option .eq. 0 ) then
c
               cmin = cmplx( 1.0 , 0.0 )
c
            else
c
               cmin = cov12 / rcov1
c
            endif
c
c   Compute the sum and difference images
c
            call csumdif ( c1 , c2 , csum , cdif , cmin , nx , ny ,
     .                     border , ilag , jlag )
c
c-----------------------------------------------------------------------
c
c   Compute the elements of the covariance matrix of the sum and
c   difference images
c
            rcovs(ir,jr) = cabs( cov( csum , csum , nx * ny ) )
c
            rcovd(ir,jr) = cabs( cov( cdif , cdif , nx * ny ) )
c
            covsd        = cov( csum , cdif , nx * ny )
c
c-----------------------------------------------------------------------
c
c   Output covariance information scaled by the variance of the sum
c   image
c
            write ( 7 , * )
c
            write ( 7 , * ) ' ilag, jlag:       ' , ilag , jlag
c
            write ( 7 , * )
     .         ' Variances for image_1, image_2:           ' ,
     .           rcov1 / rcovs(ir,jr) , rcov2 / rcovs(ir,jr)
c
            write ( 7 , * )
     .         ' Variances for sum and diff. images:       ' ,
     .           rcovs(ir,jr) / rcovs(ir,jr) ,
     .           rcovd(ir,jr) / rcovs(ir,jr)
c
            write ( 7 , * )
     .         ' Covariance amp. and phase for 2 images:   ' ,
     .           cabs( cov12 ) / rcovs(ir,jr) ,
     .           raddeg * atan2( aimag( cov12 ) , real( cov12 ) )
c
            write ( 7 , * )
     .         ' Covariance amp. and phase for s/d images: ' ,
     .           cabs( covsd ) / rcovs(ir,jr) ,
     .           raddeg * atan2( aimag( covsd ) , real( covsd ) )
c
c-----------------------------------------------------------------------
c
c   Determine the best lag to use and store so that the best sum and
c   difference images can be returned
c
            if ( ( ilag .eq. - maxlag .and. jlag .eq. - maxlag ) .or.
     .           ( rcovd(ir,jr) / rcovs(ir,jr) ) .lt. covbest ) then
c
               covbest  = rcovd(ir,jr) / rcovs(ir,jr)
c
               ilagbest = ilag
               jlagbest = jlag
c
               cminbest = cmin
c
            endif
c
c-------------------  Calculation for each lag  ------------------------
c-----------------------------------------------------------------------
c
         enddo         !  ilag
c
      enddo            !  jlag
c
      write ( 7 , * )
c
      write ( 7 , * ) ' ilagbest, jlagbest: ' , ilagbest , jlagbest
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Make sure the sum and difference images correspond to the best lag
c
      if ( ( ilagbest .ne. maxlag ) .or. ( jlagbest .ne. maxlag ) ) then
c
         call csumdif ( c1 , c2 , csum , cdif , cminbest , nx , ny ,
     .                  border , ilagbest , jlagbest )
c
      endif
c
c   Normalize both images by the standard deviation of the sum image
c
      jr        = jlagbest + 1 + maxlag
c
      ir        = ilagbest + 1 + maxlag
c
      scale     = 1.0 / sqrt( rcovs(ir,jr) )
c
      csum(:,:) = scale * csum(:,:)
c
      cdif(:,:) = scale * cdif(:,:)
c
      return
      end
c**
c***********************************************************************
c**
      subroutine equal2imgs ( c1 , c2 , nx , ny , blocksize )
c**
c***********************************************************************
c**
c   Equalize the power between two complex images using local weights
c
c      c1      :  First image
c
c      c2      :  Second image
c
c      nx, ny  :  Dimensions
c
c      block   :  Block size for local weighting
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      integer nx , ny , blocksize , i , j , ib , jb , iblock , jblock ,
     .        ii , jj
c
      complex c1(nx,ny) , c2(nx,ny) , csum1 , csum2
c
      real    scale
c
      do jb = 1 , ny / blocksize
c
         if ( jb .eq. ny / blocksize ) then
c
            jblock = ny - ( jb - 1 ) * blocksize
c
         else
c
            jblock = blocksize
c
         endif
c
         do ib = 1 , nx / blocksize
c
            if ( ib .eq. nx / blocksize ) then
c
               iblock = nx - ( ib - 1 ) * blocksize
c
            else
c
               iblock = blocksize
c
            endif
c
            csum1 = cmplx( 0.0 , 0.0 )
c
            csum2 = cmplx( 0.0 , 0.0 )
c
            do j = 1 , jblock
c
               jj = ( jb - 1 ) * blocksize + j
c
               do i = 1 , iblock
c
                  ii    = ( ib - 1 ) * blocksize + i
c
                  csum1 = csum1 + c1(ii,jj) * conjg( c1(ii,jj) )
c
                  csum2 = csum2 + c2(ii,jj) * conjg( c2(ii,jj) )
c
               enddo
c
            enddo
c
            scale = sqrt( abs( csum1 / csum2 ) )
c
            write ( 7 , * ) ' IBlock, JBlock, Scale: ' , ib , jb , scale
c
            do j = 1 , jblock
c
               jj = ( jb - 1 ) * blocksize + j
c
               do i = 1 , iblock
c
                  ii        = ( ib - 1 ) * blocksize + i
c
                  c2(ii,jj) = scale * c2(ii,jj)
c
               enddo
c
            enddo
c
         enddo
c
      enddo
c
      return
      end
c**
c***********************************************************************
c**
      subroutine csumdif ( c1 , c2 , csum , cdif , cmin , nx , ny ,
     .                     border , ilag , jlag )
c**
c***********************************************************************
c**
c   Compute the weighted sum and difference of two complex images
c
c      c1      :  First image
c
c      c2      :  Second image
c
c      csum    :  Sum image          ( c2 + cmin * c1 )
c
c      cdif    :  Difference image   ( c2 - cmin * c1 )
c
c      cmin    :  Complex weight for the first image
c
c      nx      :  x-dimension of images
c
c      ny      :  y-dimension of images
c
c      border  :  Width of zeroes at edge of sum and difference images
c
c      ilag    :  x-lag between images
c
c      jlag    :  y-lag between images
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      integer nx , ny , border , ilag , jlag , im , ip , jm , jp ,
     .        ilagm , ilagp , jlagm , jlagp , border_use , maxlag_use
c
      complex c1(nx,ny) , c2(nx,ny) , csum(nx,ny) , cdif(nx,ny) , cmin
c
c-----------------------------------------------------------------------
c
c   Ensure that the array limits are not violated by using a larger
c   border if necessary
c
      maxlag_use = max( iabs( ilag ) , iabs( jlag ) )
c
      border_use = max( border , ( maxlag_use + 1 ) / 2 )
c
c   Split the lag up into positive and negative terms to center the
c   result
c
      jlagp      = jlag / 2
      jlagm      = jlagp - jlag
c
      ilagp      = ilag / 2
      ilagm      = ilagp - ilag
c
c   Compute the sum and difference images
c
c   First, zero the entire images
c
      csum(:,:)  = cmplx( 0.0 , 0.0 )
c
      cdif(:,:)  = cmplx( 0.0 , 0.0 )
c
c   Next, scan the interior of the images computing the weighted sum
c   and differences
c
      do j = 1 + border_use , ny - border_use
c
         jp = j + jlagp
         jm = j + jlagm
c
         do i = 1 + border_use , nx - border_use
c
            ip        = i + ilagp
            im        = i + ilagm
c
            csum(i,j) = 0.5 * ( c2(im,jm) + cmin * c1(ip,jp) )
c
            cdif(i,j) = 0.5 * ( c2(im,jm) - cmin * c1(ip,jp) )
c
         enddo
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine addimg ( ci , cibig , nxci , nyci , nxcbig , nycbig ,
     .                    nxoff , nyoff , option , scale )
C**
C***********************************************************************
C**
c   Adds complex image (ci) to a larger array (cibig)
c
c   Aside from the two arrays and their dimensions, the important
c   arguments are:
c
c             nxoff, nyoff :  Offsets to the overlay position
c
c             option       :  0 for replace, 1 for addition
c
c             scale        :  Scale factor
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nxci , nyci , nxcbig , nycbig , option , i , j ,
     .        nxoff , nyoff
c
      complex ci(nxci,nyci) , cibig(nxcbig,nycbig)
c
      real    scale
c
      if ( option .eq. 0 ) then
c
c   Replace section of big image with small image
c
         do j = 1 , nyci
c
            do i = 1 , nxci
c
               if ( i + nxoff .le. nxcbig .and. i + nxoff .gt. 0 .and.
     .              j + nyoff .le. nycbig .and. j + nyoff .gt. 0 )
     .              cibig(i+nxoff,j+nyoff) = scale * ci(i,j)
c
            enddo
c
         enddo
c
      else if ( option .eq. 1 ) then
c
c   Add small image to big image
c
         do j = 1 , nyci
c
            do i = 1 , nxci
c
               if ( i + nxoff .le. nxcbig .and. i + nxoff .gt. 0 .and.
     .              j + nyoff .le. nycbig .and. j + nyoff .gt. 0 )
     .              cibig(i+nxoff,j+nyoff) = cibig(i+nxoff,j+nyoff) +
     .                                       scale *  ci(i,j)
c
            enddo
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine unwrap ( c , p , nx )
C**
C***********************************************************************
C**
      implicit none
c
      integer nx , i
c
      real    p(nx) , unwrap_p , atan2p
c
      complex c(nx) , d_cov
c
c   Compute un-wrapped phase, starting at the middle
c
      p(1+nx/2) = atan2p( aimag( c(1+nx/2) ) , real( c(1+nx/2) ) )
c
c   Do the right half forwards from the middle
c
      do i = 2 + nx / 2 , nx
c
         d_cov    = c(i) * conjg( c(i-1) )
c
         unwrap_p = atan2p( aimag( d_cov ) , real( d_cov ) )
c
         p(i)     = p(i-1) + unwrap_p
c
      enddo
c
c   Do the left half backwards from the middle
c
      do i = nx / 2 , 1 , - 1
c
         d_cov    = c(i+1) * conjg( c(i) )
c
         unwrap_p = atan2p( aimag( d_cov ) , real( d_cov ) )
c
         p(i)     = p(i+1) - unwrap_p
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine accel ( z , n , b , dt , a )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i
c
      complex z(n) , zsum1 , zsum2
c
      real    b , dt , pi , arg1 , arg2 , argdif , a
c
      pi = atan2( 0.0 , - 1.0 )
c
c  
c  Estimates signal acceleration and amplitude rate for time series
c  from velocities and amplitudes from the first and second half
c  of time series.
c
c     z(n)  :  Complex time series
c     n     :  Number of samples
c     b     :  Acceleration of signal ( Hz/s )
c     dt    :  Time separation between samples
c     a     :  Amplitude rate of signal
c
c  Compute signal covariance for first and second half of time series
c
      zsum1 = cmplx( 0.0 , 0.0 )
      zsum2 = cmplx( 0.0 , 0.0 )
c
      do i = 1 , n / 2 - 1
c
         zsum1 = zsum1 + z(i) * conjg( z(i+1) )
c
      enddo
c
      do i = 1 + n / 2 , n - 1
c
         zsum2 = zsum2 + z(i) * conjg( z(i+1) )
c
      enddo
c
c  Compute phase of signal covariance (proportional to velocity)
c  for first and second half of time series 
c
      if ( zsum1 .eq. cmplx( 0.0 , 0.0 ) ) then
c
         arg1 = 0.0
c
      else
c
         arg1 = atan2( aimag( zsum1 ) , real( zsum1 ) )
c
      endif
c
      if ( zsum2 .eq. cmplx( 0.0 , 0.0 ) ) then
c
         arg2 = 0.0
c
      else
c
         arg2 = atan2( aimag( zsum2 ) , real( zsum2 ) )
c
      endif
c
c  Compute phase difference (proportional to velocity difference) from
c  first to second half of time series; if necessary, dealias phase diff
c
      argdif = arg2 - arg1
      if ( argdif .gt.   pi ) argdif = argdif - 2.0 * pi
      if ( argdif .lt. - pi ) argdif = argdif + 2.0 * pi
c
c  Estimate fractional amplitude rate (1/amp) * d[amp]/dt from the
c  first to second half of the time series. Use the equivalent
c  expression d[ln(amp)]/dt and the ratio amp2/amp1 to compute the rate.
c  The 0.1 terms were added to save the logarithm if the amp in the
c  first or second half (but not both) of the time series was zero.
c
      if ( cabs( zsum1 ) .eq. 0.0 .or. cabs( zsum2 ) .eq. 0.0 ) then
c
         a = 0.0
c
      else
c
         a = 0.5 * alog( ( 0.1 * cabs( zsum1 ) + cabs( zsum2 ) ) /
     .                   ( 0.1 * cabs( zsum2 ) + cabs( zsum1 ) ) ) /
     .                       ( float( n ) * dt / 2.0 )
c
      endif
c
c  Estimate acceleration in Hz/s as d[phi]/dt / (2*pi*dt)
c
      b = argdif / ( 2.0 * pi * dt * dt * float( n / 2 ) )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine cshift ( cin , n , nshift , cout )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , nshift , i , ishift , jshift
c
      complex cin(n) , cout(n)
c
c   cout is cin shifted by nshift points
c
      if ( nshift .ge. 0 ) then
c
         ishift = nshift
c
      else
c
         ishift = nshift + n * ( 1 + abs( nshift ) / n )
c
      endif
c
      do i = 1 , n
c
         jshift       = 1 + mod( i + ishift - 1 , n )
         cout(jshift) = cin(i)
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine dchirp ( z , n , bavg , vel , sg0 , d0 , dd , dt ,
     .                    nterms , aavg , lambda , dc , dw )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , nterms , i , nlag , m , m2 , m22 , m222 , mlim
c
      complex z(n) , cj
c
      real    b1 , b2 , b3 , bavg , dt , a1 , a2 , a3 , aavg , lambda ,
     .        vel , d0 , dd , sg0 , pi , arg , ph , dc , dw , dcnew ,
     .        dwnew , sg0new , d0new , argnew , snr , snrold , snrnew ,
     .        bavg0 , dtm2
c
      logical done
c
c   Computes estimates of target acceleration, signal + noise power,
c   centroid, and Doppler width from a time series for a single range
c   cell.  Recursively improves estimates using filtered versions of
c   time series centered on the target.  All processing is done in the
c   time domain.
c
c    z(n)   :  Complex time series
c    n      :  Number of samples
c    bavg   :  Estimated target acceleration ( Hz/s )
c    vel    :  Estimated target velocity ( m/s )
c    sg0    :  Estimated target + noise power
c    d0     :  Estimated square of Doppler width ( non-dimensional )
c    dd     :  Initial square of Doppler width ( non-dimensional )
c    dt     :  Time separation between samples
c
c    nterms :  abs(nterms) =0,1,2 means despin time series with nothing,
c              current velocity estimate, current velocity and
c              acceleration estimates before making new
c              parameter estimates
c
c              nterms < 0 means perform recursive estimate improvement
c
c    aavg   :  Estimated target amplitude time rate of change
c    lambda :  Wavelength ( m )
c    dc     :  Estimated target Doppler centroid ( Hz )
c    dw     :  Estimated target Doppler width ( Hz )
c  
      pi   = atan2( 0.0 , - 1.0 )
      cj   = cmplx( 0.0 , 1.0 )
c
      nlag = 1
c
c   Estimate spectral parameters dc, dw, d0, sg0, and arg from signal
c   covariance using entire time series
c
      call onelag ( z , n , dt , nlag , dc , dw , d0 , sg0 , arg )
c
      vel  = dc * lambda / 2.0
      dd   = d0
c
c   If desired, remove amplitude variation, velocity, and possibly 
c   acceleration from time series to improve estimates of spectral
c   parameters.  Removing velocity by itself does not improve spectral
c   estimates, but having the target at dc will be necessary for the
c   next step
c
      if ( nterms .eq. 0 ) then ! don't estimate/remove amp, vel, accel
c
         aavg = 0.0
         bavg = 0.0
c
      else
c
c        Estimate accel and amp variation using three subsequences
c
         call accel ( z(1)     , n / 2 , b1 , dt , a1 )
         call accel ( z(1+n/4) , n / 2 , b2 , dt , a2 )
         call accel ( z(1+n/2) , n / 2 , b3 , dt , a3 )
c
         aavg = ( a1 + a2 + a3 ) / 3.0 ! ave amp variation
         bavg = ( b1 + b2 + b3 ) / 3.0 ! ave accel over 3 subsequences
c
c        Remove vel, amp variation, and possibly acceleration
c
         do i = 1 , n
c
            if ( iabs( nterms ) .eq. 1 ) ph = arg * float( i ) ! vel
c
            if ( iabs( nterms ) .eq. 2 ) ph = arg * float( i )
     .             + bavg * pi * ( dt * float( i - n / 2 ) ) ** 2
c
            z(i) = z(i) * cexp(
     .             cj * ph - 0.5 * aavg * float( i - n / 2 ) * dt )
c
         enddo
c
c        Compute improved spectral estimates
c
         call onelag ( z , n , dt , nlag , dcnew , dw , dd , sg0 , arg )
c
         dc  = dc + dcnew
         vel = dc * lambda / 2.0
c
      endif
c
c   Compute how many times the filtered time series centered on the
c   target can be cut in half and still end up with 16 points in the
c   smoothed time series.  16 points represent the minimum desired
c   length time series that can be cleanly broken into four disjoint
c   sections allowing three independent acceleration estimates
c
      m    = n
      mlim = 0
c
      do while ( m .ge. 32 )
c
        m    = m / 2
        mlim = mlim + 1
c
      enddo
c
c   Recursively improve spectral estimates using filtered versions
c   of time series centered on the target.  All processing is done in
c   the time domain.  Nyquist is effectively halved on each iteration.
c
      if ( nterms .lt. 0 .and. mlim .gt. 0 ) then
c
c        Make more passes using filtered versions of the time series
c
         m      = 0
         snrold = snr( sg0 , dw , dt , 1.0 )
         bavg0  = bavg
         done   = .false.
c
         do while ( m .lt. mlim .and. ( .not. done ) )
c
            m    = m + 1
            m2   = 2 ** m
            m22  = 2 * m2
            m222 = 2 * m22
c
            dtm2 = float( m2 ) * dt
c
c           Decimate 1/2 and smooth time series
c
            do i = 1 , n / m2
c
               z(i) = 0.5 * ( z(2*i-1) + z(2*i) )
c
            enddo
c
c           Compute improved spectral estimates
c
            call onelag ( z , n / m2 , dtm2 , nlag , dcnew ,
     .                    dwnew , d0new , sg0new , argnew )
c
            snrnew = snr( float( m2 ) * sg0new , dwnew , dt , 1.0 )
c          
c           If snr increases and Doppler width decreases, update
c           spectral estimates, despin the target to dc, and iterate
c
            if ( snrnew .gt. snrold .and. dwnew .lt. dw ) then
c
               snrold = snrnew
               sg0    = float( m2 ) * sg0new
               dw     = dwnew
               d0     = d0new
               arg    = argnew
c
               call accel ( z(1)        , n / m22 , b1 , dtm2 , a1 )
               call accel ( z(1+n/m222) , n / m22 , b2 , dtm2 , a2 )
               call accel ( z(1+n/m22)  , n / m22 , b3 , dtm2 , a3 )
c
               if ( nterms .eq. - 2 ) then
c
                  bavg = bavg0 + ( b1 + b2 + b3 ) / 3.0
c
               else
c
                  bavg = ( b1 + b2 + b3 ) / 3.0
c
               endif
c
               dc  = dc + dcnew
               vel = dc * lambda / 2.0
c
               do i = 1 , n / m2
c
                  z(i) = z(i) * cexp( cj * arg * float( i ) )
c
               enddo
c
            else
c
               done = .true.
c
            endif
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine hwindo ( w , n , nt , nf , work )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , nt , i , nf
c
      real    w(nf) , pi , work(nf)
c
c   Compute n-pt (odd) Hanning window with nt-pt tapers.
c
      pi = 4.0 * atan( 1. )
c
      do i = 1 , nt
c
         work(i) = 1.0 -
     .   cos( ( pi * float( i - 1 ) ) / float ( 2 * ( nt - 1 ) ) ) ** 2
c
      enddo
c
      do i = nt + 1 , n / 2 + 1
c
         work(i) = 1.
c
      enddo
c
      do i = n / 2 + 2 , n
c
         work(i) = work(n-i+1)
c
      enddo
c
c   Zero pad to length nf
c
      do i = n + 1 , nf
c
         work(i) = 0.0
c
      enddo
c
c   Center about bin 1
c
      call rshift ( work , nf , - n / 2 , w )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine onelag ( c , n , dt , nlag , dc , dw , d , pwr , phi )
C**
C***********************************************************************
C**
c
c   Estimate the spectral parameters from the signal covariance
c   calculated at a single lag
c
c      c(n)  :  Complex time series
c      n     :  Number of samples
c      dt    :  Time separation between samples
c      nlag  :  Lag at which covariance is to be estimated
c      dc    :  Estimated Doppler centroid ( Hz )
c      dw    :  Estimated Doppler width ( Hz )
c      d     :  Square of Doppler width ( non-dimensional )
c      pwr   :  Total signal + noise power
c      phi   :  Phase of covariance ( radians )
c
      implicit none
c
      integer n , nlag , i
c
      complex c(n) , cov , dcov
c
      real    dt , dc , dw , d , cova , pwr , phi , pi
c
      pi   = atan2( 0.0 , - 1.0 )
      cov  = cmplx( 0.0 , 0.0 )
      cova = 0.0
      pwr  = 0.0
c
      do i = 1 , n - nlag
c
         dcov = c(i) * conjg( c(i+nlag) )          ! lagged product
         cov  = cov + dcov                         ! signal covariance
         cova = cova + cabs( dcov )                ! signal power
         pwr  = pwr + cabs( c(i) * conjg( c(i) ) ) ! signal + noise
c
      enddo
c
c   Compute total signal + noise power
c
      pwr  = pwr + 0.5 * ( cabs( c(n) * conjg( c(n) ) ) -
     .                     cabs( c(1) * conjg( c(1) ) ) )
      pwr  = amax1( 0.0 , pwr ) ! possibly trap above minus sign
c
c   Compute phase of signal covariance
c
      if ( cov .eq. cmplx( 0.0 , 0.0 ) ) then
c
         phi = 0.0
c
      else
c
         phi = atan2( aimag( cov ) , real( cov ) )
c
      endif
c
c   Estimate Doppler centroid from covariance phase
c
      dc   = phi / ( 2.0 * pi * float( nlag ) * dt )
c
c   Estimate square of Doppler width from the signal covariance
c   magnitude and signal power
c
      if ( cova .eq. 0.0 ) then
c
         d = 1.0
c
      else
c
         d = amax1( 0.0 , 1.0 - cabs( cov ) / cova )
c
      endif
c
c   Compute Doppler width in Hz
c
      dw   = sqrt( 2.0 * d ) / ( 2.0 * pi * float( nlag ) * dt )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine rshift ( rin , n , nshift , rout )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , nshift , i , ishift , jshift
c
      real    rin(n) , rout(n)
c
c   real rout is real rin shifted by nshift points
c
      if ( nshift .ge. 0 ) then
c
         ishift = nshift
c
      else
c
         ishift = nshift + n * ( 1 + abs( nshift ) / n )
c
      endif
c
      do  i = 1 , n
c
         jshift       = 1 + mod( i + ishift - 1 , n )
         rout(jshift) = rin(i)
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine sort ( ra , n , jindex )
C**
C***********************************************************************
C**
      implicit none
c
      integer n , i , ir , j , k , iindex , jindex(n)
c
      real    ra(n) , rra
c
c   Replaces ra in sorted ascending order with associated index
c   using Heapsort.
c
      if ( n .lt. 2 ) return
c
      k  = n / 2 + 1
      ir = n
c
      do i = 1 , n
c
         jindex(i) = i
c
      enddo
c
   10 continue
c
      if ( k .gt. 1 ) then
c
         k          = k - 1
         rra        = ra(k)
         iindex     = jindex(k)
c
      else
c
         rra        = ra(ir)
         ra(ir)     = ra(1)
         iindex     = jindex(ir)
         jindex(ir) = jindex(1)
         ir         = ir - 1
c
         if ( ir .eq. 1 ) then
c
            ra(1)     = rra
            jindex(1) = iindex
c
            return
c
         endif
c
      endif
c
      i = k
      j = k + k
c
   20 if ( j .le. ir ) then
c
         if ( j .lt. ir ) then
c
            if ( ra(j) .lt. ra(j+1) ) j = j + 1
c
         end if
c
         if ( rra. lt. ra(j) ) then
c
            ra(i)     = ra(j)
            jindex(i) = jindex(j)
            i         = j
            j         = j + j
c
         else
c
            j        = ir + 1
c
         end if
c
         go to 20
c
      endif
c
      ra(i)     = rra
      jindex(i) = iindex
c
      go to 10
c
      end
C**
C***********************************************************************
C**
      subroutine taylor ( nx , sll_db , weights )
C**
C***********************************************************************
C**
c
c   This subroutine calculates the Taylor illumination function as
c   defined on pg. 249 of Barton, HANDBOOK OF RADAR MEASUREMENT.
c
c   Input parameters:
c
c        nx     = # of data points
c        sll_db = sidelobe level in dB (This should negative)
c
c   Output array:
c
c        weights = Taylor illumination function . This is normalized 
c                  to be one at the center.         
c
      implicit none
c
      integer     nmax , nx , nbar , m , n , ix
c
      parameter ( nmax = 30 )
c
      real        weights(nx) , f(nmax) , pi , sll_db , x , a , dwt ,
     .            sigma2 , fact , df , denom , xnum
c
      data pi / 3.14159265 /
c
      x    = 10.0 ** ( abs( sll_db ) / 20.0 )
      x    = x + sqrt( x ** 2 - 1.0 )
      a    = alog( x ) / pi
      nbar = 2 * a ** 2 + 1.5
c
      if ( nbar .gt. ( nmax + 1 ) ) then
         stop ' nbar is too large '
      endif
c
      sigma2 = float( nbar ) ** 2 /
     .         ( a ** 2 + ( float( nbar ) - 0.5 ) ** 2 )
      fact = 1.0
c
      do m = 1 , nbar - 1
c
         dwt = 1.0
c
         do n = 1 , nbar - 1
c
            df = ( 1.0 - float( m ) ** 2 / sigma2 /
     .           ( a ** 2 + ( float( n ) - 0.5 ) ** 2 ) )
c
            if ( n .ne. m ) df = df / ( 1.0 - float( m ) ** 2 /
     .                                        float( n ) ** 2 )
c
            dwt = dwt * df
c
         enddo
c
         f(m) = fact * dwt / 2.0
         fact = - fact
c
      enddo
c
      do ix = 1 , nx
c
         a     = float( ix - 1 ) / float( nx - 1 ) - 0.5
         xnum  = 0.0
         denom = 0.0
c
         do m = 1 , nbar - 1
c
            xnum  = xnum + f(m) * cos( 2.0 * pi * float( m ) * a )
            denom = denom + f(m)
c
         enddo
c
         weights(ix) = ( 1.0 + 2.0 * xnum ) / ( 1.0 + 2.0 * denom )
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine rcomp ( br , clight , dr0 , rc , crc , pt , ntr , ncr ,
     .                   wtr , dtr , mode , work , nwork )
C**
C***********************************************************************
C**
c
c   Purpose:  Range compresses a single pulse given the desired values
c             of the mo-comp center, dr0, and the mo-comp velocity, vr0.
c
c***********************************************************************
c
      implicit none
c
      integer ntr , ncr , mode , ir , nwork
c
      real    br , clight , dr0 , pt , phase , wtr(ntr) , work(nwork) ,
     .        const , dtr
c
      complex cpfast
c
      real    rc(2,ncr)
      complex crc(ncr)
c
      const = ( 2.0 * 1.0E+12 * br * dr0 ) / clight
c
      do ir = 1 , ntr
c
         phase   = pt + const * dtr * float( ir - 1 - ntr / 2 )
c
         crc(ir) = crc(ir) * cpfast( phase )
c
      enddo
c
c   Compress in range and store in range-compressed array.  Allow for
c   up-chirp (br>0) or down-chirp (br<0)
c
      if ( mode .le. 100 ) then
c
         if ( br .gt. 0.0 ) then
c
            call wtfft ( rc , wtr , rc , work , + 1 , - ntr / 2 ,
     .                   ntr , ncr )
c
         else
c
            call wtfft ( rc , wtr , rc , work , - 1 , - ntr / 2 ,
     .                   ntr , ncr )
c
         endif
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      subroutine chirp_compress ( rawdat , compressed , chirp , ntr ,
     .                            br , sampr , bw_frac , work , first ,
     .                            wtr , taywtr , signal , noise )
c**
c***********************************************************************
c**
      implicit none
c
      integer ntr , k , frac , nn(1)
c
      real    rawdat(2,ntr) , chirp(2,2*ntr) , compressed(2,2*ntr) ,
     .        work(2,2*ntr) , pi , br , sampr , bw_frac , dt , f ,
     .        fftscale , taywtr , wtr(2*ntr) , noise , signal
c
      complex chirptmp , ctmp
c
      logical first
c
      save
c
c-----------------------------------------------------------------------
c
      if ( first ) then
c
         pi         = atan2( 0.0 , - 1.0 )
c
         dt         = 1.0 / sampr
c
         chirp(:,:) = 0.0
c
         frac       = int( ( 1.0 - bw_frac ) * float( ntr ) )
c
         call taylor ( 2 * ntr - 2 * frac , taywtr , wtr(1+frac) )
c
         do k = 1 + frac , 2 * ntr - frac
c
c   Frequency in Hz
c
            f          = float( k - 1 - ntr ) /
     .                 ( float( 2 * ntr ) * dt )
c
c   The positive sign to the quadratic term corresponds to the match
c   filter for a positive quadratic phase (concave upwards) for the
c   signal in the time domain.  Given a positive quadratic signal, the
c   match filter to it is a negative quadratic term; the Fourier
c   transform of this is a positive quadratic term in the frequency
c   domain.
c
            chirptmp   = cexp( cmplx( 0.0 , ( pi / br ) * f ** 2 ) )
c
            if ( mod( k , 2 ) .eq. 0 ) chirptmp = - chirptmp
c
            chirp(1,k) =  real( chirptmp ) * wtr(k)
            chirp(2,k) = aimag( chirptmp ) * wtr(k)
c
         enddo
c
c        write ( 6 , * ) ' Chirp function formed'
c
         first = .false.
c
      endif
c
c-----------------------------------------------------------------------
c
c             Load data array
c
c   Zero-fill to double the original size and FFT the pulse, flipping
c   signs to put DC at center
c
         compressed(:,:) = 0.0
c
         fftscale = 1.0 / sqrt( float( ntr ) )  ! Preserve Std. Dev.
c
         do k = 2 , ntr , 2
c
            compressed(:,k-1+ntr/2) = rawdat(:,k-1) * fftscale
c
            compressed(:,k+ntr/2)   = - rawdat(:,k) * fftscale
c
         enddo
c
c-----------------------------------------------------------------------
c
c             FFT the data
c
         nn(1) = 2 * ntr
c
         call fourt ( compressed , nn , 1 , - 1 , 1 , work , 4 * ntr )
c
         noise  = 0.0
c
         signal = 0.0
c
         if ( frac .gt. 2 ) then
c
            do k = 1 , frac / 2
c
               noise  = noise  + compressed(1,k) ** 2
     .                         + compressed(2,k) ** 2
c
               noise  = noise  + compressed(1,2*ntr+1-k) ** 2
     .                         + compressed(2,2*ntr+1-k) ** 2
c
               signal = signal + compressed(1,ntr+k) ** 2
     .                         + compressed(2,ntr+k) ** 2
c
               signal = signal + compressed(1,ntr+1-k) ** 2
     .                         + compressed(2,ntr+1-k) ** 2
c
            enddo
c
         else
c
            noise  = 1.0
c
            signal = 1.0
c
         endif
c
c-----------------------------------------------------------------------
c
c             Multiply by FFT of the chirp
c
         do k = 1 , 2 * ntr
c
            ctmp              = compressed(1,k)
c
            compressed(1,k)   = compressed(1,k) * chirp(1,k) -
     .                          compressed(2,k) * chirp(2,k)
c
            compressed(2,k)   = ctmp            * chirp(2,k) +
     .                          compressed(2,k) * chirp(1,k)
c
         enddo
c
c-----------------------------------------------------------------------
c
c             FFT the pulse to get range-compressed output
c
         fftscale = 1.0 / sqrt( float( 2 * ntr ) )  ! Preserve Std. Dev.
c
         do k = 2 , 2 * ntr , 2
c
            compressed(:,k-1) = + compressed(:,k-1) * fftscale
c
            compressed(:,k)   = - compressed(:,k)   * fftscale
c
         enddo
c
         nn(1) = 2 * ntr
c
         call fourt ( compressed , nn , 1 , 1 , 1 , work , 4 * ntr )
c
c   Flip signs to correct for the fact that DC is in the middle of the
c   frequency domain
c
         do k = 2 , 2 * ntr , 2
c
            compressed(:,k)   = - compressed(:,k)
c
         enddo
c
c   There are now 2 * ntr range cells, but only the middle ones are
c   fully compressed.
c
c-----------------------------------------------------------------------
c
      return
      end
c**
c***********************************************************************
c**
      subroutine iqstats ( rdata , ntr , imean , qmean , iqcor ,
     .                     iscale , qscale , a , iqrms , imin , imax ,
     .                     qmin , qmax )
c**
c***********************************************************************
c**
      implicit none
c
      integer ntr , k
c
      real    rdata(2,ntr) , ibar , qbar , imean , qmean , iibar ,
     .        qqbar , iscale , qscale , iqbar , iqcor , a , iqrms ,
     .        imin , imax , qmin , qmax
c
c   Compute the mean I/Q values for this pulse so that mean I/Q values
c   and correlation between the channels can be adaptively removed
c
      ibar  = 0.0
      qbar  = 0.0
c
      imin  = rdata(1,1)
      qmin  = rdata(2,1)
c
      imax  = rdata(1,1)
      qmax  = rdata(2,1)
c
      iibar = 0.0
      qqbar = 0.0
c
      iqbar = 0.0
c
      do k = 1 , ntr
c
         imin  = min( imin , rdata(1,k) )
c
         imax  = max( imax , rdata(1,k) )
c
         qmin  = min( qmin , rdata(2,k) )
c
         qmax  = max( qmax , rdata(2,k) )
c
         ibar  = ibar  + rdata(1,k)
c
         iibar = iibar + rdata(1,k) ** 2
c
         qbar  = qbar  + rdata(2,k)
c
         qqbar = qqbar + rdata(2,k) ** 2
c
         iqbar = iqbar + rdata(1,k) * rdata(2,k)
c
      enddo
c
      imean  = ibar / float( ntr )
      qmean  = qbar / float( ntr )
c
      iibar  = amax1( 0.0 , iibar / float( ntr ) - imean ** 2 )
      qqbar  = amax1( 0.0 , qqbar / float( ntr ) - qmean ** 2 )
c
      iqbar  = iqbar / float( ntr ) - imean * qmean
c
      if ( iibar .gt. 0.0 .and. qqbar .gt. 0.0 ) then
c
         iscale = sqrt( 2.0 * iibar / ( iibar + qqbar ) )
         qscale = sqrt( 2.0 * qqbar / ( iibar + qqbar ) )
c
         iqcor  = iqbar / sqrt( iibar * qqbar )
c
         a      = - iqbar / qqbar
c
         iqrms  = sqrt( 0.5 * ( iibar + qqbar ) )
c
      else
c
         iscale = 0.0
         qscale = 0.0
         iqcor  = 0.0
         a      = 0.0
         iqrms  = 0.0
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      real function ctrast ( ci , nx , ny )
C**
C***********************************************************************
C**
c   Compute contrast function for a complex image, ci(nx,ny)
c
      implicit none
c
      integer           nx , ny , i , j
c
      double precision  delint , ibar , isqbar
c
      complex ci(nx,ny)
c
      ibar   = 0.0
      isqbar = 0.0
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            delint = cabs( ci(i,j) ) ** 2
            ibar   = ibar + delint
            isqbar = isqbar + delint ** 2
c
         enddo
c
      enddo
c
      ibar   = ibar / float( nx * ny )
      isqbar = isqbar / float( nx * ny )
c
      if ( ibar .ne. 0.0 ) then
c
         ctrast = isqbar / ibar ** 2
c
      else
c
         ctrast = 0.0
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      real function snr ( sg0 , dw , dt , n0 )
C**
C***********************************************************************
C**
      implicit none
c
      real sg0 , dw , dt , n0
c
c  Computes signal-to-noise ratio de-weighted by Doppler width
c
c    snr    = weight * power( signal + noise ) / power( noise )
c    weight = ( Nyquist width / Doppler width )
c
      if ( dw .lt. 0.001 / dt ) then  ! trap potential singularity
c
         snr = ( sg0 / n0 ) * 500.0
c
      else
c
         snr = ( sg0 / n0 ) * ( 1.0 / ( dt * 2.0 * dw ) )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      complex function cpfast0 ( p )
C**
C***********************************************************************
C**
c
c   Compute CEXP(J*2*PI*P) fast by pre-computing a table of values and
c   then using nearest neighbor interpolation
c
      implicit none
c
      real         p
c
      integer      ndelta , kdelta
c
      parameter  ( ndelta = 16385 )
c
      complex      cdelta(ndelta)
c
      real         ddelc , deltac , pi
c
      logical      first
c
      save         first , cdelta , ddelc
c
      data         first / .true. /
c
c***********************************************************************
c
c            Pre-compute complex exponentials
c
      if ( first ) then
c
         first = .false.
         pi    = atan2( 0.0 , - 1.0 )         !  Pi
c
         ddelc = 1.0 / float( ndelta - 1 )
c
         do kdelta = 1 , ndelta
c
            deltac         = ddelc * float( kdelta - 1 )
            cdelta(kdelta) = cexp( cmplx( 0.0 , 2.0 * pi * deltac ) )
c
         enddo
c
      endif
c
      deltac  = amod( p , 1.0 )
      if ( deltac .lt. 0.0 ) deltac = deltac + 1.0
c
      kdelta  = 1 + nint( deltac / ddelc )
c
      cpfast0 = cdelta(kdelta)
c
      return
      end
C**
C***********************************************************************
C**
      complex function cpfast ( p )
C**
C***********************************************************************
C**
c
c   Compute CEXP(J*2*PI*P) fast by pre-computing a table of values and
c   then using nearest neighbor interpolation
c
      implicit none
c
      logical first / .true. /
      real    pi
      save    first , pi
c      
      real    p
c
c***********************************************************************
c
      if ( first ) then
         pi    = atan2( 0.0 , - 1.0 )         !  Pi
         first = .false.
      endif
c
      cpfast = cexp( cmplx( 0.0 , 2.0 * pi * p ) )
c
      return
      end
C**
C***********************************************************************
C**
      integer function twopwr ( k )
C**
C***********************************************************************
C**
c
c   Computes the power of two greater or equal to the argument
c
      implicit none
c
      integer k
c
      twopwr = 1
c
      do while ( twopwr .lt. k )
c
         twopwr = 2 * twopwr
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      real function mflops ( work , nft )
C**
C***********************************************************************
C**
c   This routine estimates the Mega-Flop rate of a computer based on
c   the speed at which it does an FFT
c
      implicit none
c
      integer nft , k , i , nn(1) , nflops , nits
c
      real    work(2,2*nft)
c
      real    dtime , timer(2) , time0 , fttime , scale , alog2
c
      time0  = dtime( timer )
      fttime = 0.0
c
      nn(1)  = nft
c
      nits   = 10000
c
      do k = 1 , nits
c
         time0  = dtime( timer )
c
         if ( mod( k , 2 ) .eq. 0 ) then
c
            call fourt ( work , nn , 1 , - 1 , 1 , work(1,1+nft) ,
     .                   2 * nft )
c
         else
c
            call fourt ( work , nn , 1 , 1 , 1 , work(1,1+nft) ,
     .                   2 * nft )
c
         endif
c
         fttime = fttime + dtime( timer )
c
         scale = 1.0 / sqrt( float( nft ) )
c
         do i = 1 , nft
c
            work(1,i) = work(1,i) * scale
            work(2,i) = work(2,i) * scale
c
         enddo
c
      enddo
c
c   Apply the formula 5*N*log2(N) to the floating point operations
c   for an N-point complex FFT; multiply by no. of FFTs.
c
      nflops = 5 * nft * ifix( alog2( float( nft ) ) ) * nits
c
      mflops = 1.0E-6 * float( nflops ) / amax1( fttime , 1.0E-6 )
c
      return
      end
C**
C***********************************************************************
C**
      real function alog2 ( x )
C**
C***********************************************************************
C**
      implicit none
c
      real x
c
      alog2 = alog10( x ) / alog10( 2.0 )
c
      return
      end
C**
C***********************************************************************
C**
      real function rtrast ( ri , nx , ny )
C**
C***********************************************************************
C**
c   Compute contrast function for a real image, ri(nx,ny)
c
      implicit none
c
      integer          nx , ny , i , j
c
      double precision ibar , isqbar
c
      real ri(nx,ny)
c
      ibar   = 0.0
      isqbar = 0.0
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            ibar   = ibar + abs(ri(i,j))
            isqbar = isqbar + ri(i,j) ** 2
c
         enddo
c
      enddo
c
      ibar   = ibar / float( nx * ny )
      isqbar = isqbar / float( nx * ny )
c
      if ( ibar .ne. 0.0 ) then
c
         rtrast = isqbar / ibar ** 2
c
      else
c
         rtrast = 0.0
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine readci ( ci , nfa , nfrb , ifile , ftype , xoffset ,
     .                    joffset , ifirst , nfdim , nrdim , ctmp ,
     .                    i2tmp , iq8bit , platdir )
C**
C***********************************************************************
C**
c   Purpose:
c
c             To read a portion of a complex image from disk
c
      implicit none
c
      integer   nfa , nfrb , xoffset , joffset , ifirst , nfdim ,
     .          nrdim , i , j , k , reclda , ibmi2 
c
      character ifile*80 , ftype*80
c
      complex   ctmp(nfdim) , ci(nfa,nfrb)
c
      integer*2 i2tmp(2,nfdim)
c
      character iq8bit(2,2,nfdim)*1
c
      real      platdir
c
      logical   reopen
c
      character lastfile*80
c
      save      lastfile
c
c-----------------------------------------------------------------------
c
      reopen   = ( ifirst .eq. 0 .or. ifile .ne. lastfile )
c
      lastfile = ifile
c
c   Read in full resolution complex image to be re-focused using
c   either Floating Point or SANDIA IQ4
c
      if ( ftype(1:1) .eq. 'f' .or. ftype(1:1) .eq. 'F' ) then
c
         if ( reopen ) open ( 9 , file = ifile ,
     .                     form   = 'UNFORMATTED' , status = 'UNKNOWN' ,
     .                     access = 'DIRECT'      ,
     .                     recl   = reclda( 0 , 2 * nfdim ) )
c
         do j = 1 , nfrb
c
            if ( j + joffset .gt. 0 .and. j + joffset .le. nrdim ) then
c
               read ( 9 , rec = j + joffset )
     .              ( ctmp(i) , i = 1 , nfdim )
c
               do i = 1 , nfa
c
                  if ( i + xoffset .gt. 0 .and. i + xoffset .le. nfdim )
     .                                    then
c
                     ci(i,j) = ctmp(i+xoffset)
c
                  else
c
                     ci(i,j) = cmplx( 0.0 , 0.0 )
c
                  endif
c
               enddo
c
            else
c
               do i = 1 , nfa
c
                  ci(i,j) = cmplx( 0.0 , 0.0 )
c
               enddo
c
            endif
c
         enddo
c
      else if ( ftype(1:1) .eq. 's' .or. ftype(1:1) .eq. 'S' ) then
c
         if ( reopen ) open ( 9 , file = ifile ,
     .        form    = 'UNFORMATTED' , status = 'UNKNOWN' ,
     .        convert = 'big_endian'  , access = 'DIRECT' ,
     .        recl    = reclda( 0 , 2 * nfdim ) )
c
         do j = 1 , nfrb
c
            if ( j + joffset .gt. 0 .and. j + joffset .le. nrdim ) then
c
               read ( 9 , rec = j + joffset )
     .              ( ctmp(i) , i = 1 , nfdim )
c
               do i = 1 , nfa
c
                  if ( i + xoffset .gt. 0 .and. i + xoffset .le. nfdim )
     .                                    then
c
                     ci(i,j) = ctmp(i+xoffset)
c
                  else
c
                     ci(i,j) = cmplx( 0.0 , 0.0 )
c
                  endif
c
               enddo
c
            else
c
               do i = 1 , nfa
c
                  ci(i,j) = cmplx( 0.0 , 0.0 )
c
               enddo
c
            endif
c
         enddo
c
      else
c
         if ( reopen ) open ( 9 , file = ifile ,
     .        form = 'UNFORMATTED' , status = 'UNKNOWN' ,
     .        access = 'DIRECT' , recl = reclda( 4 * nfdim , 0 ) )
c
         do j = 1 , nfrb
c
            if ( j + joffset .gt. 0 .and. j + joffset .le. nrdim ) then
c
               read ( 9 , rec = j + joffset )
     .              ( ( i2tmp(k,i) , k = 1 , 2 ) , i = 1 , nfdim )
c
               if ( ibmi2() .eq. 1  ) call bswap ( 2 * nfdim , iq8bit )
c
               do i = 1 , nfa
c
                  if ( i + xoffset .gt. 0 .and. i + xoffset .le. nfdim )
     .                                    then
c
                     ci(i,j) = cmplx( float( i2tmp(1,i+xoffset) ) ,
     .                                float( i2tmp(2,i+xoffset) ) )
c
                  else
c
                     ci(i,j) = cmplx( 0.0 , 0.0 )
c
                  endif
c
               enddo
c
            else
c
               do i = 1 , nfa
c
                  ci(i,j) = cmplx( 0.0 , 0.0 )
c
               enddo
c
            endif
c
         enddo
c
      endif
c
c   For left look, conjugate the data so that it conforms to range and
c   doppler coordinates in the same convention as the RDRTec ISAR code.
c
      if ( platdir .gt. 0.0 ) then
c
         do j = 1 , nfrb
c
            do i = 1 , nfa
c
               ci(i,j) = conjg( ci(i,j) )
c
            enddo
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine iocimg ( lunit , irec0 , nfa , nfr , option , ci , ne )
C**
C***********************************************************************
C**
      implicit none
c
      integer lunit , irec0 , nfa , nfr , i , j , option , ne
c
      complex ci(nfa,nfr)
c
c-----------------------------------------------------------------------
c
c   Read or write a section of a complex strip-map image
c
      do i = 1 , nfa
c
         if ( option .eq. 0 ) then
c
            read ( lunit , rec = irec0 + i , err = 200 )
     .                         ( ci(i,j) , j = 1 , nfr )
c
         else if ( option .eq. 1 ) then
c
            write ( lunit , rec = irec0 + i , err = 300 )
     .                         ( ci(i,j) , j = 1 , nfr )
c
         else if ( option .eq. 2 ) then
c
            read ( lunit , rec = irec0 + nfa + 1 - i ,
     .             err = 200 ) ( ci(i,j) , j = 1 , nfr )
c
         else if ( option .eq. 3 ) then
c
            write ( lunit , rec = irec0 + nfa + 1 - i , err = 300 )
     .                         ( ci(i,j) , j = 1 , nfr )
c
         endif
c
      enddo
c
c-----------------------------------------------------------------------
c
c   Normal return
c
      ne = 0
      return
c
c   Read Error
c
  200 continue
      ne = 2
      return
c
c   Write Error
c
  300 continue
      ne = 3
      return
c
c-----------------------------------------------------------------------
c
      end
C**
C***********************************************************************
C**
      subroutine getfrm ( iframe , pframe , ci , nfa , nfr , mode ,
     .                    iunit , ifile , ne )
C**
C***********************************************************************
C**
      implicit none
c
      integer   iframe , pframe , nfa , nfr , mode , iunit , ne ,
     .          offset , i , j , jj , reclda
c
      complex   ci(nfa,nfr)
c
      character ifile*80 , fname*80
c
      if ( mode .eq. 1 ) then
c
         if ( iframe .gt. 0 ) then
c
            call numfile ( ifile , fname , iframe )
c
            open ( 8 , file = fname , form = 'UNFORMATTED' ,
     .             access = 'DIRECT' , recl = reclda( 0 , 2 * nfa ) )
c
            do j = 1 , nfr
c
               read ( 8 , rec = j ) ( ci(i,j) , i = 1 , nfa )
c
            enddo
c
            close ( 8 )
c
         else
c
            open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .             access = 'DIRECT' , recl = reclda( 0 , 2 * nfa ) )
c
            do j = 1 , nfr
c
               jj = j + ( - iframe - 1 ) * nfr
c
               read ( 8 , rec = jj ) ( ci(i,j) , i = 1 , nfa )
c
            enddo
c
            close ( 8 )
c
         endif
c
         ne = 0
c
      else
c
         offset = ( iframe - 1 ) * nfa
         call iocimg ( iunit , offset , nfa , nfr , 2 , ci , ne )
c
      endif
c
      if ( ne .eq. 0 ) then
c
         pframe = abs( iframe )  !  Frame in memory now
c
      else
c
         write ( 6 , * ) ' Frame no. ' , abs( iframe ) , ' not found'
         pframe = 0
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine oframe ( ci , mi , bi , pframe , ifile , nx , ny ,
     .                    nxmid , nymid , dbinc , nlocal , pixbar ,
     .                    color , imgslk , taywta , rlook , flook )
C**
C***********************************************************************
C**
      implicit none
c
      integer   pframe , nx , ny , reclda , i , j , nlocal , pixbar ,
     .          ndot , color , nxmid , nxp , nymid , nyp , ilast , jj ,
     .          nxmid_use , lastdot , imgslk , rlook , flook
c
      complex   ci(nx,ny)
c
      real      mi(nx,ny,2) , dbinc , taywta
c
      character bi(nx,ny)*1
c
      character ifile*80 , ofile*80 , ufile*80
c
      logical   transp
c
      transp = ( ( ifile(1:2) .eq. 't_' ) .or.
     .           ( ifile(3:4) .eq. 't_' ) )
c
      call multiLook ( mi , ci , rlook , flook , nx , ny )
c     
c     do j = 1 , ny
c
c        do i = 1 , nx
c
c           mi(i,j) = cabs( c(i,j,1) ) ** 2
c
c        enddo
c
c     enddo
c
      ndot = lastdot( ifile ) - 1
c
c   Don't number files if frame number is negative
c
      if ( pframe .ge. 0 ) then
c
         call numfile ( ifile , ofile , pframe )
c
         ufile = ofile
c
      else
c
         if ( color .ne. 0 ) then
c
            ofile  = ifile(1:ndot) // '.rmc'
c
         else
c
            ofile  = ifile(1:ndot) // '.rmd'
c
         endif
c
         ufile  = ' '
c
      endif
c
c   Make 8-bit image
c
      call outdbs ( ufile , mi , bi , nx , ny , dbinc , nlocal ,
     .              pixbar , color , imgslk , taywta )
c
c   Write auxiliary information into 8-bit image in the lowest numbers
c
      if ( color .ne. 0 ) call tcolor ( bi , nx , ny )
c
      if ( pframe .lt. 0 ) then
c
         if ( pframe .eq. - 1 ) ilast = 0
c
         close ( 8 )
c
         if ( transp ) then
c
            nxmid_use = min( nx , nxmid )
c
            nxp       = ( nx - nxmid_use ) / 2
c
            open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .                 status = 'UNKNOWN' , access = 'DIRECT' ,
     .                 recl = reclda( ny , 0 ) )
c
            do i = 1 , nxmid_use
c
               write ( 8 , rec = ilast + i )
     .               ( bi(nx+1-i-nxp,j) , j = 1 , ny )
c
            enddo
c
            ilast = ilast + nxmid_use
c
         else
c
            open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .                 status = 'UNKNOWN' , access = 'DIRECT' ,
     .                 recl = reclda( nx * ny , 0 ) )
c
            write ( 8 , rec = - pframe ) bi
c
         endif
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Write complex data out to file if file name starts with 'c_'
c
      if ( index(ifile,'c_') .eq. 1 ) then
c
         if ( pframe .ge. 0 ) then
c
c   Standard floating point file
c
            ofile = 'ci' // ufile(3:)
c
            close ( 9 )
c
            if ( transp ) then
c
               open ( 9 , file = ofile , form = 'UNFORMATTED' ,
     .                    status = 'UNKNOWN' , access = 'DIRECT' ,
     .                    recl = reclda( 0 , 2 * ny ) )
c
               do i = 1 , nx
c
                  write ( 9 , rec = i ) ( ci(i,j) , j = 1 , ny )
c
               enddo
c
            else
c
               open ( 9 , file = ofile , form = 'UNFORMATTED' ,
     .                    status = 'UNKNOWN' , access = 'DIRECT' ,
     .                    recl = reclda( 0 , 2 * nx ) )
c
               do j = 1 , ny
c
                  write ( 9 , rec = j ) ( ci(i,j) , i = 1 , nx )
c
               enddo
c
            endif
c
            close ( 9 )
c
         else
c
            ofile = 'ci' // ifile(3:)
c
            ofile = ofile(1:index(ofile,' ')-1) // '.flt'
c
            open ( 9 , file = ofile , form = 'UNFORMATTED' ,
     .                 status = 'UNKNOWN' , access = 'DIRECT' ,
     .                 recl = reclda( 0 , 2 * nx ) )
c
            do j = 1 , ny
c
               jj = j + ( abs( pframe ) - 1 ) * ny
c
               write ( 9 , rec = jj ) ( ci(i,j) , i = 1 , nx )
c
            enddo
c
c   SANDIA IQ4 file
c
c           ofile = ofile(1:index(ofile,'.')) // 'iq4'
c
c           call putiq4 ( ci , nx , ny , bi , ofile , pframe )
c
         endif
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine outdbs ( ofile , image , imageb , mxpixx , mxpixy ,
     .                    dbinc , nlocal , pixbar , color , imgslk ,
     .                    taywta )
C**
C***********************************************************************
C**
c   This routine converts an intensity image to an eight-bit dB scaled
c   image.
c
c   If dbinc > 0, the image is scaled downwards from the peak value
c
c   If dbinc < 0, the image is scaled upwards from the noise level.  The
c   noise is determined as the smallest average intensity of the four
c   borders.  This value is set to 16, a value which allows 60 dB of
c   dynamic range if dbinc = -0.25.
c
      implicit none
c
      character ofile*80
c
      integer   mxpixx , mxpixy , i , j , ix , nlocal , pixbar ,
     .          color , imgslk
c
      character imageb(mxpixx,mxpixy)*1
c
      real      image(mxpixx,mxpixy) , zmin , zmax , zbar , zsdv ,
     .          db , dbinc , dbinc2 , noise , noisdn , noisup ,
     .          noislt , noisrt , floor , taywta
c
      integer   ii , jj , imin , imax , jmin , jmax , imin2 , imax2 ,
     .          jmin2 , jmax2 , nloclp , nloc2 , nbigx , nbigy ,
     .          pixbp , nerror
c
      real      zlavg , zlmax , zref , zrefm , bigbox
c
c   Calculate global min, max, mean, std. dev.
c
      call stats ( image , mxpixx * mxpixy , zmin , zmax , zbar , zsdv )
c
      if ( zbar .eq. 0.0 .and. zsdv .eq. 0 ) then
c
         floor = 1.0
c
      else
c
         floor = 0.1 * ( zbar ** 2 ) / amax1( zsdv , abs( zbar ) )
c
      endif
c
      if ( dbinc .lt. 0.0 .or. imgslk .ne. 0 ) then
c
         noisdn = 0.0
         noisup = 0.0
         noislt = 0.0
         noisrt = 0.0
c
c   Define the noise as the logarithmic (geometric) mean of the
c   intensity around the boundary.  However, ignore very small values
c   by defining a floor value derived from the large scale stats.
c
         do i = 1 , mxpixx
c
            noisdn = noisdn + alog( amax1( floor , image(i,1) ) )
            noisup = noisup + alog( amax1( floor , image(i,mxpixy) ) )
c
         enddo
c
         do j = 1 , mxpixy
c
            noislt = noislt + alog( amax1( floor , image(1,j) ) )
            noisrt = noisrt + alog( amax1( floor , image(mxpixx,j) ) )
c
         enddo
c
         noisup = exp( noisup / float( mxpixx ) ) ** 0.25
         noisdn = exp( noisdn / float( mxpixx ) ) ** 0.25
         noislt = exp( noislt / float( mxpixy ) ) ** 0.25
         noisrt = exp( noisrt / float( mxpixy ) ) ** 0.25
c
         noise  = noisdn * noisup * noisrt * noislt
         noise  = 0.5 * ( noise + floor )
c
      endif
c
      if ( imgslk .ne. 0 ) then
         call reduce_sidelobes( image , noise , mxpixx , mxpixy ,
     .                          taywta )
      endif
      
      if ( abs( dbinc ) .lt. 99.0 ) then
c
c   Use global stats
c
         dbinc2 = abs( dbinc )
c
         do j = 1 , mxpixy
c
            do i = 1 , mxpixx
c
               if ( image(i,j) .gt. 0.0 ) then
c
                  if ( dbinc .gt. 0 ) then
c
                     ix = 255 + nint( db( image(i,j) / zmax ) / dbinc2 )
c
                  else
c
                     ix = 32 + nint( db( image(i,j) / noise ) / dbinc2 )
c
                  endif
c
               else
c
                  ix = 1
c
               endif
c
               ix = min( ix , 255 )
               ix = max( ix ,   0 )
c
               imageb(i,j) = char( ix )
c
            enddo
c
         enddo
c
      else
c
c   Use local stats
c
         if ( nlocal .lt. 2 ) then
c
            nloclp = nint( 0.5 * float( mxpixx * mxpixy ) ** 0.25 )
            nloclp = min( nloclp , mxpixx )
            nloclp = min( nloclp , mxpixy )
c
         else
c
            nloclp = nlocal
c
         endif
c
c   Limit mean pixel intensity to a reasonable range
c
         pixbp = pixbar
         if ( pixbar .lt. 40  ) pixbp = 40
         if ( pixbar .gt. 216 ) pixbp = 216
c
         nloc2  = nloclp / 2
c
         nbigx  = mxpixx / nloclp
         if ( mod( mxpixx , nloclp ) .ne. 0 ) nbigx = nbigx + 1
c
         nbigy  = mxpixy / nloclp
         if ( mod( mxpixy , nloclp ) .ne. 0 ) nbigy = nbigy + 1
c
c   Loop over boxes of size nloclp*nloclp.  Define a larger box 4 times
c   the size to calculate the local intensity for scaling.
c
c-----------------------------------------------
c
c        Define loop
c
         jmin  = 1 - nloclp
         jmax  = 0
c
         do 203 jj = 1 , nbigy
c
         jmin  = jmin + nloclp
         jmax  = min( jmax + nloclp , mxpixy )
c
         jmin2 = jmin - nloc2
         jmax2 = jmax + nloc2
c
         if ( jmin2 .lt. 1 ) then
c
            jmin2 = 1
            jmax2 = nloclp + 2 * nloc2
c
         endif
c
         if ( jmax2 .gt. mxpixy ) then
c
            jmax2 = mxpixy
            jmin2 = mxpixy + 1 - nloclp - 2 * nloc2
c
         endif
c
         imin  = 1 - nloclp
         imax  = 0
c
         do 203 ii = 1 , nbigx
c
            imin  = imin + nloclp
            imax  = min( imax + nloclp , mxpixx )
c
            imin2 = imin - nloc2
            imax2 = imax + nloc2
c
            if ( imin2 .lt. 1 ) then
c
               imin2 = 1
               imax2 = nloclp + 2 * nloc2
c
            endif
c
            if ( imax2 .gt. mxpixx ) then
c
               imax2 = mxpixx
               imin2 = mxpixx + 1 - nloclp - 2 * nloc2
c
            endif
c
c        Define loop
c
c-----------------------------------------------
c
c   Compute local average
c
            zlavg = 0.0
            zlmax = 0.0
c
            do j = jmin2 , jmax2
c
               do i = imin2 , imax2
c
                  zlavg = zlavg + alog( amax1( floor , image(i,j) ) )
                  zlmax = amax1( zlmax , image(i,j) )
c
               enddo
c
            enddo
c
c   Subtract the doughnut hole
c
            do j = jmin , jmax
c
               do i = imin , imax
c
                  zlavg = zlavg - alog( amax1( floor , image(i,j) ) )
c
               enddo
c
            enddo
c
            bigbox = float( nloclp + 2 * nloc2 ) ** 2 -
     .               float( ( imax + 1 - imin ) * ( jmax + 1 - jmin ) )
c
            zlavg = exp( zlavg / bigbox )
c
c   Define the local reference as a combination of the local and global
c   values
c
            zref  = sqrt( zlavg * zbar )
            zrefm = sqrt( zlmax * zmax )
c
            if ( dbinc .gt. 0.0 ) then
c
               dbinc2 = db( zrefm / zref ) / float( 255 - pixbp )
c
            else
c
               dbinc2 = db( zref / noise ) / float( pixbp - 32 )
c
            endif
c
            dbinc2 = amax1( dbinc2 , 1.0 / 256.0 )
c
            do j = jmin , jmax
c
               do i = imin , imax
c
                  if ( image(i,j) .gt. 0.0 ) then
c 
                     ix = pixbp +
     .                    nint( db( image(i,j) / zref ) / dbinc2 )
c
                  else
c
                     ix = 1
c
                  endif
c
                  ix = min( ix , 255 )
                  ix = max( ix ,   0 )
c
                  if ( color .ne. 0 ) then
c
                     imageb(i,j) = char( 4 * ( ix / 4 ) )
c
                  else
c
                     imageb(i,j) = char( ix )
c
                  endif
c
               enddo
c
            enddo
c
c        End of loop
c
  203    continue  !  Loop over boxes
c
      endif
c
      if ( ofile .ne. ' ' ) call dabyte ( 'w' , imageb , mxpixx ,
     .                                     mxpixy , ofile , nerror )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine dabyte ( option , array , ni , nj , ofile , nerror )
C**
C***********************************************************************
C**
c
c   Outputs a direct access file of a byte array
c
c        array  : data array contain values to be read or written
c
c        ofile  : name of file
c
c        ni     : fast varying index
c
c        nj     : slow varying index
c
c        option : allows user to select a read or write option
c
      implicit none
c
      character ofile*80
      character option*1
      integer   ni , nj , nerror
c
c   data array
c
      character array(ni,nj)*1
c
      integer   i , j , reclda
c
      logical   transp
c
      open ( unit = 8 , file = ofile , form = 'UNFORMATTED' ,
     .       status = 'UNKNOWN' , access = 'DIRECT' ,
     .       recl = reclda( ni , 0 ) , err = 2000 )
c
      transp = ( ofile(1:2) .eq. 't_' )
c
      if ( option .eq. 'r' .or. option .eq. 'R' ) then
c
         if ( transp ) then
c
            do i = 1 , ni
c
               read ( 8 , rec = i , err = 1000 )
     .              ( array(i,j) , j = 1 , nj )
c
            enddo
c
         else
c
            do j = 1 , nj
c
               read ( 8 , rec = j , err = 1000 )
     .              ( array(i,j) , i = 1 , ni )
c
            enddo
c
         endif
c
      elseif ( option .eq. 'w' .or. option .eq. 'W' ) then
c
         if ( transp ) then
c
            do i = 1 , ni
c
               write ( 8 , rec = i , err = 1000 )
     .               ( array(i,j) , j = 1 , nj )
c
            enddo
c
         else
c
            do j = 1 , nj
c
               write ( 8 , rec = j , err = 1000 )
     .               ( array(i,j) , i = 1 , ni )
c
            enddo
c
         endif
c
      endif
c
      close ( 8 )
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      nerror = 0  !  Normal return
c
      return
c
 1000 nerror = 1  !  Error return - end of file encountered
c
      write ( 6 , * ) ' End of file encountered in DABYTE'
      write ( 7 , * ) ' End of file encountered in DABYTE'
c
      return
c
 2000 nerror = 2  !  Error return - file I/O error on open or read
c
      write ( 6 , * ) ' File I/O error encountered in DABYTE'
      write ( 7 , * ) ' File I/O error encountered in DABYTE'
c
      return
      end
C**
C***********************************************************************
C**
      subroutine bmovie ( echo , bytes , bframe , ofile , nf , nr )
C**
C***********************************************************************
C**
      implicit none
c
      integer   bframe , nf , nr , i , j , ndot , lastdot , reclda
c
      character ofile*80 , ofile2*80
c
      real      echo(nf,nr) , emax
c
      character bytes(nf,nr)*1
c
      emax = 0.0
c
      do j = 1 , nr
c
         do i = 1 , nf
c
            emax = amax1( emax , abs( echo(i,j) ) )
c
         enddo
c
      enddo
c
      do j = 1 , nr
c
         do i = 1 , nf
c
            bytes(i,j) = char( 128 +
     .                         nint( ( 127.0 / emax ) * echo(i,j) ) )
c
         enddo
c
      enddo
c
      ndot   = lastdot( ofile ) - 1
      ofile2 = ofile(1:ndot) // '.rmd'
c
      close ( 8 )
c
      open ( 8 , file   = ofile2    , form   = 'unformatted' ,
     .           status = 'unknown' , access = 'direct' ,
     .           recl   = reclda( nf * nr , 0 ) )
c
      write ( 8 , rec = bframe ) bytes
c
      close ( 8 )
c
      return
      end
C**
C***********************************************************************
C**
      subroutine apack6 ( rdwrt , r , npts , iofile , nerror ,
     .                    rmin , rmax , rbar , rsdv , scale )
C**
C***********************************************************************
C**
      implicit none
c
      character rdwrt*1
c
      character iofile*80
c
      integer   npts , i , nerror , nptsin
c
      real      r(npts) , rmin , rmax , rbar , rsdv , scale
c
      if ( rdwrt .eq. 'r' .or. rdwrt .eq. 'R' ) then
c
         open ( 8 , file = iofile , form = 'FORMATTED' ,
     .              status = 'OLD' , err = 1000 )
c
         read ( 8 , * ) nptsin , rmin , rmax , rbar , rsdv , scale
c
         read ( 8 , '(8f8.0)' , end = 1000 , err = 1000 ) r
c
         do i = 1 , npts
c
            r(i) = r(i) / scale
c
         enddo
c
      else
c
         open ( 8 , file = iofile , form = 'FORMATTED' ,
     .              status = 'UNKNOWN' , err = 1000 )
c
         call stats ( r , npts , rmin , rmax , rbar , rsdv )
c
         scale = 999999.0 / amax1( abs( rmax ) , abs( rmin ) )
c
         write ( 8 , * ) npts , rmin , rmax , rbar , rsdv , scale
c
         write ( 8 , '(8i8)' ) ( nint( scale * r(i) ) , i = 1 , npts )
c
      endif
c
      nerror = 0
c
 1000 nerror = 1
c
      return
      end
C**
C***********************************************************************
C**
      subroutine putiq4 ( c , nx , ny , i2_b , ofile , kframe )
C**
C***********************************************************************
C**
c   Write a Sandia IQ4 file from a complex image
c
      implicit none
c
      integer   nx , ny , i , j , k , reclda , kframe , offset ,
     .          iq , lowbyte , highbyte , itemp
c
      complex   c(nx,ny)       !  Complex image
c
      character i2_b(2,2,nx)*1    !  8-bit array
c
      character ofile*80          !  Output IQ4 file
c
      real      iqmax , scale  !  Max I/Q and scale factor
c
      real      rtemp
c
c-----------------------------------------------------------------------
c
c   Determine the scale factor from the max real or imaginary part of
c   the complex image
c
      iqmax = 0.0
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            iqmax = amax1( iqmax , abs(  real( c(i,j) ) ) )
            iqmax = amax1( iqmax , abs( aimag( c(i,j) ) ) )
c
         enddo
c
      enddo
c
      scale = 32767.0 / iqmax
c
c-----------------------------------------------------------------------
c
c   Output the 16-bit image values
c
      open ( 8 , file = ofile , form = 'UNFORMATTED' ,
     .           status = 'UNKNOWN' , access = 'DIRECT' ,
     .           recl = reclda( 4 * nx , 0 ) )
c
      if ( kframe .eq. 0 ) then
c
         offset = 0
c
      else
c
         offset = ( abs( kframe ) - 1 ) * ny
c
      endif
c
      do j = 1 , ny
c
         do i = 1 , nx
c
            do iq = 1 , 2
c
               if ( iq .eq. 1 ) then
                  rtemp =  real( c(i,j) )
               else
                  rtemp = aimag( c(i,j) )
               endif
c
               rtemp        = scale * rtemp
c
               rtemp        = amin1( + 32767.0 , rtemp )
               rtemp        = amax1( - 32767.0 , rtemp )
c
               itemp        = nint( rtemp )
c
               if ( itemp .lt. 0 ) itemp = 65536 + itemp
c
               lowbyte      = itemp / 256
c
               highbyte     = itemp - 256 * lowbyte
c
               i2_b(1,iq,i) = char( lowbyte )
               i2_b(2,iq,i) = char( highbyte )
c
            enddo
c
         enddo
c
         write ( 8 , rec = j + offset) ( ( ( i2_b(k,iq,i) , k = 1 , 2 )
     .                                     , iq = 1 , 2 ) , i = 1 , nx )
c
      enddo
c
      close ( 8 )
c
c-----------------------------------------------------------------------
c
      return
      end
C**
C***********************************************************************
C**
      subroutine packto8 ( ri , nx , ny , outfil , bytes , iopt )
C**
C***********************************************************************
C**
      implicit none
c
      integer   nx , ny , iopt
c
      real      ri(2,nx,ny)                !  Complex, stored as real
c
      character outfil*80
c
      byte      bytes(nx*ny)
c
      real      xmin , xmax , xbar , xsdv , scale 
c
      integer   jr , ji , i , j , outlen , reclda
c
      call stats ( ri , 2 * nx * ny , xmin , xmax , xbar , xsdv )

      write ( 6 , '(2x,4f13.3)' ) xmin , xmax , xbar , xsdv
      write ( 7 , '(2x,4f13.3)' ) xmin , xmax , xbar , xsdv
c
      scale  = 127.0 /  amax1( abs( xmax ) , abs( xmin ) )
c
      if ( iopt .eq. 0 ) then
c
         outlen = 2 * nx          !  Write out in normal order
c
      else
c
         outlen = 2 * ny          !  Write out in transpose order
c
      endif
c
      open ( 9 , file = outfil , form = 'UNFORMATTED' ,
     .           status = 'UNKNOWN' , access = 'DIRECT' ,
     .           recl = reclda( outlen , 0 ) )
c
      if ( iopt .eq. 0 ) then
c
         do j = 1 , ny
c
            do i = 1 , nx
c
               jr = 2 * i - 1
               ji = 2 * i
c
               bytes(jr) = min( 127 , max( - 127 ,
     .                             nint( scale * ri(1,i,j) ) ) )
               bytes(ji) = min( 127 , max( - 127 ,
     .                             nint( scale * ri(2,i,j) ) ) )
c
            enddo
c
            write ( 9 , rec = j ) ( bytes(i) , i = 1 , 2 * nx )
c
         enddo
c
      else
c
         do i = 1 , nx
c
            do j = 1 , ny
c
               jr = 2 * j - 1
               ji = 2 * j
c
               bytes(jr) = min( 127 , max( - 127 ,
     .                                     nint( scale * ri(1,i,j) ) ) )
               bytes(ji) = min( 127 , max( - 127 ,
     .                                     nint( scale * ri(2,i,j) ) ) )
c
            enddo
c
            write ( 9 , rec = i ) ( bytes(j) , j = 1 , 2 * ny )
c
         enddo
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine put16bit ( r , nx , i2_b , scale , outunit , recnum )
C**
C***********************************************************************
C**
c   Write 16-bit records of real numbers on a computer without integer*2
c
c   r(nx)      :  Array of real numbers to be written to file
c
c   nx         :  Size of array
c
c   i2_b(2,nx) :  Character*1 array for 16-bit data
c
c   scale      :  Scale factor for conversion to 16 bits
c
c   outunit    :  Unit number for output to direct access file
c
c   recnum     :  Record number for output
c
c-----------------------------------------------------------------------
c
c   Notes:
c
c      1.  Unit 'outunit' must have been opened with the following
c          statement:
c
c            open ( outunit , file = yourfile , form = 'unformatted' ,
c                   access = 'direct' , recl = 2 * nx )
c
c      2.  The form of the output is Sun/SGI integer*2.  If the binary]
c          file is transferred to an IBM-PC the low and high order bytes
c          must be swapped.
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer   nx , i , outunit , recnum , itemp , lowbyte , highbyte
c
      real      r(nx) , scale , rtemp
c
      character i2_b(2,nx)*1             !  8-bit array
c
c-----------------------------------------------------------------------
c
      do i = 1 , nx
c
         rtemp     = 32767.0 * scale * r(i)
c
         rtemp     = amin1( + 32767.0 , rtemp )
         rtemp     = amax1( - 32767.0 , rtemp )
c
         itemp     = nint( rtemp )
c
         if ( itemp .lt. 0 ) itemp = 65536 + itemp
c
         lowbyte   = itemp / 256
c
         highbyte  = itemp - 256 * lowbyte
c
         i2_b(1,i) = char( lowbyte )
         i2_b(2,i) = char( highbyte )
c
      enddo
c
c   So that IQ4 files are defined according to the Sun integer*2 byte
c   order, swap the bytes if we have detected the IBM-PC order
c
      write ( outunit , rec = recnum ) i2_b
c
c-----------------------------------------------------------------------
c
      return
      end
c**
c***********************************************************************
c**
      subroutine numfile ( infile , ifile , k )
c**
c***********************************************************************
c**
      implicit none
c
      character infile*80 , ifile*80
c
      integer   k , nblank , ntest
c
      character num*4
c
c   Make frame number string to append to file name
c
      num = '0000'
c 
      if ( k .lt. 10                      ) write ( num(4:4) , '(i1)' )k
      if ( k .ge. 10   .and. k .lt. 100   ) write ( num(3:4) , '(i2)' )k
      if ( k .ge. 100  .and. k .lt. 1000  ) write ( num(2:4) , '(i3)' )k
      if ( k .ge. 1000 .and. k .lt. 10000 ) write ( num(1:4) , '(i4)' )k
c
      ntest  = 80
c
      nblank = 0
c
      do while ( ntest .ne. 0 )
c
         if ( infile(ntest:ntest) .eq. ' ' ) then
c
            ntest = ntest - 1
c
         else
c
            nblank = ntest + 1
c
            ntest  = 0
c
         endif
c
      enddo
c
      if ( infile(nblank-1:nblank-1) .eq. '.' ) then
c
         ifile  = infile(1:nblank-1) // num
c
         if ( k .eq. 0 ) ifile = ifile(1:nblank-1+4) // '.raw'
c
      else
c
         ifile  = infile(1:nblank-1) // '.' // num
c
         if ( k .eq. 0 ) ifile = ifile(1:nblank-1+5) // '.raw'
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      subroutine sun_flt ( nx , ny , work , infile , outfile )
c**
c***********************************************************************
c**
c   Convert a file from Sun floating point to PC floating point format
c
      implicit none
c
      character infile*80 , outfile*80
c
      integer   nx , ny , j , reclda
c
      real      work(nx) , iqmin , iqmax , iqbar , iqsdv
c
      open ( 51 , file    = infile   , form = 'unformatted' ,
     .            access  = 'direct' , recl = reclda( 0 , nx ) ,
     .            convert = 'big_endian' )
c
      open ( 52 , file    = outfile  , form = 'unformatted' ,
     .            access  = 'direct' , recl = reclda( 0 , nx ) )
c
      do j = 1 , ny
c
         read  ( 51 , rec = j ) work
c
         write ( 52 , rec = j ) work
c
         call stats ( work , nx , iqmin , iqmax , iqbar , iqsdv )
c
         write (  6 , '(i6,9f8.2)' ) j , iqmin , iqmax , iqbar , iqsdv
c
         write ( 16 , '(i6,9f8.2)' ) j , iqmin , iqmax , iqbar , iqsdv
c
      enddo
c
      close ( 51 )
c
      close ( 52 )
c
      return
      end
c**
c***********************************************************************
c**
      subroutine decode_double ( raw , float8 )
c**
c***********************************************************************
c**
c   Convert between Sun IEEE double precision floating point and
c   PC double precision floating point
c
      implicit none
c
      integer       j
c
      byte          raw(8)
c
      real*8        float8
c
      real*8        temp
c
      byte          temp_b(8) 
c
      equivalence ( temp , temp_b )
c
c   Read the bytes backwards into the temporary variable
c
      do j = 1 , 8
c
         temp_b(j)   = raw(9-j)
c
      enddo
c
c   Transfer to output argument
c
      float8 = temp
c
      return
      end
c**
c***********************************************************************
c**
      subroutine sun2flt ( ibyt , obyt , n )
c**
c***********************************************************************
c**
c   Convert between Sun IEEE floating point and PC floating point
c
c      ibyt  :  Bytes of the input floating point data
c
c      obyt  :  Bytes of the output floating point data
c
c      n     :  Number of 4-byte floating point values
c
c   Notes:
c
c        1.  The input and output arrays can be the same
c
c        2.  This routine is its own inverse - it also converts
c            back to the Sun format
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer n , k
c
      byte    ibyt(4,n) , obyt(4,n) , btmp(4)
c
      do k = 1 , n
c
c   Use temporary variable in case the two arrays are the same
c
         btmp(:)   = ibyt(:,k)
c
c   Big-Endian to Little-Endian conversion
c
         obyt(1,k) = btmp(4)
c
         obyt(2,k) = btmp(3)
c
         obyt(3,k) = btmp(2)
c
         obyt(4,k) = btmp(1)
c
      enddo
c
      return
      end
c**
c***********************************************************************
c**
      subroutine hex2int ( hex , nhex , a )
c**
c***********************************************************************
c**
      implicit none
c
      character hex*80
c
      integer   nhex , a(nhex) , j , temp
c
      do j = 1 , nhex
c
         temp = ichar( hex(j:j) )
c
         if ( temp .gt. 47 .and. temp .lt. 58 ) then
c
c                 Digits 0-9
c
            a(j) = temp - 48
c
         else if ( temp .gt. 64 .and. temp .lt. 71 ) then
c
c                 Digits A-F
c
            a(j) = temp - 55
c
         else
c
            write ( 6 , * ) ' Non-HEX digit found!' , j ,
     .                        temp , ' , ' , hex
            read ( 5 , * )
c
            a(j) = 0
c
         endif
c
      enddo
c
      return
      end
c**
c***********************************************************************
c**
      integer function ibmi2 ()
c**
c***********************************************************************
c**
      implicit none
c
      integer*2     i2
c
      character     b(2)*1
c
      equivalence ( i2 , b )
c
      logical       first / .true. /
c
      integer       ibmpc / 0      /
c
      save
c
c-----------------------------------------------------------------------
c
c   Determine the order of bytes in an INTEGER*2 variable.  The control
c   integer 'ibmpc' is set to one for a PC; 0 otherwise.  This value is
c   used in reading images in SANDIA's IQ4 format.  It is assumed that
c   the files were written as INTEGER*2 (Short) variables on a SUN
c   workstation or equivalent computer.  If ibmpc=1, then the routine
c   'readci' swaps the bytes before converting the integers to floats.
c
      if ( first ) then
c
         i2 = 1
c
         write ( 6 , * ) i2 , ichar( b(1) ) , ichar( b(2) )
c
         if ( ichar( b(1) ) .eq. 1 ) then
c
            ibmpc = 1
c
            write ( 6 , * ) ' IBM-PC Integer format detected'
c
         else
c
            ibmpc = 0
c
            write ( 6 , * ) ' SUN Integer format detected'
c
         endif
c
         first = .false.
c
      endif
c
      ibmi2 = ibmpc
c
      return
      end
c
c***********************************************************************
c
      integer function lastdot ( ch80var )
c
c***********************************************************************
c
c   Purpose:  Find the last dot in a file name
c
      implicit none
c
      integer   ndot , ntest , nblank
c
      character ch80var*80
c
      ndot    = index( ch80var , '.' )
      lastdot = ndot
c
      do while ( ndot .ne. 0 )
c
         ndot = index( ch80var(lastdot+1:80) , '.' )
c
         if ( ndot .gt. 0 ) lastdot = lastdot + ndot
c
      enddo
c
      if ( lastdot .le. 0 ) then
c
         ntest  = 80
c
         nblank = 0
c
         do while ( ntest .ne. 0 )
c
            if ( ch80var(ntest:ntest) .eq. ' ' ) then
c
               ntest = ntest - 1
c
            else
c
               nblank = ntest + 1
c
               ntest  = 0
c
            endif
c
         enddo
c
         lastdot = max( 1 , nblank )
c
      endif
c
      return
      end
