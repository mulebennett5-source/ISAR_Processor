c
c***********************************************************************
c
c            File 'mdiags.f' - Routines for conveying diagnostic
c                              information from subroutine 'movie'
c                              to the console or to the log file
c
c***********************************************************************
c
C**
C***********************************************************************
C**
      subroutine movie2 ( drc , rsmin , rsmax , drf , rfmin , rfmax ,
     .                    mprat , dfc , dff , dfmp )
C**
C***********************************************************************
C**
c   Write out parameters of the finite difference grid to be used in
c   the image formation process
c
      implicit none
c
      real    drc , rsmin , rsmax , drf , rfmin , rfmax , dfc , dff ,
     .        dfmp
c
      integer mprat
c
c   Write to screen
c
      write ( 6 , '(/,6x,a32,f12.4)' )
     .        ' Sub-image range sampling (m) ' , drc
      write ( 6 , '(6x,a32,f12.2)' )
     .        ' Sub-image minimum range (m)  ' , rsmin
      write ( 6 , '(6x,a32,f12.2)' )
     .        ' Sub-image maximum range (m)  ' , rsmax
c
      write ( 6 , '(/,6x,a32,f12.4)' )
     .        ' Final range sampling (m)     ' , drf
      write ( 6 , '(6x,a32,f12.2)' )
     .        ' Final image minimum range (m)' , rfmin
      write ( 6 , '(6x,a32,f12.2)' )
     .        ' Final image maximum range (m)' , rfmax
c
      write ( 6 , '(/,6x,a32,i12)' )
     .        ' Master Particle ratio        ' , mprat
      write ( 6 , '(6x,a32,f12.3)' )
     .        ' Coarse Doppler sampling (Hz) ' , dfc
      write ( 6 , '(6x,a32,f12.3)' )
     .        ' Fine Doppler sampling (Hz)   ' , dff
      write ( 6 , '(6x,a32,f12.3,//)' )
     .        ' Master Particle sampling (Hz)' , dfmp
c
c   Write to log file
c
      write ( 7 , '(/,6x,a32,f12.4)' )
     .        ' Sub-image range sampling (m) ' , drc
      write ( 7 , '(6x,a32,f12.2)' )
     .        ' Sub-image minimum range (m)  ' , rsmin
      write ( 7 , '(6x,a32,f12.2)' )
     .        ' Sub-image maximum range (m)  ' , rsmax
c
      write ( 7 , '(/,6x,a32,f12.4)' )
     .        ' Final range sampling (m)     ' , drf
      write ( 7 , '(6x,a32,f12.2)' )
     .        ' Final image minimum range (m)' , rfmin
      write ( 7 , '(6x,a32,f12.2)' )
     .        ' Final image maximum range (m)' , rfmax
c
      write ( 7 , '(/,6x,a32,i12)' )
     .        ' Master Particle ratio        ' , mprat
      write ( 7 , '(6x,a32,f12.3)' )
     .        ' Coarse Doppler sampling (Hz) ' , dfc
      write ( 7 , '(6x,a32,f12.3)' )
     .        ' Fine Doppler sampling (Hz)   ' , dff
      write ( 7 , '(6x,a32,f12.3,//)' )
     .        ' Master Particle sampling (Hz)' , dfmp
c
      return
      end
C**
C***********************************************************************
C**
      subroutine timing ( isub , iframe , work , nwork , frrate , opt )
C**
C***********************************************************************
C**
c   Write out information of timing of processes in the ISAR movie
c   formation
c
      implicit none
c
      include 'motime.h'
c
      include 'realtime.h'
c
      integer  isub , iframe , nwork , nft , opt
c
      real     totime , inctim , a , g , h , frrate ,
     .         work(nwork) , speed , mflops , time0
c
      if ( opt .eq. 0 ) then
c
c   Initialize timing variables
c
         time0  = dtime( timer )    !  Initialize timer
         actime = 0.0               !  Acquisition
         rctime = 0.0               !  Range compression
         sitime = 0.0               !  Sub-image formation
         frtime = 0.0               !  Frame generation ,
         pgtime = 0.0               !  Phase Gradient Autofocus
         tgtime = 0.0               !  Target Detection
         iqtime = 0.0               !  Data Fetching or Generating
         ohtime = 0.0               !  Overhead functions
c
      else
c
         write ( 6 , '(//,a42,/)' )
     .        '     --------  Movie is Complete  --------'
         write ( 7 , '(//,a42,/)' )
     .        '     --------  Movie is Complete  --------'
c
         write ( 6 , '(/,2x,a34,i12)' )
     .      ' Sub-Images Computed:             ' , isub - 1 ,
     .      ' Fine Resolution Frames Computed: ' , iframe
         write ( 7 , '(/,2x,a34,i12)' )
     .      ' Sub-Images Computed:             ' , isub - 1 ,
     .      ' Fine Resolution Frames Computed: ' , iframe
c
         totime = ohtime + iqtime + actime + rctime + tgtime + sitime +
     .            frtime + pgtime
         totime = amax1( totime , 1.0E-6 )
         a      = 100.0 / totime
c
         write ( 7 , '(//,a43,/)' )
     .        '     --------  Timing ( sec and % )  --------'
c
         write ( 7 , '(2x,a29,2f12.2)' )
     .      ' Overhead:                   ' , ohtime , a * ohtime ,
     .      ' Data Fetching:              ' , iqtime , a * iqtime ,
     .      ' Acquisition:                ' , actime , a * actime ,
     .      ' Range Compression:          ' , rctime , a * rctime ,
     .      ' Target Detection:           ' , tgtime , a * tgtime ,
     .      ' Sub-image Formation:        ' , sitime , a * sitime ,
     .      ' Frame Generation:           ' , frtime , a * frtime ,
     .      ' Phase Gradient Autofocus:   ' , pgtime , a * pgtime ,
     .      ' Total:                      ' , totime , 100.0
c
         g      = 1.0 / float( max(1,iframe) )
         inctim = g * ( rctime + tgtime + sitime + frtime + pgtime )
         inctim = amax1( inctim , 1.0E-6 )
         h      = 100.0 * g / inctim
c
         write ( 7 , '(//,a51,/)' )
     .        '     -----  Incremental Timing ( sec and % )  -----'
c
         write ( 7 , '(2x,a29,2f12.2)' )
     .      ' Range Compression:          ' , g * rctime , h * rctime ,
     .      ' Target Detection:           ' , g * tgtime , h * tgtime ,
     .      ' Sub-image Formation:        ' , g * sitime , h * sitime ,
     .      ' Frame Generation:           ' , g * frtime , h * frtime ,
     .      ' Phase Gradient Autofocus:   ' , g * pgtime , h * pgtime ,
     .      ' Total:                      ' , inctim     , 100.0
c
         nft   = 256
c
         speed = mflops( work , nft )
c
         write ( 7 , '(//,a44,f11.2)' )
     .        '   CPU Speed (MFLOPS):                       ' , speed
c
         speed = frrate * inctim * speed
c
         if ( rt_img .ne. 0 ) speed = 2.0 * speed
c
         write ( 6 , '(/,a44,f11.2,//)' )
     .        '   Speed Required for real-time (MFLOPS):    ' , speed
c
         write ( 7 , '(/,a44,f11.2,//)' )
     .        '   Speed Required for real-time (MFLOPS):    ' , speed
c
      endif
c
      return
      end
