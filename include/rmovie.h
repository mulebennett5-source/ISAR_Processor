c
c-----------------------------------------------------------------------
c
c                   Rmovie variables and constants
c
      integer       maxmaxx , maxmaxy , maxx , maxy	, movie2_maxx ,
     .              movie_maxx , movie_maxy , ascan_minheight
c
c   Biggest movie frame allowed
c
      parameter   ( maxmaxx         = 2048 )
      parameter   ( maxmaxy         = 1080 )
c
c   Set Maximum window sizes
c
      parameter   ( movie_maxx      = maxmaxx + 256 )
      parameter   ( movie_maxy      = maxmaxy + 576 )
c
      parameter   ( ascan_minheight = 50  )
      parameter   ( movie2_maxx     = movie_maxx / 2 )
c
      logical       do_sct
c
      common      / screenlimits / maxx , maxy , do_sct
c
c-----------------------------------------------------------------------
c
c   Structure for a window definition
c      
      structure / win_def / 
c
c   Pixel coordinates
c
         integer*2 x1_phys
         integer*2 x2_phys
         integer*2 y1_phys
         integer*2 y2_phys
c
c   Floating point counterparts of the pixel coordinates
c
         real*8    x1_win
         real*8    x2_win
         real*8    y1_win
         real*8    y2_win
c
c   Invert the vertical axis or not
c
         logical*2 finvert 
c
      end structure     
c
c-----------------------------------------------------------------------
c
c   Color table offset
c
      integer                 ncoff
      parameter             ( ncoff = 20 )
c
c   Size of IPR plot window
c
      integer                 ipr_total
      parameter             ( ipr_total = 41 )
c
      integer                 single , naptime , bottom , top ,
     .                        ctable(0:127+ncoff) , ctable4(0:15)
c
      logical                 newall , plotipr , bmpmovie , rootbmp ,
     .                        tcolor
c
      real                    quality(13,100000)
c
      common / rmovie_parms / single , naptime , bottom , top , newall ,
     .                        plotipr , bmpmovie , rootbmp , tcolor ,
     .                        ctable , ctable4 , quality
c
c-----------------------------------------------------------------------
c
