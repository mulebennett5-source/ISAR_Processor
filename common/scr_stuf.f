C**
C***********************************************************************
C**
      subroutine set_cur_win ( cur_win )
C**
C***********************************************************************
C**
      use dflib
c
      implicit none
c
      include           'rmovie.h'
c
      integer*2          status
c
      record / win_def / cur_win
c
c-----------------------------------------------------------------------
c
c   CALL SETVIEWPORT (x1, y1, x2, y2) 
c
c                     x1, y1 (Input) INTEGER(2). Physical coordinates
c                                    for upper-left corner of viewport. 
c
c                     x2, y2 (Input) INTEGER(2). Physical coordinates
c                                    for lower-right corner of viewport. 
c
c-----------------------------------------------------------------------
c
      if ( cur_win.finvert ) then
c
         call setviewport ( cur_win.x1_phys , cur_win.y2_phys ,
     .                      cur_win.x2_phys , cur_win.y1_phys )
c
      else
c
         call setviewport ( cur_win.x1_phys , cur_win.y1_phys ,
     .                      cur_win.x2_phys , cur_win.y2_phys )
c
      endif
c
c-----------------------------------------------------------------------
c
c   result = SETWINDOW (finvert, wx1, wy1, wx2, wy2) 
c
c                       finvert (Input) LOGICAL(2). Direction of
c                               increase of the y-axis. If finvert is
c                               .TRUE., the y-axis increases from the
c                               window bottom to the window top (as
c                               Cartesian coordinates). If finvert is
c                               .FALSE., the y-axis increases from the
c                               window top to the window bottom (as
c                               pixel coordinates). 
c
c
c                       wx1, wy1 (Input) REAL(8). Window coordinates
c                                for upper-left corner of window. 
c
c                       wx2, wy2 (Input) REAL(8). Window coordinates
c                                for lower-right corner of window.
c
c-----------------------------------------------------------------------
c
      status = setwindow( cur_win.finvert ,
     .                    cur_win.x1_win , cur_win.y1_win ,
     .                    cur_win.x2_win , cur_win.y2_win )
c
      if ( status .eq. 0 ) stop ' Cannot set window'
c
c-----------------------------------------------------------------------
c
      return
      end
C**
C***********************************************************************
C**
      subroutine setup_video ( myscreen , w , h )
C**
C***********************************************************************
C**
      use dflib
c
      implicit none
c
      include                'rmovie.h'
c
      record / windowconfig / myscreen
c
      integer                 colidx
c
      integer                 r , g , b
c
      byte                    col148(4,128+ncoff) , col16(4,16)
c
      integer                 status , col , nstart
c
      logical                 wstatus
c
      integer*2               w , h
c
c-----------------------------------------------------------------------
c
      myscreen.numxpixels  = maxx
      myscreen.numypixels  = maxy
      myscreen.numtextcols = w
      myscreen.numtextrows = h
      myscreen.numcolors   = -1
      myscreen.title       =
     .   " ISAR - John R. Bennett "C
      myscreen.fontsize    = - 1
c
      wstatus = setwindowconfig( myscreen )
c
      if ( .not. wstatus ) wstatus = setwindowconfig( myscreen )
c
c-----------------------------------------------------------------------
c
      if ( myscreen.numcolors .eq. 256 ) then
c
         if ( tcolor ) then
c
            col148(1,ncoff+1) = 0
            col148(2,ncoff+1) = 0
            col148(3,ncoff+1) = 0
            col148(4,ncoff+1) = 0
c
            col148(3,ncoff+2) = 0
            col148(2,ncoff+2) = 254     !  Green for good targets
            col148(1,ncoff+2) = 0
c
            col148(3,ncoff+3) = 0       !  Blue for box outline
            col148(2,ncoff+3) = 0
            col148(1,ncoff+3) = 254
c
            col148(3,ncoff+4) = 254     !  Red for targets rejected by target box
            col148(2,ncoff+4) = 0
            col148(1,ncoff+4) = 0
c
            nstart = 5
c
         else
c
            nstart = 2
c
         endif
c
         do col = ncoff + nstart , ncoff + 128
c
c           col148(1,col) = 2 * ( 128 - col )  !  Blue ---> Yellow
c
            col148(1,col) = 2 * ( col - ncoff - 1 )
            col148(2,col) = 2 * ( col - ncoff - 1 )
            col148(3,col) = 2 * ( col - ncoff - 1 )
c
         enddo
c
         do colidx = ncoff , 127 + ncoff
c
            r              = col148(3,colidx+1)
            g              = col148(2,colidx+1)
            b              = col148(1,colidx+1)
c
            ctable(colidx) = rgbtointeger( r , g , b )
c
         enddo
c
c   Index Color Index Color 
c
c   0  $BLACK     8  $GRAY  
c   1  $BLUE      9  $LIGHTBLUE  
c   2  $GREEN    10  $LIGHTGREEN  
c   3  $CYAN     11  $LIGHTCYAN  
c   4  $RED      12  $LIGHTRED  
c   5  $MAGENTA  13  $LIGHTMAGENTA  
c   6  $BROWN    14  $YELLOW  
c   7  $WHITE    15  $BRIGHTWHITE  
c
         ctable(0)  = $black
         ctable(1)  = $blue
         ctable(2)  = $green
         ctable(3)  = $cyan
         ctable(4)  = $red
         ctable(5)  = $magenta
         ctable(6)  = $brown
         ctable(7)  = $white
         ctable(8)  = $gray
         ctable(9)  = $lightblue
         ctable(10) = $lightgreen
         ctable(11) = $lightcyan
         ctable(12) = $lightred
         ctable(13) = $lightmagenta
         ctable(14) = $yellow
         ctable(15) = $brightwhite
c
         status     = remapallpalettergb( ctable )
c
         status     = setcolorrgb( ctable(127+ncoff) )
c
         status     = settextcolorrgb( ctable(127+ncoff) )
c
         status     = setbkcolorrgb( ctable(ncoff) )
c
      else if ( myscreen.numcolors .eq. 16 ) then
c
         col16(1,1) = 0
         col16(2,1) = 0
         col16(3,1) = 0
         col16(4,1) = 0
c
         do col = 2 , 16
c
c           col16(1,col) = 16 * ( 16 - col )
c
            col16(1,col) = 16 * ( col - 1 )
            col16(2,col) = 16 * ( col - 1 )
            col16(3,col) = 16 * ( col - 1 )
            col16(4,col) = 0
c
         enddo
c
         do colidx = 0 , 15
c
            r = col16(3,colidx+1)
            g = col16(2,colidx+1)
            b = col16(1,colidx+1)
c
            ctable4(colidx) = rgbtointeger( r , g , b )
c
         enddo
c
         status = remapallpalettergb( ctable4 )
c
         status = setcolorrgb( ctable(15) )
c
         status = settextcolorrgb( ctable4(15) )
c
         status = setbkcolorrgb( ctable4(0) )
c
      else
c
         stop ' Cannot set 16 or 256 colors mode'
c
      endif
c
c-----------------------------------------------------------------------
c
c     write ( 6 , * ) ' Number of colors   = ' , myscreen.numcolors
c     write ( 6 , * ) ' Number of xpixels  = ' , myscreen.numxpixels
c     write ( 6 , * ) ' Number of ypixels  = ' , myscreen.numypixels
c     write ( 6 , * ) ' Number of textcols = ' , myscreen.numtextcols
c     write ( 6 , * ) ' Number of textrows = ' , myscreen.numtextrows
c     write ( 6 , * ) ' Fontsize           = ' , myscreen.fontsize
c     write ( 6 , * ) ' Bitsperpixel       = ' , myscreen.bitsperpixel
c
c     read ( 5 , * )
c
c     call clearscreen ( $gclearscreen )
c
c-----------------------------------------------------------------------
c
      return
      end
