C**
C***********************************************************************
C**
      subroutine set_windows ( nfr , nfa , myscreen , secondary , 
     .                         transpose , ascan , root_win ,
     .                         movie_win , movie2_win , ascan_win  ,
     .                         ascant_win , x_ipr_win , y_ipr_win )
C**
C***********************************************************************
C**
      use                     dflib
c
      implicit none
c
      include                'rmovie.h'
c
      record / windowconfig / myscreen
c
      integer                 nfr , nfa , nx , ny , win_maxx ,
     .                        win_maxy , y1phys , y2phys
c
      logical                 secondary , transpose , ascan
c
      record / win_def      / root_win , movie_win , movie2_win ,
     .                        ascan_win , ascant_win , x_ipr_win ,
     .                        y_ipr_win
c
c   Root window is entire drawing area
c
      root_win.x1_phys = 0
      root_win.x2_phys = myscreen.numxpixels - 1
      root_win.y1_phys = 0
      root_win.y2_phys = myscreen.numypixels - 1
      root_win.x1_win  = 0
      root_win.x2_win  = myscreen.numxpixels - 1
      root_win.y1_win  = 0
      root_win.y2_win  = myscreen.numypixels - 1
      root_win.finvert = .false.
c
      if ( transpose ) then
c
         if ( secondary ) stop 'Transpose and Secondary'
c
         nx = nfr
         ny = nfa
c
      else
c
         nx = nfa
         ny = nfr
c
      endif
c
c   Use only standard (non-inverted) windows
c
      movie_win.finvert  = .false.
      movie2_win.finvert = .false.
      ascan_win.finvert  = .false.
      ascant_win.finvert = .false.
c
      win_maxx           = min( nx , root_win.x2_phys )
c
      win_maxy           = min( movie_maxy - ascan_minheight , ny )
c
      ascant_win.y2_phys = root_win.y2_phys - 3
      ascant_win.y1_phys = ascant_win.y2_phys - ascan_minheight
c
      ascant_win.x1_win  = ( nx - win_maxx ) / 2 + 1
      ascant_win.x2_win  = ascant_win.x1_win + win_maxx
c
      movie_win.y2_phys  = min( ascant_win.y1_phys - 4 ,
     .                     ( root_win.y2_phys + win_maxy ) / 2 )
      movie_win.y1_phys  = movie_win.y2_phys - win_maxy + 1
c
      movie_win.x1_win   = ( nx - win_maxx ) / 2 + 1
      movie_win.x2_win   = movie_win.x1_win + win_maxx - 1
      movie_win.y1_win   = ( ny - win_maxy ) / 2 + 1
      movie_win.y2_win   = movie_win.y1_win + win_maxy - 1
c
      win_maxy           = min( movie_maxy , ny )
c
      y1phys             = ( root_win.y2_phys - win_maxy ) / 2
      y2phys             = ( root_win.y2_phys + win_maxy ) / 2 - 1
c
      if ( secondary ) then
c
c   Position secondary window against left edge,
c   Primary and ascan against right edge.
c
         win_maxx           = min( movie2_maxx , nx )
c
         ascan_win.x2_phys  = root_win.x2_phys - 3
         ascan_win.x1_phys  = ascan_win.x2_phys - ascan_minheight + 1
c
         ascan_win.y1_win   = ( ny - win_maxy ) / 2 + 1
         ascan_win.y2_win   = ascan_win.y1_win + win_maxy - 1
c
         movie_win.x2_phys  = root_win.x2_phys - ascan_minheight - 6
         movie_win.x1_phys  = movie_win.x2_phys - win_maxx + 1
         movie_win.y1_phys  = y1phys
         movie_win.y2_phys  = y2phys
c
         movie2_win.x1_phys = 3
         movie2_win.x2_phys = movie2_win.x1_phys + win_maxx - 1
         movie2_win.y1_phys = y1phys
         movie2_win.y2_phys = y2phys
c
         movie_win.x1_win   = ( nx - win_maxx ) / 2 + 1
         movie_win.x2_win   = movie_win.x1_win + win_maxx - 1
         movie_win.y1_win   = ( ny - win_maxy ) / 2 + 1
         movie_win.y2_win   = movie_win.y1_win + win_maxy - 1
c
         movie2_win.x1_win  = ( nx - win_maxx ) / 2 + 1
         movie2_win.x2_win  = movie2_win.x1_win + win_maxx - 1
         movie2_win.y1_win  = ( ny - win_maxy ) / 2 + 1
         movie2_win.y2_win  = movie2_win.y1_win + win_maxy - 1
c
      else
c
         win_maxx           = min( root_win.x2_phys - ascan_minheight ,
     .                             nx )
c
c                         Ascan window
c
         ascan_win.x2_phys  = root_win.x2_phys - 3
         ascan_win.x1_phys  = ascan_win.x2_phys - ascan_minheight + 1
c
         ascan_win.y1_win   = ( ny - win_maxy ) / 2 + 1
         ascan_win.y2_win   = ascan_win.y1_win + win_maxy - 1
c
c                         Movie window
c
         movie_win.x2_phys  = min( ascan_win.x1_phys - 4 ,
     .                        ( root_win.x2_phys + win_maxx ) / 2 )
         movie_win.x1_phys  = movie_win.x2_phys - nx + 1
         movie_win.y1_phys  = y1phys
         movie_win.y2_phys  = y2phys
c
         movie_win.x1_win   = ( nx - win_maxx ) / 2 + 1
         movie_win.x2_win   = movie_win.x1_win + win_maxx - 1
         movie_win.y1_win   = ( ny - win_maxy ) / 2 + 1
         movie_win.y2_win   = movie_win.y1_win + win_maxy - 1
c
      endif
c
c   Align the A-Scan/IPR windows with the primary movie window
c
      ascant_win.x1_phys = movie_win.x1_phys
      ascant_win.x2_phys = movie_win.x2_phys
c
      ascan_win.y1_phys  = movie_win.y1_phys
      ascan_win.y2_phys  = movie_win.y2_phys
c
      x_ipr_win          = ascant_win
      y_ipr_win          = ascan_win
c
      return
      end
C**
C***********************************************************************
C**
      subroutine frames ( root_win , movie_win , movie2_win ,
     .                    ascan_win , ascant_win , transpose ,
     .                    secondary , ascan , x_ipr_win , y_ipr_win ,
     .                    log_data , logfile )
C**
C***********************************************************************
C**
      use                dflib
c
      implicit none
c
      include           'rmovie.h'
c
      include           'parselog.h'
c
      record / win_def / root_win , movie_win , movie2_win , ascan_win ,
     .                   ascant_win , x_ipr_win , y_ipr_win
c
      logical            transpose , secondary , ascan , logfile
c
      record / log_dat / log_data
c
      integer            nxpix , nypix , numfonts
c
      integer*2          status
c
      record / xycoord / s
c
      character          rngtxt*40, ftxt*40
c
      real               fmin , fmax
c
      if ( secondary .and. transpose ) stop 'Transpose and side-by-side'
c
      status = rectangle( $GBORDER , root_win.x1_phys ,
     .                    root_win.y1_phys , root_win.x2_phys ,
     .                    root_win.y2_phys )
c
      status = rectangle( $GBORDER , movie_win.x1_phys - 1 ,
     .                               movie_win.y1_phys - 1 ,
     .                               movie_win.x2_phys + 1 ,
     .                               movie_win.y2_phys + 1 )
c
      if ( secondary ) status = rectangle( $GBORDER ,
     .                 movie2_win.x1_phys - 1 , movie2_win.y1_phys - 1 ,
     .                 movie2_win.x2_phys + 1 , movie2_win.y2_phys + 1 )
c
      if ( ascan ) then
c
c        if ( transpose ) then
c
                       status = rectangle( $GBORDER ,
     .                 ascant_win.x1_phys - 1 , ascant_win.y1_phys - 1 ,
     .                 ascant_win.x2_phys + 1 , ascant_win.y2_phys + 1 )
c
c        else
c
                       status = rectangle( $GBORDER ,
     .                 ascan_win.x1_phys - 1 , ascan_win.y1_phys - 1 ,
     .                 ascan_win.x2_phys + 1 , ascan_win.y2_phys + 1 )
c
c        endif
c
      endif
c
      if ( plotipr )   status = rectangle( $GBORDER ,
     .                 x_ipr_win.x1_phys - 1 , x_ipr_win.y1_phys - 1 ,
     .                 x_ipr_win.x2_phys + 1 , x_ipr_win.y2_phys + 1 )
c
      if ( plotipr )   status = rectangle( $GBORDER ,
     .                 y_ipr_win.x1_phys - 1 , y_ipr_win.y1_phys - 1 ,
     .                 y_ipr_win.x2_phys + 1 , y_ipr_win.y2_phys + 1 )
c
      if ( logfile ) then
c
         numfonts = initializefonts()
c
          if ( setfont( 't''Arial''h12w8i' ) .LT. 0 )
     .         stop 'Error: can not execute routine setfont'
c
         write ( rngtxt , '(F8.1,'' /'',F7.1,'' /'',F5.1)' )
     .        log_data.rmin , log_data.rmax , log_data.rinc
c
         nxpix = movie_win.x2_phys - movie_win.x1_phys
         nypix = movie_win.y2_phys - movie_win.y1_phys
c
         if ( transpose ) then
c
            fmin = - 0.5 * log_data.dfine * float( nypix )
            fmax = fmin + log_data.dfine * float( nypix - 1 )
c
            write ( ftxt , '(F7.1,'' /'',F6.1,'' /'',F4.1)' )
     .            fmin , fmax , log_data.dfine
c
            call setgtextrotation ( 0 )
            call moveto ( int2( movie_win.x1_phys + nxpix / 2 ) ,
     .                    int2( movie_win.y1_phys - 26 ) , s )
            call outgtext ( 'R(m)' )
            call moveto (int2( movie_win.x1_phys + nxpix / 4 ) ,
     .                   int2( movie_win.y1_phys - 13) , s )
            call outgtext ( rngtxt )
c
            call setgtextrotation ( 900 )
            call moveto ( int2( movie_win.x1_phys - 26 ) ,
     .                    int2( movie_win.y1_phys + nypix / 2 ) , s )
            call outgtext ( 'F(Hz)' )
            call moveto ( int2( movie_win.x1_phys - 13 ) ,
     .               int2 (movie_win.y1_phys + ( nypix / 4 ) * 3 ) , s )
            call outgtext ( ftxt )
c
         else
c
            fmin = - 0.5 * log_data.dfine * float( nxpix )
            fmax = fmin + log_data.dfine * float( nxpix - 1 )
c
            write ( ftxt , '(F7.1,'' /'',F6.1,'' /'',F4.1)' )
     .              fmin , fmax , log_data.dfine
c
            call setgtextrotation ( 0 )
            call moveto ( int2( movie_win.x1_phys + nxpix / 2 ) ,
     .                    int2( movie_win.y1_phys - 32 ) , s )
            call outgtext ( 'F(Hz)' )
            call moveto ( int2( movie_win.x1_phys + nxpix / 4 ) ,
     .                    int2( movie_win.y1_phys - 18 ) , s )
            call outgtext ( ftxt )
c
            call setgtextrotation ( 900 )
            call moveto ( int2( movie_win.x1_phys - 32 ) ,
     .                    int2( movie_win.y1_phys + nypix / 2 ) , s )
            call outgtext ( 'R(m)' )
            call moveto ( int2( movie_win.x1_phys - 18 ) ,
     .               int2( movie_win.y1_phys + ( nypix / 4 ) * 3 ) , s )
            call outgtext ( rngtxt )
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
      subroutine text_info ( framenum , ascan , ascan_win , ascant_win ,
     .                       produce , log_data , transpose , logfile ,
     .                       infile )
C**
C***********************************************************************
C**
      use                dflib
c
      implicit none
c
      include           'rmovie.h'
c
      include           'parselog.h'
c
      integer            framenum , nfirst , nblank , i
c
      character          strasmin*8 , strasmax*8
c
      character          qualstring*70
c
      character          frame*5
c
      character          infile*80
c
      logical            ascan , transpose , logfile , produce
c
      record / win_def / ascan_win , ascant_win
c
      record / log_dat / log_data
c
      real               sasmin , sasmax
c
      record / rccoord / p
c
      if ( produce ) then
c
         call settextposition ( 2 , 3 , p )
         call outtext ( 'RMOVIE' )
c
      endif
c
      write ( frame , '(I5)' ) framenum
      call settextposition ( 2 , 11 , p )
      call outtext ( 'Frame: ' // frame )
c
      call settextposition ( 2 , 25 , p )
c
      if ( transpose ) then
c
         call outtext ( '  Transpose  ' )
c
      else
c
         call outtext ( '  Normal  ' )
c
      endif
c
c-----------------------------------------------------------------------
c
c   Write ascan limits on the screen
c
      if ( ascan ) then
c
         if ( transpose ) then
c
            sasmin = ascant_win.y2_win
            sasmax = ascant_win.y1_win
c
         else
c
            sasmin = - ascan_win.x2_win
            sasmax = - ascan_win.x1_win
c
         endif
c
         write ( strasmin , '(F8.1)' ) - 10.0 * sasmin
         write ( strasmax , '(F8.1)' ) 10.0 * sasmax
c
         call settextposition ( 2 , 41 , p )
         call outtext ( 'Ascan (dB) min: ' // strasmin
     .                         // ' max: ' // strasmax )
c
      endif
c
      if ( logfile ) then
c
c-----------------------------------------------------------------------
c
c   Write contrast and time information on the screen
c
         call settextposition ( 3 , 2 , p )
c
         nblank = index( infile , ' ' ) - 1
         nfirst = max( 1 , nblank - 15 )
c
         call outtext ( log_data.strcontrast(framenum)(1:60)
     .                  // ' ' // infile(nfirst:nblank) )
c
c-----------------------------------------------------------------------
c
c   Write quality data on the screen
c
         if ( do_sct ) then
c
            write ( qualstring , '(7x,4f14.3)' )
     .            ( quality(i,framenum) , i = 10 , 13 )
c
            call settextposition ( 6 , 2 , p )
c
            call outtext ( qualstring )
c
         endif
c
c-----------------------------------------------------------------------
c
      endif
c
      return
      end
