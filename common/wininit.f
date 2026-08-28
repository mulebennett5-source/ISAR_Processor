c
c***********************************************************************
c
c   This file contains routines which are specific to MicroSoft Windows
c   programs compiled by the DEC Visual Fortran compiler
c
c   The first routine, wininit, initializes the screen and sets the
c   text for the About Box.
c
c   The second routine, set_cur_win, sets the current active window to
c   a given window description.
c
c   The third routine, setup_video, sets up the color table.
c
c   The fourth routine, INITIALSETTINGS, is a special name which, when
c   linked with a DEC Fortran program, is called during the program
c   initialization process.  The initial window sizes are set here.
c**
c***********************************************************************
c**
      subroutine wininit ( w , h )
c**
c***********************************************************************
c**
      use                     dflib
c
      implicit none
c
      integer                 iresult
c
      logical                 result
c
      record / windowconfig / myscreen
c
      TYPE ( qwinfo )         qwi
c
      integer*2               w , h
c
c-----------------------------------------------------------------------
c
      iresult  = setexitqq( qwin$exitnopersist )
c
c   Set the text for the ABOUT-Box on the help menu
c
      result   = aboutboxqq(
     .         '   RDRTec ISAR\n   Version 3.31\n Jan. 24, 2005'C )
c
      qwi.type = qwin$set
      qwi.w    = w
      qwi.h    = h
c
      result   = setwsizeqq( 0 , qwi )
c
      call clearscreen ( $gclearscreen )
c
c   Set the video display mode.
c
      call setup_video ( myscreen , w , h )
c
c-----------------------------------------------------------------------
c
c   The folowing set of function calls appears to be necessary to make
c   the window maximum - I don't know why.
c
      result     = getwsizeqq( 0 , qwin$sizecurr , qwi )
c
      write ( 6 , * ) qwi.w , qwi.h
c
      result     = getwsizeqq( 0 , qwin$sizemax , qwi )
c
      qwi.type   = qwin$set
c
      result     = setwsizeqq( 0 , qwi )
c
      call clearscreen ( $gclearscreen )
c
c-----------------------------------------------------------------------
c
      return
      end
c**
c***********************************************************************
c**
      include 'scr_stuf.f'
c**
c***********************************************************************
c**
c**
c***********************************************************************
c**
      LOGICAL(4) FUNCTION INITIALSETTINGS()
c**
c***********************************************************************
c**
      use              dflib
!
      implicit none
!
      integer          i
!
      TYPE ( qwinfo )  qwi
!
! Set window frame size.
!
      qwi.x           = 0
      qwi.y           = 0
      qwi.w           = 680
      qwi.h           = 600
!
      qwi.type        = QWIN$SET
!
      i               = SetWSizeQQ( QWIN$FRAMEWINDOW , qwi )
!      
      INITIALSETTINGS = .true.
!
      END FUNCTION INITIALSETTINGS
