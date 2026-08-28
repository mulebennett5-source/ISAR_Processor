c
c***********************************************************************
c
c   This code implements the ship simulator.  It was written by
c   Chris Campo
c
c***********************************************************************
c
c     Subroutine getuparms
c
c     Purpose: Get User-Defined Parameters for Ship Rotation from a
c              Text File
c
c     inyaw     = Initial (time=0) yaw angle(radians)
c     inpitch   = Initial (time=0) pitch angle(radians)
c     inroll    = Initial (time=0) roll angle(radians)
c     avelyaw   = Constant yaw angular velocity
c     avelpitch = Constant pitch angular velocity
c     avelroll  = Constant roll angular velocity
c     aacyaw    = Constant yaw angular acceleration
c     aacpitch  = Constant pitch angular acceleration 
c     aacroll   = Constant roll angular acceleration
c
c     infile    = name of input text file with user input parameters
c     inunit    = unit on which input file will be opened.
c
c***********************************************************************
c
      subroutine getuparms ( inyaw , inpitch , inroll ,
     .                       avelyaw , avelpitch , avelroll ,
     .                       aacyaw , aacpitch , aacroll , pyawa ,
     .                       pyawp , pyawt , ppitcha , ppitchp ,
     .                       ppitcht , prolla , prollp , prollt ,
     .                       infile , inunit , nerror )
c
      implicit none
c
      real      inyaw , inpitch , inroll , avelyaw , avelpitch ,
     .          avelroll , aacyaw , aacpitch , aacroll , pyawa ,
     .          pyawp , pyawt , ppitcha , ppitchp , ppitcht ,
     .          prolla , prollp , prollt
c
      character infile*80

      integer   inunit , nerror
c
      open ( file = infile , form = 'formatted' , status = 'old' ,
     .       unit = inunit , err = 2000 )
c
      read ( inunit , 100 , end = 1000 , err = 2000 ) inyaw 
      read ( inunit , 100 , end = 1000 , err = 2000 ) inpitch
      read ( inunit , 100 , end = 1000 , err = 2000 ) inroll
      read ( inunit , 100 , end = 1000 , err = 2000 ) avelyaw
      read ( inunit , 100 , end = 1000 , err = 2000 ) avelpitch
      read ( inunit , 100 , end = 1000 , err = 2000 ) avelroll
      read ( inunit , 100 , end = 1000 , err = 2000 ) aacyaw
      read ( inunit , 100 , end = 1000 , err = 2000 ) aacpitch 
      read ( inunit , 100 , end = 1000 , err = 2000 ) aacroll
c
c   Optionally, read periodic components of yaw, pitch and roll
c
      read ( inunit , * , end = 3000 , err = 3000 )
     .                                 pyawa   , pyawp   , pyawt  
      read ( inunit , * , end = 3000 , err = 3000 )
     .                                 ppitcha , ppitchp , ppitcht  
      read ( inunit , * , end = 3000 , err = 3000 )
     .                                 prolla  , prollp  , prollt  
c
      close ( inunit )
c
      write ( 6 , * )
      write ( 6 , * )
      write ( 6 , '(3x,a30,f10.4)' ) 'INITIAL YAW ANGLE(radians):   ' ,
     .                                inyaw 
      write ( 6 , '(3x,a30,f10.4)' ) 'INITIAL PITCH ANGLE(radians): ' ,
     .                                inpitch 
      write ( 6 , '(3x,a30,f10.4)' ) 'INITIAL ROLL ANGLE(radians):  ' ,
     .                                inroll 
      write ( 6 , '(3x,a30,f10.4)' ) 'YAW ANGULAR VELOCITY:         ' ,
     .                                avelyaw 
      write ( 6 , '(3x,a30,f10.4)' ) 'PITCH ANGULAR VELOCITY:       ' ,
     .                                avelpitch 
      write ( 6 , '(3x,a30,f10.4)' ) 'ROLL ANGULAR VELOCITY:        ' ,
     .                                avelroll 
      write ( 6 , '(3x,a30,f10.4)' ) 'YAW ANGULAR ACCELERATION:     ' ,
     .                                aacyaw 
      write ( 6 , '(3x,a30,f10.4)' ) 'PITCH ANGULAR ACCELERATION:   ' ,
     .                                aacpitch 
      write ( 6 , '(3x,a30,f10.4)' ) 'ROLL ANGULAR ACCELERATION:    ' ,
     .                                aacroll  
      write ( 6 , * )
      write ( 6 , * ) pyawa   , pyawp   , pyawt
      write ( 6 , * ) ppitcha , ppitchp , ppitcht
      write ( 6 , * ) prolla  , prollp  , prollt
      write ( 6 , * )
c
      write ( 7 , * )
      write ( 7 , * )
      write ( 7 , '(3x,a30,f10.4)' ) 'INITIAL YAW ANGLE(radians):   ' ,
     .                                inyaw 
      write ( 7 , '(3x,a30,f10.4)' ) 'INITIAL PITCH ANGLE(radians): ' ,
     .                                inpitch 
      write ( 7 , '(3x,a30,f10.4)' ) 'INITIAL ROLL ANGLE(radians):  ' ,
     .                                inroll 
      write ( 7 , '(3x,a30,f10.4)' ) 'YAW ANGULAR VELOCITY:         ' ,
     .                                avelyaw 
      write ( 7 , '(3x,a30,f10.4)' ) 'PITCH ANGULAR VELOCITY:       ' ,
     .                                avelpitch 
      write ( 7 , '(3x,a30,f10.4)' ) 'ROLL ANGULAR VELOCITY:        ' ,
     .                                avelroll 
      write ( 7 , '(3x,a30,f10.4)' ) 'YAW ANGULAR ACCELERATION:     ' ,
     .                                aacyaw 
      write ( 7 , '(3x,a30,f10.4)' ) 'PITCH ANGULAR ACCELERATION:   ' ,
     .                                aacpitch 
      write ( 7 , '(3x,a30,f10.4)' ) 'ROLL ANGULAR ACCELERATION:    ' ,
     .                                aacroll  
      write ( 7 , * )
      write ( 7 , * ) pyawa   , pyawp   , pyawt
      write ( 7 , * ) ppitcha , ppitchp , ppitcht
      write ( 7 , * ) prolla  , prollp  , prollt
      write ( 7 , * )
c
  100 format( f10.4 )
  200 format( i4 )
c
      nerror = 0
c
      return
c
 3000 pyawa   = 0.0
      pyawp   = 0.0
      pyawt   = 0.0
      ppitcha = 0.0
      ppitchp = 0.0
      ppitcht = 0.0
      prolla  = 0.0
      prollp  = 0.0
      prollt  = 0.0
c
      nerror  = 0
c
      return
c
 1000 nerror = 1
c
      write ( 6 , * ) ' End of file encountered in GETUPARMS'
      write ( 7 , * ) ' End of file encountered in GETUPARMS'
c
      return
c
 2000 nerror = 2
c
      write ( 6 , * ) ' File I/O error encountered in GETUPARMS'
      write ( 7 , * ) ' File I/O error encountered in GETUPARMS'
c
      return 
      end
c
c      *************************************************************
c      * SUBROUTINE GETSHPCOORD                                    *       
c      *************************************************************
c
      subroutine getshpcoord ( tgt_type , x_scat , y_scat , z_scat ,
     .                         rcs , svnois , svtime , svfreq ,
     .                         svamp , spfreq , spdur, spamp ,
     .                         spstart , maxscats , numscats ,
     .                         infile , inunit , nerror )
c
      implicit none
c
      integer   maxscats , numscats , tgt_type(maxscats)
c
      real      x_scat(maxscats) , y_scat(maxscats) ,
     .          z_scat(maxscats) , svnois(maxscats) ,
     .          svtime(maxscats) , svfreq(maxscats) ,
     .          svamp(maxscats) , rcs(maxscats) ,
     .          spfreq(maxscats) , spdur(maxscats) , 
     .          spamp(maxscats) , spstart(maxscats)
c
      character infile*80
c
      integer   inunit , iloop , nerror 
c
c      OPEN INPUT FILE CONTAINING SCATTERER COORDS
c
      open( file = infile , form = 'formatted' , status = 'old' ,
     .      unit = inunit , err = 2000 )
c
c       READ IN NUMBER OF SCATTERERS
c
      read( inunit , 100 , end = 1000 , err = 2000 ) numscats
  100 format( i4 )
c
c      -------------------------------------------------------------
c      CHECK TO MAKE SURE WE ARE DEALING WITH THE CORRECT NUMBER OF
c      SCATTERERS
c      -------------------------------------------------------------
c
      if ( numscats .gt. maxscats ) then 
c
         write ( 6 , * ) '---------------------------------------------'
         write ( 6 , * ) ' ERROR IN SUBROUTINE GET_SHP_COORD:'
         write ( 6 , * ) ' '
         write ( 6 , * ) ' No. Of Scatterers (Hardwired Parameter) = ' ,
     .                     maxscats
         write ( 6 , * ) ' No. Of Scatterers (Passed Variable)     = ' ,
     .                     numscats
         write ( 6 , * ) ' '
         write ( 6 , * ) ' YOU SHOULD INCREASE NTESTM IN THE FILE'
         write ( 6 , * ) ' SARPRM.H TO SOLVE THIS PROBLEM !!!!!!'
         write ( 6 , * )
         write ( 6 , * ) ' Hardwired should be greater than or equal to'
         write ( 6 , * ) ' the passed variable. Bailing..... '
         write ( 6 , * ) '---------------------------------------------'
c
         stop
c
      endif
c
c       LOOP OVER READING IN ALL SCATTERER COORDS
c
      do iloop = 1 , numscats
c
c     Target type 1 is ship particle         
c                 2 is ship particle, independent vibration
c                 3 is ship particle, periodic rcs
c
c
         read( inunit , * , end = 1000 , err = 2000 ) 
     .      tgt_type(iloop) ,
     .      x_scat(iloop) , y_scat(iloop) , z_scat(iloop) , rcs(iloop)
c
         if ( tgt_type(iloop) .eq. 2 ) then
c
            read( inunit , * , end = 1000 , err = 2000 ) 
     .      svnois(iloop) , svtime(iloop) , svfreq(iloop) , svamp(iloop)
c
         else if ( tgt_type(iloop) .eq. 3 ) then
c
            read( inunit , * , end = 1000 , err = 2000 ) 
     .      spfreq(iloop) , spdur(iloop) , spamp(iloop) , spstart(iloop)
c
         endif
c
      end do
c
c      CLOSE INPUT FILE CONTAINING SCATTERER COORDS
c
      close( inunit )
c
      nerror = 0
      return
c
 1000 nerror = 1
c
      write ( 6 , * ) ' End of file encountered in GETSHPCOORD'
      write ( 7 , * ) ' End of file encountered in GETSHPCOORD'
c
      return
c
 2000 nerror = 2
c
      write ( 6 , * ) ' File I/O error encountered in GETSHPCOORD'
      write ( 7 , * ) ' File I/O error encountered in GETSHPCOORD'
c
      return
      end
c
c***********************************************************************
c
c     Subroutine rot_ship
c   
c     Purpose: rotate a ship as time goes on according to specification
c              of initial ship orientation (yaw,pitch,roll),
c              constant angular velocities, and constant angular 
c              accelerations (all input). Based on Greg Medlin's
c              MATLAB rotation routine .
c   
c     INPUT VARS IN ARGUMENT LIST
c     numscats                   = Number of scatterers on the body.
c     inyaw,inpitch,inroll       = ORIGINAL angles of orientation of
c                                  body (i.e.: at time = 0).
c     x_scat,y_scat,z_scat       = 3d coords of scatterers on body.
c     avelyaw,avelpitch,avelroll = angular velocities.
c     aacyaw,aacpitch,aacroll    = angular accelerations.
c     time                       = time elapsed from initial time.
c   
c     OUTPUT VARS IN ARGUMENT LIST
c     newxsc,newysc,newzsc       = new 3d coords after rotation.
c    
c***********************************************************************
c
      subroutine rot_ship ( numscats , inyaw , inpitch , inroll ,
     .                      avelyaw , avelpitch , avelroll ,
     .                      aacyaw , aacpitch , aacroll , pyawa ,
     .                      pyawp , pyawt , ppitcha , ppitchp ,
     .                      ppitcht , prolla , prollp , prollt ,
     .                      x_scat , y_scat , z_scat ,
     .                      newxsc , newysc , newzsc , time )
c
      implicit none
c
c     PASSED VARIABLES
c
      integer numscats
c
      real    inyaw , inpitch , inroll , avelyaw , avelpitch ,
     .        avelroll , aacyaw , aacpitch , aacroll , pyawa , pyawp ,
     .        pyawt , ppitcha , ppitchp , ppitcht , prolla , prollp ,
     .        prollt , time
c
      real    x_scat(numscats) , y_scat(numscats) , z_scat(numscats)
c
      real    newxsc(numscats) , newysc(numscats) , newzsc(numscats)
c
c     LOCAL VARIABLES
c
      real    yaw_rot(3,3) , pitch_rot(3,3) , roll_rot(3,3) ,
     .        tmp_rot(3,3) , rotmat(3,3) , colvec(3) , outcol(3)
c
      real    outyaw , outpitch , outroll , twopi
c
      integer iloop
c      
c-----------------------------------------------------------------------
c
c   CALCULATE CURRENT YAW, PITCH, ROLL ANGLES ACCORDING TO:
c
c             Ang = AngVel * Time + 1/2 AngAcc * Time^2 + OrigAng
c                 + amp * cos( ( twopi / period ) * ( t - t0 ) )
c
c-----------------------------------------------------------------------
c
c   Previous value + Linear + Quadratic
c
      outyaw   = inyaw   + avelyaw   * time + 0.5 * aacyaw   * time ** 2
      outpitch = inpitch + avelpitch * time + 0.5 * aacpitch * time ** 2
      outroll  = inroll  + avelroll  * time + 0.5 * aacroll  * time ** 2
c
c   Add time periodic components
c
      twopi    = 2.0 * atan2( 0.0 , - 1.0 )
c
      outyaw   = outyaw   + pyawa   *
     .           cos( ( twopi / pyawp )   * ( time - pyawt   ) )
c
      outpitch = outpitch + ppitcha *
     .           cos( ( twopi / ppitchp ) * ( time - ppitcht ) )
c
      outroll  = outroll  + prolla  *
     .           cos( ( twopi / prollp )  * ( time - prollt  ) )
c
      !---------------------------------
      ! CONSTRUCT YAW ROTATION MATRIX
      !---------------------------------
c
      ! FIRST ROW
      yaw_rot(1,1)   = cos( outyaw )
      yaw_rot(1,2)   = - sin( outyaw )
      yaw_rot(1,3)   = 0.0
c
      ! SECOND ROW
      yaw_rot(2,1)   = sin( outyaw )
      yaw_rot(2,2)   = cos( outyaw )
      yaw_rot(2,3)   = 0.0
c      
      ! THIRD ROW
      yaw_rot(3,1)   = 0.0
      yaw_rot(3,2)   = 0.0
      yaw_rot(3,3)   = 1.0
c
      !---------------------------------
      ! CONSTRUCT PITCH ROTATION MATRIX
      !---------------------------------
c
      ! FIRST ROW
      pitch_rot(1,1) = cos( outpitch )
      pitch_rot(1,2) = 0.0
      pitch_rot(1,3) = - sin( outpitch )
c
      ! SECOND ROW
      pitch_rot(2,1) = 0.0
      pitch_rot(2,2) = 1.0
      pitch_rot(2,3) = 0.0
c
      ! THIRD ROW
      pitch_rot(3,1) = sin( outpitch )
      pitch_rot(3,2) = 0.0
      pitch_rot(3,3) = cos( outpitch )
c
      !---------------------------------
      ! CONSTRUCT ROLL ROTATION MATRIX
      !---------------------------------
c
      ! FIRST ROW
      roll_rot(1,1)  = 1.0
      roll_rot(1,2)  = 0.0
      roll_rot(1,3)  = 0.0
c
      ! SECOND ROW
      roll_rot(2,1)  = 0.0
      roll_rot(2,2)  = cos( outroll )
      roll_rot(2,3)  = sin( outroll )
c
      ! THIRD ROW
      roll_rot(3,1)  = 0.0
      roll_rot(3,2)  = - sin( outroll )
      roll_rot(3,3)  = cos( outroll )
c
      !-----------------------
      ! ROTATE THE SHIP
      !-----------------------
c
      ! CREATE COMPOSITE ROTATION MATRIX
      call mat_mult ( pitch_rot , yaw_rot , tmp_rot )
      call mat_mult ( tmp_rot , roll_rot , rotmat )
      
      ! ROTATE THE COORDS OF ALL SCATTERERS ON THE SHIP
      do iloop = 1 , numscats
c
      ! FILL COLUMN VECTOR WITH X,Y,Z COORDS OF CURRENT SCATTERER
         colvec(1)     = x_scat(iloop)
         colvec(2)     = y_scat(iloop)
         colvec(3)     = z_scat(iloop)
c
      ! MULT COLVEC OF OLD COORDS BY COMPOSITE ROTATION MATRIX
         call matvec_mul ( rotmat , colvec , outcol )

      ! GET NEW ROTATED COORDS OF SCATTERERS
         newxsc(iloop) = outcol(1)
         newysc(iloop) = outcol(2)
         newzsc(iloop) = outcol(3)
c
      enddo
c
      return
      end
c
c***********************************************************************
c
c         Version for compatibility with the 3D-ISAR algorithm
c
c     Subroutine rot_ship2: Version using extrinsic definition for
c                           the Euler rotation angles. Also, changes
c                           the signs of the periodic terms.
c   
c     Purpose: rotate a ship as time goes on according to specification
c              of initial ship orientation (yaw,pitch,roll),
c              constant angular velocities, and constant angular 
c              accelerations (all input). Based on Greg Medlin's
c              MATLAB rotation routine .
c   
c     INPUT VARS IN ARGUMENT LIST
c     numscats                   = Number of scatterers on the body.
c     inyaw,inpitch,inroll       = ORIGINAL angles of orientation of
c                                  body (i.e.: at time = 0).
c     x_scat,y_scat,z_scat       = 3d coords of scatterers on body.
c     avelyaw,avelpitch,avelroll = angular velocities.
c     aacyaw,aacpitch,aacroll    = angular accelerations.
c     time                       = time elapsed from initial time.
c   
c     OUTPUT VARS IN ARGUMENT LIST
c     newxsc,newysc,newzsc       = new 3d coords after rotation.
c    
c***********************************************************************
c
      subroutine rot_ship2 ( numscats , inyaw , inpitch , inroll ,
     .                       avelyaw , avelpitch , avelroll ,
     .                       aacyaw , aacpitch , aacroll , pyawa ,
     .                       pyawp , pyawt , ppitcha , ppitchp ,
     .                       ppitcht , prolla , prollp , prollt ,
     .                       x_scat , y_scat , z_scat ,
     .                       newxsc , newysc , newzsc , time )
c
      implicit none
c
c     PASSED VARIABLES
c
      integer numscats
c
      real    inyaw , inpitch , inroll , avelyaw , avelpitch ,
     .        avelroll , aacyaw , aacpitch , aacroll , pyawa , pyawp ,
     .        pyawt , ppitcha , ppitchp , ppitcht , prolla , prollp ,
     .        prollt , time , tilt
c
      real    x_scat(numscats) , y_scat(numscats) , z_scat(numscats)
c
      real    newxsc(numscats) , newysc(numscats) , newzsc(numscats)
c
c     LOCAL VARIABLES
c
      real    yaw_rot(3,3) , pitch_rot(3,3) , roll_rot(3,3) ,
     .        tmp_rot(3,3) , rotmat(3,3) , colvec(3) , outcol(3)
c
      real    outyaw , outpitch , outroll , twopi
c
      integer iloop
c      
c-----------------------------------------------------------------------
c
c   CALCULATE CURRENT YAW, PITCH, ROLL ANGLES ACCORDING TO:
c
c             Ang = AngVel * Time + 1/2 AngAcc * Time^2 + OrigAng
c                 + amp * cos( ( twopi / period ) * ( t - t0 ) )
c
c-----------------------------------------------------------------------
c
c   Previous value + Linear + Quadratic
c
      outyaw   = inyaw   + avelyaw   * time + 0.5 * aacyaw   * time ** 2
      outpitch =         + avelpitch * time + 0.5 * aacpitch * time ** 2
      outroll  =         + avelroll  * time + 0.5 * aacroll  * time ** 2
c
c   Add time periodic components
c
      twopi    = 2.0 * atan2( 0.0 , - 1.0 )
c
      outyaw   = outyaw   - pyawa   *
     .           cos( ( twopi / pyawp )   * ( time - pyawt   ) )
c
      outpitch = outpitch - ppitcha *
     .           cos( ( twopi / ppitchp ) * ( time - ppitcht ) )
c
      outroll  = outroll  - prolla  *
     .           cos( ( twopi / prollp )  * ( time - prollt  ) )
c
      !---------------------------------
      ! CONSTRUCT YAW ROTATION MATRIX
      !---------------------------------
c
      ! FIRST ROW
      yaw_rot(1,1)   = cos( outyaw )
      yaw_rot(1,2)   = - sin( outyaw )
      yaw_rot(1,3)   = 0.0
c
      ! SECOND ROW
      yaw_rot(2,1)   = sin( outyaw )
      yaw_rot(2,2)   = cos( outyaw )
      yaw_rot(2,3)   = 0.0
c      
      ! THIRD ROW
      yaw_rot(3,1)   = 0.0
      yaw_rot(3,2)   = 0.0
      yaw_rot(3,3)   = 1.0
c
      !---------------------------------
      ! CONSTRUCT PITCH ROTATION MATRIX
      !---------------------------------
c
      tilt           = inpitch + cos( outyaw ) * outpitch +
     .                           sin( outyaw ) * outroll
c
      ! FIRST ROW
      pitch_rot(1,1) = cos( tilt )
      pitch_rot(1,2) = 0.0
      pitch_rot(1,3) = - sin( tilt )
c
      ! SECOND ROW
      pitch_rot(2,1) = 0.0
      pitch_rot(2,2) = 1.0
      pitch_rot(2,3) = 0.0
c
      ! THIRD ROW
      pitch_rot(3,1) = sin( tilt )
      pitch_rot(3,2) = 0.0
      pitch_rot(3,3) = cos( tilt )
c
      !---------------------------------
      ! CONSTRUCT ROLL ROTATION MATRIX
      !---------------------------------
c
      ! FIRST ROW
      roll_rot(1,1)  = 1.0
      roll_rot(1,2)  = 0.0
      roll_rot(1,3)  = 0.0
c
      ! SECOND ROW
      roll_rot(2,1)  = 0.0
      roll_rot(2,2)  = 1.0
      roll_rot(2,3)  = 0.0
c
      ! THIRD ROW
      roll_rot(3,1)  = 0.0
      roll_rot(3,2)  = 0.0
      roll_rot(3,3)  = 1.0
c
      !-----------------------
      ! ROTATE THE SHIP
      !-----------------------
c
      ! CREATE COMPOSITE ROTATION MATRIX
      call mat_mult ( pitch_rot , yaw_rot , tmp_rot )
      call mat_mult ( tmp_rot , roll_rot , rotmat )
      
      ! ROTATE THE COORDS OF ALL SCATTERERS ON THE SHIP
      do iloop = 1 , numscats
c
      ! FILL COLUMN VECTOR WITH X,Y,Z COORDS OF CURRENT SCATTERER
         colvec(1)     = x_scat(iloop)
         colvec(2)     = y_scat(iloop)
         colvec(3)     = z_scat(iloop)
c
      ! MULT COLVEC OF OLD COORDS BY COMPOSITE ROTATION MATRIX
         call matvec_mul ( rotmat , colvec , outcol )

      ! GET NEW ROTATED COORDS OF SCATTERERS
         newxsc(iloop) = outcol(1)
         newysc(iloop) = outcol(2)
         newzsc(iloop) = outcol(3)
c
      enddo
c
      return
      end
c
c***********************************************************************
c
c     SUBROUTINE MATVEC_MUL
c     Routine to perform matrix multiplication on a 
c     3 x 3 matrix with a column vector.
c
c***********************************************************************
c
      subroutine matvec_mul ( mat1 , mat2 , matprod )
c
      implicit none
c
      real mat1(3,3) , mat2(3) , matprod(3)
c
      matprod(1) = mat1(1,1)*mat2(1) + mat1(1,2)*mat2(2) +
     .             mat1(1,3)*mat2(3)
      matprod(2) = mat1(2,1)*mat2(1) + mat1(2,2)*mat2(2) +
     .             mat1(2,3)*mat2(3)
      matprod(3) = mat1(3,1)*mat2(1) + mat1(3,2)*mat2(2) +
     .             mat1(3,3)*mat2(3)
c
      return 
      end
c
c***********************************************************************
c
c     SUBROUTINE MAT_MULT
c     Routine to perform matrix multiplication on two 
c     3 x 3 matrices.
c
c***********************************************************************
c
      subroutine mat_mult ( mat1 , mat2 , matprod )
c
      implicit none
c
      real mat1(3,3) , mat2(3,3) , matprod(3,3)
c 
      matprod(1,1) = mat1(1,1) * mat2(1,1) + mat1(1,2) * mat2(2,1) +
     .               mat1(1,3) * mat2(3,1)
      matprod(1,2) = mat1(1,1) * mat2(1,2) + mat1(1,2) * mat2(2,2) +
     .               mat1(1,3) * mat2(3,2)
      matprod(1,3) = mat1(1,1) * mat2(1,3) + mat1(1,2) * mat2(2,3) +
     .               mat1(1,3) * mat2(3,3)
c
      matprod(2,1) = mat1(2,1) * mat2(1,1) + mat1(2,2) * mat2(2,1) +
     .               mat1(2,3) * mat2(3,1)
      matprod(2,2) = mat1(2,1) * mat2(1,2) + mat1(2,2) * mat2(2,2) +
     .               mat1(2,3) * mat2(3,2)
      matprod(2,3) = mat1(2,1) * mat2(1,3) + mat1(2,2) * mat2(2,3) +
     .               mat1(2,3) * mat2(3,3)
c
      matprod(3,1) = mat1(3,1) * mat2(1,1) + mat1(3,2) * mat2(2,1) +
     .               mat1(3,3) * mat2(3,1)
      matprod(3,2) = mat1(3,1) * mat2(1,2) + mat1(3,2) * mat2(2,2) +
     .               mat1(3,3) * mat2(3,2)
      matprod(3,3) = mat1(3,1) * mat2(1,3) + mat1(3,2) * mat2(2,3) +
     .               mat1(3,3) * mat2(3,3)
c
      return
      end
