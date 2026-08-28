C**
C***********************************************************************
C**
      REAL FUNCTION UNI ( JD )
C**
C***********************************************************************
C**
C***BEGIN PROLOGUE  UNI
C***DATE WRITTEN   810915
C***REVISION DATE  830805
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    BLUE, JAMES, SCIENTIFIC COMPUTING DIVISION, NBS
C             KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, COMPUTER SCIENCE DEPT., WASH STATE UNIV
C
C***PURPOSE  THIS ROUTINE GENERATES QUASI UNIFORM RANDOM NUMBERS ON [0,1
C             AND CAN BE USED ON ANY COMPUTER WITH WHICH ALLOWS INTEGERS
C             AT LEAST AS LARGE AS 32767.
C***DESCRIPTION
C
C       THIS ROUTINE GENERATES QUASI UNIFORM RANDOM NUMBERS ON THE INTER
C       [0,1).  IT CAN BE USED WITH ANY COMPUTER WHICH ALLOWS
C       INTEGERS AT LEAST AS LARGE AS 32767.
C
C
C   USE
C       FIRST TIME....
C                   Z = UNI(JD)
C                     HERE JD IS ANY  N O N - Z E R O  INTEGER.
C                     THIS CAUSES INITIALIZATION OF THE PROGRAM
C                     AND THE FIRST RANDOM NUMBER TO BE RETURNED AS Z.
C       SUBSEQUENT TIMES...
C                   Z = UNI(0)
C                     CAUSES THE NEXT RANDOM NUMBER TO BE RETURNED AS Z.
C
C
C..................................................................
C   NOTE: USERS WHO WISH TO TRANSPORT THIS PROGRAM FROM ONE COMPUTER
C         TO ANOTHER SHOULD READ THE FOLLOWING INFORMATION.....
C
C   MACHINE DEPENDENCIES...
C      MDIG = A LOWER BOUND ON THE NUMBER OF BINARY DIGITS AVAILABLE
C              FOR REPRESENTING INTEGERS, INCLUDING THE SIGN BIT.
C              THIS VALUE MUST BE AT LEAST 16, BUT MAY BE INCREASED
C              IN LINE WITH REMARK A BELOW.
C
C   REMARKS...
C     A. THIS PROGRAM CAN BE USED IN TWO WAYS:
C        (1) TO OBTAIN REPEATABLE RESULTS ON DIFFERENT COMPUTERS,
C            SET 'MDIG' TO THE SMALLEST OF ITS VALUES ON EACH, OR,
C        (2) TO ALLOW THE LONGEST SEQUENCE OF RANDOM NUMBERS TO BE
C            GENERATED WITHOUT CYCLING (REPEATING) SET 'MDIG' TO THE
C            LARGEST POSSIBLE VALUE.
C     B. THE SEQUENCE OF NUMBERS GENERATED DEPENDS ON THE INITIAL
C          INPUT 'JD' AS WELL AS THE VALUE OF 'MDIG'.
C          IF MDIG=16 ONE SHOULD FIND THAT
C            THE FIRST EVALUATION
C              Z=UNI(305) GIVES Z=.027832881...
C            THE SECOND EVALUATION
C              Z=UNI(0) GIVES   Z=.56102176...
C            THE THIRD EVALUATION
C              Z=UNI(0) GIVES   Z=.41456343...
C            THE THOUSANDTH EVALUATION
C              Z=UNI(0) GIVES   Z=.19797357...
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***END PROLOGUE  UNI
C
      IMPLICIT NONE
C
      INTEGER I , J , K , M(17) , M1 , M2 , MDIG , J0 , J1 ,
     .        K0 , K1 , JSEED , JD
C
      SAVE
C
C***FIRST EXECUTABLE STATEMENT  UNI
C
      IF ( JD .NE. 0 ) THEN
C
         MDIG = 32
C
C   BE SURE THAT MDIG AT LEAST 16...
C
         IF ( MDIG .LT. 16 ) WRITE ( * , * ) 'UNI--MDIG LESS THAN 16'
C
         M1    = 2 ** ( MDIG - 2 ) + ( 2 ** ( MDIG - 2 ) - 1 )
         M2    = 2 ** ( MDIG / 2 )
         JSEED = MIN0( IABS( JD ) , M1 )
         IF ( MOD( JSEED , 2 ) .EQ. 0 ) JSEED = JSEED - 1
         K0    = MOD( 9069 , M2 )
         K1    = 9069 / M2
         J0    = MOD( JSEED , M2 )
         J1    = JSEED / M2
C
         DO I = 1 , 17
C
            JSEED = J0 * K0
            J1    = MOD( JSEED / M2 + J0 * K1 + J1 * K0 , M2 / 2 )
            J0    = MOD( JSEED , M2 )
            M(I)  = J0 + M2 * J1
C
         ENDDO
C
         I = 5
         J = 17
C
      ENDIF
C
C  BEGIN MAIN LOOP HERE
C
      K = M(I) - M(J)
      IF ( K .LE. 0 ) K = K + M1
C
      M(J) = K
C
      I = I - 1
      IF ( I .EQ. 0 ) I = 17
C
      J = J - 1
      IF ( J .EQ. 0 ) J = 17
C
      UNI = FLOAT( K ) / FLOAT( M1 )
c
      if ( uni .eq. 0.0 ) then
c
         write ( 6 , * ) ' UNI = 0 '
         read ( 5 , * )
c
      endif
C
      RETURN
      END
C**
C***********************************************************************
C**
      REAL FUNCTION RNOR ( JD )
C**
C***********************************************************************
C**
C***BEGIN PROLOGUE  RNOR
C***DATE WRITTEN   810915
C***REVISION DATE  830805
C***CATEGORY NO.  L6A14
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, COMPUTER SCIENCE DEPT., WASH STATE UNIV
C
C***PURPOSE  GENERATES QUASI NORMAL RANDOM NUMBERS, WITH MEAN ZERO AND
C             UNIT STANDARD DEVIATION, AND CAN BE USED WITH ANY COMPUTER
C             WITH INTEGERS AT LEAST AS LARGE AS 32767.
C***DESCRIPTION
C
C       RNOR GENERATES QUASI NORMAL RANDOM NUMBERS WITH ZERO MEAN AND
C       UNIT STANDARD DEVIATION.
C       IT CAN BE USED WITH ANY COMPUTER WITH INTEGERS AT LEAST AS
C       LARGE AS 32767.
C
C
C   USE
C       FIRST TIME....
C                   Z = RNOR(JD)
C                     HERE JD IS ANY  N O N - Z E R O  INTEGER.
C                     THIS CAUSES INITIALIZATION OF THE PROGRAM
C                     AND THE FIRST RANDOM NUMBER TO BE RETURNED AS Z.
C       SUBSEQUENT TIMES...
C                   Z = RNOR(0)
C                     CAUSES THE NEXT RANDOM NUMBER TO BE RETURNED AS Z.
C
C.....................................................................
C
C    NOTE: USERS WHO WISH TO TRANSPORT THIS PROGRAM TO OTHER
C           COMPUTERS SHOULD READ THE FOLLOWING ....
C
C   MACHINE DEPENDENCIES...
C      MDIG = A LOWER BOUND ON THE NUMBER OF BINARY DIGITS AVAILABLE
C              FOR REPRESENTING INTEGERS, INCLUDING THE SIGN BIT.
C              THIS MUST BE AT LEAST 16, BUT CAN BE INCREASED IN
C              LINE WITH REMARK A BELOW.
C
C   REMARKS...
C     A. THIS PROGRAM CAN BE USED IN TWO WAYS:
C        (1) TO OBTAIN REPEATABLE RESULTS ON DIFFERENT COMPUTERS,
C            SET 'MDIG' TO THE SMALLEST OF ITS VALUES ON EACH, OR,
C        (2) TO ALLOW THE LONGEST SEQUENCE OF RANDOM NUMBERS TO BE
C            GENERATED WITHOUT CYCLING (REPEATING) SET 'MDIG' TO THE
C            LARGEST POSSIBLE VALUE.
C     B. THE SEQUENCE OF NUMBERS GENERATED DEPENDS ON THE INITIAL
C          INPUT 'JD' AS WELL AS THE VALUE OF 'MDIG'.
C          IF MDIG=16 ONE SHOULD FIND THAT
C            THE FIRST EVALUATION
C              Z=RNOR(87) GIVES  Z=-.40079207...
C            THE SECOND EVALUATION
C              Z=RNOR(0) GIVES   Z=-1.8728870...
C            THE THIRD EVALUATION
C              Z=RNOR(0) GIVES   Z=1.8216004...
C            THE FOURTH EVALUATION
C              Z=RNOR(0) GIVES   Z=.69410355...
C            THE THOUSANDTH EVALUATION
C              Z=RNOR(0) GIVES   Z=.96782424...
C
C***REFERENCES  MARSAGLIA & TSANG, "A FAST, EASILY IMPLEMENTED
C                 METHOD FOR SAMPLING FROM DECREASING OR
C                 SYMMETRIC UNIMODAL DENSITY FUNCTIONS", TO BE
C                 PUBLISHED IN SIAM J SISC 1983.
C***END PROLOGUE  RNOR
C
      IMPLICIT NONE
C
      REAL    V(65) , W(65) , RMAX , AA , B , C , C1 , C2 , PC , XN ,
     .        X , Y , S , UNI
C
      INTEGER M(17) , I , J , I1 , J1 , M1 , M2 , MDIG , JD , JSEED ,
     .        J0 , K0 , K1
C
      LOGICAL FLIP11 , FLIP12
C
      SAVE
C
      DATA    AA , B , C / 12.37586 , 0.4878992 , 12.67706 /
C
      DATA    C1 , C2 , PC , XN / 0.9689279 , 1.301198 , 0.1958303E-1 ,
     .                            2.776994 /
C
      DATA    V /  .3409450,  .4573146,  .5397793,  .6062427,  .6631691,
     .  .7136975,  .7596125,  .8020356,  .8417227,  .8792102,  .9148948,
     .  .9490791,  .9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917,
     . 1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635,
     . 1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929,
     . 1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454,
     . 1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422,
     . 1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166,
     . 1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713,
     . 2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117,
     . 2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943/
C
C***FIRST EXECUTABLE STATEMENT  RNOR
C
      IF ( JD .NE. 0 ) THEN
C
         MDIG = 32
C
C   BE SURE THAT MDIG AT LEAST 16...
C
         IF ( MDIG .LT. 16 ) WRITE ( * , * ) 'RNOR--MDIG LESS THAN 16'
C
         M1 = 2 ** ( MDIG - 2 ) + ( 2 ** ( MDIG - 2 ) - 1 )
         M2 = 2 ** ( MDIG / 2 )
C
         JSEED = MIN0( IABS( JD ) , M1 )
C
         IF ( MOD( JSEED , 2 ) .EQ. 0 ) JSEED = JSEED - 1
C
         K0 = MOD( 9069 , M2 )
         K1 = 9069 / M2
         J0 = MOD( JSEED , M2 )
         J1 = JSEED / M2
C
         DO I = 1 , 17
C
            JSEED = J0 * K0
            J1    = MOD( JSEED / M2 + J0 * K1 + J1 * K0 , M2 / 2 )
            J0    = MOD( JSEED , M2 )
            M(I)  = J0 + M2 * J1
C
         ENDDO
C
         J1 = 17
         I1 = 5
C
         RMAX = 1. / FLOAT( M1 )
C
         DO I = 1 , 65
C
            W(I) = RMAX * V(I)
C
         ENDDO
C
C   SEED UNIFORM (0,1) GENERATOR
C
         RNOR = UNI( JD )
C
      ENDIF  !  JD .NE. 0
C
C   FAST PART...
C
      I = M( I1 ) - M( J1 )
      IF ( I .LT. 0 ) I = I + M1
C
      M(J1) = I
C
      I1 = I1 - 1
      IF ( I1 .EQ. 0 ) I1 = 17
C
      J1 = J1 - 1
      IF ( J1 .EQ. 0 ) J1 = 17
C
      J    = MOD( I , 64 ) + 1
C
      RNOR = FLOAT( I ) * W( J+1 )
C
      IF ( MOD( I / M2 , 2 ) .EQ. 0 ) RNOR = - RNOR
C
C   Fast part complete - test if slow part necessary
C
C   In the original code, there were two ways to exit - a normal return
C   and an execution of 'GO TO 11', which does a flip operation and
C   then returns.
C
      FLIP11 = .FALSE.
      FLIP12 = .FALSE.
C
      IF ( ABS( RNOR ) .GT. V(J) ) THEN
C
C   SLOW PART; AA IS A*F(0)
C
         X = ( ABS( RNOR ) - V(J) ) / ( V(J+1) - V(J) )
         Y = UNI( 0 )
         S = X + Y
C
         IF ( S .GT. C2 ) FLIP11 = .TRUE.  !  Accomplishes 'GO TO 11'
C
         IF ( .NOT. FLIP11 ) THEN
C
            IF ( S .GT. C1 ) THEN
C
               FLIP12 = .FALSE.
C
               IF ( Y .GT. C - AA * EXP( -.5 * ( B - B * X ) ** 2 ) )
     .              FLIP12 = .TRUE.        !  Accomplishes 'GO TO 11'
C
               IF ( .NOT. FLIP12 ) THEN
C
                  IF ( EXP( -.5 * V(J+1) ** 2 ) + Y * PC / V(J+1) .GT.
     .                         EXP( - 0.5 * RNOR ** 2 ) ) THEN
C
C   TAIL PART; 3.855849 IS 0.5 * XN ** 2
C
   22                S = XN - ALOG( UNI( 0 ) ) / XN
C
                     IF ( ( 3.855849 + ALOG( UNI( 0 ) ) - XN * S ) .GT.
     .                  ( -.5 * S ** 2 ) )                   GO TO 22
C
                     RNOR = SIGN( S , RNOR )
C
                  ENDIF
C
               ENDIF
C
            ENDIF
C
         ENDIF
C
      ENDIF
C
   11 IF ( FLIP11 .OR. FLIP12 ) RNOR = SIGN( B - B * X , RNOR )
C
      RETURN
      END
