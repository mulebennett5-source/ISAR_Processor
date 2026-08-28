C***********************************************************************
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                     SCPLib.for,  Version 3.51                     **
C**                                                                   **
C**                                                                   **
C**                Simple Command Parser Routines for the             **
C**                                                                   **
C**                         ERIM Ocean Model                          **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                           Release Date                            **
C**                           May 15, 1991                            **
C**                                                                   **
C**                           Last Update                             **
C**                           Sept. 6, 2002                           **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                       Public Domain Software                      **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**    Routines by John R. Bennett and James T. Clinthorne, ERIM      **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**   MODULE    TYPE            DESCRIPTION           LAST MODIFIED   **
C**   -------------------------------------------------------------   **
C**                                                                   **
C**   COMCHK     INT. FUNCTION   Test Command vs list  Release Date   **
C**   GETCMD     SUBROUTINE      Get Next Command      Nov. 20, 1999  **
C**   GETLIN     INT. FUNCTION   Get line from input   Feb. 23, 2000  **
C**   GETI4      INT. FUNCTION   String to Integer     Release Date   **
C**   GETNI4     INT. FUNCTION   String to N Integers  Release Date   **
C**   GETR4      INT. FUNCTION   String to Real        Release Date   **
C**   GETNR4     INT. FUNCTION   String to N Reals     Release Date   **
C**   LFTSTR     INT. FUNCTION   Left Justify          Release Date   **
C**   RHTSTR     INT. FUNCTION   Right Justify         Release Date   **
C**   LOCASE     SUBROUTINE      To Lower Case         Release Date   **
C**   RCDERR     SUBROUTINE      Output error message  Release Date   **
C**   REAL_AFTER SUBROUTINE      Real after string     May 1, 1997    **
C**   INT_AFTER  SUBROUTINE      Integer after string  May 1, 1997    **
C**   READ_C80   SUBROUTINE      Read w/o comments     Feb. 23 2000   **
C**   LEFTCH     SUBROUTINE      Puts integer to char  Sept. 6 2002   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C***********************************************************************
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION COMCHK ( COMBUF , COMS , COMNUM , NCMD )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   This routine accepts a character string as input, and
C   determines if the string matches any command of a given set of
C   legal commands.
C
C  The parameters to COMCHK are defined as follows:
C
C      COMBUF :  Command string to be checked
C
C      COMS   :  Array of command list elements (see below).
C
C      COMNUM :  The maximum number of commands in the command
C                list.
C
C      NCMD   :  The number of the command matched
C
C      COMCHK :  Describes the success or failure of the call
C                to COMCHK
C
C                  =  0 => COMMAND MATCHED THE NCMD'TH
C                          COMMAND IN THE LIST
C                  =  1 => DID NOT MATCH
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      INTEGER   COMNUM , NCMD , K
C
      CHARACTER COMBUF*80 , COMS(COMNUM)*80
C
      DO K = 1 , COMNUM
C        
         IF ( COMBUF .EQ. COMS(K) ) GO TO 60
C
      ENDDO
C
      GO TO 70
C
C  The K'th element of the command list matched the string
C
   60 NCMD   = K
      COMCHK = 0
      RETURN
C
C   No command matched the string
C
   70 COMCHK = 1
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE GETCMD ( PROMPT , EXT , ISTR , TMPSTR , LHS , RHS ,
     .                    EFLAG , LUOUT )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C  AUTHORS: John Bennett, Jim Clinthorne
C     DATE: April 28, 1988
C
C   Inputs:
C
C PROMPT ==> Character*6 variable containing desired prompt
C EXT    ==> Character*4 variable containing .extension for command file
C            (i.e., .EOM)
C ISTR   ==> must be initially set to 0
C
C   Output:
C
C TMPSTR ==> Character*80 result containing single command word
C    LHS ==> Character*80 result containing left hand side of =
C    RHS ==> Character*80 result containing right hand side of =
C  EFLAG ==> Logical flag , true if = in string
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER ARRSTR(20)*80 , TMPSTR*80 , LHS*80 , RHS*80 , CMDSTR*80
      CHARACTER TEMP*6
      CHARACTER EXT*3
      CHARACTER PROMPT*3
C
      LOGICAL   EFLAG
C
      INTEGER   COMFLG , LUN , EPOS , J , LENGTH , ISTR , LUOUT ,
     .          NSTR , GETLIN , GETSTT
C
      SAVE
C
      DATA NSTR / 0 / , COMFLG / 0 /, LUN / 5 /
C
      TMPSTR = ' '
      EFLAG  = .FALSE.
C
      IF ( ISTR .EQ. 0 ) LUN = 5
C
      IF ( COMFLG .EQ. 2 ) THEN
         TEMP = ' ' // PROMPT // '} '
      ELSE
         TEMP = ' ' // PROMPT // '> '
      ENDIF
C
C   GETLIN returned NSTR command strings from the line - interpret them
C
      IF ( ISTR .GT. NSTR .OR. ISTR .EQ. 0 ) THEN
C
         ISTR = 1
         IF ( COMFLG .NE. 2 ) WRITE ( 6 , '(A6,$)' ) TEMP
         GETSTT = GETLIN( LUN , ARRSTR , NSTR )
C
      ENDIF
C
C   If the command line is the last line of a command file,
C   set the input to come from the keyboard again
C
      IF ( GETSTT .EQ. 2 ) THEN
C
         CLOSE ( UNIT = 3 )
         LUN    = 5
         COMFLG = 0
         NSTR   = 0
C
      ENDIF
C
C   If there are no commands on the line, get another line
C
      IF ( GETSTT .NE. 0 ) RETURN
C
C   If the commands are coming from a command file, echo them
C
      IF ( COMFLG .EQ. 2  ) WRITE ( 6 , '(A6,A40)' ) TEMP , ARRSTR(ISTR)
C
C   GETLIN returned NSTR command strings from the line - send them back in
C          in TMPSTR
C
      TMPSTR = ARRSTR(ISTR)
C
C   If the first character is '@', execute a command file
C
      IF ( TMPSTR (1:1) .EQ. '@' ) THEN
C
C  Record the command to the log file
C
        IF ( LUOUT .NE. 0 ) WRITE ( LUOUT , '(A)' ) TMPSTR
C
         CMDSTR = TMPSTR(2:80)
         TMPSTR = '@'
C
      ENDIF
C
      LENGTH = INDEX(TMPSTR , ' ') - 1
C
      EPOS = INDEX(TMPSTR(1:LENGTH),'=')
C
      IF ( EPOS .NE. 0 ) THEN
C
         IF ( ( EPOS .EQ. 1 ) .OR. ( EPOS .EQ. LENGTH ) ) THEN
C
            CALL RCDERR ( '%GETCMD-E-BADEQ, Not an equation.' )
            CALL RCDERR ( '%GETCMD, Remainder of line ignored.' )
C
         ELSE
C
            EFLAG = .TRUE.
            LHS   = TMPSTR(1:EPOS-1)
            RHS   = TMPSTR(EPOS+1:LENGTH)
C
         ENDIF
C
      ENDIF
C
C   Read input from a command file if requested look for a '@'
C
      IF ( TMPSTR(1:1) .EQ. '@' ) THEN
C
         COMFLG = 1
         J      = INDEX ( CMDSTR , '.' )
C
C   COMFLG ==> a flag to flag errors on attaching input command file
C
         IF ( J .EQ. 0 ) THEN
C
            J = INDEX( CMDSTR , ' ' ) - 1
            IF ( J .NE. 0 ) CMDSTR = CMDSTR(1:J) // '.' // EXT
C
         ENDIF
C
         OPEN( UNIT = 3 , FILE = CMDSTR , STATUS='OLD' , ERR = 100 )
         TMPSTR = ' '
         COMFLG = 2
         LUN    = 3
C
      ENDIF
C
  100 IF ( COMFLG .EQ. 1 ) THEN
C
        CALL RCDERR ( '%GETCMD-E-COMFERR, Command file open error.' )
        COMFLG = 0
C
      ENDIF
C
C   Count the number of commands serviced
C
      ISTR = ISTR + 1
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION GETLIN ( LUN , ARRSTR , NSTR )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   Purpose:
C
C        To parse an input line into strings, assuming the separator
C   is a space character and that anything following a comment character
c   is ignored.
C**                                                                   **
C***********************************************************************
C**                                                                   **
C
C   Purpose:
C
C        To parse an input line into strings, assuming the separator
C   is a space character and that all folowing a semi-colon, a pound
C   sign or an exclamation point is ignored.
C
C        Modifications Jan. 8, 1998
C
C        1.  Convert tabs to spaces
C
C        2.  Allow spaces before and after equal signs
C
C        Modifications Feb. 23, 2000
C
C        1.  Convert all control characters to spaces
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER CBUF*80 , CTMP*80
C
      INTEGER   ISTR , POINTR , REMAIN , LENGTH , LUN , NSTR , LEN ,
     .          BRKPOS , BRK2 , BRK3 , BRK4 , ILET , NEWPTR
C
      CHARACTER ARRSTR(20)*80
C
      LOGICAL   FOUND_SPACE
C
C-----------------------------------------------------------------------
C
      LENGTH = LEN( CBUF )
C
      IF ( LUN .EQ. 5 ) THEN
C
         READ ( LUN , '(A80)' ) CBUF
C
      ELSE
C
         READ ( LUN , '(A80)' , END = 200 ) CBUF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   Strip off comments - assume that anything beyond a semi-colon,
C   a pound sign, an exclamation point or a percent sign is comment
C   information.
C
      BRKPOS = INDEX( CBUF , ';' )
      BRK2   = INDEX( CBUF , '#' )
      BRK3   = INDEX( CBUF , '!' )
      BRK4   = INDEX( CBUF , '%' )
C
      IF ( BRK2 .NE. 0 ) THEN
C
         IF ( BRKPOS .EQ. 0 ) THEN
            BRKPOS = BRK2
         ELSE
            BRKPOS = MIN( BRKPOS , BRK2 )
         ENDIF
C
      ENDIF
C
      IF ( BRK3 .NE. 0 ) THEN
C
         IF ( BRKPOS .EQ. 0 ) THEN
            BRKPOS = BRK3
         ELSE
            BRKPOS = MIN( BRKPOS , BRK3 )
         ENDIF
C
      ENDIF
C
      IF ( BRK4 .NE. 0 ) THEN
C
         IF ( BRKPOS .EQ. 0 ) THEN
            BRKPOS = BRK4
         ELSE
            BRKPOS = MIN( BRKPOS , BRK4 )
         ENDIF
C
      ENDIF
C
      IF ( BRKPOS .NE. 0 ) THEN
C
         IF ( BRKPOS .EQ. 1 ) THEN
            CBUF = ' '
         ELSE
            CBUF = CBUF(1:BRKPOS-1)
         ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   Convert upper case letter to lower case
C
      DO ISTR = 1 , 80
C
         ILET = ICHAR( CBUF ( ISTR : ISTR ) )
C
         IF ( ILET .GE. 65 .AND. ILET .LE. 90 ) THEN
C
            ILET = ILET + 32
            CBUF(ISTR:ISTR) = CHAR( ILET )
C
         ENDIF
C
C   Set tabs or other control characters to spaces
C
         IF ( ILET .LE. 26 ) CBUF(ISTR:ISTR) = CHAR( 32 )
C
      ENDDO
C
C-----------------------------------------------------------------------
C
C   Remove excess spaces - leave one to separate strings and remove all
C   around equal signs
C
      CTMP        = CBUF
      CBUF        = ' '
      FOUND_SPACE = .FALSE.
      NEWPTR      = 0
C
      DO ISTR = 1 , 80
C
         ILET = ICHAR( CTMP ( ISTR : ISTR ) )
C
         IF ( ( ILET .NE. 32 ) .OR. ( .NOT. FOUND_SPACE ) ) THEN
C
            NEWPTR = NEWPTR + 1
C
C   If this is an equal sign and there has already been a space, then
C   cancel bumping up the pointer to put the equal right after the
C   previous non-space character
C
            IF ( FOUND_SPACE .AND. CTMP ( ISTR : ISTR ) .EQ. '=' )
     .           NEWPTR = NEWPTR - 1
C
            CBUF( NEWPTR : NEWPTR ) = CHAR( ILET )
C
         ENDIF
C
C   Do not add a space after either a space or an equal sign
C
         FOUND_SPACE = ( ILET .EQ. 32 ) .OR.
     .                 ( CTMP ( ISTR : ISTR ) .EQ. '=' )
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      NSTR   = 0
C
      POINTR = 1
      REMAIN = LENGTH
C
      DO ISTR = 1 , 80
C
         IF ( REMAIN .LE. 0 ) GO TO 100
C
         BRKPOS = INDEX( CBUF(POINTR:LENGTH) , ' ' ) + POINTR - 1
C
         IF ( BRKPOS .NE. 0 ) THEN
C
            IF ( BRKPOS .GT. POINTR ) THEN
C
               NSTR         = NSTR + 1
               ARRSTR(NSTR) = CBUF(POINTR:BRKPOS-1)
C
            ENDIF
C
            REMAIN = REMAIN - BRKPOS + POINTR - 1
            POINTR = BRKPOS + 1
C
         ELSE
C
            NSTR         = NSTR + 1
            ARRSTR(NSTR) = CBUF(POINTR:LENGTH)
            GO TO 100
C
         ENDIF
C
      ENDDO
C
  100 CONTINUE
C
      IF ( NSTR .GT. 0 ) THEN
C
         GETLIN = 0
C
      ELSE
C
         GETLIN = 1
C
      ENDIF
C
      RETURN
C
  200 GETLIN = 2
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION GETI4 ( STRI4 , LSTR , I4VAL )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   THIS ROUTINE ACCEPTS A STRING AS INPUT AND SEARCHES FOR
C   VALID ASCII VALUES. THE RETURNED VALUE IS A TWO BYTE INTEGER.
C   IF VALID CHARACTERS ARE PRESENT A ZERO STATUS IS RETURNED.
C
C   THE PARAMETERS ARE:
C
C        IASC  : USED IN COMPARING LEGAL ASCII VALUES FOR (DIGITS 1 - 9)
C
C        ICHK  : LOOP 20 RUNS THROUGH ASCII CHARACTERS 48 - 57
C
C        ISTR  : POINTER INTO INPUT STRING "STRI4"
C
C        I4VAL : FOUR BYTE INTEGER VALUE; RETURNED
C
C        LSTR  : LENGTH OF STRING IN BYTES; INPUT
C
C        STRI4 : CHARACTER STRING CONTAINING NUMBER; INPUT
C
C        RETURN STATUS:
C
C              = 0 => GOOD STATUS
C              = 1 => BAD STATUS
C
C        USER NOTE:
C
C              IF A FLOATING POINT VALUE IS CONTAINED IN THE INPUT
C              STRING IT WILL BE EVALUATED BY TRUNCATION.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER STRI4*80
      CHARACTER TMPSTR*8 , TMPSR1*8
C
      INTEGER   LSTR , ICHK , ISTR , IASC , I4VAL , ICHAR , ISTAT ,
     .          RHTSTR
C
C   CHECK EACH CHARACTER IN STRR4 AGAINST A LIST OF LEGAL ASCII VALUES
C   (DIGITS 0 - 9, '.' AND '-', FOR NEGATIVE VALUES). IF AN ILLEGAL
C   DIGIT OCCURS EXIT WITH ERROR STATUS; I.E., GETI4 = 1. IF A LEGAL
C   CHARACTER IS MATCHED FOR THE CURRENT CHARACTER GO BACK TO LABEL 10
C   AND CHECK THE NEXT CHARACTER AGAINST THE LIST OF POSSIBLE ASCII
C   VALUES.
C
      DO 10 ISTR = 1 , LSTR
C
      IASC = ICHAR ( STRI4 ( ISTR : ISTR ) )
      IF ( IASC .EQ. 32 .OR. IASC .EQ. 46 ) GO TO 25
C
      DO 20 ICHK = 48 , 57
      IF ( IASC .NE. ICHK ) GETI4 = 1
      IF ( IASC .EQ. ICHK .OR. (IASC .EQ. 45
     .         .AND. ISTR .EQ. 1 ) ) THEN
         GETI4 = 0
         GO TO 10
      ENDIF
C
   20 CONTINUE
C
      IF ( GETI4 .EQ. 1 ) RETURN
C
   10 CONTINUE
C
   25 GETI4 = 0
      ISTR = ISTR - 1
C
      TMPSTR = STRI4 ( 1 : ISTR )
      ISTAT  = RHTSTR(TMPSTR,8,TMPSR1)
      READ ( TMPSR1 , '(I8)' ) I4VAL
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION GETNI4 ( STRI4 , LSTR , I4ARR , NARGS )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   THIS ROUTINE ACCEPTS A STRING AS INPUT AND SEARCHES FOR
C   VALID ASCII VALUES. THE RETURNED VALUE IS A FLOATING POINT
C   NUMBER. IF VALID CHARACTERS ARE PRESENT, A ZERO STATUS IS
C   RETURNED.
C
C   THE PARAMETERS ARE:
C
C        I4ARR : INTEGER ARRAY RETURNED BY FUNCTION
C
C        LSTR  : LENGTH OF STRING IN BYTES; INPUT
C
C        STRI4 : CHARACTER STRING CONTAINING NUMBER; INPUT
C
C        NARGS : NUMBER OF NUMERICAL ARGUMENTS IN LIST
C
C        RETURN STATUS:
C
C              = 0 => GOOD STATUS
C              = 1 => BAD STATUS
C
C        USER NOTE:
C
C              FLOATING POINT VALUES WILL BE TRUNCATED
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER STRI4*80
C
      INTEGER   LSTR , LENGTH , EPOS1 , EPOS2 , EPOS3 , EPOS ,
     .          GETI4 , NARGS , I , LSTRTMP
C
      INTEGER   I4ARR(*)
C
      DO I = 1 , LSTR
         IF ( STRI4( (LSTR+1) - I:(LSTR+1) - I ) .NE. ' ' ) GO TO 7
      ENDDO
C
C   Fix March 12, 1999 - JRB
C
C   Use temporary variable for LSTR since the calling routine may send
C   LSTR as a hardwired constant.
C
    7 LSTRTMP = ( LSTR + 1 ) - I
C
      DO 10 NARGS = 1 , 21
C
         EPOS1 = INDEX( STRI4(1:LSTRTMP) , ',' )
         EPOS2 = INDEX( STRI4(1:LSTRTMP) , ':' )
         EPOS3 = INDEX( STRI4(1:80) , ' ' )
         IF ( ( EPOS1 .NE. 0 ) .OR. ( EPOS2 .NE. 0 ) .AND.
     .              ( NARGS .LE. 20 ) ) THEN
C
            EPOS = MAX( EPOS1 , EPOS2 )
            LENGTH = EPOS - 1
C
            IF ( GETI4( STRI4 , LENGTH , I4ARR(NARGS) )
     .           .EQ. 0 ) THEN
               GETNI4 = 0
            ELSE
               GETNI4 = 1
            ENDIF
C
         ELSE IF ( NARGS .GE. 21 ) THEN
C
            GETNI4 = 1
C
         ELSE IF ( ( ( EPOS3 .NE. 0 ) .AND. ( EPOS1 .EQ. 0 )
     .         .AND. ( EPOS2 .EQ. 0 ) ) .AND. ( NARGS .LE. 20 ) ) THEN
C
            GETNI4 = 1
            LENGTH = EPOS3 - 1
            IF ( STRI4(LENGTH:LENGTH) .NE. ',' )
     .         GETNI4 = GETI4( STRI4 , LENGTH , I4ARR(NARGS) )
            RETURN
C
         ENDIF
C
         STRI4(1:80) = STRI4(EPOS+1:80)
         IF ( STRI4 .EQ. ' ' .OR. GETNI4 .EQ. 1 ) THEN
            GETNI4 = 1
            RETURN
         ENDIF
C
   10 CONTINUE
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION GETR4 ( STRR4 , LSTR , R4VAL )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   THIS ROUTINE ACCEPTS A STRING AS INPUT AND SEARCHES FOR
C   VALID ASCII VALUES. THE RETURNED VALUE IS A FLOATING POINT
C   NUMBER. IF VALID CHARACTERS ARE PRESENT, A ZERO STATUS IS
C   RETURNED.
C
C   THE PARAMETERS ARE:
C
C        CHKPNT : CHECK FOR A "." IF MORE THAN ONE '.', EXIT WITH ERROR
C
C        ICHK   : LOOP 20 RUNS THROUGH ASCII CHARACTERS 48 - 57
C
C        ISTR   : POINTER INTO INPUT STRING "STRR4"
C
C        R4VAL  : FOUR BYTE FLOATING POINT VALUE; RETURNED
C
C        IASC   : USED IN COMPARING LEGAL ASCII VALUES FOR (DIGITS 1 - 9)
C
C        LSTR   : LENGTH OF STRING IN BYTES; INPUT
C
C        STRR4  : CHARACTER STRING CONTAINING NUMBER; INPUT
C
C        RETURN STATUS:
C
C              = 0 => GOOD STATUS
C              = 1 => BAD STATUS
C
C        USER NOTE:
C
C              INTEGER DATA VALUES WILL BE CONVERTED TO FLOATING POINT.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER STRR4*80
      CHARACTER TMPSTR*12
      INTEGER   LSTR , ISTR , IASC , CHKPNT , ICHAR
      REAL      R4VAL
C
C   CHECK EACH CHARACTER IN STRR4 AGAINST A LIST OF LEGAL ASCII VALUES
C   (DIGITS 0 - 9, '.' AND '-', FOR NEGATIVE VALUES). IF AN ILLEGAL
C   DIGIT OCCURS EXIT WITH ERROR STATUS; I.E., GETR4 = 1. IF A LEGAL
C   CHARACTER IS MATCHED FOR THE CURRENT CHARACTER GO BACK TO LABEL 10
C   AND CHECK THE NEXT CHARACTER AGAINST THE LIST OF POSSIBLE ASCII
C   VALUES.
C
      CHKPNT = 0
C
      DO 10 ISTR = 1 , LSTR
C
      IASC = ICHAR( STRR4(ISTR:ISTR) )
C
      GETR4 = 1
C
      IF ( IASC .EQ. 46 ) THEN
C
         GETR4  = 0
         CHKPNT = CHKPNT + 1
C
      ELSE IF ( IASC .EQ. 32 ) THEN
C
         GETR4 = 0
         GO TO 25
C
      ELSE IF ( IASC .LE. 57 .AND. IASC .GE. 48 ) THEN
C
         GETR4 = 0
C
      ELSE IF ( IASC .EQ. 45 .AND. ISTR .EQ. 1 ) THEN
C
         GETR4 = 0
C
      ENDIF
C
      IF ( GETR4 .EQ. 1 ) RETURN
C
   10 CONTINUE
C
   25 CONTINUE
C
C   IF MORE THAN ONE '.' IS DETECTED, RETURN WITH ERROR STATUS,
C   ELSE IF NO '.' IS DETECTED APPEND A '.' TO THE END OF THE STRING.
C
      IF ( CHKPNT .GT. 1 ) THEN
C
         GETR4 = 1
         RETURN
C
      ELSE IF ( CHKPNT .EQ. 0 ) THEN
C
         STRR4(ISTR:ISTR) = '.'
         ISTR = ISTR + 1
C
      ENDIF
C
C   DECREMENT POINTER INTO STRING STRR4 AND CONVERT TO A 4 BYTE FLOATING
C   POINT NUMBER, RETURN WITH A 0 ERROR STATUS.
C
      ISTR   = ISTR - 1
      TMPSTR = STRR4(1:ISTR)
      READ ( TMPSTR , '(F12.0)' ) R4VAL
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION GETNR4 ( STRR4 , LSTR , R4ARR , NARGS )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   THIS ROUTINE ACCEPTS A STRING AS INPUT AND SEARCHES FOR
C   VALID ASCII VALUES. THE RETURNED VALUE IS A FLOATING POINT
C   NUMBER. IF VALID CHARACTERS ARE PRESENT, A ZERO STATUS IS
C   RETURNED.
C
C   THE PARAMETERS ARE:
C
C        R4ARR : FLOATING POINT ARRAY RETURNED BY FUNCTION
C
C        LSTR  : LENGTH OF STRING IN BYTES; INPUT
C
C        STRR4 : CHARACTER STRING CONTAINING NUMBER; INPUT
C
C        NARGS : NUMBER OF NUMERICAL ARGUMENTS IN LIST
C
C        RETURN STATUS:
C
C              = 0 => GOOD STATUS
C              = 1 => BAD STATUS
C
C        USER NOTE:
C
C              INTEGER DATA VALUES WILL BE CONVERTED TO FLOATING POINT.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER STRR4*80
      INTEGER   LSTR , LENGTH , EPOS1 , EPOS2 , EPOS3 , EPOS ,
     .          GETR4 , NARGS , I
      REAL      R4ARR(*)
C
C
      DO 5 I = 1 , LSTR
         IF ( STRR4( (LSTR+1) - I:(LSTR+1) - I ) .NE. ' ' ) GO TO 7
    5 CONTINUE
    7 LSTR = (LSTR+1) - I
C
      DO 10 NARGS = 1 , 21
C
         EPOS1 = INDEX( STRR4(1:LSTR) , ',' )
         EPOS2 = INDEX( STRR4(1:LSTR) , ':' )
         EPOS3 = INDEX( STRR4(1:80) , ' ' )
C
         IF ( ( ( EPOS1 .NE. 0 ) .OR. ( EPOS2 .NE. 0 ) ) .AND.
     .        ( NARGS .LE. 20 ) ) THEN
C
            EPOS = MAX( EPOS1 , EPOS2 )
            LENGTH = EPOS - 1
C
            IF ( GETR4( STRR4 , LENGTH , R4ARR(NARGS) ) .EQ. 0 ) THEN
               GETNR4 = 0
            ELSE
               GETNR4 = 1
            ENDIF
C
         ELSE IF ( NARGS .GE. 21 ) THEN
C
            GETNR4 = 1
C
         ELSE IF ( ( (EPOS3 .NE. 0) .AND. (EPOS1 .EQ. 0)
     .           .AND. (EPOS2 .EQ. 0)) .AND. ( NARGS .LE. 20 ) ) THEN
C
            GETNR4 = 1
            LENGTH = EPOS3 - 1
            IF ( STRR4(LENGTH:LENGTH) .NE. ',' )
     .         GETNR4 = GETR4( STRR4 , LENGTH , R4ARR(NARGS) )
            RETURN
C
         ENDIF
C
         STRR4(1:80) = STRR4(EPOS+1:80)
         IF ( STRR4 .EQ. ' ' .OR. GETNR4 .EQ. 1 ) THEN
            GETNR4 = 1
            RETURN
         ENDIF
C
   10 CONTINUE
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION LFTSTR ( STRIN , LSTR , STROUT )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   Written by: James Clinthorne, APRIL, 1988
C
C   Function to left justify a given input character string
C
C   STRIN  : Input character string
C   STROUT : Left justified output string
C
C   STATUS :
C           0 = Successfull left justify
C           1 = String was null
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER*(*) STRIN , STROUT
      CHARACTER     TMPSTR*1
      INTEGER       I , J , LSTR
C
      STROUT = ' '
      DO 20 J = 1 , LSTR
C
         TMPSTR = STRIN(J:J)
C
         IF(TMPSTR .NE. ' ') THEN
C
           DO 10 I = 1 , (LSTR-(J-1))
              STROUT(I:I) = STRIN(J+I-1:J+I-1)
   10      CONTINUE
           LFTSTR = 0
           RETURN
C
         ENDIF
C
   20 CONTINUE
C
      LFTSTR = 1
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      INTEGER FUNCTION RHTSTR ( STRIN , LSTR , STROUT )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C   Written by: James Clinthorne, APRIL, 1988
C
C   Function to right justify a given input character string
C
C   STRIN  : Input character string
C   STROUT : right justified output string
C
C   STATUS :
C           0 = Successfull right justify
C           1 = String was null
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      INTEGER       I , J , LSTR
      CHARACTER*(*) STRIN , STROUT
      CHARACTER     TMPSTR*1
C
      STROUT = ' '
      DO 20 J = 1 , LSTR
C
         TMPSTR = STRIN(LSTR-J+1:LSTR-J+1)
C
         IF(TMPSTR .NE. ' ') THEN
C
           DO 10 I = 1 , (LSTR-(J-1))
              STROUT(LSTR-I+1:LSTR-I+1) = STRIN(LSTR-J+2-I:LSTR-J+2-I)
   10      CONTINUE
           RHTSTR = 0
           RETURN
C
         ENDIF
C
   20 CONTINUE
C
      RHTSTR = 1
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE LOCASE ( STR , LENSTR )
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER STR*80 , NEWSTR*80
C
      INTEGER   LENSTR , I , IL
C
      LENSTR = 80
C
      NEWSTR = ' '
C
      DO 10 I = 1 , LENSTR
C
      IL = ICHAR ( STR(I:I) )
      IF ( IL .GE. 65 .AND. IL .LE. 90 ) IL = IL + 32
C
      IF ( I .EQ. 1 ) THEN
         NEWSTR = CHAR(IL)
      ELSE
         NEWSTR = NEWSTR(1:I-1) // CHAR(IL)
      ENDIF
C
   10 CONTINUE
C
      STR = NEWSTR
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE RCDERR ( ERRSTR )
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER*(*) ERRSTR
C
      INTEGER       LENGTH
C
      LENGTH = INDEX( ERRSTR , '.' ) - 1
      IF ( LENGTH .LE. 0 ) LENGTH = 1
C
      WRITE ( 6 , * )  ERRSTR(1:LENGTH)
C
      RETURN
      END
C**
C***********************************************************************
C**
      subroutine real_after ( in_line , str , val , control )
C**
C***********************************************************************
C**
      implicit none
c
      character*(*) in_line , str
c
      real          val
c
      integer       idx , idxend , dummy , lftstr , control
c      
      character     tmpstr*80
c      
      idx    = index( in_line , str )
      idx    = idx + len( str )
      idxend = len( in_line )
c
      if ( control .eq. 0 ) then
c
         dummy  = lftstr( in_line(idx:idxend) ,
     .                    len( in_line(idx:idxend) ) , tmpstr )
c
      else
c
         tmpstr = in_line(idx:idx+control)
c
      endif   
c
      read( tmpstr , '(f7.3)' ) val
c      
      return 
      end
C**
C***********************************************************************
C**
      subroutine int_after ( in_line , str , val , control )
C**
C***********************************************************************
C**
      implicit none
c
      character*(*) in_line , str
c
      integer       val , dummy , lftstr , idx , idxend , control
c      
      character     tmpstr*80
c      
      idx    = index( in_line , str )
      idx    = idx + len( str )
      idxend = len( in_line )
c
      if ( control .eq. 0 ) then
c
         dummy  = lftstr( in_line(idx:idxend) ,
     .                    len( in_line(idx:idxend) ) , tmpstr )
c
      else
c
         tmpstr = in_line(idx:idx+control)
c
      endif
c
      read ( tmpstr , '(I7)' ) val
c      
      return 
      end
c**
c***********************************************************************
c**
      subroutine read_c80  ( luin , luout , string80 )
c**
c***********************************************************************
c**
      implicit none
c
c   Read a line of text from the logical unit luin, ignoring comments
c
c-----------------------------------------------------------------------
c
      integer      luin , luout , ipound , iexcla , isemic , iperct ,
     .             icomnt , idummy , lftstr
c
      character    string80*80 , temp80*80
c
      logical      comment
c
      comment = .true.
c
      do while ( comment )
c
         read ( luin , '(a)' ) temp80
c
c   Find first comment character in line
c
         ipound = 81
c
         iexcla = 81
c
         isemic = 81
c
         iperct = 81
c
         if ( index(temp80,'#') .ne. 0 ) ipound = index(temp80,'#')
c
         if ( index(temp80,'!') .ne. 0 ) iexcla = index(temp80,'!')
c
         if ( index(temp80,';') .ne. 0 ) isemic = index(temp80,';')
c
         if ( index(temp80,'%') .ne. 0 ) iperct = index(temp80,';')
c
         icomnt = min( ipound , iexcla )
c
         icomnt = min( icomnt , isemic )
c
         icomnt = min( icomnt , iperct )
c
c   If the line is all comment then write out the line to the output
c   file
c
         if ( icomnt .eq. 1 .or. temp80(1:icomnt-1) .eq. ' ' ) then
c
            comment = .true.
c
c           write ( luout , * ) icomnt , ipound , iexcla , isemic
c
            write ( luout , '(a79)' ) temp80(1:79)
c
         else
c
            idummy  = lftstr ( temp80(1:icomnt-1) , icomnt-1 ,
     .                         string80 )
c
            comment = .false.
c
         endif
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      subroutine leftch ( number , chrs , nchrs )
C**
C***********************************************************************
C**
c
c   Writes the integer 'number' into the left side of the character
c   string 'chrs' and returns the number of characters required (nchrs)
c
      implicit none
c
      integer   number , nchrs
c
      character chrs*9
c
      if ( number .lt. 10 ) then
c
         write( chrs , '(i1)' ) number
         nchrs = 1
c
      else if ( number .lt. 100 ) then
c
         write( chrs , '(i2)' ) number
         nchrs = 2
c
      else if ( number .lt. 1000 ) then
c
         write( chrs , '(i3)' ) number
         nchrs = 3
c
      else if ( number .lt. 10000 ) then
c
         write( chrs , '(i4)' ) number
         nchrs = 4
c
      else if ( number .lt. 100000 ) then
c
         write( chrs , '(i5)' ) number
         nchrs = 5
c
      else if ( number .lt. 1000000 ) then
c
         write( chrs , '(i6)' ) number
         nchrs = 6
c
      else if ( number .lt. 10000000 ) then
c
         write( chrs , '(i7)' ) number
         nchrs = 7
c
      else if ( number .lt. 100000000 ) then
c
         write( chrs , '(i8)' ) number
         nchrs = 8
c
      else if ( number .lt. 1000000000 ) then
c
         write( chrs , '(i9)' ) number
         nchrs = 9
c
      endif
c
      return
      end
C***********************************************************************
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                     HLPLib.for,  Version 3.50                     **
C**                                                                   **
C**                                                                   **
C**                       Help Routines for the                       **
C**                                                                   **
C**                         ERIM Ocean Model                          **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                           Release Date                            **
C**                           May 15, 1991                            **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                       Public Domain Software                      **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**   MODULE    TYPE            DESCRIPTION           LAST MODIFIED   **
C**   -------------------------------------------------------------   **
C**                                                                   **
C**                                                                   **
C**              Original routines by R. Darin Dean, ERIM             **
C**                                                                   **
C**   HLP        SUBROUTINE      VAX-like Help         Release Date   **
C**   CMDGET     SUBROUTINE      Used by HLP           Release Date   **
C**   CAPS       SUBROUTINE      Used by HLP           Release Date   **
C**   HELP       SUBROUTINE      Used by HLP           Release Date   **
C**   CMDPRC     SUBROUTINE      Used by HLP           Release Date   **
C**   DISSUB     SUBROUTINE      Used by HLP           Release Date   **
C**   ABREV      SUBROUTINE      Used by HLP           Release Date   **
C**   HELPER     SUBROUTINE      Basic calling routine May 1, 1997    **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C**                                                                   **
C***********************************************************************
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE HLP ( COMAND , HLPFIL )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C        This routine is designed to allow the user to issue commands
C     at the help level.  It must first read the entire help text file
C     to determine where within it the command words and their
C     corresponding text are.  There are 2 help levels, a TOPIC and
C     a SUBTOPIC level.  Not all commands have a SUBTOPIC level.
C     Those that do will display their SUBTOPICS and give the user a
C     chance to interrogate these SUBTOPIC commands.  Abbreviated
C     command words at the help level are allowed, but if more than
C     one command word is found that matches all matches will be
C     displayed and the topic level will not be changed.  SUBTOPICS
C     can only be accessed when at the SUBTOPIC level.  Matches will
C     be based on the last word of the user's command sequence.  Input
C     commands can be any mixture of upper and lower case.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER TEMSTR*80 , COMAND*80 , HLPFIL*80
      CHARACTER VALIDC(100)*20
C
      INTEGER   INDX(100,2) , TOPICL , COMMN , NCOMDS , RETCNT ,
     .          ITOPIC , RECN , NOSPC , RECLDA
C
      OPEN ( UNIT = 11 , FILE = HLPFIL , ACCESS = 'DIRECT' ,
     .       FORM = 'FORMATTED' , STATUS = 'OLD' ,
     .       RECL = RECLDA( 80 , 0 ) , ERR = 250 )
C
C      REWIND ( UNIT = 11 )
C
      RECN  = 1
      COMMN = 1
C
C
C      This section of code reads in the file opened above one record
C   at a time and keeps track of where it finds command words and the
C   topic level value of each command.  The file it reads must be a
C   direct access file with 80 byte records.  The contents of the file
C   must contain a "1" in the first column followed by at least one
C   space and then the command word in capitals for each main topic, and
C   a "2" in the first column followed by at least one space and the
C   command word in capitals for each sub-topic.  The sub-topics for a
C   main topic must also follow it and its text within the file.
C
C
 10   FORMAT ( A80 )
C
 20   READ ( 11 , REC = RECN , FMT = 10 , ERR = 100 ) TEMSTR
C
      IF ( TEMSTR .EQ. 'EOF' ) GO TO 100
C
      IF ( TEMSTR(1:1) .NE. ' ' ) THEN
C
         INDX(COMMN,1) = RECN
C
         IF ( TEMSTR(1:1) .EQ. '1' ) THEN
            INDX(COMMN,2) = 1
         ELSE
            INDX(COMMN,2) = 2
         ENDIF
C
         CALL CMDGET ( TEMSTR )
         VALIDC(COMMN) = TEMSTR
         COMMN         = COMMN + 1
C
      ENDIF
C
      RECN = RECN + 1
      GO TO 20
C
 100  INDX(COMMN,1) = RECN
      TOPICL        = 1
      NCOMDS        = COMMN - 1
      COMMN         = 1
      RETCNT        = 0
      ITOPIC        = 0
C
C      For each command issued at the help level, call the CAPS routine.
C   This routine will return the last word of the command string most
C   recently issued in capital letter form.
C
C
 200  CALL CAPS ( COMAND )

      WRITE ( 6 , * ) ' '
C
C
C      If the last word of the command string was HELP then call CMDPRC
C   to display the text for the HELP command and then call the special
C   HELP routine to display all of the main topics there is text for.
C   If the command was anything else besides a carriage return then just
C   call the CMDPRC routine to display its textual information as well
C   as any sub-topics it may have.
C
C
      IF ( COMAND .EQ. 'HELP' .AND. TOPICL .EQ. 1 ) THEN
C
         CALL CMDPRC ( COMAND , INDX , NCOMDS , VALIDC , TOPICL ,
     .                 COMMN )
         CALL HELP ( INDX , NCOMDS , VALIDC )
C
      ELSE IF ( COMAND .NE. ' ' ) THEN
C
         CALL CMDPRC ( COMAND , INDX , NCOMDS , VALIDC , TOPICL ,
     .                 COMMN )
C
      ENDIF
C
C
C      Display the correct prompt for the present help level, and wait
C   for user input.
C
C
      WRITE ( 6 , * )  ' '
      IF ( TOPICL .EQ. 1 ) THEN
         WRITE ( 6 , '(A9,$)' ) ' Topic?  '
      ELSE
         WRITE ( 6 , '(A12,$)' ) ' Subtopic?  '
      ENDIF
C
      READ ( 5 , '(A80)' ) COMAND
C
C      Keep track of how many times in a row the carriage return is hit
C   without any other keys being hit, and update the help level.
C
C
      IF ( COMAND .EQ. ' ' ) THEN
C
         RETCNT = RETCNT + 1
         IF ( RETCNT .EQ. 2 .OR.
     .        ( RETCNT .EQ. 1 .AND. TOPICL .EQ. 1 ) ) GO TO 300
C
         IF ( TOPICL .EQ. 2 ) THEN
           TOPICL = 1
           ITOPIC = 1
         ENDIF
C
      ELSE
C
         RETCNT = 0
C
      ENDIF
C
      GO TO 200
C
  250 NOSPC = INDEX( HLPFIL , ' ' )
      WRITE ( 6 , * ) ' Error opening help file: ' , HLPFIL(1:NOSPC-1)
C
  300 CLOSE ( UNIT = 11 )
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE CMDGET ( TEMSTR )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C      This routine simply takes the command record from the help file
C   and strips the level indicator off leaving just the command word for
C   storage in the valid command word array.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER TEMSTR*80
C
      INTEGER   X
C
      DO 50 X = 1 , 80
         IF ( TEMSTR(X:X) .NE. '1' .AND. TEMSTR(X:X) .NE. '2' .AND.
     .        TEMSTR(X:X) .NE. ' ' ) THEN
            TEMSTR = TEMSTR(X:80)
            RETURN
         ENDIF
 50   CONTINUE
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE CAPS ( COMAND )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C      This routine takes the incoming command string and returns just
C   the last word in capital letter form in the same string variable.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER COMAND*80
C
      INTEGER   X , I
C
      LOGICAL   FOUNDC
C
      FOUNDC = .FALSE.
C
      DO 10 X = 80 , 1 , -1
C
         IF ( COMAND(X:X) .NE. ' ' ) THEN
C
            I = ICHAR( COMAND(X:X) )
            IF ( I .GE. 97 .AND. I .LE. 122 ) THEN
               I = I - 32
               COMAND(X:X) = CHAR(I)
            ENDIF
            FOUNDC = .TRUE.
C
         ELSE IF ( ( COMAND(X:X) .EQ. ' ' ) .AND. ( FOUNDC ) ) THEN
C
            COMAND = COMAND(X+1:80)
            RETURN
C
         ENDIF
C
   10 CONTINUE
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE HELP ( INDX , NCOMDS , VALIDC )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C      This routine displays all of the command words that there
C   is help information available for.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER TEMSTR*80
      CHARACTER VALIDC(100)*20
C
      INTEGER   INDX(100,2) , NCOMDS , X , Y
C
      WRITE ( 6 , '(A25)' ) '  Information available :'
      WRITE ( 6 , '(A1)' ) ' '
C
      Y = 0
C
      DO 10 X = 1 , NCOMDS
C
         IF ( INDX(X,2) .EQ. 1 ) THEN
C
            TEMSTR(1+Y*20:(Y+1)*20) = VALIDC(X)
            Y = Y + 1
C
            IF ( Y .EQ. 4 .OR. X .EQ. NCOMDS ) THEN
C
                WRITE ( 6 , '(A80)' ) ' '//TEMSTR(1:79)
                TEMSTR = ' '
                Y      = 0
C
            ENDIF
C
         ENDIF
C
   10 CONTINUE
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE CMDPRC ( TCMD , INDX , NCOMDS , VALIDC , TOPICL ,
     .                    COMMN )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C      This routine actually finds all matches to the user input
C   command and displays the textual help information.  It is also
C   responsible for setting the correct topic level.  It will also
C   display all sub-topics of a command as well as the help
C   information for these sub-topics if requested by the user.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER TCMD*80 , TEMSTR*80
      CHARACTER VALIDC(100)*20 , COMAND*20
C
      INTEGER   INDX(100,2) , NCOMDS , TOPICL , COMMN , X , Y ,
     .          NUMFND
C
      LOGICAL   VALID , FOUND , ABREV
C
      X      = 1
      VALID  = .TRUE.
      FOUND  = .FALSE.
      COMAND = TCMD
C
C      If at sub-topic level and not at the end of the help file
C   commands then check next command of the help file for a match.
C
 5    IF ( ( TOPICL .EQ. 2 ) .AND. ( VALID ) .AND.
     .                             ( COMMN+X .LE. NCOMDS ) ) THEN
C
C      If next sub-topic command matches then display its help text and
C   go on to next command word of help file.
C
         IF ( ( INDX(COMMN+X,2) .EQ. 2 ) .AND.
     .        ( ABREV( COMAND , VALIDC(COMMN+X) ) ) ) THEN
            DO 50 Y = INDX(COMMN+X,1) , INDX(COMMN+X+1,1) - 1
               READ ( 11 , REC = Y , FMT = 10 , ERR = 1000 ) TEMSTR
 10            FORMAT ( A80 )
               WRITE ( 6 , '(A80)' ) TEMSTR(2:80)
 50         CONTINUE
            WRITE ( 6 , '(A1)' ) ' '
            FOUND = .TRUE.
C
C      If no matches are found at the sub-topic level then display
C   appropriate error message and the attempted match.
C
         ELSE IF ( INDX(COMMN+X,2) .EQ. 1 ) THEN
            IF ( .NOT. FOUND ) THEN
               TEMSTR = ' No information available on : ' //
     .                  VALIDC(COMMN) // COMAND
               WRITE ( 6 , '(A80)' ) TEMSTR
                    WRITE ( 6 , '(A1)' ) ' '
                    CALL DISSUB ( COMMN , INDX , VALIDC , NCOMDS ,
     .                            TOPICL )
            ENDIF
            VALID = .FALSE.
         ENDIF
         X = X + 1
      ELSE
         GO TO 100
      ENDIF
C
      GO TO 5
C
C     If at the main topic level then check the main topics for matches.
C
 100  IF ( TOPICL .EQ. 1 ) THEN
         NUMFND = 0
C
C     For each command of the help file check for a match.
C
         DO 150 Y = 1 , NCOMDS
C
C      If a match is found then update the number of matches found,
C   display the help text for the corresponding command and display the
C   subtopics for the corresponding command.
C
            IF ( ( INDX(Y,2) .EQ. 1 ) .AND.
     .           ( ABREV( COMAND , VALIDC(Y) ) ) ) THEN
               NUMFND = NUMFND + 1
C
               DO 160 X = INDX(Y,1) , INDX(Y+1,1) - 1
                  READ ( 11 , REC = X , FMT = 10 , ERR = 1000 ) TEMSTR
                  WRITE ( 6 , '(A80)' ) TEMSTR(2:80)
 160           CONTINUE
C
               COMMN = Y
               CALL DISSUB ( COMMN , INDX , VALIDC , NCOMDS , TOPICL )
            ENDIF
 150     CONTINUE
C
C      If more than one match was found then do not change topic level
C   from level one regardless of last matches sub-topic status.  If no
C   matches are found then display appropriate error message and re-
C   display the topics that there is help text for.
C
         IF ( NUMFND .GT. 1 ) THEN
            TOPICL = 1
         ELSE IF ( NUMFND .LT. 1 ) THEN
            TEMSTR = ' No information available on : ' // COMAND
            WRITE ( 6 , '(A80)' ) TEMSTR
            WRITE ( 6 , '(A1)' ) ' '
            CALL HELP ( INDX , NCOMDS , VALIDC )
         ENDIF
C
      ENDIF
C
      RETURN
C
 1000 WRITE ( 6 , * ) ' Error reading help file.'
C
      RETURN
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      SUBROUTINE DISSUB ( COMMN , INDX , VALIDC , NCOMDS , TOPICL )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C      This routine is responsible for displaying the list of sub-topic
C   commands that there is help text for.  It will also change the topic
C   level if it finds sub-topics for the incoming command.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER TEMSTR*80
      CHARACTER VALIDC(100)*20
C
      INTEGER   INDX(100,2) , COMMN , NCOMDS , X , TOPICL
C
      X = COMMN + 1
C
 10   IF ( X .GT. NCOMDS .OR. INDX(X,2) .EQ. 1 ) THEN
         GO TO 100
      ELSE
         IF ( X .EQ. COMMN + 1 ) THEN
            WRITE ( 6 , '(A1)' ) ' '
            WRITE ( 6 , '(A25)' ) '  Information available :'
            WRITE ( 6 , '(A1)' )
            TOPICL = 2
         ENDIF
         IF ( INDX(X+3,2) .EQ. 2 ) THEN
            TEMSTR = VALIDC(X) // VALIDC(X+1) // VALIDC(X+2) //
     .               VALIDC(X+3)
         ELSE IF ( INDX(X+2,2) .EQ. 2 ) THEN
            TEMSTR = VALIDC(X) // VALIDC(X+1) // VALIDC(X+2)
         ELSE IF ( INDX(X+1,2) .EQ. 2 ) THEN
            TEMSTR = VALIDC(X) // VALIDC(X+1)
         ELSE
            TEMSTR = VALIDC(X)
         ENDIF
         WRITE ( 6 , '(A80)' ) ' ' // TEMSTR(1:79)
         X = X + 4
         TEMSTR = ' '
      ENDIF
C
      GO TO 10
C
  100 RETURN
C
      END
C**                                                                   **
C***********************************************************************
C**                                                                   **
      LOGICAL FUNCTION ABREV ( COMAND , VALCOM )
C**                                                                   **
C***********************************************************************
C**                                                                   **
C        This routine will return the truth value in the function's
C     name if the possibly abbreviated command string COMAND matches
C     any part of the valid command string VALCOM.
C**                                                                   **
C***********************************************************************
C**                                                                   **
      IMPLICIT NONE
C
      CHARACTER COMAND*80
      CHARACTER VALCOM*20
C
      INTEGER   X , CHDIFF
C
      ABREV = .TRUE.
C
      DO 10 X = 1 , 20
C
         IF ( COMAND(X:X) .EQ. ' ' ) THEN
            RETURN
         ELSE
            CHDIFF = IABS( ICHAR( COMAND(X:X) ) - ICHAR( VALCOM(X:X) ) )
            IF ( CHDIFF .NE. 0 .AND. CHDIFF .NE. 32 ) THEN
               ABREV = .FALSE.
               RETURN
            ENDIF
         ENDIF
C
   10 CONTINUE
C
      RETURN
      END
C**
C***********************************************************************
C**
      subroutine helper ( comand , hlpfil )
C**
C***********************************************************************
C**
c   This routine reads a formatted help file, converts it to a direct
c   access file, and then calls a help utility to present it to the user
c   in the manner of the VAX help system.
c
      implicit none
c
      character comand*80 , hlpfil*80 , hlpfilda*80 , tmpstr*80
c
      integer   i , reclda
c
      logical   first
c
      save      first
c
      data      first / .true. /
c
      hlpfilda = 'help.da'
c
      if ( first ) then
c
         first = .false.
c
         open ( unit = 1 , file = hlpfil , access = 'SEQUENTIAL' ,
     .          status = 'old' )
c
         open ( unit = 2 , file = hlpfilda , access = 'DIRECT' ,
     .          form = 'FORMATTED' , status = 'UNKNOWN' ,
     .          recl = reclda( 80 , 0 ) , err = 1000 )
c
         i = 0
         do while ( .true. )
c
            read ( 1 , '(a80)' , end = 100 ) tmpstr
c
            i = i + 1
            write ( 2 , '(a80)' , rec = i ) tmpstr
c
         enddo
c
  100    close ( 1 )
         close ( 2 )
c
      endif
c
      call hlp ( comand , hlpfilda )
c
      return
c
 1000 continue
c
      call rcderr ( ' Error in opening help files.' )
c
      return
      end
