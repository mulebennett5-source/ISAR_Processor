C*-----------------------------------------------------------------------
C
C   Information derived from the pulse header of the Advanced Profile
C   Radar operated by Texas Instruments for the Naval Research Laboratory
C
      BYTE         DTYPE    !data type (image,iqpre,iqpost768,iqpost512)
      REAL         ACALT    !aircraft altitude in meters
      BYTE         SPARR    !sparse array threshold (0:255)
      INTEGER*2    DATE(3)  !date (year,month,day)
      REAL         ACLAT    !aircraft latitude in degrees north
      REAL         ACLON    !aircraft longitude in degrees
      LOGICAL      HIRES    !high resolution, boolean
      REAL         RANGE    !slant range in meters
      INTEGER*2    ACHDG    !aircraft heading in bams
      BYTE         RMODE    !radar mode (0:15)
      INTEGER*2    ANTAZ    !antenna pointing angle in bams
      BYTE         DQUAL    !data quality (0:7)
      LOGICAL      ASVAL    !aspect valid; boolean
      INTEGER*2    TAGNO    !tws target index (0:511)
      BYTE         SSTAT    !sea state (0:7)
      BYTE         FPATH    !flight path (,fcircle,hcircle,flyby,flyout,flyin,,)

      REAL         NSVEL    !north/south velocity in m/s north
      LOGICAL      VERTP    !vertical polarization; boolean
      REAL         EWVEL    !east/west velocity in m/s east
      INTEGER*2    ASANG    !aspect angle in degrees, 511 - unknown
      INTEGER*2    ATILT    !antenna tilt
      BYTE         APRTR    !aperture
      INTEGER*2    ASERR    !aspect error in degrees
      BYTE         ECODE    !event code (0:7)
      INTEGER*2    TIME(4)  !time (hrs,mins,sec,.001sec)
      LOGICAL      RNGAV    !reduced resolution (range averaging);boolean
      LOGICAL      ITEST    !imaging on test pulse; boolean
      INTEGER*2    PRFNO    !pulse repetition frequency
      BYTE         PULSE    !pulse characteristics
      INTEGER*4    SWATH    !swath counter 
      LOGICAL      TLOCK    !target locked; boolean
      REAL         PHASE    !phase in degrees
      INTEGER*2    TKCEL    !track cell number (0:512)
      REAL         AGCGN    !agc gain in decibels
      REAL         MANGN    !manual gain in decibels
C
      COMMON / HEADER / DTYPE , ACALT , SPARR , DATE , ACLAT , ACLON ,
     .                  HIRES , RANGE , ACHDG , RMODE , ANTAZ , DQUAL ,
     .                  ASVAL , TAGNO , SSTAT , FPATH , NSVEL , VERTP ,
     .                  EWVEL , ASANG , ATILT , APRTR , ASERR , ECODE , 
     .                  TIME , RNGAV , ITEST , PULSE , PRFNO , SWATH , 
     .                  TLOCK , PHASE , TKCEL , AGCGN , MANGN  
C
C*---------------------------------------------------------------------
