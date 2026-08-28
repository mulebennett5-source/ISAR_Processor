c
c                     File units for Rmovie
c
      integer         asunit , rmmunit , rmmunit2 , isartunit ,
     .                movieunit , rmdunit1 , rmdunit2 , childunit ,
     .                iprunit , pltunit
c
      parameter     ( asunit    = 11 )  !  A-Scan data
      parameter     ( rmmunit   = 12 )  !  Movie replay data
      parameter     ( rmmunit2  = 13 )  !  Secondary movie replay data
      parameter     ( isartunit = 14 )  !  Log file
      parameter     ( movieunit = 15 )  !  Replay instructions
      parameter     ( rmdunit1  = 16 )  !  Primary raw data file
      parameter     ( rmdunit2  = 17 )  !  Secondary raw data file
      parameter     ( childunit = 18 )  !  Child window I/O
      parameter     ( iprunit   = 19 )  !  Impulse Response file
      parameter     ( pltunit   = 20 )  !  IPR plot file
c      
      