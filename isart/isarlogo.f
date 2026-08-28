C**
C***********************************************************************
C**
      subroutine isarlogo ( lun )
C**
C***********************************************************************
C**
      implicit none
c
      integer  lun
c
c   Print a listing of those responsible for the development of ISART
c
      write ( lun , 5 )
    5 format ( //
     .  '                             Program ISAR-T'//,
     .  '                     An Inverse-SAR Movie Processor'//,
     .  '                             VERSION 3.95'/,
     .  '                           January 17, 2018'//,
     .  '                          COPYRIGHT 1994-2017'//,
     .  '                            John R. Bennett'/,
     .  '                         Escondido, California'//,
     .  '                               AUTHORS'/,
     .  '                           John R. Bennett'/ ,
     .  '                         Kenneth A. Melendez'/ ,
     .  '                           David S. Brown'/ ,
     .  '                        Christopher J. Campo'/ ,
     .  '                          Barton P. Schade'/ ,
     .  '                          Julie M. Jauregui'// )
c
      return
      end
