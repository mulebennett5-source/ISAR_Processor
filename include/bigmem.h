c
c***********************************************************************
c
c                       File: 'bigmem.h'
c
c***********************************************************************
c
c
c   Define all array memory by pointing into a single large array
c
c   Memory sizes:
c
c         nbig   :  The number of words of memory used for real,
c                   complex or integer arrays.
c
c         nbigc  :  Number of complex elements, = nbig / 2
c
c         nbigb  :  Number of bytes used for 8-bit variables
c
c         nbigi2 :  Number of integer*2 variables. = nbigb / 2
c
      integer       nbig , nbigc , nbigi2 , nbigb
c
c   Fundamental memory size (words)
c
      parameter   ( nbig   = 120000 * 1024 )
c
c   Consider all real, complex, and integer memory as one array
c
      parameter   ( nbigc  = nbig / 2     )
c
      real          big(nbig)
c
      complex       bigc(nbigc)
c
      integer       bigi(nbig)
c
      equivalence ( big , bigc ) , ( big , bigi )
c
c   Consider all character*1 and integer*2 memory as one array
c
      parameter   ( nbigb  = nbig         )
c
      parameter   ( nbigi2 = nbigb / 2    )
c
      character     bigb(nbigb)*1
c
      byte          big8(nbigb)
c
      integer*2     bigi2(nbigi2)
c
      equivalence ( bigb , bigi2 ) , ( bigb , big8 )
c
c   Put both arrays in a common block so that all routines can use it
c
      common / bigmem / big , bigb
c
c***********************************************************************
c