C**
C***********************************************************************
C**
      subroutine image3 ( ifile , ofile , nx , ny , cin , csbimg ,
     .                    work , bytes , mprat , foption , clutterc )
C**
C***********************************************************************
C**
      implicit none
c
      character     ifile*80 , ofile*80
c
      integer       nx , ny , foption , reclda , i , j , it , mprat ,
     .              clutterc , mframes , mmovie , kmovie
c
      complex       cin(nx,ny)
c
      complex       csbimg(1+nx/mprat,ny,2*mprat)
c
      character     bytes(nx,ny)*1
c
      real          work(nx,ny) , pi , inmin , inmax , inbar , insdv ,
     .              thresh , topsg0 , tgtsg0 , oneside 
c
      complex       cw(1024)
c
      real          rw(2,1024)
c
      equivalence ( cw , rw )
c
      logical       double
c
c-----------------------------------------------------------------------
c
      if ( foption .eq. 1 ) then
c
         open ( 8 , file = ifile , form = 'UNFORMATTED' ,
     .              access = 'DIRECT' , recl = reclda( 0 , 2 * nx ) )
c
         do j = 1 , ny
c
            read ( 8 , rec = j ) ( cin(i,j) , i = 1 , nx )
c
         enddo
c
         close ( 8 )
c
      endif
c
c-----------------------------------------------------------------------
c
      mmovie  = 2 * mprat - 1   !  Size of a minimovie
      mframes = mmovie          !  At least one mini-movie
c
c   Check how many optional mini-movies were selected
c
c     Bit 0  :  Spatially Flattened image
c
c     Bit 1  :  Temporally Flattened image
c
c     Bit 2  :  Time-Cancelled image
c
      if ( mod( clutterc , 2 )     .eq. 1 ) mframes = mframes + mmovie
      if ( mod( clutterc / 2 , 2 ) .eq. 1 ) mframes = mframes + mmovie
      if ( mod( clutterc / 4 , 2 ) .eq. 1 ) mframes = mframes + mmovie
c
      call make_rmg ( 1 + nx / mprat , ny , mframes , 1 , ofile )
c
      write ( 6 , * ) ' Movie X-Pixels: ' , 1 + nx / mprat
      write ( 6 , * ) ' Movie Y-Pixels: ' , ny
      write ( 6 , * ) ' Movie Frames:   ' , mframes
c
      pi = atan2( 0.0 , - 1.0 )
c
c   Get the clutter cancellation parameters
c
      call vparms ( thresh , topsg0 , tgtsg0 , oneside )
c
c   Convert original complex image from frequency-range to
c   frequency-range-time
c
      double = .false.
c
      call imgeni ( cin , csbimg , nx , ny , mprat , cw , rw ,
     .              work , ny , double )
c
c   Default movie - no clutter reduction
c
      kmovie = 1
c
      call oframe ( csbimg(1,1,2) , work , bytes , - kmovie , ofile ,
     .              1 + nx / mprat , ny * ( 2 * mprat - 1 ) , 0 , 0 ,
     .              0.25 , 0 , 0 , 0 )
c
c   Flatten the cross section field so that bright targets don't spill
c   over and cause false alarms when blurred
c
      if ( mod( clutterc , 2 ) .eq. 1 ) then
c
         call flatten ( cin , work , nx , ny , topsg0 ,
     .                  inmin , inmax , inbar , insdv )
c
         write ( 6 , '(/,a48)' )
     .      ' Stats for image intensity (min, max, avg, sdv)'
         write ( 6 , '(4e12.4)' ) inmin , inmax , inbar , insdv
         write ( 6 , * )
c
         call imgeni ( cin , csbimg , nx , ny , mprat , cw , rw ,
     .                 work , ny , double )
c
c   Spatially flattened clutter movie
c
         kmovie = kmovie + 1
c
         call oframe ( csbimg(1,1,2) , work , bytes , - kmovie , ofile ,
     .                 1 + nx / mprat , ny * ( 2 * mprat - 1 ) , 0 , 0 ,
     .                 0.2 , 0 , 0 , 0 )
      endif
c
c
c   Flatten in time - removes effects of glints from flat walls, etc.
c
      if ( mod( clutterc / 2 , 2 ) .eq. 1 ) then
c
         call tflatn ( csbimg , nx , ny , mprat , 32 , topsg0 , work ,
     .                 1.0 )
c
c   Temporally flattened clutter movie
c
         kmovie = kmovie + 1
c
         call oframe ( csbimg(1,1,2) , work , bytes , - kmovie , ofile ,
     .                 1 + nx / mprat , ny * ( 2 * mprat - 1 ) , 0 , 0 ,
     .                 0.2 , 0 , 0 , 0 )
c
      endif
c
c   Clutter cancel in time - removes stationary targets
c
      if ( mod( clutterc / 4 , 2 ) .eq. 1 ) then
c
         do it = 1 , mprat
c
            do i = 1 , 1 + nx / mprat
c
               do j = 1 , ny
c
                  csbimg(i,j,it+1)         = 0.5 * ( csbimg(i,j,it+1) -
     .                                       csbimg(i,j,2*mprat-it+1) )
c
                  csbimg(i,j,2*mprat-it+1) = csbimg(i,j,it+1)
c
               enddo
c
            enddo
c
         enddo
c
         kmovie = kmovie + 1
c
         call oframe ( csbimg(1,1,2) , work , bytes , - kmovie , ofile ,
     .                 1 + nx / mprat , ny * ( 2 * mprat - 1 ) , 0 , 0 ,
     .                 0.2 , 0 , 0 , 0 )
c
      endif
c
      return
      end
c**
c***********************************************************************
c**
      subroutine make_rmg ( nx , ny , mframes , nskip , ofile )
c**
c***********************************************************************
c**
      implicit none
c
c   Create a simple text file which allows RMOVIE to automatically make
c   a movie which does not come from the SAIC ISAR processor
c
      integer      nx , ny , mframes , nskip , ndot , lastdot
c
      character    ofile*80 , ofile_g*80
c
c-----------------------------------------------------------------------
c
      ndot    = lastdot( ofile ) - 1
c
      ofile_g = ofile(1:ndot) // '.rmg'
c
      write ( 7 , * ) ofile(1:40)
      write ( 7 , * ) ofile_g(1:40)
      write ( 7 , * ) ndot
c
      open ( 9 , file = ofile_g , form = 'formatted' ,
     .           status = 'unknown' )
c
      write ( 9 , '(i10,a32)' ) nx      , ' Movie X-Pixels'
      write ( 9 , '(i10,a32)' ) ny      , ' Movie Y-Pixels'
      write ( 9 , '(i10,a32)' ) mframes , ' Movie Frames'
      write ( 9 , '(i10,a32)' ) nskip   , ' Skip'
c
      return
      end
