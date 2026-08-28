C**
C***********************************************************************
C**
      subroutine tcolor ( bi , nx , ny )
C**
C***********************************************************************
C**
c   This routine writes auxiliary information into an 8-bit image using
c   higher order bits for color-coding pixels
c
      implicit none
c
      integer   nx , ny , iloop , ip , jp , k , kk , kdots , goodtg ,
     .          badtg , box , kbox , box_c
c
      character bi(nx,ny)*1
c
      real      f0 , r0 , f1 , f2 , r1 , r2 , fp , rp , kdots1 ,
     .          range_t0 , freq_t0
c
      include  'sarprm.h'
c
      include  'updates.h'
c
      include  'tglist.h'
c
c-----------------------------------------------------------------------
c
c   Remove the two lowest order bits so that they can be used for color
c   coding of killed targets, good targets and the box outline.
c
      do jp = 1 , ny
c
         do ip = 1 , nx
c
            bi(ip,jp) = char( 4 * ( ichar( bi(ip,jp) ) / 4 ) )
c
         enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
      badtg = 3  !  Flag value for killed targets
c
      if ( otype .eq. 2 .and. curve .ne. 0 ) then
c
c   Lower left corner range and frequency
c
         f0 = fcenuse - float( nx / 2 ) * dff
c
         r0 = rcenuse - float( ny / 2 ) * drf
c
c-----------------------------------------------------------------------
c
c   Make accepted targets green or blue
c
         if ( omega_valid .eq. 2 ) then
c
            goodtg = 1
            box    = 2
c
         else
c
            goodtg = 2
            box    = 1
c
         endif
c
         do iloop = 1 , nlist
c
            if ( iflag(iloop) .eq. editac ) then
c
c   Correct range and frequency to center of aperture time
c
               range_t0      =   range(iloop) + 0.5 * lambda *
     .                           freq(iloop) * ( tbar - time(iloop) )
c
               freq_t0       =   freq(iloop)  + accel(iloop) *
     .                                          ( tbar - time(iloop) )
c
               ipixtg(iloop) = 1 + nint( ( freq_t0  - f0 ) / dff )
               jpixtg(iloop) = 1 + nint( ( range_t0 - r0 ) / drf )
c
               ip            = ipixtg(iloop)
               jp            = jpixtg(iloop)
c
               ip            = max( 1 , min( ip , nx ) )
               jp            = max( 1 , min( jp , ny ) )
c
               if ( mod( ichar( bi(ip,jp) ) , 4 ) .eq. 0 )
     .            bi(ip,jp) = char( ichar( bi(ip,jp) ) + goodtg )
c
            endif
c
         enddo
c
c-----------------------------------------------------------------------
c
c   Make rejected targets red
c
         do iloop = 1 , nlist
c
            if ( iflag(iloop) .eq. 101 ) then
c
               range_t0      =   range(iloop) + 0.5 * lambda *
     .                           freq(iloop) * ( tbar - time(iloop) )
c
               freq_t0       =   freq(iloop)  + accel(iloop) *
     .                                          ( tbar - time(iloop) )
c
               ipixtg(iloop) = 1 + nint( ( freq_t0  - f0 ) / dff )
               jpixtg(iloop) = 1 + nint( ( range_t0 - r0 ) / drf )
c
               ip            = ipixtg(iloop)
               jp            = jpixtg(iloop)
c
               ip            = max( 1 , min( ip , nx ) )
               jp            = max( 1 , min( jp , ny ) )
c
               if ( mod( ichar( bi(ip,jp) ) , 4 ) .eq. 0 )
     .            bi(ip,jp) = char( ichar( bi(ip,jp) ) + badtg )
c
            endif
c
         enddo
c
c-----------------------------------------------------------------------
c
c   Make box outline green (if omega_valid) blue (if not omega_valid)
c
c   Also, if color < 0 make the outline of the a priori box using the
c   bad target color
c
         do kbox = 1 , 2
c
            if ( kbox .eq. 1 ) then
c
               box_c = box
c
            else
c
               if ( color .lt. 0 ) then
c
                  box_c = badtg
c
               else
c
                  box_c = 0
c
               endif
c
            endif
c
            do k = 1 , 4
c
               ip = 1 + nint( ( corner(k,1,kbox) - f0 ) / dff )
               jp = 1 + nint( ( corner(k,2,kbox) - r0 ) / drf )
c
               ip = max( 2 , min( ip , nx - 1 ) )
               jp = max( 2 , min( jp , ny - 1 ) )
c
               if ( mod( ichar( bi(ip,jp) ) , 4 ) .eq. 0 )
     .              bi(ip,jp) = char( ichar( bi(ip,jp) ) + box_c )
c
               if ( mod( ichar( bi(ip+1,jp) ) , 4 ) .eq. 0 )
     .              bi(ip+1,jp) = char( ichar( bi(ip+1,jp) ) + box_c )
c
               if ( mod( ichar( bi(ip,jp+1) ) , 4 ) .eq. 0 )
     .              bi(ip,jp+1) = char( ichar( bi(ip,jp+1) ) + box_c )
c
               if ( mod( ichar( bi(ip-1,jp) ) , 4 ) .eq. 0 )
     .              bi(ip-1,jp) = char( ichar( bi(ip-1,jp) ) + box_c )
c
               if ( mod( ichar( bi(ip,jp-1) ) , 4 ) .eq. 0 )
     .              bi(ip,jp-1) = char( ichar( bi(ip,jp-1) ) + box_c )
c
c   Plot intermediate points
c
               f1 = corner(k,1,kbox) - f0
               r1 = corner(k,2,kbox) - r0
c
               if ( k .eq. 1 ) then
c
                  f2 = corner(2,1,kbox) - f0
                  r2 = corner(2,2,kbox) - r0
c
               else if ( k .eq. 2 ) then
c
                  f2 = corner(4,1,kbox) - f0
                  r2 = corner(4,2,kbox) - r0
c
               else if ( k .eq. 3 ) then
c
                  f2 = corner(1,1,kbox) - f0
                  r2 = corner(1,2,kbox) - r0
c
               else if ( k .eq. 4 ) then
c
                  f2 = corner(3,1,kbox) - f0
                  r2 = corner(3,2,kbox) - r0
c
               endif
c
               kdots = ifix( ( sqrt( ( abs( f2 - f1 ) / dff ) ** 2 +
     .                               ( abs( r2 - r1 ) / drf ) ** 2 ) )
     .                         / 2.0 )
c
               kdots1 = float( kdots + 1 )
c
               do kk = 1 , kdots
c
                  fp = f1 + ( float( kk ) / kdots1 ) * ( f2 - f1 )
                  rp = r1 + ( float( kk ) / kdots1 ) * ( r2 - r1 )
c
                  ip = 1 + nint( fp / dff )
                  jp = 1 + nint( rp / drf )
c
                  ip = max( 1 , min( ip , nx ) )
                  jp = max( 1 , min( jp , ny ) )
c
                  if ( mod( ichar( bi(ip,jp) ) , 4 ) .eq. 0 )
     .                 bi(ip,jp) = char( ichar( bi(ip,jp) ) + box_c )
c
               enddo
c
            enddo
c
         enddo
c
      endif
c
      return
      end
