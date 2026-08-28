C**
C***********************************************************************
C**
      subroutine adjust ( isub , dr0 , vr0 , ntp , ncr , nsr , nspf ,
     .                    nakeep , curtime , dts , r0min , r0max ,
     .                    r0bar , rfmmin , rfmmax , cstr , csbimg ,
     .                    omega0 , vslip , vrbar , vrsave , pt , ntg ,
     .                    isub0 , afocus , nabuff , tgrbar , tgrmin ,
     .                    tgrmax , nr1 , nre , nr1fm , nrefm )
C**
C***********************************************************************
C**
c
c   Report information to the radar control system so that it can
c   modify the way the data is being collected and so it can adjust
c   the range center of the image.
c
c   Also, this routine makes the initial acquisition mode decisions on
c   the scene center and motion compensation velocity and reads in
c   information for some special research modes.
c
      implicit none
c
c-----------------------------------------------------------------------
c
c   Use the global target list and standard parameters stored in common
c
      include 'tglist.h'
c
      include 'sarprm.h'
c
      include 'updates.h'
c
c-----------------------------------------------------------------------
c
      integer isub , ntp , ncr , nsr , nspf , nakeep , nshift , isbold ,
     .        ntg , isub0 , nabuff , nr1 , nre , nr1fm , nrefm , ncf ,
     .        ncfm , noff
c
      complex cstr(ntp,ncr)               !  Range-compressed signal
c                                            history
c
      complex csbimg(nakeep,nsr,nabuff)   !  3-D buffer of sub-images
c
      real    dr0 , vr0 , curtime , dts , r0min , r0max ,
     .        r0bar , rfmmin , rfmmax , rcent , omega0 , drfmin ,
     .        drfmax , rwidth , dvr0dt , vslip , eslip , vadd , vrbar ,
     .        vrsave , pt , afocus , tgrbar , tgrmin , tgrmax 
c
      logical no_rnotch
c
c-----------------------------------------------------------------------
c
c   The integer flag 'alias' is used to control the decisions of the
c   acquisition mode:
c
c   alias = 2  : A special mode for the NRL Advanced Profile Radar.  In
c                this mode there are two mo-comp velocities, one for the
c                phase (vr0) and vslip, a correction to this value which
c                moves the radar cross section without affecting the
c                phase.
c
c   alias = 1  : Normal mode - let the program decide where to put the
c                mo-comp point even to the point of deciding which alias
c                velocity is correct.
c
c   alias = 0  : Form the image at the center of the range swath of the
c                raw input data and ignore the information from the
c                target detector.
c
c--------------- Manual Mo-Comp Modes ----------------------------------
c
c   alias = -1 : Program prompts for the initial mo-comp center and
c                velocity and for the acceleration.  The velocity then
c                increases linearly at this rate; the range therefore
c                varies quadratically with time.
c
c   alias = -2 : Program prompts for the initial mo-comp range and
c                velocity and then prompts every sub-image for a new
c                velocity.
c
c   alias = -3 : This special mode is designed for the NRL Advanced
c                profile radar, which has range errors which are not
c                consistent with the phase information.  Every sub-image
c                the program reads a value of the 'slip' velocity,
c                vslip.  This value is added to the velocity of the
c                mo-comp point but is not used to correct the doppler
c                frequency.
c
c   alias = -4 : Same as alias = -2, except the program also prompts
c                for a rotation rate of the object or scene.  This value
c                is then used to compute the acceleration field, over-
c                riding the information normally calculated by GETACC.
c                This gives a result approximately equivalent to that
c                of the polar format algorithm applied to a uniformly
c                rotating scene.
c
c***********************************************************************
c
         no_rnotch = .not. ( ( notch .eq. 2 ) .or. ( notch .ge. 3 ) )
c
c-----------------------------------------------------------------------
c
c   If there are user-specified range or velocity notches, use them
c   to set the initial mo-comp parameters
c
         if ( isub .eq. isub0 ) then
c
            tgrbar = 0.0
            tgrmin = 0.0
            tgrmax = 0.0
c
            if ( notch .eq. 2 .or. notch .ge. 3 ) then
c
               dr0    = rnotch
               rfmmin = - 0.5 * drntch
               rfmmax = + 0.5 * drntch
c
            endif
c
            if ( notch .eq. 1 .or. notch .ge. 3 ) then
c
               vr0    = vnotch
               vrsave = vr0
c
            endif
c
         endif
c
c-----------------------------------------------------------------------
c   Set the windows for target search and fine mo-comp
c
         if ( notch .gt. 3 .or.
     .        ( no_rnotch .and.
     .        ( ( isub .eq. isub0 ) .or.
     .          ( isub .le. 0 .and. ntg .eq. 0 ) ) ) ) then
c
c   If in wide search mode, search for targets over all range cells
c   (except for buffer of 4 cells at each end)
c
            nr1    = 1   + 4
            nre    = ncr - 4
            nr1fm  = nr1
            nrefm  = nre
c
            tgrbar = 0.0
c
            rfmmin = rfmin
            rfmmax = rfmax
c
         else
c
c   If in narrow search mode, search over those range cells which fall
c   within the fine resolution image
c
            ncf    = int( float( nsr ) * ( ( rfmax - rfmin ) /
     .                                     ( rsmax - rsmin ) ) )
            nr1    = max( 1 + 4 , 1 + ( ncr - ncf ) / 2 )
            nre    = min( ncr - 4 , nr1 + ncf - 1 )
c
c   Apply the fine motion compensation algorithm over the range cells
c   where targets have been identified.  Allow four range resolution
c   cells on each side and limit the region to between 8 cells and the
c   size of the fine image.
c
            ncfm   = 8 + ifix( ( rfmmax - rfmmin ) / drc )
            noff   = ifix( ( rfmmax + rfmmin ) / ( 2.0 * drc ) )
            nr1fm  = max( nr1 , noff + 1 + ( ncr - ncfm ) / 2 )
            nrefm  = min( nre , nr1fm + ncfm - 1 )
c
         endif
c
         if ( quiet .gt. 1 )
     .   write ( 6 , * ) ' nr1, nre, nr1fm, nrefm' ,
     .                     nr1 , nre , nr1fm , nrefm ,
     .                     r0min , r0max , rfmmin , rfmmax
         write ( 7 , * ) ' nr1, nre, nr1fm, nrefm' ,
     .                     nr1 , nre , nr1fm , nrefm ,
     .                     r0min , r0max , rfmmin , rfmmax
c
c***********************************************************************
c**************  Acquisition and Mo-Comp range limits  *****************
c
c   This section is done every time the target detector is called
c
      if ( isub .le. 2 .or. mod( isub , 2 ) .eq. 0 ) then
c
c   If no targets were detected, set range limits to the center of the
c   desired range notch (notch=2,3) or to a box the size of the fine
c   resolution image centered at zero
c
         if ( ntg .eq. 0 ) then
c
            if ( notch .eq. 2 .or. notch .eq. 3 ) then
c
               r0min = - 0.5 * drntch
               r0max = + 0.5 * drntch
c
            else
c
               r0min = rfmin
               r0max = rfmax
c
            endif
c
         else
c
            r0min = tgrmin
            r0max = tgrmax
c
         endif
c
      endif
c
c-----------------------------------------------------------------------
c
c   Update the mo-comp parameters to reflect the estimate of the
c   position and velocity of the dominant scatterers
c
      if ( alias .gt. 0 .and. isub .le. 1  ) then    !  Acquisition mode
c
         write ( 6  , * )
         write ( 6  , * )
     .    ' Acquisition Mode: Setting mean Doppler and Mocomp position'
         write ( 7 , * )
         write ( 7 , * )
     .    ' Acquisition Mode: Setting mean Doppler and Mocomp position'
c
c   Mo-comp range
c
         r0bar  = tgrbar
c
         r0min  = r0min - r0bar 
         r0max  = r0max - r0bar
c
         dr0    = dr0 + r0bar - 4.0 * ( vr0 + vrbar ) * dts
c
         rfmmin = r0min - 2.0 * drc   !  Min. range for fine mo-comp
         rfmmax = r0max + 2.0 * drc   !  Max. range for fine mo-comp
c
         vr0    = vrbar + vrsave      !  Update from last velocity
         vrsave = vr0                 !  Store last velocity
         pt     = 0.0                 !  Integrated phase
c
      else                            !  Narrow search mode
c
c   Adjust the region used for fine Mo-Comp to the detected targets
c
         rwidth = ( rfmmax - rfmmin ) + 8.0 * drc
c
         drfmin = 0.125 * ( r0min - rfmmin )
         drfmin = amin1( drfmin , + 0.125 * rwidth )
         drfmin = amax1( drfmin , - 0.125 * rwidth )
c
         drfmax = 0.125 * ( r0max - rfmmax )
         drfmax = amin1( drfmax , + 0.125 * rwidth )
         drfmax = amax1( drfmax , - 0.125 * rwidth )
c
         rfmmin = rfmmin + drfmin
         rfmmax = rfmmax + drfmax
c
      endif
c
c**********  END of Acquisition and Mo-Comp range limits  **************
c***********************************************************************
c
c------------------------  Manual Mo-comp  -----------------------------
c
      if ( alias .lt. 0 ) then            ! Manual mode
c
         if ( isub .eq. 0 ) then
c 
            if ( alias .eq. - 4 .or. alias.eq. - 8 ) then
c
               write ( 6 , * )
     .         ' Manual Polar Format Option - ENTER Omega (Rad/Sec):'
               read ( 5 , * ) omega0
c
            endif
c
            if ( alias .eq. - 6 ) then
c
               dr0 = 0.0
c
               if ( finemc .eq. 0 ) vr0 = 0.0
c
               write ( 6 , * )
     .         ' Manual Focus Option - ENTER Afocus (Hz/Sec):'
c
               read ( 5 , * ) afocus
c
            else if ( alias .eq. - 4 ) then
c
               write ( 6 , * )
     .         ' Manual Acquisition Mode - ENTER DR0, VR0:'
c
               read ( 5 , * ) dr0 , vr0
c
            endif
c
            if ( finemc .eq. 0 .and. alias .eq. - 1 ) then
c
               write ( 6 , * ) ' ENTER DVR0DT:'
c
               read ( 5 , * ) dvr0dt
c
               isbold = 0
c
            endif
c
            r0min  = rfmin
            r0max  = rfmax
            rfmmin = r0min
            rfmmax = r0max
c
         endif
c
         if ( isub .gt. 0 .and. finemc .eq. 0 .and.
     .                           alias .eq. - 1 ) then
c
            vr0    = vr0 + dvr0dt * dts * float( isub - isbold )
            isbold = isub
c
         endif
c
         if ( isub .gt. 0 .and. finemc .eq. 0 .and.
     .        ( alias .eq. - 2 .or. alias .eq. - 8 ) ) then
c
            read ( 5 , * ) vr0
c
         endif
c
         if ( isub .gt. 0 .and. alias .eq. - 3 ) then
c
            read ( 5 , * ) vslip
c
         endif
c
      endif
c
c----------------------  END of Manual Mo-comp  ------------------------
c
      if ( quiet .gt. 1 )
     .write ( 6 , '(/,1x,a27,5f10.4,/)' ) ' t, r0bar, vrbar, dr0, vr0:'
     .                          , curtime , r0bar , vrbar , dr0 , vr0
      write ( 7 , '(/,1x,a27,5f10.4,/)' ) ' t, r0bar, vrbar, dr0, vr0:'
     .                          , curtime , r0bar , vrbar , dr0 , vr0
c
c   If the image is too far off-center, shift the motion compensation
c   point and then rewrite history so that all data is consistent with
c   this value
c
      rcent = 0.5 * ( rfmmin + rfmmax )
c              
      if ( alias .gt. 0 .and. isub .gt. 0 .and.
     .     abs( rcent ) .gt. float( nsr / 4 ) * drc ) then
c
c   Shift the range by an integer number of coarse resolution range 
c   cells
c
         nshift = nint( rcent / drc )
c
         call moverc ( isub , drc , ntp , ncr , nsr , nabuff , 
     .                 nakeep , dr0 , r0min , r0max , r0bar ,
     .                 rfmmin , rfmmax , cstr , csbimg , nshift )
c
      endif
c
      return
      end
C**
C***********************************************************************
C**
      subroutine moverc ( isub , drc , ntp , ncr , nsr , nabuff , 
     .                    nakeep , dr0 , r0min , r0max , r0bar ,
     .                    rfmmin , rfmmax , cstr , csbimg , nshift )
C**
C***********************************************************************
C**
c
c   Purpose:  To shift the range center by an integer number of coarse
c             resolution range cells.  This requires redefining the
c             variables which describe the range and it requires the
c             rewriting of history so that all arrays of data and
c             tables are consistent with this new center.
c
      implicit none
c
c-----------------------------------------------------------------------
c
c   Use the global target list stored in common
c
      include 'tglist.h'
c
c-----------------------------------------------------------------------
c
      integer isub , ntp , ncr , nsr , nabuff , nakeep , nshift ,
     .        nsubs , isb , itr , isr , ifa , iloop , ir , mshift , itgt
c
      complex cstr(ntp,ncr)               !  Range-compressed signal
c                                            history
c
      complex csbimg(nakeep,nsr,nabuff)   !  3-D buffer of sub-images
c
      real    dr0 , drc , r0min , r0max , r0bar , rfmmin , rfmmax ,
     .        rshift
c
c-----------------------------------------------------------------------
c
c   Shift the range an integral number of coarse range cells
c
      rshift = drc * float( nshift )
c
      write ( 6 , '(/,a24,i6,a6,f8.2,a7)' )
     .      ' Range origin shifted by' ,
     .        nshift , ' cells' , rshift , ' meters'
      write ( 7 , '(/,a24,i6,a6,f8.2,a7)' )
     .      ' Range origin shifted by' ,
     .        nshift , ' cells' , rshift , ' meters'
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Add this shift to the mo-comp range - the range relative to the
c   raw data stream, and subtract it from the local estimates of the
c   object relative to this point.
c
      dr0    = dr0    + rshift
c
      r0bar  = r0bar  - rshift
      r0min  = r0min  - rshift
      r0max  = r0max  - rshift
      rfmmin = rfmmin - rshift
      rfmmax = rfmmax - rshift
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c   Re-write history: Modify the sub-image array, the range-compressed
c   signal history array and the target list as if the new center had
c   always been used.
c
c-----------------------------------------------------------------------
c
c   Adjust the sub-image array and range-compressed signal history to
c   the new center
c
      nsubs = min( isub , nabuff )
c
      if ( nshift .gt. 0 ) then
c
c   Do all the sub-images computed so far
c
         do isb = 1 , nsubs
c
c   Move sub-image data down
c
            do isr = 1 , nsr - nshift
c
               do ifa = 1 , nakeep
c
                  csbimg(ifa,isr,isb) = csbimg(ifa,isr+nshift,isb)
c
               enddo
c
            enddo
c
c   Fill sub-image top with zeroes
c
            do isr = nsr - nshift + 1 , nsr
c
               do ifa = 1 , nakeep
c
                  csbimg(ifa,isr,isb) = cmplx( 0.0 , 0.0 )
c
               enddo
c
            enddo
c
         enddo
c
c   Move signal history down
c
         do ir = 1 , ncr - nshift
c
            do itr = 1 , ntp
c
               cstr(itr,ir) = cstr(itr,ir+nshift)
c
            enddo
c
         enddo
c
c   Fill signal history top with zeroes
c
         do ir = ncr - nshift + 1 , ncr
c
            do itr = 1 , ntp
c
               cstr(itr,ir) = cmplx( 0.0 , 0.0 )
c
            enddo
c
         enddo
c
      else
c
         mshift = - nshift
c
         do isb = 1 , nsubs
c
c   Move sub-image data up
c
            do iloop = 1 , nsr - mshift
c
               isr = nsr + 1 - iloop
c
               do ifa = 1 , nakeep
c
                  csbimg(ifa,isr,isb) = csbimg(ifa,isr-mshift,isb)
c
               enddo
c
            enddo
c
c   Fill sub-image bottom  with zeroes
c
            do isr = 1 , mshift
c
               do ifa = 1 , nakeep
c
                  csbimg(ifa,isr,isb) = cmplx( 0.0 , 0.0 )
c
               enddo
c
            enddo
c
         enddo
c
c   Move signal history up
c
         do iloop = 1 , ncr - mshift
c
            ir = ncr + 1 - iloop
c
            do itr = 1 , ntp
c
               cstr(itr,ir) = cstr(itr,ir-mshift)
c
            enddo
c
         enddo
c
c   Fill signal history bottom with zeroes
c
         do ir = 1 , mshift
c
            do  itr = 1 , ntp
c
               cstr(itr,ir) = cmplx( 0.0 , 0.0 )
c
            enddo
c
         enddo
c
      endif
c
c-----------------------------------------------------------------------
c
c   Adjust the range of all targets presently in the target list
c
      do itgt = 1 , nlist
c
         if ( time(itgt) .ge. 0.0 ) range(itgt) = range(itgt) - rshift
c
      enddo
c
c-----------------------------------------------------------------------
c
      return
      end
