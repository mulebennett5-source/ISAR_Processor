C**
C***********************************************************************
C**
      subroutine fixem ( sariq , nr , np , ib , qb , iib , qqb , iqb ,
     .                   nedit , npdone )
C**
C***********************************************************************
C**
      implicit none
c
c   The routine performs the following tasks:
c
c        A. Editing out anomolous points apparently due to
c           interference from another radar.
c
c        B. Removing the mean i and q values.
c
c        C. Rotating the i and q values to remove correlation between
c           them.
c
c        D. Rescaling i and q so that they have the same variance.
c
c   The routine uses estimates of the statistics averaged over all
c   pulses and smoothed over range times using a weight function
c   with nwts values
c
      integer     nwts
      parameter ( nwts = 129 )
c
      real        wts(nwts)
c 
c   The input data is a complex array with the following dimensions:
c
c         nr     :  Fast-time samples
c
c         np     :  Total pulses
c
c         npdone :  Number of pulses already corrected
c
      integer     nr , np , npdone , nedit
c
      real        sariq(2,nr,np)
c
c   Local variables
c
      integer     irange , ipulse , ipass , npass1 , npass2 , npass3 ,
     .            npass4 , iwt 
c
      real        ibar , qbar , iibar , qqbar , iqbar , ib(nr) ,
     .            qb(nr) , iib(nr) , qqb(nr) , iqb(nr) , pi , theta ,
     .            ctheta , eps , stheta , scalei , scaleq , scale ,
     .            sariq1 , icor , qcor , iicor , qqcor , iqcor , awt ,
     .            smth , thresh , ibar_s , qbar_s , iibar_s , qqbar_s ,
     .            iqbar_s
c
      save        pi ,  wts , ibar_s , qbar_s , iibar_s , qqbar_s ,
     .            iqbar_s , thresh
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      pi     = atan2( 0.0 , - 1.0 )
c
      if ( npdone .eq. 0 ) then
c
         nedit = 0
c
c   Calculate initial statistics of raw data
c
         call iqstat ( sariq(1,1,npdone+1) , nr , np-npdone , 
     .                 ib , qb , iib , qqb , iqb ,
     .                 ibar , qbar , iibar , qqbar , iqbar )
c
         ibar_s  = ibar
         qbar_s  = qbar
         iibar_s = iibar
         qqbar_s = qqbar
         iqbar_s = iqbar
c
         thresh = 3.8 * sqrt( 0.5 *
     .            ( iibar - ibar ** 2 + qqbar - qbar ** 2 ) )
c
         write ( 7 , '(///,a47,/)' )
     .        ' Statistics at stages in correction of raw data'
         write ( 7 , '(a32,2e14.6)' )
     .                           ' Mean of i, q for raw data    :  ' ,
     .                             ibar , qbar
         write ( 7 , '(a32,3e14.6)' )
     .                           ' Mean of ii, qq, iq for raw data:' ,
     .                             iibar , qqbar , iqbar
c
c   Generate weights for smoothing calibration statistics
c
         do iwt = 1 , nwts
c
            awt      = float( iwt - nwts / 2 ) / 8.0
            wts(iwt) = exp( - pi * awt ** 2 )
c
         enddo
c
      endif
c
c   Numbers of passes for each correction stage
c
      if ( npdone .eq. 0 ) then
c
         npass1 = 3
c
      else
c
         npass1 = 1
c
      endif
c
      npass2 = 1
      npass3 = 1
      npass4 = 1
c
c***********************************************************************
c
c                     Stage 1 - remove bad points
c
c   The flag is equal to 4 times the standard deviation of i or q.
c
      do ipass = 1 , npass1
c
         do ipulse = npdone + 1 , np
c
            do irange = 1 , nr
c
c   If edit flag is set then replace both i and q with zeros
c
               if (
     .         ( abs( sariq(1,irange,ipulse) - ibar_s ) .gt. thresh )
     .                            .or.
     .         ( abs( sariq(2,irange,ipulse) - qbar_s ) .gt. thresh ) )
     .                            then
c
                  nedit                  = nedit + 1
c
                  sariq(1,irange,ipulse) = ibar_s
                  sariq(2,irange,ipulse) = qbar_s
c
               endif
c
            enddo
c
         enddo
c
         call iqstat ( sariq(1,1,npdone+1) , nr , np-npdone , 
     .                 ib , qb , iib , qqb , iqb ,
     .                 ibar , qbar , iibar , qqbar , iqbar )
c
c   Calculate the editing threshold using the latest values, but use
c   time smoothed values for mean removal, rotation, and balancing.
c
         eps     = float( np - npdone ) / float( max( np , npdone ) )
c
         ibar_s  = eps * ibar  + ( 1.0 - eps ) * ibar_s 
         qbar_s  = eps * qbar  + ( 1.0 - eps ) * qbar_s 
         iibar_s = eps * iibar + ( 1.0 - eps ) * iibar_s 
         qqbar_s = eps * qqbar + ( 1.0 - eps ) * qqbar_s 
         iqbar_s = eps * iqbar + ( 1.0 - eps ) * iqbar_s 
c
         thresh  = 3.8 * sqrt( 0.5 *
     .             ( iibar - ibar ** 2 + qqbar - qbar ** 2 ) )
c
      enddo
c
      write ( 7 , '(a33,i14)' )    ' Number of points edited:        ' ,
     .                               nedit
      write ( 7 , '(a33,2e14.6)' ) ' Mean of i, q after edits     :  ' ,
     .                               ibar , qbar
      write ( 7 , '(a33,3e14.6)' ) ' Mean of ii, qq, iq after edits: ' ,
     .                               iibar , qqbar , iqbar
c
c***********************************************************************
c
c                  Stage 2 - Subtract mean I and Q
c
      do ipass = 1 , npass2
c
         do irange = 1 , nr
c
            icor = smth( ib , wts , nwts , nr , irange )
            qcor = smth( qb , wts , nwts , nr , irange )
c
            do ipulse = npdone + 1 , np
c
               sariq(1,irange,ipulse) = sariq(1,irange,ipulse) - icor
               sariq(2,irange,ipulse) = sariq(2,irange,ipulse) - qcor
c
            enddo
c
         enddo
c
         call iqstat ( sariq(1,1,npdone+1) , nr , np-npdone , 
     .                 ib , qb , iib , qqb , iqb ,
     .                 ibar , qbar , iibar , qqbar , iqbar )
c
      enddo
c
      write ( 7 , '(a33,2e14.6)' ) ' Mean of i, q after mean rem.:   ' ,
     .                               ibar , qbar
      write ( 7 , '(a33,3e14.6)' ) ' Mean of ii, qq, iq for rotation:' ,
     .                               iibar , qqbar , iqbar
c
c**********************************************************************
c
c                    Stage 3 - Rotate I and Q
c
      do ipass = 1 , npass3
c
         do irange = 1 , nr
c
            iicor = smth( iib , wts , nwts , nr , irange )
            qqcor = smth( qqb , wts , nwts , nr , irange )
            iqcor = smth( iqb , wts , nwts , nr , irange )
c
            if ( iicor .eq. qqcor ) then
c
               theta = pi / 4.0
c
            else
c
               theta = 0.5 * atan( ( 2.0 * iqcor ) / ( iicor - qqcor ) )
c
            endif
c
            ctheta = cos( theta )
            stheta = sin( theta )
c
            do ipulse = npdone + 1 , np
c
               sariq1 = sariq(1,irange,ipulse)

               sariq(1,irange,ipulse) = ctheta * sariq(1,irange,ipulse)
     .                                + stheta * sariq(2,irange,ipulse)
               sariq(2,irange,ipulse) = - stheta * sariq1
     .                                + ctheta * sariq(2,irange,ipulse)
c
            enddo
c
         enddo
c
         call iqstat ( sariq(1,1,npdone+1) , nr , np-npdone , 
     .                 ib , qb , iib , qqb , iqb ,
     .                 ibar , qbar , iibar , qqbar , iqbar )
c
      enddo
c
      write ( 7 , '(a33,3e14.6)' ) ' Mean of ii, qq, iq for balance: ' ,
     .                               iibar , qqbar , iqbar
c
c***********************************************************************
c
c                     Stage 4 - Scale I and Q
c
      do ipass = 1 , npass4
c
         do irange = 1 , nr
c
            iicor  = smth( iib , wts , nwts , nr , irange )
            qqcor  = smth( qqb , wts , nwts , nr , irange )
            scale  = sqrt( iicor + qqcor )
c
            if ( iicor .ne. 0.0 ) then
               scalei = scale / sqrt( 2.0 * iicor )
            else
               scalei = 1.
            endif
c
            if ( qqcor .ne. 0 ) then
               scaleq = scale / sqrt( 2.0 * qqcor )
            else
               scaleq = 1.
            endif
c
            do ipulse = npdone + 1 , np
c
               sariq(1,irange,ipulse) = scalei * sariq(1,irange,ipulse)
               sariq(2,irange,ipulse) = scaleq * sariq(2,irange,ipulse)
c
            enddo
c
         enddo
c
         call iqstat ( sariq(1,1,npdone+1) , nr , np-npdone , ib , 
     .                 qb , iib , qqb , iqb ,
     .                 ibar , qbar , iibar , qqbar , iqbar )
c
      enddo
c
      write ( 7 , '(a33,3e14.6)' ) ' Final means of ii, qq, iq:      ' ,
     .                               iibar , qqbar , iqbar
c
      return
      end
C**
C***********************************************************************
C**
      subroutine iqstat ( sariq , nr , np , ib , qb , iib , qqb , iqb ,
     .                    ibar , qbar , iibar , qqbar , iqbar )
C**
C***********************************************************************
C**
      implicit none
c
      integer  nr , np , irange , ipulse
c
      real     sariq(2,nr,np) , ib(nr) , qb(nr) , iib(nr) , qqb(nr) ,
     .         iqb(nr) , ibar , qbar , iibar , qqbar , iqbar
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      ibar  = 0.0
      qbar  = 0.0
      iibar = 0.0
      qqbar = 0.0
      iqbar = 0.0
c
      do irange = 1 , nr
c
         ib(irange)  = 0.0
         qb(irange)  = 0.0
         iib(irange) = 0.0
         qqb(irange) = 0.0
         iqb(irange) = 0.0
c
         do ipulse = 1 , np
c
            ib(irange)  = ib(irange)  + sariq(1,irange,ipulse)
            qb(irange)  = qb(irange)  + sariq(2,irange,ipulse)
            iib(irange) = iib(irange) + sariq(1,irange,ipulse) ** 2
            qqb(irange) = qqb(irange) + sariq(2,irange,ipulse) ** 2
            iqb(irange) = iqb(irange) + sariq(1,irange,ipulse) *
     .                                  sariq(2,irange,ipulse)
c
         enddo
c
         ib(irange)  = ib(irange)  / float( np )
         qb(irange)  = qb(irange)  / float( np )
         iib(irange) = iib(irange) / float( np )
         qqb(irange) = qqb(irange) / float( np )
         iqb(irange) = iqb(irange) / float( np )
c
         ibar        = ibar  + ib(irange)  / float( nr )
         qbar        = qbar  + qb(irange)  / float( nr )
         iibar       = iibar + iib(irange) / float( nr )
         qqbar       = qqbar + qqb(irange) / float( nr )
         iqbar       = iqbar + iqb(irange) / float( nr )
c
      enddo
c
      return
      end
C**
C***********************************************************************
C**
      real function smth ( y , wts , nwts , nr , irange )
C**
C***********************************************************************
C**
      implicit none
c
      integer iwt , nwts , nr , irange , nwhalf , irmin , ir
c
      real    y(nr) , wts(nwts) , ytot , wtot
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      nwhalf = nwts / 2
c
      ytot   = 0.0
      wtot   = 0.0
c
      irmin  = min( irange , nr - irange + 1 )
      irmin  = min( irmin , nwhalf )
c
      do iwt = ( nwhalf + 2 ) - irmin , nwhalf + irmin 
c
         ir = irange - 1 - nwhalf + iwt
c
         if ( ir .ge. 1 .and. ir .le. nr ) then
c
            ytot = ytot + wts(iwt) * y(ir)
            wtot = wtot + wts(iwt)
c
         endif
c
      enddo
c
      smth = ytot / wtot
c
      return
      end
