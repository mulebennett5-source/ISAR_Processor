C**
C***********************************************************************
C**
      subroutine getAccTen(nf,nr,nt,dotdot,fStart,rStart,tStart,deltaF,
     .                     deltaR,deltaT,lambda,drtotal,timetxt,phitxt,
     .                     thetatxt,rcen,fcen,npulse,arbar,afbar)
C**
C***********************************************************************
C**
      implicit none
c
c Input:
c
c    nf         - number of frequency points
c    nr         - number of range points
c    nt         - number of time points
c    fStart     - (Hz) initial frequency coordinate of dotDot array
c    rStart     - (m) initial range coordinate of dotDot array
c    tStart     - (s) initial time coordinate of dotDot array
c    deltaF     - (Hz) frequency sampling of dotDot array
c    deltaR     - (m) range sampling of dotDot array
c    deltaT     - (s) time sampling of dotDot array
c    lambda     - (m) Wavelength
c    drtotal    - (m) Total mocomp range from input file
c    timetxt    - Time from text input file
c    phitxt     - Platform azimuth angle from text input file
c    thetatxt   - Platform elevation angle from text input file
c
c Output:
c
c    dotDot     - 3-D master particle acceleration array (Hz/second)
c    rCen       - (m) range centroid of distributed target
c    fCen       - (Hz) Doppler centroid of distributed target
c 

      integer nr , nf , nt , npulse , neart , neartp , ndt_use ,
     .        kt , kr , kf , option
c
      real    Ar(nt) , Af(nt) , time(nt) , fStart , rStart , tStart ,
     .        deltaF , deltaR , deltaT , lambda , dtime , delta ,
     .        timetxt(npulse) , phitxt(npulse) , thetatxt(npulse) ,
     .        dt_use , theta , phidotm , phidotp , phidot , phidotdotm ,
     .        phidotdotp , phidotdot , thetadotm , thetadotp ,
     .        thetadot , thetadotdotm , thetadotdotp , thetadotdot ,
     .        rCen , fCen , tantheta , r , f , arbar , afbar ,
     .        drtotal(npulse)
      
      real    dotdot(nf,nr,nt)

      dtime   = (timetxt(npulse)-timetxt(1))/float(npulse-1)
      
      do kt = 1 , nt

         time(kt)     = tStart + ( kt - 1 ) * deltaT;

      enddo

      do kt = 1 , nt

         ndt_use      = int(0.5+deltaT/dtime)
         dt_use       = ndt_use*dtime

         neart        = min(1+int((time(kt)-timetxt(1))/dtime),
     .                      npulse)
         neartp       = min(neart+1,npulse)
         delta        = 1+((time(kt)-timetxt(1))/dtime)-neart

         theta        = (1-delta) * thetatxt(neart) +
     .                      delta * thetatxt(neartp)

         phidotm      = (phitxt(neart+ndt_use)-phitxt(neart-ndt_use))
     .                  /(2*dt_use)
         phidotp      = (phitxt(neartp+ndt_use)-
     .                   phitxt(neartp-ndt_use))/(2*dt_use)
         phidot       = (1-delta) * phidotm + delta * phidotp

         thetadotm    = (thetatxt(neart+ndt_use)  -
     .                   thetatxt(neart-ndt_use))  /(2*dt_use)
         thetadotp    = (thetatxt(neartp+ndt_use) -
     .                   thetatxt(neartp-ndt_use)) /(2*dt_use)
         thetadot     = (1-delta) * thetadotm + delta * thetadotp

         phidotdotm   = (phitxt(neart+ndt_use) +
     .                   phitxt(neart-ndt_use) - 2*phitxt(neart)) /
     .                   (dt_use**2)
         phidotdotp   = (phitxt(neartp+ndt_use) +
     .                   phitxt(neartp-ndt_use)
     .                  - 2*phitxt(neartp)) / (dt_use**2)
         phidotdot    = (1-delta) * phidotdotm  + delta * phidotdotp

         thetadotdotm = (thetatxt(neart+ndt_use) +
     .                   thetatxt(neart-ndt_use) - 2*thetatxt(neart))
     .                   /(dt_use**2);
         thetadotdotp = (thetatxt(neartp+ndt_use)+
     .                   thetatxt(neartp-ndt_use)-2*thetatxt(neartp))
     .                   /(dt_use**2);
         thetadotdot  = (1-delta)*thetadotdotm+delta*thetadotdotp

         tantheta     = tan(theta)

         Ar(kt)       = -(phidot**2+thetadot**2)+(thetadotdot/tantheta)
     .                  -thetadot*phidotdot/(tantheta*phidot)
     .                  -2*(thetadot/tantheta)**2

         Af(kt)       = (phidotdot/phidot)+2*thetadot/tantheta

      enddo

      arbar = 0.0
      afbar = 0.0
      do kt = 1 , nt
         arbar = arbar+Ar(kt)
         afbar = afbar+Af(kt)
      enddo
      arbar = arbar/float(nt)
      afbar = afbar/float(nt)
c
      Ar(:)  = arbar
      Af(:)  = afbar

c Fill acceleration array

      rCen   = 0
      fCen   = 0

      do kr = 1 , nr

         r = rStart+(kr-1)*deltaR

         do kf = 1 , nf

            f               = fStart+(kf-1)*deltaF

            do kt = 1 , nt

               dotDot(kf,kr,kt) = (2/lambda)*r*Ar(kt)+f*Af(kt)

            enddo
               
         enddo
    
      enddo

      return
      end
