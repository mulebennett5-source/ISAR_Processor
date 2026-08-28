c
c Source:  http://www.mech.utah.edu/~brannon/public/rotation.pdf
c
c---.----1----.----2----.----3----.----4----.----5----.----6----.----7--
      SUBROUTINE EU2DC(EU,DC,ierr)
c
c INPUT
c -----
c EU: The Euler angles
c
c OUTPUT
c ------
c DC: The direction cosine matrix.
c
************************************************************************
c
      implicit none
c      integer ierr
      real*8  eu(3),dc(3,3)
      real*8  cphi,sphi,cth,sth,cpsi,spsi
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      ierr    = 0
      cphi    = COS(eu(1))
      sphi    = SIN(eu(1))
      cth     = COS(eu(2))
      sth     = SIN(eu(2))
      cpsi    = COS(eu(3))
      spsi    = SIN(eu(3))
      dc(1,1) = cphi*cpsi - cth*sphi*spsi
      dc(2,1) = cpsi*sphi + cphi*cth*spsi
      dc(3,1) = spsi*sth
      dc(1,2) = -(cpsi*cth*sphi) - cphi*spsi
      dc(2,2) = cphi*cpsi*cth - sphi*spsi
      dc(3,2) = cpsi*sth
      dc(1,3) = sphi*sth
      dc(2,3) = -(cphi*sth)
      dc(3,3) = cth
      RETURN
      END
c---.----1----.----2----.----3----.----4----.----5----.----6----.----7--
      SUBROUTINE DC2EU(DC,EULER,ierr)
c
c This routine converts an orthogonal proper rotation matrix to euler
c angles.
c
c INPUT
c -----
c DC: direction cosine matrix for the rotation TENSOR
C DC(i,j) is defined as the inner product between the ith lab
C base vector with the jth rotated base vector.
c
c OUTPUT
c ------
C EULER: Euler angles {phi,theta,psi}, which describe the rotation
C of the lab triad by the following procedure...First rotate
C the triad an angle phi about its z-axis. Then rotate the new
C triad an angle theta about its own new x-axis. Then rotate
C the newer triad an angle psi about its own z-axis.
c IERR: =0 if success, =1 otherwise
c
c***********************************************************************
c
      implicit none
c
      real*8     pzero,pone,puny,small
      parameter (pzero=0.0d0,pone=0.1d1,small=0.1d-8,puny=0.1d-20)
c
      integer ierr
      real*8     dc(3,3), euler(3)
c
      real*8     dum,dum11,dum12,dum21,dum22,cth,sthsq,theta,psi,phi
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      cth   = dc(3,3)
      sthsq = pone-cth*cth
c
      if ( sthsq .lt. puny ) then
c
         theta = acos(cth)
         psi   = pzero
         phi   = atan2(dc(2,1),dc(1,1))
         return
c
      endif
c
      dum11 = dc(1,1)+((dc(2,3)*dc(3,2) + dc(1,3)*dc(3,1)*dc(3,3))/sthsq)
      dum12 = dc(1,2)-(dc(2,3)*dc(3,1) - dc(1,3)*dc(3,2)*dc(3,3))/sthsq
      dum21 = dc(2,1)-(dc(1,3)*dc(3,2) - dc(2,3)*dc(3,1)*dc(3,3))/sthsq
      dum22 = dc(2,2)+(dc(1,3)*dc(3,1) + dc(2,3)*dc(3,2)*dc(3,3))/sthsq
      dum   = sqrt(dum11**2+dum22**2+dum12**2+dum21**2)
c
      if ( dum .gt. small ) then
c
c The DC matrix is not proper orthogonal
         ierr = 1
         return
c
      endif
c
      theta    = acos(cth)
      phi      = atan2(dc(1,3),-dc(2,3))
      psi      = atan2(dc(3,1), dc(3,2))
      euler(1) = phi
      euler(2) = theta
      euler(3) = psi
      ierr     = 0
c
      return
      end
