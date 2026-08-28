c**
c***********************************************************************
c**
      double precision function clight_avg ( alt_km )
c**
c***********************************************************************
c**
c   Purpose: To compute the speed of light averaged over the path to
c            the ground from altitude alt_km
c
c            alt_km      :  Altitude in km
c
c            clight_avg  :  Speed of light in m/s averaged to surface
c
c-----------------------------------------------------------------------
c
      implicit none
c
      real             alt_km
c
      double precision clight_vacuum , refractivity_avg , alt_km7 ,
     .                 refrac_surface
c
      clight_vacuum  = 2.99792458D+8  !  Speed of light in a vacuum
c
      refrac_surface = 313.0D0        !  Refractivity at surface
c
      alt_km7        = dble( alt_km ) / 7.0D0
c
      if ( alt_km .lt. 0.001 ) then
c
         refractivity_avg = refrac_surface
c
      else
c
         refractivity_avg = refrac_surface * ( 1.0D0 / alt_km7 ) *
     .                      ( 1.0D0 - dexp( - alt_km7 ) )
c
      endif
c
      clight_avg = clight_vacuum / ( 1.0D0 + 1.0D-6 * refractivity_avg )
c
      return
      end
