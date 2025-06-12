 /*@@
   @file      GRHydro_ShockTubeM.F90
   @date      Sep 23, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   Initial data of the shock tube type - MHD version.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)
#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define Bconsx(i,j,k) Bcons(i,j,k,1)
#define Bconsy(i,j,k) Bcons(i,j,k,2)
#define Bconsz(i,j,k) Bcons(i,j,k,3)

#define OOSQRT2 (0.7071067811865475244008442)
#define OOSQRT3 (0.5773502691896257645091489)
#define OOSQRT6 (0.4082482904638630163662140)


 /*@@
   @routine    GRHydro_shocktubeM
   @date       Sat Jan 26 02:53:49 2002
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   Initial data for shock tubes, parallel to
   a coordinate axis. Either Sods problem or the standard shock tube.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Expansion and alteration of the test code from GRAstro_Hydro, 
   written by Mark Miller.
   @endhistory 

@@*/

subroutine GRHydro_shocktubeM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det, DIRECTION_TINY
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr
  CCTK_REAL :: bvcxl,bvcyl,bvczl,bvcxr,bvcyr,bvczr
  CCTK_REAL :: ux,uy,uz,ut,tmp,tmp2,tmp3
  
  
  bvcxl = Bx_init
  bvcyl = By_init
  bvczl = Bz_init
  bvcxr = Bx_init
  bvcyr = By_init
  bvczr = Bz_init

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  do i=1,nx
    do j=1,ny
      do k=1,nz
        if (CCTK_EQUALS(shocktube_type,"diagshock")) then
!!$ The diagshock choice yields a shock plane perpendicular to the fixed vector (1,1,1)
!!$ This could be changed, but would require 3 new params containing the new shock direction
           direction = x(i,j,k) - shock_xpos + &
                y(i,j,k) - shock_ypos + z(i,j,k) - shock_zpos

        else if (CCTK_EQUALS(shocktube_type,"diagshock2d")) then
!!$ The diagshock choice yields a shock plane perpendicular to the fixed vector (1,1,0), with similarity in the z-dir.
!!$ This could be changed, but would require 2 new params containing the new shock direction
           direction = x(i,j,k) - shock_xpos + &
                y(i,j,k) - shock_ypos

        else if (CCTK_EQUALS(shocktube_type,"xshock")) then
           direction = x(i,j,k) - shock_xpos
        else if (CCTK_EQUALS(shocktube_type,"yshock")) then
           direction = y(i,j,k) - shock_ypos
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
           direction = z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"sphere")) then
           direction = sqrt((x(i,j,k)-shock_xpos)**2+&
                (y(i,j,k)-shock_ypos)**2+&
                (z(i,j,k)-shock_zpos)**2)-shock_radius
        end if
        if (CCTK_EQUALS(shock_case,"Simple")) then
          rhol = 10.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 2.d0
          epsr = 1.d-6
        else if (CCTK_EQUALS(shock_case,"Sod")) then
          rhol = 1.d0
          rhor = 0.125d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 1.5d0
          epsr = 1.2d0
!!$This line only for polytrope, k=1
!!$          epsr = 0.375d0
        else if (CCTK_EQUALS(shock_case,"Blast")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 1500.d0
          epsr = 1.5d-2

!!$ The following shocktubes are from Balsara, ApJSS 132 p83 (2001),
!!$ http://adsabs.harvard.edu/abs/2001ApJS..132...83B.
!!$   All use n=1600 cells, over domain x=[-0.5,0.5]
!!$   All assume ideal-gas or gamma-law EOS, the first test uses GAMMA=2. while the rest GAMMA=5./3.

!!$ Unmagnetized Test 1 (rel. Brio & Wu 1988 by van Putten 1993) of Balsara 2001 -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara0")) then
          rhol = 1.0d0
          rhor = 0.125d0
          velxl = 0.0d0
          velxr = 0.0d0
          velyl = 0.0d0
          velyr = 0.0d0
          velzl = 0.0d0
          velzr = 0.0d0
          bvcxl=0.0d0
          bvcxr=0.0d0
          bvcyl=0.0d0
          bvcyr=0.0d0
          bvczl=0.0d0
          bvczr=0.0d0
          epsl = 1.0d0/rhol
          epsr = 0.1d0/rhor

!!$ Test 1 (rel. Brio & Wu 1988 by van Putten 1993) of Balsara 2001 -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara1")) then
          rhol = 1.0d0
          rhor = 0.125d0
          velxl = 0.0d0
          velxr = 0.0d0
          velyl = 0.0d0
          velyr = 0.0d0
          velzl = 0.0d0
          velzr = 0.0d0
          bvcxl=0.5d0
          bvcxr=0.5d0
          bvcyl=1.0d0
          bvcyr=-1.0d0
          bvczl=0.0d0
          bvczr=0.0d0
          epsl = 1.0d0/rhol
          epsr = 0.1d0/rhor

!!$ Test 2 (blast wave) of Balsara 2001  -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara2")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          bvcxl=5.0d0
          bvcxr=5.0d0
          bvcyl=6.d0
          bvcyr=0.7d0
          bvczl=6.d0
          bvczr=0.7d0
          epsl = 1.5d0*30.0d0/rhol
          epsr = 1.5d0*1.0d0/rhor

!!$ Test 3 (blast wave) of Balsara 2001  -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara3")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          bvcxl=10.0d0
          bvcxr=10.0d0
          bvcyl=7.d0
          bvcyr=0.7d0
          bvczl=7.d0
          bvczr=0.7d0
          epsl = 1.5d0*1000.0d0/rhol
          epsr = 1.5d0*0.1d0/rhor

!!$ Test 4 (rel. version of Noh 1987) of Balsara 2001  -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara4")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.999d0
          velxr = -0.999d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          bvcxl=10.0d0
          bvcxr=10.0d0
          bvcyl=7.d0
          bvcyr=-7.d0
          bvczl=7.d0
          bvczr=-7.d0
          epsl = 1.5d0*0.1d0/rhol
          epsr = 1.5d0*0.1d0/rhor

!!$ Test 5 (non-coplanar set of waves) of Balsara 2001  -- compare at t=0.55
       else if (CCTK_EQUALS(shock_case,"Balsara5")) then
          rhol = 1.08d0
          rhor = 1.d0
          velxl = 0.4d0
          velxr = -0.45d0
          velyl = 0.3d0
          velyr = -0.2d0
          velzl = 0.2d0
          velzr = 0.2d0
          bvcxl=2.0d0
          bvcxr=2.0d0
          bvcyl=0.3d0
          bvcyr=-0.7d0
          bvczl=0.3d0
          bvczr=0.5d0
          epsl = 1.5d0*0.95d0/rhol
          epsr = 1.5d0*1.0d0/rhor

!!$  "Generic Alfven Test of  Giacomazzo and Rezzolla J.Comp.Phys (2006)
       else if (CCTK_EQUALS(shock_case,"Alfven")) then
          rhol = 1.d0
          rhor = 0.9d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.3d0
          velyr = 0.d0
          velzl = 0.4d0
          velzr = 0.d0
          bvcxl=1.0d0
          bvcxr=1.0d0
          bvcyl=6.d0
          bvcyr=5.d0
          bvczl=2.d0
          bvczr=2.d0
          epsl = 1.5d0*5.d0/rhol
          epsr = 1.5d0*5.3d0/rhor

!!$ The following 9 tests are from Komissarov MNRAS 303 p343 (1999),
!!$ http://adsabs.harvard.edu/abs/1999MNRAS.303..343K. 
!!$   Note that the data is specified in terms of the 4-velocity, so some conversion is necessary. 
!!$   All assume ideal-gas or gamma-law EOS with GAMMA=4./3. 

!!$ Fast Shock test of Komissarov 1999 -- compare  at t=2.5 , n=40, x=[-1,1]
       else if (CCTK_EQUALS(shock_case,"Komissarov1")) then
          rhol = 1.d0
          rhor = 25.48d0
          bvcxl=20.d0
          bvcxr=20.d0
          bvcyl=25.02d0
          bvcyr=49.d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d0  /rhol
          epsr = 3.d0 * 367.5d0 /rhor

          ux=25.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=1.091d0
          uy=0.3923d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

!!$ Slow Shock test of Komissarov 1999 -- compare at t=2. , n=200,  x=[-0.5,1.5]
       else if (CCTK_EQUALS(shock_case,"Komissarov2")) then
          rhol = 1.d0
          rhor = 3.323d0
          bvcxl=10.d0
          bvcxr=10.d0
          bvcyl=18.28d0
          bvcyr=14.49d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d0  /rhol
          epsr = 3.d0 * 55.36d0 /rhor

          ux=1.53d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=0.9571d0
          uy=-0.6822d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

!!$ Switch-off Fast test of Komissarov 1999 -- compare  at t=1. , n=150, x=[-1,1] 
       else if (CCTK_EQUALS(shock_case,"Komissarov3")) then
          rhol = 0.1d0
          rhor = 0.562d0
          bvcxl=2.d0
          bvcxr=2.d0
          bvcyl=0.d0
          bvcyr=4.710d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d0  /rhol
          epsr = 3.d0 * 10.d0 /rhor

          ux=-2.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=-0.212d0
          uy=-0.590d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

!!$ Switch-on Slow test of Komissarov 1999 -- compare  at t=2. , n=150, x=[-1,1.5]
       else if (CCTK_EQUALS(shock_case,"Komissarov4")) then
          rhol = 1.78d-2
          rhor = 1.d-2
          bvcxl=1.d0
          bvcxr=1.d0
          bvcyl=1.022d0
          bvcyr=0.d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 0.1d0  /rhol
          epsr = 3.d0 * 1.d0 /rhor

          ux=-0.765d0
          uy=-1.386d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=0.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

!!$ Alfven wave test of Komissarov 1999 -- compare  at t=2. , n=200, x=[-1,1.5]
!!$   Needs special setup -- FIX
       else if (CCTK_EQUALS(shock_case,"Komissarov5")) then
          rhol = 1.d0
          rhor = 1.d0
          bvcxl=3.d0
          bvcxr=3.d0
          bvcyl=3.d0
          bvcyr=-6.857d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d0  /rhol
          epsr = 3.d0 * 1.d0 /rhor

          ux=0.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=3.70d0
          uy=5.76d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

!!$ Compound wave test of Komissarov 1999 -- compare  at t=0.1,0.75,1.5 , n=200, x=[-0.5,1.5]
!!$   Needs special setup -- FIX
       else if (CCTK_EQUALS(shock_case,"Komissarov6")) then
          rhol = 1.d0
          rhor = 1.d0
          bvcxl=3.d0
          bvcxr=3.d0
          bvcyl=3.d0
          bvcyr=-6.857d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d0  /rhol
          epsr = 3.d0 * 1.d0 /rhor

          ux=0.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=3.70d0
          uy=5.76d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

!!$ Shock Tube 1  test of Komissarov 1999 -- compare  at t=1. , n=400, x=[-1,1.5]
       else if (CCTK_EQUALS(shock_case,"Komissarov7")) then
          rhol = 1.d0
          rhor = 1.d0
          bvcxl=1.d0
          bvcxr=1.d0
          bvcyl=0.d0
          bvcyr=0.d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d3  /rhol
          epsr = 3.d0 * 1.d0 /rhor

          velxl = 0.
          velyl = 0.
          velzl = 0.

          velxr = 0.
          velyr = 0.
          velzr = 0.

!!$ Shock Tube 2 test of Komissarov 1999 -- compare  at t=1. , n=500, x=[-1.25,1.25]
       else if (CCTK_EQUALS(shock_case,"Komissarov8")) then
          rhol = 1.d0
          rhor = 1.d0
          bvcxl=0.d0
          bvcxr=0.d0
          bvcyl=20.d0
          bvcyr=0.d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 30.d0  /rhol
          epsr = 3.d0 * 1.d0 /rhor

          velxl = 0.
          velyl = 0.
          velzl = 0.

          velxr = 0.
          velyr = 0.
          velzr = 0.

!!$ Collision test of Komissarov 1999 -- compare  at t=1.22 , n=200, x=[-1,1]
       else if (CCTK_EQUALS(shock_case,"Komissarov9")) then
          rhol = 1.d0
          rhor = 1.d0
          bvcxl=10.d0
          bvcxr=10.d0
          bvcyl=10.d0
          bvcyr=-10.d0
          bvczl=0.d0
          bvczr=0.d0
          epsl = 3.d0 * 1.d0  /rhol
          epsr = 3.d0 * 1.d0 /rhor

          ux=5.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut

          ux=-5.d0
          uy=0.d0
          uz=0.d0
          ut = sqrt(1. + ux*ux + uy*uy + uz*uz)
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut

       else
          call CCTK_WARN(0,"Shock case not recognized")
        end if

        DIRECTION_TINY = 0.0001*CCTK_DELTA_SPACE(1)

        if ( ((change_shock_direction==0).and.(direction .lt. DIRECTION_TINY)).or.& 
             ((change_shock_direction==1).and.(direction .gt. DIRECTION_TINY)) ) then

!!$ Left state

          rho(i,j,k) = rhol
          velx(i,j,k) = velxl
          vely(i,j,k) = velyl
          velz(i,j,k) = velzl
          eps(i,j,k) = epsl
          Bvecx(i,j,k)=bvcxl
          Bvecy(i,j,k)=bvcyl
          Bvecz(i,j,k)=bvczl
        else

!!$ Right state

          rho(i,j,k) = rhor
          velx(i,j,k) = velxr
          vely(i,j,k) = velyr
          velz(i,j,k) = velzr
          eps(i,j,k) = epsr
          Bvecx(i,j,k)=bvcxr
          Bvecy(i,j,k)=bvcyr
          Bvecz(i,j,k)=bvczr
        end if

        
        if (CCTK_EQUALS(shocktube_type,"yshock")) then
!!$ Cycle x,y,z forward           
           tmp=velx(i,j,k)
           velx(i,j,k)=velz(i,j,k)
           velz(i,j,k)=vely(i,j,k)
           vely(i,j,k)=tmp
           tmp=Bvecx(i,j,k)
           Bvecx(i,j,k)=Bvecz(i,j,k)
           Bvecz(i,j,k)=Bvecy(i,j,k)
           Bvecy(i,j,k)=tmp
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
!!$ Cycle x,y,z backward           
           tmp=velx(i,j,k)
           velx(i,j,k)=vely(i,j,k)
           vely(i,j,k)=velz(i,j,k)
           velz(i,j,k)=tmp
           tmp=Bvecx(i,j,k)
           Bvecx(i,j,k)=Bvecy(i,j,k)
           Bvecy(i,j,k)=Bvecz(i,j,k)
           Bvecz(i,j,k)=tmp
        else if (CCTK_EQUALS(shocktube_type,"diagshock")) then
!!$ Rotated basis vectors necessary to evaluate the orthogonal matrix elements:
!!$ xhat = 1/sqrt(3)[1,1,1], yhat = 1/sqrt(2)[-1,1,0]; zhat = 1/sqrt(6)[-1,-1,2]
!!$ Orthogonal matrix constructed from the tensor product between the original
!!$ cartesian basis x's and new basis vectors xhat's, rotated towards the diagonal 
!!$ shock normal and tangent directions.
           tmp  = OOSQRT3*velx(i,j,k) - OOSQRT2*vely(i,j,k) - OOSQRT6*velz(i,j,k)
           tmp2 = OOSQRT3*velx(i,j,k) + OOSQRT2*vely(i,j,k) - OOSQRT6*velz(i,j,k)
           tmp3 = OOSQRT3*velx(i,j,k)                  + 2.d0*OOSQRT6*velz(i,j,k)
           velx(i,j,k)=tmp
           vely(i,j,k)=tmp2
           velz(i,j,k)=tmp3
           tmp  = OOSQRT3*Bvecx(i,j,k) - OOSQRT2*Bvecy(i,j,k) - OOSQRT6*Bvecz(i,j,k)
           tmp2 = OOSQRT3*Bvecx(i,j,k) + OOSQRT2*Bvecy(i,j,k) - OOSQRT6*Bvecz(i,j,k)
           tmp3 = OOSQRT3*Bvecx(i,j,k)                   + 2.d0*OOSQRT6*Bvecz(i,j,k)
           Bvecx(i,j,k)=tmp
           Bvecy(i,j,k)=tmp2
           Bvecz(i,j,k)=tmp3
        else if (CCTK_EQUALS(shocktube_type,"diagshock2d")) then
!!$ New basis:
!!$ xhat = 1/sqrt(2)[1,1,0], yhat = 1/sqrt(2)[-1,1,0]; zhat = [0,0,1]
           tmp  = OOSQRT2*velx(i,j,k) - OOSQRT2*vely(i,j,k)
           tmp2 = OOSQRT2*velx(i,j,k) + OOSQRT2*vely(i,j,k)
           velx(i,j,k)=tmp
           vely(i,j,k)=tmp2
           tmp  = OOSQRT2*Bvecx(i,j,k) - OOSQRT2*Bvecy(i,j,k)
           tmp2 = OOSQRT2*Bvecx(i,j,k) + OOSQRT2*Bvecy(i,j,k)
           Bvecx(i,j,k)=tmp
           Bvecy(i,j,k)=tmp2
        endif
           
        det=SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

        if (CCTK_EQUALS(GRHydro_eos_type,"Polytype")) then
          call Prim2ConPolyM(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),&
               w_lorentz(i,j,k))
        else
          call Prim2ConGenM(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),&
               w_lorentz(i,j,k))
        end if
    enddo
    enddo
  enddo

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0
  Bconsrhs = 0.d0


  return
  
end subroutine GRHydro_shocktubeM

subroutine GRHydro_Diagshock_BoundaryM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz, sten, stenp1, minsum, maxsum
  CCTK_INT :: inew, jnew, knew
  CCTK_INT :: xoff,yoff,zoff,indsum
  CCTK_REAL :: det
  CCTK_INT :: status
  character(len=200) :: error_message
  integer :: gindex_tau, gindex_dens, gindex_scon, gindex_Bcons
  CCTK_INT, parameter :: num_groups = 4
  CCTK_INT :: groups(num_groups)


  sten=GRHydro_stencil
  stenp1=GRHydro_stencil + 1

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  minsum = 3*stenp1
  maxsum = nx+ny+nz-3*stenp1 + 3

  xoff=0
  yoff=0
  zoff=0
        
!!$ Loop over the cubical domain 6 faces:

  !!$ 1) xmin:
  do k = stenp1, nz-sten
    do j = stenp1, ny-sten
      do i = 1, sten

        indsum = i+j+k
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          knew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - stenp1 + 1
          jnew = ny - stenp1 + 1
          knew = nz - stenp1 + 1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
           call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
          dens(i,j,k) = dens(inew,jnew,knew)
        end if
           
        sx(i,j,k) = sx(inew,jnew,knew)
        sy(i,j,k) = sy(inew,jnew,knew)
        sz(i,j,k) = sz(inew,jnew,knew)
        tau(i,j,k) = tau(inew,jnew,knew)
        Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
        Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
        Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,knew)
        endif

      end do 
    end do
  end do 

  !!$ 2) xmax:
  do k = stenp1, nz-sten
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx

        indsum = i+j+k
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          knew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - stenp1 + 1
          jnew = ny - stenp1 + 1
          knew = nz - stenp1 + 1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
           call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
           dens(i,j,k) = dens(inew,jnew,knew)
        end if
           
        sx(i,j,k) = sx(inew,jnew,knew)
        sy(i,j,k) = sy(inew,jnew,knew)
        sz(i,j,k) = sz(inew,jnew,knew)
        tau(i,j,k) = tau(inew,jnew,knew)
        Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
        Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
        Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,knew)
        endif

      end do 
    end do
  end do 

  !!$ 3) ymin:
  do k = stenp1, nz-sten
    do j = 1, sten
      do i =  stenp1, nx-sten

        indsum = i+j+k
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          knew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - stenp1 + 1
          jnew = ny - stenp1 + 1
          knew = nz - stenp1 + 1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
           call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
           dens(i,j,k) = dens(inew,jnew,knew)
        end if
           
        sx(i,j,k) = sx(inew,jnew,knew)
        sy(i,j,k) = sy(inew,jnew,knew)
        sz(i,j,k) = sz(inew,jnew,knew)
        tau(i,j,k) = tau(inew,jnew,knew)
        Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
        Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
        Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,knew)
        endif

      end do 
    end do
  end do 

  !!$ 4) ymax:
  do k = stenp1, nz-sten
    do j = ny-sten+1, ny
      do i = stenp1, nx-sten

        indsum = i+j+k
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          knew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - stenp1 + 1
          jnew = ny - stenp1 + 1
          knew = nz - stenp1 + 1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
           call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
           dens(i,j,k) = dens(inew,jnew,knew)
        end if
           
        sx(i,j,k) = sx(inew,jnew,knew)
        sy(i,j,k) = sy(inew,jnew,knew)
        sz(i,j,k) = sz(inew,jnew,knew)
        tau(i,j,k) = tau(inew,jnew,knew)
        Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
        Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
        Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,knew)
        endif

      end do 
    end do
  end do 

!!$  the following two cases are different for a 3d diagonal than a 2d diagonal

  !!$ 5) zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten
      do i =  stenp1, nx-sten
 
        indsum = i+j+k
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          knew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - stenp1 + 1
          jnew = ny - stenp1 + 1
          knew = nz - stenp1 + 1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
           call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
           dens(i,j,k) = dens(inew,jnew,knew)
        end if
           
        sx(i,j,k) = sx(inew,jnew,knew)
        sy(i,j,k) = sy(inew,jnew,knew)
        sz(i,j,k) = sz(inew,jnew,knew)
        tau(i,j,k) = tau(inew,jnew,knew)
        Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
        Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
        Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,knew)
        endif

      end do 
    end do
  end do 

  !!$ 6) zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = stenp1, nx-sten
 
        indsum = i+j+k
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          knew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - stenp1 + 1
          jnew = ny - stenp1 + 1
          knew = nz - stenp1 + 1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
           call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
           dens(i,j,k) = dens(inew,jnew,knew)
        end if
           
        sx(i,j,k) = sx(inew,jnew,knew)
        sy(i,j,k) = sy(inew,jnew,knew)
        sz(i,j,k) = sz(inew,jnew,knew)
        tau(i,j,k) = tau(inew,jnew,knew)
        Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
        Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
        Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,knew)
        endif

      end do 
    end do
  end do 

!!$ Loop over the cubical domain 12 edges:

  !!$ 1) xmin,ymin:
  do k = stenp1, nz-sten
    do j = 1, sten
      do i = 1, sten

         indsum = i+j+k
!!$ indsum > maxsum impossible!
         if(indsum.lt.minsum) then
            inew = stenp1
            jnew = stenp1
            knew = stenp1
         else
            call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         end if

!!$  Don't poison edges!
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 2) xmin,ymax:
  do k = stenp1, nz-sten
    do j =  ny-sten+1, ny
      do i = 1, sten

         indsum = i+j+k
!!$ indsum > maxsum and indsum < minsum both impossible!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 3) xmax,ymin:
  do k = stenp1, nz-sten
    do j = 1, sten
      do i = nx-sten+1, nx 

         indsum = i+j+k
!!$ indsum > maxsum and indsum < minsum both impossible!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 4) xmax,ymax:
  do k = stenp1, nz-sten
    do j = ny-sten+1, ny 
      do i = nx-sten+1, nx 

         indsum = i+j+k
!!$ indsum < minsum impossible!
         if(indsum.gt.maxsum) then
            inew = nx - stenp1 + 1
            jnew = ny - stenp1 + 1
            knew = nz - stenp1 + 1
         else
            call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         end if
         
!!$  Don't poison edges!
         dens(i,j,k) = dens(inew,jnew,knew)

         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 5) ymin,zmin:
  do k = 1, sten 
    do j = 1, sten 
      do i =  stenp1, nx-sten

         indsum = i+j+k
!!$ indsum > maxsum impossible!
         if(indsum.lt.minsum) then
            inew = stenp1
            jnew = stenp1
            knew = stenp1
         else
            call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         end if
         
!!$  Don't poison edges!
         dens(i,j,k) = dens(inew,jnew,knew)

         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif
 
      end do 
    end do
  end do 

  !!$ 6) ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten 
      do i =  stenp1, nx-sten
 
         indsum = i+j+k
!!$ indsum > maxsum and indsum < minsum both impossible!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 7) ymax,zmin:
  do k = 1, sten 
    do j = ny-sten+1, ny
      do i =  stenp1, nx-sten
 
         indsum = i+j+k
!!$ indsum > maxsum and indsum < minsum both impossible!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 8) ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny
      do i =  stenp1, nx-sten
 
         indsum = i+j+k
!!$ indsum < minsum impossible!
         if(indsum.gt.maxsum) then
            inew = nx - stenp1 + 1
            jnew = ny - stenp1 + 1
            knew = nz - stenp1 + 1
         else
            call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         end if

!!$  Don't poison edges!
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 9) xmin,zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten 
      do i = 1, sten 
 
         indsum = i+j+k
!!$ indsum > maxsum impossible!
         if(indsum.lt.minsum) then
            inew = stenp1
            jnew = stenp1
            knew = stenp1
         else
            call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         end if
         
!!$  Don't poison edges!
         dens(i,j,k) = dens(inew,jnew,knew)

         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 10) xmin,zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = 1, sten 

         indsum = i+j+k
!!$ indsum > maxsum and indsum < minsum both impossible!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 11) xmax,zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx 
 
         indsum = i+j+k
!!$ indsum > maxsum and indsum < minsum both impossible!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 12) xmax,zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx 

         indsum = i+j+k
!!$ indsum < minsum impossible!
         if(indsum.gt.maxsum) then
            inew = nx - stenp1 + 1
            jnew = ny - stenp1 + 1
            knew = nz - stenp1 + 1
         else
            call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         end if
                  
!!$  Don't poison edges!
         dens(i,j,k) = dens(inew,jnew,knew)

         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

!!$ Loop over the cubical domain 8 corners:

  !!$ 1) xmin,ymin,zmin:
  do k = 1, sten
    do j = 1, sten
      do i = 1, sten
 
!!$ indsum < minsum guaranteed!

         inew = stenp1
         jnew = stenp1
         knew = stenp1

!!$  Don't poison corners!
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 2) xmin,ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten
      do i = 1, sten
 
         indsum = i+j+k
!!$ indsum < maxsum and indsum > minsum guaranteed!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 3) xmin,ymax,zmin:
  do k = 1, sten
    do j = ny-sten+1, ny 
      do i = 1, sten
 
         indsum = i+j+k
!!$ indsum < maxsum and indsum > minsum guaranteed!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 4) xmin,ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny
      do i = 1, sten
 
         indsum = i+j+k
!!$ indsum < maxsum and indsum > minsum guaranteed!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 5) xmax,ymin,zmin:
  do k = 1, sten
    do j = 1, sten
      do i = nx-sten+1, nx
 
         indsum = i+j+k
!!$ indsum < maxsum and indsum > minsum guaranteed!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 6) xmax,ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten
      do i = nx-sten+1, nx
 
         indsum = i+j+k
!!$ indsum < maxsum and indsum > minsum guaranteed!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 7) xmax,ymax,zmin:
  do k = 1, sten
    do j = ny-sten+1, ny 
      do i = nx-sten+1, nx
 
         indsum = i+j+k
!!$ indsum < maxsum and indsum > minsum guaranteed!
         call find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)
         dens(i,j,k) = dens(inew,jnew,knew)         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

  !!$ 8) xmax,ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny 
      do i = nx-sten+1, nx
 
!!$ indsum > maxsum guaranteed!

         inew = nx - stenp1 + 1
         jnew = ny - stenp1 + 1
         knew = nz - stenp1 + 1

!!$  Don't poison corners!
         dens(i,j,k) = dens(inew,jnew,knew)
         
         sx(i,j,k) = sx(inew,jnew,knew)
         sy(i,j,k) = sy(inew,jnew,knew)
         sz(i,j,k) = sz(inew,jnew,knew)
         tau(i,j,k) = tau(inew,jnew,knew)
         Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
         Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
         Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
         if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,knew)
         endif

      end do 
    end do
  end do 

!!$ Synchronize ghost zones (first time):

  call CCTK_GroupIndex(gindex_tau,"GRHydro::tau")
  call CCTK_GroupIndex(gindex_dens,"GRHydro::dens")
  call CCTK_GroupIndex(gindex_scon,"GRHydro::scon")
  call CCTK_GroupIndex(gindex_Bcons,"GRHydro::Bcons")

  groups(1) = gindex_tau
  groups(2) = gindex_dens
  groups(3) = gindex_scon
  groups(4) = gindex_Bcons

  call CCTK_SyncGroupsI(status,cctkGH,num_groups,groups)
  if(status.lt.0) then
    write(error_message,'(a33,i3,a5,i3,a7)' ) &
    'CCTK_SyncGroupsI returned status ',status,' for ',num_groups,'groups.'
    call CCTK_WARN(CCTK_WARN_ABORT,error_message)
  endif
  
  if(clean_divergence.ne.0) then
    call CCTK_SyncGroup(status,cctkGH,"GRHydro::psidc")
    if(status.lt.0) then
      write(error_message,'(a31,i3,a19)' ) &
      'CCTK_SyncGroup returned status ',status,' for GRHydro::psidc'
      call CCTK_WARN(CCTK_WARN_ABORT,error_message)
    endif
  endif 


!!$ Fix the indsum < minsum area in one go?
  do k=1,minsum-3
     do j=1,minsum-k-2
        do i=1,minsum-k-j-1

           indsum=i+j+k
           if(dens(i,j,k).lt.0.0d0) then
              if(dens(j,i,k).gt.0.0d0) then
                 inew=j
                 jnew=i
                 knew=k
              else if(dens(k,j,i).gt.0.0d0) then
                 inew=k
                 jnew=j
                 knew=i
              else if(dens(i,k,j).gt.0.0d0) then
                 inew=i
                 jnew=k
                 knew=j
              else
                 inew=stenp1
                 jnew=stenp1
                 knew=stenp1
              endif

              dens(i,j,k) = dens(inew,jnew,knew)
              sx(i,j,k) = sx(inew,jnew,knew)
              sy(i,j,k) = sy(inew,jnew,knew)
              sz(i,j,k) = sz(inew,jnew,knew)
              tau(i,j,k) = tau(inew,jnew,knew)
              Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
              Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
              Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
              if(clean_divergence.ne.0) then 
                 psidc(i,j,k)=psidc(inew,jnew,knew)
              endif
           end if

        end do
     end do
  end do

!!$ Fix the indsum > maxsum area in one go?
  do k=nz-3*stenp1+4,nz
     do j=ny-3*stenp1+4+(nz-k),ny
        do i=nx-3*stenp1+4+(nz-k)+(ny-j),nx

           indsum=i+j+k
           xoff=nx-i
           yoff=ny-j
           zoff=nz-k
           if(dens(i,j,k).lt.0.0d0) then

              if(dens(nx-yoff,ny-xoff,k).gt.0.0d0) then
                 inew=nx-yoff
                 jnew=ny-xoff
                 knew=k
              else if(dens(nx-zoff,j,nz-xoff).gt.0.0d0) then
                 inew=nx-zoff
                 jnew=j
                 knew=nz-xoff
              else if(dens(i,ny-zoff,nz-yoff).gt.0.0d0) then
                 inew=i
                 jnew=ny-zoff
                 knew=nz-xoff
              else
                 inew=nx - stenp1 + 1
                 jnew=ny - stenp1 + 1
                 knew=nz - stenp1 + 1
              endif
              
              dens(i,j,k) = dens(inew,jnew,knew)
              sx(i,j,k) = sx(inew,jnew,knew)
              sy(i,j,k) = sy(inew,jnew,knew)
              sz(i,j,k) = sz(inew,jnew,knew)
              tau(i,j,k) = tau(inew,jnew,knew)
              Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
              Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
              Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
              if(clean_divergence.ne.0) then 
                 psidc(i,j,k)=psidc(inew,jnew,knew)
              endif
           end if
           
        end do
     end do
  end do
  

!!$ Fix the ymin face in the region indsum < minsum
!!$  do k = stenp1, nz-sten
!!$    do j = 1, sten
!!$      do i =  stenp1, nx-sten
!!$
!!$        indsum = i+j
!!$        if(indsum.lt.minsum.and.dens(i,j,k).lt.0.0d0) then
!!$          inew = j
!!$          jnew = i
!!$          if(dens(inew,jnew,k).lt.0.0d0) then
!!$            inew = stenp1
!!$            jnew = stenp1
!!$          endif 
!!$           
!!$          dens(i,j,k) = dens(inew,jnew,k)
!!$          sx(i,j,k) = sx(inew,jnew,k)
!!$          sy(i,j,k) = sy(inew,jnew,k)
!!$          sz(i,j,k) = sz(inew,jnew,k)
!!$          tau(i,j,k) = tau(inew,jnew,k)
!!$          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
!!$          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
!!$          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
!!$          if(clean_divergence.ne.0) then 
!!$            psidc(i,j,k)=psidc(inew,jnew,k)
!!$          endif
!!$        end if
!!$
!!$      end do 
!!$    end do
!!$  end do 

!!$ Fix the xmin face in the region indsum < minsum
!!$  do k = stenp1, nz-sten
!!$    do j = stenp1, ny-sten
!!$      do i = 1, sten
!!$
!!$        indsum = i+j
!!$        if(indsum.lt.minsum.and.dens(i,j,k).lt.0.0d0) then
!!$          inew = j
!!$          jnew = i
!!$          if(dens(inew,jnew,k).lt.0.0d0) then
!!$            inew = stenp1
!!$            jnew = stenp1
!!$          endif 
!!$           
!!$          dens(i,j,k) = dens(inew,jnew,k)
!!$          sx(i,j,k) = sx(inew,jnew,k)
!!$          sy(i,j,k) = sy(inew,jnew,k)
!!$          sz(i,j,k) = sz(inew,jnew,k)
!!$          tau(i,j,k) = tau(inew,jnew,k)
!!$          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
!!$          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
!!$          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
!!$          if(clean_divergence.ne.0) then 
!!$            psidc(i,j,k)=psidc(inew,jnew,k)
!!$          endif
!!$        end if
!!$
!!$
!!$      end do 
!!$    end do
!!$  end do 

!!$ Fix the ymax face in the region indsum > maxsum
!!$  do k = stenp1, nz-sten
!!$    do j = ny-sten+1, ny
!!$      do i = stenp1, nx-sten
!!$
!!$        indsum = i+j
!!$        if(indsum.gt.maxsum.and.dens(i,j,k).lt.0.0d0) then
!!$          inew = j-ny+nx
!!$          jnew = i-nx+ny
!!$          if(dens(inew,jnew,k).lt.0.0d0) then
!!$            inew = nx - stenp1
!!$            jnew = ny - stenp1
!!$          endif 
!!$           
!!$          dens(i,j,k) = dens(inew,jnew,k)
!!$          sx(i,j,k) = sx(inew,jnew,k)
!!$          sy(i,j,k) = sy(inew,jnew,k)
!!$          sz(i,j,k) = sz(inew,jnew,k)
!!$          tau(i,j,k) = tau(inew,jnew,k)
!!$          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
!!$          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
!!$          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
!!$          if(clean_divergence.ne.0) then 
!!$            psidc(i,j,k)=psidc(inew,jnew,k)
!!$          endif
!!$        endif
!!$
!!$      end do 
!!$    end do
!!$  end do 

!!$ Fix the xmax face in the region indsum > maxsum
!!$  do k = stenp1, nz-sten
!!$    do j = stenp1, ny-sten
!!$      do i = nx-sten+1, nx
!!$
!!$        indsum = i+j
!!$        if(indsum.gt.maxsum.and.dens(i,j,k).lt.0.0d0) then
!!$          inew = j-ny+nx
!!$          jnew = i-nx+ny
!!$          if(dens(inew,jnew,k).lt.0.0d0) then
!!$            inew = nx - stenp1
!!$            jnew = ny - stenp1
!!$          endif 
!!$           
!!$          dens(i,j,k) = dens(inew,jnew,k)
!!$          sx(i,j,k) = sx(inew,jnew,k)
!!$          sy(i,j,k) = sy(inew,jnew,k)
!!$          sz(i,j,k) = sz(inew,jnew,k)
!!$          tau(i,j,k) = tau(inew,jnew,k)
!!$          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
!!$          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
!!$          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
!!$          if(clean_divergence.ne.0) then 
!!$            psidc(i,j,k)=psidc(inew,jnew,k)
!!$          endif
!!$        endif
!!$
!!$      end do 
!!$    end do
!!$  end do 

!!$ Synchronize ghost zones (second time):

  call CCTK_SyncGroupsI(status,cctkGH,num_groups,groups)
  if(status.lt.0) then
    write(error_message,'(a33,i3,a5,i3,a7)' ) &
    'CCTK_SyncGroupsI returned status ',status,' for ',num_groups,'groups.'
    call CCTK_WARN(CCTK_WARN_ABORT,error_message)
  endif
  
  if(clean_divergence.ne.0) then
    call CCTK_SyncGroup(status,cctkGH,"GRHydro::psidc")
    if(status.lt.0) then
      write(error_message,'(a31,i3,a19)' ) &
      'CCTK_SyncGroup returned status ',status,' for GRHydro::psidc'
      call CCTK_WARN(CCTK_WARN_ABORT,error_message)
    endif
  endif 

!!$  do k=1,nz
!!$     if(k.lt.stenp1)then
!!$              zoff=k-stenp1
!!$           else if(k.gt.nz-stenp1+1) then
!!$              zoff=k-(nz-stenp1+1)
!!$           else
!!$              zoff=0
!!$           endif
!!$
!!$     do j=1,ny
!!$        if(j.lt.stenp1)then
!!$           yoff=j-stenp1
!!$        else if(j.gt.ny-stenp1+1) then
!!$           yoff=j-(ny-stenp1+1)
!!$        else
!!$           yoff=0
!!$        endif
!!$
!!$        do i=1,nx
!!$           if(i.lt.stenp1) then
!!$              xoff=i-stenp1
!!$           else if(i.gt.nx-stenp1+1) then
!!$              xoff=i-(nx-stenp1+1)
!!$           else
!!$              xoff=0
!!$           endif
!!$
!!$           indsum = i+j+k
!!$           
!!$           if( (xoff.ne.0.or.yoff.ne.0.or.zoff.ne.0) .and. &
!!$                indsum.ge.minsum.and.indsum.le.maxsum) then
!!$    We can map the point to the interior diagonal, orthogonal to the shock.
!!$
!!$              inew=indsum/3
!!$              jnew=(indsum-inew)/2
!!$              knew=indsum-inew-jnew
!!$              
!!$              dens(i,j,k) = dens(inew,jnew,knew)
!!$              sx(i,j,k) = sx(inew,jnew,knew)
!!$              sy(i,j,k) = sy(inew,jnew,knew)
!!$              sz(i,j,k) = sz(inew,jnew,knew)
!!$              tau(i,j,k) = tau(inew,jnew,knew)
!!$              Bconsx(i,j,k)=Bconsx(inew,jnew,knew)
!!$              Bconsy(i,j,k)=Bconsy(inew,jnew,knew)
!!$              Bconsz(i,j,k)=Bconsz(inew,jnew,knew)
!!$             if(clean_divergence.ne.0) then 
!!$               psidc(i,j,k)=psidc(inew,jnew,knew)
!!$             endif
!!$
!!$           endif
!!$           
!!$        enddo
!!$     enddo
!!$  enddo
 
end subroutine GRHydro_Diagshock_BoundaryM

subroutine GRHydro_Diagshock2D_BoundaryM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, kc, nx, ny, nz, sten, stenp1, minsum, maxsum, inew, jnew, knew
  CCTK_INT :: xoff,yoff,zoff,indsum,numoff
  CCTK_REAL :: det
  CCTK_INT :: status
  character(len=200) :: error_message
  integer :: gindex_tau, gindex_dens, gindex_scon, gindex_Bcons
  CCTK_INT, parameter :: num_groups = 4
  CCTK_INT :: groups(num_groups)
  
  sten=GRHydro_stencil
  stenp1=GRHydro_stencil + 1
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  minsum = 2*stenp1
  maxsum = nx+ny-2*sten
  
  xoff=0
  yoff=0
  zoff=0

!!$ Loop over the cubical domain 6 faces:

  !!$ 1) xmin:
  do k = stenp1, nz-sten
    do j = stenp1, ny-sten
      do i = 1, sten

        xoff=i-stenp1
        yoff=0
        zoff=0
        indsum = i+j
        inew = i - xoff
        jnew = j + xoff
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - sten
          jnew = ny - sten
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
          dens(i,j,k) = dens(inew,jnew,k)
        end if
           
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 2) xmax:
  do k = stenp1, nz-sten
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx

        xoff=i-(nx-stenp1+1)
        yoff=0
        zoff=0
        indsum = i+j
        inew = i - xoff
        jnew = j + xoff
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - sten
          jnew = ny - sten
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
          dens(i,j,k) = dens(inew,jnew,k)
        end if
           
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 3) ymin:
  do k = stenp1, nz-sten
    do j = 1, sten
      do i =  stenp1, nx-sten

        xoff=0
        yoff=j-stenp1
        zoff=0
        indsum = i+j
        inew = i + yoff
        jnew = j - yoff
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - sten
          jnew = ny - sten
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
          dens(i,j,k) = dens(inew,jnew,k)
        end if
           
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 4) ymax:
  do k = stenp1, nz-sten
    do j = ny-sten+1, ny
      do i = stenp1, nx-sten

        xoff=0
        yoff=j-(ny-stenp1+1)
        zoff=0
        indsum = i+j
        inew = i + yoff
        jnew = j - yoff
        if(indsum.lt.minsum) then
          inew = stenp1
          jnew = stenp1
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else if(indsum.gt.maxsum) then
          inew = nx - sten
          jnew = ny - sten
          dens(i,j,k) = -1.0d0 !!$ poison at shock propagation direction corners
        else
          dens(i,j,k) = dens(inew,jnew,k)
        end if
           
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 5) zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 6) zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

!!$ Loop over the cubical domain 12 edges:

  !!$ 1) xmin,ymin:
  do k = stenp1, nz-sten
    do j = 1, sten
      do i = 1, sten

        !!$ We already know that for this edge indsum.lt.minsum.
        inew = stenp1
        jnew = stenp1
           
        dens(i,j,k) = dens(inew,jnew,k)
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 2) xmin,ymax:
  do k = stenp1, nz-sten
    do j =  ny-sten+1, ny
      do i = 1, sten

        xoff=i-stenp1
        yoff=j-(ny-stenp1+1)
        zoff=0
 
        numoff = max(abs(xoff),abs(yoff))
        if( abs(xoff).eq.numoff) then
          inew = i + numoff
          jnew = j - numoff
        else
          jnew = j - numoff
          inew = i + numoff
        endif  

        dens(i,j,k) = dens(inew,jnew,k)
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 3) xmax,ymin:
  do k = stenp1, nz-sten
    do j = 1, sten
      do i = nx-sten+1, nx 

        xoff=i-(nx-stenp1+1)
        yoff=j-stenp1
        zoff=0
        numoff = max(abs(xoff),abs(yoff))
        if( abs(xoff).eq.numoff) then
          inew = i - numoff
          jnew = j + numoff
        else
          jnew = j + numoff
          inew = i - numoff
        endif  
           
        dens(i,j,k) = dens(inew,jnew,k)
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 4) xmax,ymax:
  do k = stenp1, nz-sten
    do j = ny-sten+1, ny 
      do i = nx-sten+1, nx 

        !!$ We already know that for this edge indsum.gt.maxsum
        inew = nx - sten
        jnew = ny - sten
          
        dens(i,j,k) = dens(inew,jnew,k)
        sx(i,j,k) = sx(inew,jnew,k)
        sy(i,j,k) = sy(inew,jnew,k)
        sz(i,j,k) = sz(inew,jnew,k)
        tau(i,j,k) = tau(inew,jnew,k)
        Bconsx(i,j,k)=Bconsx(inew,jnew,k)
        Bconsy(i,j,k)=Bconsy(inew,jnew,k)
        Bconsz(i,j,k)=Bconsz(inew,jnew,k)
        if(clean_divergence.ne.0) then 
          psidc(i,j,k)=psidc(inew,jnew,k)
        endif

      end do 
    end do
  end do 

  !!$ 5) ymin,zmin:
  do k = 1, sten 
    do j = 1, sten 
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 6) ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten 
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 7) ymax,zmin:
  do k = 1, sten 
    do j = ny-sten+1, ny
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 8) ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 9) xmin,zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten 
      do i = 1, sten 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 10) xmin,zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = 1, sten 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 11) xmax,zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 12) xmax,zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

!!$ Loop over the cubical domain 8 corners:

  !!$ 1) xmin,ymin,zmin:
  do k = 1, sten
    do j = 1, sten
      do i = 1, sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 2) xmin,ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten
      do i = 1, sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 3) xmin,ymax,zmin:
  do k = 1, sten
    do j = ny-sten+1, ny 
      do i = 1, sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 4) xmin,ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny
      do i = 1, sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 5) xmax,ymin,zmin:
  do k = 1, sten
    do j = 1, sten
      do i = nx-sten+1, nx
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 6) xmax,ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten
      do i = nx-sten+1, nx
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 7) xmax,ymax,zmin:
  do k = 1, sten
    do j = ny-sten+1, ny 
      do i = nx-sten+1, nx
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 8) xmax,ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny 
      do i = nx-sten+1, nx
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

!!$ Synchronize ghost zones (first time):

  call CCTK_GroupIndex(gindex_tau,"GRHydro::tau")
  call CCTK_GroupIndex(gindex_dens,"GRHydro::dens")
  call CCTK_GroupIndex(gindex_scon,"GRHydro::scon")
  call CCTK_GroupIndex(gindex_Bcons,"GRHydro::Bcons")

  groups(1) = gindex_tau
  groups(2) = gindex_dens
  groups(3) = gindex_scon
  groups(4) = gindex_Bcons

  call CCTK_SyncGroupsI(status,cctkGH,num_groups,groups)
  if(status.lt.0) then
    write(error_message,'(a33,i3,a5,i3,a7)' ) &
    'CCTK_SyncGroupsI returned status ',status,' for ',num_groups,'groups.'
    call CCTK_WARN(CCTK_WARN_ABORT,error_message)
  endif
  
  if(clean_divergence.ne.0) then
    call CCTK_SyncGroup(status,cctkGH,"GRHydro::psidc")
    if(status.lt.0) then
      write(error_message,'(a31,i3,a19)' ) &
      'CCTK_SyncGroup returned status ',status,' for GRHydro::psidc'
      call CCTK_WARN(CCTK_WARN_ABORT,error_message)
    endif
  endif 

!!$ Fix the ymin face in the region indsum < minsum
  do k = stenp1, nz-sten
    do j = 1, sten
      do i =  stenp1, nx-sten

        indsum = i+j
        if(indsum.lt.minsum.and.dens(i,j,k).lt.0.0d0) then
          inew = j
          jnew = i
          if(dens(inew,jnew,k).lt.0.0d0) then
            inew = stenp1
            jnew = stenp1
          endif 
           
          dens(i,j,k) = dens(inew,jnew,k)
          sx(i,j,k) = sx(inew,jnew,k)
          sy(i,j,k) = sy(inew,jnew,k)
          sz(i,j,k) = sz(inew,jnew,k)
          tau(i,j,k) = tau(inew,jnew,k)
          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
          if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,k)
          endif
        end if

      end do 
    end do
  end do 

!!$ Fix the xmin face in the region indsum < minsum
  do k = stenp1, nz-sten
    do j = stenp1, ny-sten
      do i = 1, sten

        indsum = i+j
        if(indsum.lt.minsum.and.dens(i,j,k).lt.0.0d0) then
          inew = j
          jnew = i
          if(dens(inew,jnew,k).lt.0.0d0) then
            inew = stenp1
            jnew = stenp1
          endif 
           
          dens(i,j,k) = dens(inew,jnew,k)
          sx(i,j,k) = sx(inew,jnew,k)
          sy(i,j,k) = sy(inew,jnew,k)
          sz(i,j,k) = sz(inew,jnew,k)
          tau(i,j,k) = tau(inew,jnew,k)
          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
          if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,k)
          endif
        end if


      end do 
    end do
  end do 

!!$ Fix the ymax face in the region indsum > maxsum
  do k = stenp1, nz-sten
    do j = ny-sten+1, ny
      do i = stenp1, nx-sten

        indsum = i+j
        if(indsum.gt.maxsum.and.dens(i,j,k).lt.0.0d0) then
          inew = j-ny+nx
          jnew = i-nx+ny
          if(dens(inew,jnew,k).lt.0.0d0) then
            inew = nx - sten
            jnew = ny - sten
          endif 
           
          dens(i,j,k) = dens(inew,jnew,k)
          sx(i,j,k) = sx(inew,jnew,k)
          sy(i,j,k) = sy(inew,jnew,k)
          sz(i,j,k) = sz(inew,jnew,k)
          tau(i,j,k) = tau(inew,jnew,k)
          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
          if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,k)
          endif
        endif

      end do 
    end do
  end do 

!!$ Fix the xmax face in the region indsum > maxsum
  do k = stenp1, nz-sten
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx

        indsum = i+j
        if(indsum.gt.maxsum.and.dens(i,j,k).lt.0.0d0) then
          inew = j-ny+nx
          jnew = i-nx+ny
          if(dens(inew,jnew,k).lt.0.0d0) then
            inew = nx - sten
            jnew = ny - sten
          endif 
           
          dens(i,j,k) = dens(inew,jnew,k)
          sx(i,j,k) = sx(inew,jnew,k)
          sy(i,j,k) = sy(inew,jnew,k)
          sz(i,j,k) = sz(inew,jnew,k)
          tau(i,j,k) = tau(inew,jnew,k)
          Bconsx(i,j,k)=Bconsx(inew,jnew,k)
          Bconsy(i,j,k)=Bconsy(inew,jnew,k)
          Bconsz(i,j,k)=Bconsz(inew,jnew,k)
          if(clean_divergence.ne.0) then 
            psidc(i,j,k)=psidc(inew,jnew,k)
          endif
        endif

      end do 
    end do
  end do 

!!$ Fix the upper and lower (5-12 below) domain edges since 
!!$ they still contain poison got from the faces in the initial
!!$ step: 

  !!$ 5) ymin,zmin:
  do k = 1, sten 
    do j = 1, sten 
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 6) ymin,zmax:
  do k = nz-sten+1, nz
    do j = 1, sten 
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 7) ymax,zmin:
  do k = 1, sten 
    do j = ny-sten+1, ny
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 8) ymax,zmax:
  do k = nz-sten+1, nz
    do j = ny-sten+1, ny
      do i =  stenp1, nx-sten
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 9) xmin,zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten 
      do i = 1, sten 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 10) xmin,zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = 1, sten 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 11) xmax,zmin:
  do k = 1, sten 
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

  !!$ 12) xmax,zmax:
  do k = nz-sten+1, nz
    do j = stenp1, ny-sten
      do i = nx-sten+1, nx 
 
        kc = nz/2

        dens(i,j,k) = dens(i,j,kc)
        sx(i,j,k) = sx(i,j,kc)
        sy(i,j,k) = sy(i,j,kc)
        sz(i,j,k) = sz(i,j,kc)
        tau(i,j,k) = tau(i,j,kc)
        Bconsx(i,j,k)=Bconsx(i,j,kc)
        Bconsy(i,j,k)=Bconsy(i,j,kc)
        Bconsz(i,j,k)=Bconsz(i,j,kc)
        if(clean_divergence.ne.0) then
          psidc(i,j,k)=psidc(i,j,kc)
        endif

      end do 
    end do
  end do 

!!$ Synchronize ghost zones (second time):

  call CCTK_SyncGroupsI(status,cctkGH,num_groups,groups)
  if(status.lt.0) then
    write(error_message,'(a33,i3,a5,i3,a7)' ) &
    'CCTK_SyncGroupsI returned status ',status,' for ',num_groups,'groups.'
    call CCTK_WARN(CCTK_WARN_ABORT,error_message)
  endif
  
  if(clean_divergence.ne.0) then
    call CCTK_SyncGroup(status,cctkGH,"GRHydro::psidc")
    if(status.lt.0) then
      write(error_message,'(a31,i3,a19)' ) &
      'CCTK_SyncGroup returned status ',status,' for GRHydro::psidc'
      call CCTK_WARN(CCTK_WARN_ABORT,error_message)
    endif
  endif 
        
!!$  do k=1,nz
!!$     if(k.lt.stenp1)then
!!$        zoff=k-stenp1
!!$     else if(k.gt.nz-stenp1+1) then
!!$        zoff=k-(nz-stenp1+1)
!!$     else
!!$        zoff=0
!!$     endif
!!$
!!$     do j=1,ny
!!$        if(j.lt.stenp1)then
!!$           yoff=j-stenp1
!!$        else if(j.gt.ny-stenp1+1) then
!!$           yoff=j-(ny-stenp1+1)
!!$        else
!!$           yoff=0
!!$        endif
!!$
!!$        do i=1,nx
!!$           if(i.lt.stenp1) then
!!$              xoff=i-stenp1
!!$           else if(i.gt.nx-stenp1+1) then
!!$              xoff=i-(nx-stenp1+1)
!!$           else
!!$              xoff=0
!!$           endif
!!$  
!!$        indsum = i+j
!!$        
!!$        if(xoff.ne.0.or.yoff.ne.0) then
!!$
!!$          if(indsum.ge.minsum.and.indsum.le.maxsum) then
!!$
!!$           numoff = max(abs(xoff),abs(yoff))
!!$           if( abs(xoff).eq.numoff) then
!!$             if(xoff.gt.0) then
!!$               inew = i - numoff
!!$               jnew = j + numoff
!!$             else
!!$               inew = i + numoff
!!$               jnew = j - numoff
!!$             endif
!!$           else
!!$             if(yoff.gt.0) then
!!$               jnew = j - numoff
!!$               inew = i + numoff
!!$             else
!!$               jnew = j + numoff
!!$               inew = i - numoff
!!$             endif
!!$           endif  
!!$!!$ Note that the following two conditions helps but
!!$!!$ don't solve the problem when running in several
!!$!!$ processors
!!$         else if(indsum.lt.minsum) then
!!$           inew = stenp1
!!$           jnew = stenp1
!!$
!!$         else if(indsum.gt.maxsum) then
!!$           inew = nx - stenp1
!!$           jnew = ny - stenp1
!!$         end if
!!$           
!!$           dens(i,j,k) = dens(inew,jnew,k)
!!$           sx(i,j,k) = sx(inew,jnew,k)
!!$           sy(i,j,k) = sy(inew,jnew,k)
!!$           sz(i,j,k) = sz(inew,jnew,k)
!!$           tau(i,j,k) = tau(inew,jnew,k)
!!$           Bconsx(i,j,k)=Bconsx(inew,jnew,k)
!!$           Bconsy(i,j,k)=Bconsy(inew,jnew,k)
!!$           Bconsz(i,j,k)=Bconsz(inew,jnew,k)
!!$           if(clean_divergence.ne.0) then 
!!$              psidc(i,j,k)=psidc(inew,jnew,k)
!!$           endif
!!$
!!$        else if( zoff.ne.0) then
!!$ 
!!$           kc = nz/2
!!$
!!$           dens(i,j,k) = dens(i,j,kc)
!!$           sx(i,j,k) = sx(i,j,kc)
!!$           sy(i,j,k) = sy(i,j,kc)
!!$           sz(i,j,k) = sz(i,j,kc)
!!$           tau(i,j,k) = tau(i,j,kc)
!!$           Bconsx(i,j,k)=Bconsx(i,j,kc)
!!$           Bconsy(i,j,k)=Bconsy(i,j,kc)
!!$           Bconsz(i,j,k)=Bconsz(i,j,kc)
!!$           if(clean_divergence.ne.0) then
!!$              psidc(i,j,k)=psidc(i,j,kc)
!!$           endif
!!$           
!!$        endif
!!$        
!!$        enddo
!!$     enddo
!!$  enddo
  
end subroutine GRHydro_Diagshock2D_BoundaryM

subroutine find_mapped_point(i,j,k,nx,ny,nz,stenp1,inew,jnew,knew)

!!$ Find the nearest mapped point in the interior, assuming it exists

  implicit none
  CCTK_INT :: i,j,k,n,nx,ny,nz,stenp1
  CCTK_INT :: xoff,yoff,zoff,indsum,minsum,maxsum
  CCTK_INT :: inew,jnew,knew,dspace,newsum
  
  minsum = 3*stenp1
  maxsum = nx+ny+nz - 3*(stenp1-1)

  indsum = i+j+k

  if(indsum.lt.minsum.or.indsum.gt.maxsum)write(6,*)'HUGE ERROR - ILLEGAL INDEX SUM!!!!!!'

  if(i.lt.stenp1)then
     xoff = i - stenp1 
  else if(i.gt.nx-stenp1+1) then
     xoff=i - (nx-stenp1+1)
  else
     xoff=0
  endif
  
  if(j.lt.stenp1)then
     yoff = j - stenp1 
  else if(j.gt.ny-stenp1+1) then
     yoff=j - (ny-stenp1+1)
  else
     yoff=0
  endif
     
  if(k.lt.stenp1)then
     zoff = k - stenp1 
  else if(k.gt.nz-stenp1+1) then
     zoff=k-(nz-stenp1+1)
  else
     zoff=0
  endif
  
  if(xoff.eq.0.and.yoff.eq.0.and.zoff.eq.0) then
     write(6,*)'Why did you try to map an interior point?'
     inew=i
     jnew=j
     knew=k
     return
  else
     
!!$ We have something to actually do!
     
     inew = i-xoff
     jnew = j-yoff
     knew = k-zoff
     dspace = xoff + yoff + zoff
     
     do n=1,abs(dspace)
        
        if(dspace.lt.0) then
           
!!$     dspace negative; point sum corrected upward
!!$     need to move other coords down
           
           if(inew.gt.stenp1) then
              inew = inew - 1
           else if(jnew.gt.stenp1) then
              jnew = jnew - 1
           else if (knew.gt.stenp1) then
              knew = knew - 1
           else
              write(6,*)'Big error!!!!'
           endif
           
        else if (dspace.gt.0) then
           
!!$     dspace positive; point sum corrected downward
!!$     need to move other coords up
           
           if(inew.lt.nx-stenp1+1) then
              inew=inew+1
           else if(jnew.lt.ny-stenp1+1) then
              jnew=jnew+1
           else if(knew.lt.nz-stenp1+1) then
              knew=knew+1
           else
              write(6,*)'Big error2!!!!'
           endif

        end if

     end do
     
     newsum = inew+jnew+knew
     
     if(inew.lt.stenp1.or.inew.gt.(nx-(stenp1-1)).or. &
          jnew.lt.stenp1.or.jnew.gt.(ny-(stenp1-1)).or. &
          knew.lt.stenp1.or.knew.gt.(nz-(stenp1-1)).or. &
          newsum.ne.indsum)write(6,*)'Big error 3!!!!'
     
  end if

  return

end subroutine find_mapped_point

