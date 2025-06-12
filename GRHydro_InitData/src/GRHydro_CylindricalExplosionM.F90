 /*@@
   @file      GRHydro_CylindricalExplosionM.F90
   @date      Apr 22, 2011
   @author    Scott Noble, Joshua Faber, Bruno Mundim
   @desc 
   Cylindrical magnetized shocks.
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


 /*@@
   @routine    GRHydro_cylindricalexplosionM
   @date       Fri Apr 22 14:51:37 EDT 2011
   @author     Scott Noble, Joshua Faber, Bruno Mundim
   @desc 
   Initial data for cylidrically symmetric shocks with 
   magnetic fields.  Meant to recreate tests demonstrated 
   in Komissarov (1999).
   @enddesc 
   @calls     
   @calledby   
   @history 
   Using GRHydro_ShockTubeM.F90 as a template.
   @endhistory 

@@*/

subroutine GRHydro_cylindricalexplosionM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det
  CCTK_REAL :: rhol, rhor, pressl, pressr
  CCTK_REAL :: bvcxl,bvcyl,bvczl
  CCTK_REAL :: tmp,tmp2,gam,r_inner,r_outer

  CCTK_REAL :: cyl_fr

!!$ Check that user selected proper magnetic field ID. Warn if wrong one is
!!$ selected and proceed to ignore the value
  if (.not. CCTK_EQUALS(initial_Bvec, "cylexp")) then
    call CCTK_WARN(1, "When using the cyclindrical explosion initial data, please also select 'cylexp' for initial_Bvec. I will proceed as if this was set")
  end if

!!$Uses Bx_init, By_init and Bz_init to set magnetic field strength
!!$Original tests had Bx_init = 0.1, 1.0   with By_init=Bz_init = 0

  bvcxl = Bx_init
  bvcyl = By_init
  bvczl = Bz_init

!!$Inner radius and outer radius (Komissarov's defaults)
!  r_inner = 8.d-1
!  r_outer = 1.d0
  r_inner = cyl_r_inner
  r_outer = cyl_r_outer

!!$Inner values
  rhor   = cyl_rho_inner
  pressr = cyl_press_inner

!!$Outer values
  rhol   = cyl_rho_outer
  pressl = cyl_press_outer

!!$Adiabatic index for test
  gam = gl_gamma

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  do i=1,nx
    do j=1,ny
      do k=1,nz

!!$direction represents the cylindrical radius here.
!!$TODO: maybe switch these over so that shocktube_tube choses the axis of the
!!$cylinder instead of the apparently random mapping
        if (CCTK_EQUALS(shocktube_type,"xshock")) then
           direction = sqrt((x(i,j,k)-shock_xpos)**2+&
                (y(i,j,k)-shock_ypos)**2)
        else if (CCTK_EQUALS(shocktube_type,"yshock")) then
           direction = sqrt((y(i,j,k)-shock_ypos)**2+&
                (z(i,j,k)-shock_zpos)**2)
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
           direction = sqrt((x(i,j,k)-shock_xpos)**2+&
                (z(i,j,k)-shock_zpos)**2)
        end if

        tmp = cyl_fr(direction,r_inner,r_outer,rhol,rhor)
        tmp2 = cyl_fr(direction,r_inner,r_outer,pressl,pressr)/( (gam - 1.d0) * tmp ) 
        rho(i,j,k) = tmp
        eps(i,j,k) = tmp2

        velx(i,j,k) = 0.d0
        vely(i,j,k) = 0.d0
        velz(i,j,k) = 0.d0
        Bvecx(i,j,k)=bvcxl
        Bvecy(i,j,k)=bvcyl
        Bvecz(i,j,k)=bvczl

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
  
end subroutine GRHydro_CylindricalExplosionM



function cyl_fr(R,r_inner,r_outer,min_f,max_f)

  implicit none

  CCTK_REAL :: cyl_fr
  CCTK_REAL :: R,r_inner,r_outer,min_f,max_f

  if(R .gt. r_outer) then 
     cyl_fr = min_f
  else if( R .lt. r_inner) then 
     cyl_fr = max_f
  else 
     cyl_fr = exp( ( log(max_f)*(r_outer - R) + log(min_f)*(R - r_inner) ) / (r_outer - r_inner) ) 
  end if

  return 

end function cyl_fr
