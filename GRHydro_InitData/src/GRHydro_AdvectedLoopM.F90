 /*@@
   @file      GRHydro_AdvectedLoopM.F90
   @date      Aug 15, 2011
   @author    Scott Noble, Joshua Faber, Bruno Mundim
   @desc 
   Advected loop test as implemented by Beckwith and Stone Astrophys.J.Suppl. 
   193 (2011) 6,  arXiv:1101.3573. 


   Other relevant references: Devore JCP 92, 142 (1991), 
                              Toth and Odstrcil JCP 128,82 (1996),
                              Gardiner and Stone JCP 227, 4123 (2008);
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
   @routine    GRHydro_AdvectedLoopM
   @date       Aug 11, 2011
   @author     Scott Noble, Joshua Faber, Bruno Mundim
   @desc 
   Initial data for advected loop test
   @enddesc 
   @calls     
   @calledby   
   @history 
   Using GRHydro_ShockTubeM.F90 as a template.
   @endhistory 

@@*/

subroutine GRHydro_AdvectedLoopM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: det,radius,vxval,vyval,vzval,gam
  CCTK_REAL :: radius_iph, radius_imh 
  CCTK_REAL :: radius_jph, radius_jmh 
  CCTK_REAL :: radius_kph, radius_kmh 
  CCTK_REAL :: rhoval,pressval
  CCTK_REAL :: r_loop,A_loop,pi
  CCTK_REAL :: dx,dy,dz
  CCTK_REAL :: range_x,range_y,range_z,range_d
  CCTK_REAL :: cos_theta, sin_theta, tan_theta
  CCTK_REAL :: Bvecx_d, Bvecy_d, Bvecz_d
  CCTK_REAL :: x_d, y_d, z_d,diaglength
  CCTK_REAL :: dx_d, dy_d, dz_d, dx_x, dz_x

!!$Adiabatic index for test:
  gam = (5.0d0/3.0d0)
  
!!$radius of the loop: 
  r_loop = 0.3d0

!!$ stregth of the A-field
  A_loop=1.0d-3

!!$pressure and density:
  rhoval   = 1.0d0
  pressval = 3.0d0

  if (CCTK_EQUALS(advectedloop_type,"2D")) then


!!$ Vx, Vy and Vz values:
    if (CCTK_EQUALS(advectedloop_case,"V^z/=0")) then

!!$vxval=0.2d0/sqrt(6.0d0)
!!$ This new choice yields a crossing time of exactly t=24
!!$ assuming -1<x<1; -0.5<y<0.5
     vxval=1.d0/12.d0
     vyval=0.5d0*vxval
     vzval=vyval

    else if (CCTK_EQUALS(advectedloop_case,"V^z=0")) then
!!$   vxval=0.2d0/sqrt(6.0d0)
     vxval=1.d0/12.d0
     vyval=0.5d0*vxval
     vzval=0.0d0
    else
     call CCTK_WARN(0,"V^z component case not recognized!")
    end if 
    
  else if (CCTK_EQUALS(advectedloop_type,"3D")) then

     vxval=0.2d0*sqrt(2.0d0)
     vyval=0.2d0

    if (CCTK_EQUALS(advectedloop_case,"V^z/=0")) then

!!$vxval=0.2d0/sqrt(6.0d0)
!!$ This new choice yields a crossing time of exactly t=5
!!$ assuming -1<x<1; -0.5<y<0.5
     vzval=0.1

   else if (CCTK_EQUALS(advectedloop_case,"V^z=0")) then
     vzval=0.0d0
    else
     call CCTK_WARN(0,"V^z component case not recognized!")
    end if 

  endif
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

!!$ Note that the 3D test wasn't deviced to be used with AMR!
  range_x = (cctk_gsh(1)-2*cctk_nghostzones(1))*dx
  range_y = (cctk_gsh(2)-2*cctk_nghostzones(2))*dy
  range_z = (cctk_gsh(3)-2*cctk_nghostzones(3))*dz

  range_d = sqrt(range_z**2+range_x**2)
  diaglength = range_x*range_z/range_d

!!$ For 3-d case, the grid is going to be assumed to be a cube

  cos_theta = range_z/range_d
  sin_theta = range_x/range_d
  tan_theta = sin_theta/cos_theta
  
  do i=1,nx
     do j=1,ny
        do k=1,nz

           rho(i,j,k)=rhoval
           press(i,j,k)=pressval
           eps(i,j,k)=press(i,j,k)/(gam-1.0d0)/rho(i,j,k)

           if (CCTK_EQUALS(advectedloop_type,"2D")) then
             velx(i,j,k)=vxval
             vely(i,j,k)=vyval
             velz(i,j,k)=vzval
             Bvecz(i,j,k)=0.0d0

             radius = sqrt(x(i,j,k)**2+y(i,j,k)**2)
  
             if (CCTK_EQUALS(advectedloop_delA,"Exact")) then
    
               if(radius.le.r_loop) then
                  Bvecx(i,j,k)=-1.0d0*A_loop*y(i,j,k)/radius
                  Bvecy(i,j,k)=A_loop*x(i,j,k)/radius
               else
                  Bvecx(i,j,k)=0.0d0
                  Bvecy(i,j,k)=0.0d0
               endif

             else if (CCTK_EQUALS(advectedloop_delA,"Numeric")) then
  
               radius_iph = max(sqrt((x(i,j,k)+0.5d0*dx)**2+y(i,j,k)**2)-r_loop,0.d0)
               radius_imh = max(sqrt((x(i,j,k)-0.5d0*dx)**2+y(i,j,k)**2)-r_loop,0.d0)
               radius_jph = max(sqrt(x(i,j,k)**2+(y(i,j,k)+0.5d0*dy)**2)-r_loop,0.d0)
               radius_jmh = max(sqrt(x(i,j,k)**2+(y(i,j,k)-0.5d0*dy)**2)-r_loop,0.d0)

!!               if(radius.le.r_loop) then
                  Bvecx(i,j,k)=-1.0d0*A_loop*(radius_jph-radius_jmh)/dy
                  Bvecy(i,j,k)=A_loop*(radius_iph-radius_imh)/dx
!!               else
!!                  Bvecx(i,j,k)=0.0d0
!!                  Bvecy(i,j,k)=0.0d0
!!               endif

             else
               call CCTK_WARN(0,"A^b differentiation not recognized!")
             end if 

           else if (CCTK_EQUALS(advectedloop_type,"3D")) then


!!$  tangential velocity should be parallel to (1,1,1) plus a normal component (-1,0,1)

             velx(i,j,k)=cos_theta*vxval-sin_theta*vzval
             vely(i,j,k)=vyval
             velz(i,j,k)=cos_theta*vzval+sin_theta*vxval

             Bvecz_d=0.0d0

!!$ x_d = (x+z)/sqrt(2) => x_d=0 is equivalent to x+z=0

             x_d = cos_theta*x(i,j,k)+sin_theta*z(i,j,k)
             y_d = y(i,j,k)
             z_d = cos_theta*z(i,j,k)-sin_theta*x(i,j,k)
          
!!$ need to make x_d periodic!

             if(x_d.gt.1.5*diaglength) then
                x_d=x_d-2.0*diaglength
             else if (x_d.gt.0.5*diaglength .and. x_d.lt.1.5*diaglength) then
                x_d=x_d-diaglength
             else if(x_d.lt.-1.5*diaglength) then
                x_d=x_d+2.0*diaglength
             else if (x_d.lt.(-0.5*diaglength) .and. x_d.gt.(-1.5*diaglength)) then
                x_d=x_d+diaglength
             endif

             radius = sqrt(x_d**2+y_d**2)
  
             if (CCTK_EQUALS(advectedloop_delA,"Exact")) then
    
               if(radius.le.r_loop) then
                  Bvecx_d=-1.0d0*A_loop*y_d/radius
                  Bvecy_d=A_loop*x_d/radius
               else
                  Bvecx_d=0.0d0
                  Bvecy_d=0.0d0
               endif

               Bvecx(i,j,k)=cos_theta*Bvecx_d-sin_theta*Bvecz_d
               Bvecy(i,j,k)=Bvecy_d
               Bvecz(i,j,k)=cos_theta*Bvecz_d+sin_theta*Bvecx_d

             else if (CCTK_EQUALS(advectedloop_delA,"Numeric")) then

!!               dx_d = cos_theta*dx+sin_theta*dz
!!               dy_d = dy
!!               dz_d = cos_theta*dz-sin_theta*dx

!!  dx_d is the change in the rotated coords induced by a step in a direction over the Cartesian grid
                dx_x = cos_theta*dx
                dz_x = sin_theta*dz

!!  These are used for exact differencing
               radius_iph = max(sqrt((x_d+0.5d0*dx_x)**2+y_d**2)-r_loop,0.d0)
               radius_imh = max(sqrt((x_d-0.5d0*dx_x)**2+y_d**2)-r_loop,0.d0)
               radius_jph = max(sqrt(x_d**2+(y_d+0.5d0*dy)**2)-r_loop,0.d0)
               radius_jmh = max(sqrt(x_d**2+(y_d-0.5d0*dy)**2)-r_loop,0.d0)
               radius_kph = max(sqrt((x_d+0.5d0*dz_x)**2+y_d**2)-r_loop,0.d0)
               radius_kmh = max(sqrt((x_d-0.5d0*dz_x)**2+y_d**2)-r_loop,0.d0)

!! see notes
!!               if(radius.le.r_loop) then
                  Bvecx(i,j,k)=-1.0d0*A_loop*cos_theta*(radius_jph-radius_jmh)/dy
                  Bvecy(i,j,k)=A_loop*(sin_theta*(radius_kph-radius_kmh)/dz + &
                       cos_theta*(radius_iph-radius_imh)/dx)
                  Bvecz(i,j,k)=-1.0d0*A_loop*sin_theta*(radius_jph-radius_jmh)/dy
!!               else
!!                  Bvecx(i,j,k)=0.0d0
!!                  Bvecy(i,j,k)=0.0d0
!!                  Bvecz(i,j,k)=0.0d0
!!               endif

             else
               call CCTK_WARN(0,"A^b differentiation not recognized!")
             end if 

           else
             call CCTK_WARN(0,"Advected loop type not recognized!")
           end if 
           
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
  
end subroutine GRHydro_AdvectedLoopM


