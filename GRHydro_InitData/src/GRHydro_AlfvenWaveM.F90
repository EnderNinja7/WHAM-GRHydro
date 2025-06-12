 /*@@
   @file      GRHydro_AlfvenWaveM.F90
   @date      Oct 10, 2011
   @author    Bruno Mundim, Joshua Faber, Scott Noble 
   @desc 
   Circularly Polarized Alfven Wave test as implemented by 
   Beckwith and Stone Astrophys.J.Suppl. 193 (2011) 6, arXiv:1101.3573, 
   and Del Zanna et. al. A&A 473, 11 (2007), arXiv:0704.3206;

   Other relevant references: 
     Stone et. al. Astrophys.J.Suppl. 178 (2008) 137, arXiv:0804.0402;  
     Gardiner and Stone JCP 227, 4123 (2008), arXiv:0712.2634;
     Toth JCP 161, 605 (2000);
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
   @routine    GRHydro_AlfvenWaveM
   @date       Oct 10, 2011
   @author     Bruno Mundim, Joshua Faber, Scott Noble
   @desc 
   Initial data for Circularly Polarized Alfven Wave test
   @enddesc 
   @calls     
   @calledby   
   @history 
   Using GRHydro_AdvectedLoopM.F90 as a template.
   @endhistory 

@@*/

subroutine GRHydro_AlfvenWaveM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT  :: i, j, k, nx, ny, nz
  CCTK_REAL :: pi,gam,AA,wnbr
  CCTK_REAL :: vxval,vyval,vzval, valf
  CCTK_REAL :: rhoval,pressval,epsval,hval
  CCTK_REAL :: Bxval, Byval, Bzval
  CCTK_REAL :: dx,dy,dz
  CCTK_REAL :: range_x,range_y,range_z,range_d
  CCTK_REAL :: cos_theta, sin_theta
  CCTK_REAL :: diaglength,xnew,vparallel,vperp,Bparallel,Bperp
  CCTK_REAL :: Bvecx_d, Bvecz_d
  CCTK_REAL :: velx_d, velz_d
  CCTK_REAL :: det
  CCTK_REAL :: t1,t2,t3

  pi=4.0d0*atan(1.0d0)

!!$Adiabatic index for this test:
  gam = (5.0d0/3.0d0)

!!$pressure, density, B^x, specific internal energy and enthalpy:
  rhoval   = 1.0d0
  pressval = alfvenwave_pressure
  Bxval    = 1.0d0
  epsval = pressval/(gam-1.0d0)/rhoval
  hval   = 1.0d0 + epsval + pressval/rhoval

!!$ DZ: rho=P=Bxval=AA = 1
!!$ Using DZ parameters, epsval=1.5, hval=3.5
  
!!$ Alfven Wave Amplitude:
  AA=1.0d0

!!$ Alfven wave speed:
  t1 = rhoval*hval+Bxval**2*(1.0d0+AA**2)
  t2 = 2.0d0*AA*Bxval**2/t1
  t3 = 0.5d0*(1.0d0+sqrt(1.0d0-t2**2))
  valf = sqrt(Bxval**2/t1/t3)

  write(*,*)'Alfven velocity:',valf

!!$ Using DZ parameters with P=1.0, we have:
!!$ t1=5.5
!!$ t2=4/11
!!$ t3=0.96577
!!$ valf=0.43389

!!$ Vx value:
   vxval=0.0d0
    
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

!!$ Alfven wave number
  
  do i=1,nx
     do j=1,ny
        do k=1,nz

           rho(i,j,k)=rhoval
           press(i,j,k)=pressval
           eps(i,j,k)=epsval
          
           if (CCTK_EQUALS(alfvenwave_type,"1D")) then

             wnbr = 2.0d0*pi/range_x

             velx(i,j,k)=vxval
             vely(i,j,k)=-valf*AA*cos(wnbr*x(i,j,k))
             velz(i,j,k)=-valf*AA*sin(wnbr*x(i,j,k))
             Bvecx(i,j,k)=Bxval
             Bvecy(i,j,k)=AA*Bxval*cos(wnbr*x(i,j,k))
             Bvecz(i,j,k)=AA*Bxval*sin(wnbr*x(i,j,k))

           else if (CCTK_EQUALS(alfvenwave_type,"2D")) then
      
             diaglength=range_x*range_y/range_d
             range_d = sqrt(range_x**2+range_y**2)
             cos_theta = range_y/range_d
             sin_theta = range_x/range_d
             wnbr = 2.0d0*pi/diaglength

             xnew = cos_theta*x(i,j,k)+sin_theta*y(i,j,k)

             vparallel=vxval
             vperp=-valf*AA*cos(wnbr*xnew)
             velx(i,j,k)=vparallel*cos_theta-vperp*sin_theta
             vely(i,j,k)=vparallel*sin_theta+vperp*cos_theta
             velz(i,j,k)=-valf*AA*sin(wnbr*xnew)
             Bparallel=Bxval
             Bperp=AA*Bxval*cos(wnbr*xnew)
             Bvecx(i,j,k)=Bparallel*cos_theta-Bperp*sin_theta
             Bvecy(i,j,k)=Bparallel*sin_theta+Bperp*cos_theta
             Bvecz(i,j,k)=AA*Bxval*sin(wnbr*xnew)

           else
             call CCTK_WARN(0,"Alfven wave case not recognized!")
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
  
end subroutine GRHydro_AlfvenWaveM


