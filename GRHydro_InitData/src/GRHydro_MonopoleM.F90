 /*@@
   @file      GRHydro_MonopoleM.F90
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

 /*@@
   @routine    GRHydro_MonopoleM
   @date       Sat Jan 26 02:53:49 2002
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   A monopole in space
   @enddesc 
   @calls     
   @calledby   
   @history 
   Expansion and alteration of the test code from GRAstro_Hydro, 
   written by Mark Miller.
   @endhistory 

@@*/

subroutine GRHydro_MonopoleM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr
  CCTK_REAL :: bvcxl,bvcyl,bvczl,bvcxr,bvcyr,bvczr
  CCTK_REAL :: ux,uy,uz,ut,tmp,tmp2,tmp3
  CCTK_REAL :: rr2,rg2
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  do i=1,nx
     do j=1,ny
        do k=1,nz
           
           rho(i,j,k) = 1.0
           velx(i,j,k) = 0.0
           vely(i,j,k) = 0.0
           velz(i,j,k) = 0.0
           eps(i,j,k) = 0.1

           Bvecx(i,j,k)=0.0
           Bvecy(i,j,k)=0.0
           Bvecz(i,j,k)=0.0

           if(CCTK_EQUALS(monopole_type,"Point")) then
              if(i.eq.nx/2+1.and.j.eq.ny/2+1.and.k.eq.nz/2+1) then
                 Bvecx(i,j,k)=Monopole_point_Bx
              else
                 Bvecx(i,j,k)=0.0
              endif
           else if(CCTK_EQUALS(monopole_type,"Gauss")) then
              rr2=x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2
              rg2=R_Gauss*R_Gauss
              if(rr2.lt.rg2) then
                 Bvecx(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
              else
                 Bvecx(i,j,k) = 0.0
              endif
           else if(CCTK_EQUALS(monopole_type,"1dalt")) then
              rr2=x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2
              rg2=R_Gauss*R_Gauss
              if(rr2.lt.rg2) then
                 Bvecx(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
              else
                 Bvecx(i,j,k) = 0.0
              endif
              if(mod(i+j+k,2).eq.0)Bvecx(i,j,k)=-1.0*Bvecx(i,j,k)
           else if(CCTK_EQUALS(monopole_type,"2dalt")) then
              rr2=x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2
              rg2=R_Gauss*R_Gauss
              if(rr2.lt.rg2) then
                 Bvecx(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
                 Bvecy(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
              else
                 Bvecx(i,j,k) = 0.0
                 Bvecy(i,j,k) = 0.0
              endif
              if(mod(i+j+k,2).eq.0)then
                 Bvecx(i,j,k)=-1.0*Bvecx(i,j,k)
!!$  Only vary one component, for different character
!!$                 Bvecy(i,j,k)=-1.0*Bvecy(i,j,k)
              endif
           else if(CCTK_EQUALS(monopole_type,"3dalt")) then
              rr2=x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2
              rg2=R_Gauss*R_Gauss
              if(rr2.lt.rg2) then
                 Bvecx(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
                 Bvecy(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
                 Bvecz(i,j,k) = exp(-1.0*rr2/rg2)-1.0/exp(1.0)
              else
                 Bvecx(i,j,k) = 0.0
                 Bvecy(i,j,k) = 0.0
                 Bvecz(i,j,k) = 0.0
              endif
!!$  Different spatial pattern!
              if(mod(i+j,2).eq.0)then
                 Bvecx(i,j,k)=-1.0*Bvecx(i,j,k)
                 Bvecy(i,j,k)=-1.0*Bvecy(i,j,k)
                 Bvecz(i,j,k)=-1.0*Bvecz(i,j,k)
              endif
           else
              call CCTK_WARN(0,"Unrecognized monopole type!!!")
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
  
end subroutine GRHydro_MonopoleM
