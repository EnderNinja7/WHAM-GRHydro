 /*@@
   @file      GRHydro_PoloidalMagFieldM.F90
   @date      Oct 31, 2011
   @author    Bruno Mundim, Joshua Faber, Scott Noble
   @desc 
   Poloidal Magnetic field implemented as in "General relativistic 
simulations of magnetized binary neutron star mergers" - 
by Yuk Tung Liu, Stuart L. Shapiro, Zachariah B. Etienne, and 
Keisuke Taniguchi - Phys. Rev. D 78, 024012 (2008)  

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
#define Avecx(i,j,k) Avec(i,j,k,1)
#define Avecy(i,j,k) Avec(i,j,k,2)
#define Avecz(i,j,k) Avec(i,j,k,3)

 /*@@
   @routine   GRHydro_PoloidalMagFieldM
   @date      Oct 31, 2011
   @author    Bruno Mundim, Joshua Faber, Scott Noble
   @desc 
   Poloidal Magnetic field implemented as in "General relativistic 
simulations of magnetized binary neutron star mergers" - 
by Yuk Tung Liu, Stuart L. Shapiro, Zachariah B. Etienne, and 
Keisuke Taniguchi - Phys. Rev. D 78, 024012 (2008)  
   @enddesc 
   @calls     
   @calledby   
   @history 
   Using GRHydro_ShockTubeM.F90 as a template.
   @endhistory 

@@*/

subroutine GRHydro_PoloidalMagFieldM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz, set_Avec
  CCTK_REAL :: det
  CCTK_REAL :: sdet
  CCTK_REAL :: dx,dy,dz
  CCTK_REAL :: rhofac, delPcut, maxP_Pcut
  CCTK_REAL :: AphiL, Ax, Ay, Az
  CCTK_REAL :: rho_dx, rho_dy, rho_dz
  CCTK_REAL :: press_dx, press_dy, press_dz
  CCTK_REAL :: Aphi_dx, Aphi_dy, Aphi_dz
  CCTK_REAL :: Ax_dy, Ax_dz, Ay_dx, Ay_dz
  integer :: have_GetMap, have_MapIsCartesian

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  set_Avec = 0
  if ( CCTK_EQUALS(Bvec_evolution_method,"GRHydro_Avec") ) then
     set_Avec = 1
  end if

  write(*,*)'GRHydro_InitData: Setting up initial poloidal magnetic field'

  ! abort if called on curvilinear grid (no support for Jacobian)
  call CCTK_IsFunctionAliased(have_GetMap, "GetMap")
  call CCTK_IsFunctionAliased(have_MapIsCartesian, "MultiPatch_MapIsCartesian")
  if (have_GetMap .ne. 0 .and. have_MapIsCartesian .ne. 0) then
     if (.not. MultiPatch_MapIsCartesian(GetMap(cctkGH)) .ne. 0) then
        call CCTK_ERROR("This algorithm does only work on Cartesian grids!!")
     endif
  end if

! Initialize to zero
  Bvec = 0.0d0

  do i=2,nx-1
   do j=2,ny-1
    do k=2,nz-1

     rhofac = 1.0d0-rho(i,j,k)/poloidal_rho_max
     delPcut = press(i,j,k)-poloidal_P_cut
     maxP_Pcut = max(delPcut,0.0d0)
     AphiL = poloidal_A_b*rhofac**poloidal_n_p*maxP_Pcut**poloidal_P_p
     Ax = -y(i,j,k)*AphiL
     Ay =  x(i,j,k)*AphiL
     Az = 0.0
 
!!     write(*,*)'Before accessing rho(i,k,k)'
!!     write(*,*)'rho(',i,',',j,',',k,') = ', rho(i,j,k)
!!     write(*,*)'Ax, Ay, Az, Aphi = ', Ax, Ay, Az,AphiL
!!     write(*,*)'rhofac = ', rhofac
!!     write(*,*)'delPcut = ', delPcut
!!     write(*,*)'maxP_Pcut = ', maxP_Pcut 

     rho_dx = 0.5d0*(rho(i+1,j,k)-rho(i-1,j,k))/dx 
     rho_dy = 0.5d0*(rho(i,j+1,k)-rho(i,j-1,k))/dy 
     rho_dz = 0.5d0*(rho(i,j,k+1)-rho(i,j,k-1))/dz 
     press_dx = 0.5d0*(press(i+1,j,k)-press(i-1,j,k))/dx 
     press_dy = 0.5d0*(press(i,j+1,k)-press(i,j-1,k))/dy 
     press_dz = 0.5d0*(press(i,j,k+1)-press(i,j,k-1))/dz 

     if (maxP_Pcut > 0.0) then
       Aphi_dx = poloidal_A_b*( &
                   -poloidal_n_p*rho_dx/poloidal_rho_max*maxP_Pcut**poloidal_P_p &
                   +rhofac**poloidal_n_p*press_dx*poloidal_P_p*maxP_Pcut**(poloidal_P_p-1))
       Aphi_dy = poloidal_A_b*( &
                   -poloidal_n_p*rho_dy/poloidal_rho_max*maxP_Pcut**poloidal_P_p &
                   +rhofac**poloidal_n_p*press_dy*poloidal_P_p*maxP_Pcut**(poloidal_P_p-1))
       Aphi_dz = poloidal_A_b*( &
                   -poloidal_n_p*rho_dz/poloidal_rho_max*maxP_Pcut**poloidal_P_p &
                   +rhofac**poloidal_n_p*press_dz*poloidal_P_p*maxP_Pcut**(poloidal_P_p-1))
     else
       Aphi_dx = 0.0
       Aphi_dy = 0.0
       Aphi_dz = 0.0
     endif

     Ax_dy = -AphiL - y(i,j,k)*Aphi_dy
     Ax_dz = -y(i,j,k)*Aphi_dz
     Ay_dx = AphiL + x(i,j,k)*Aphi_dx
     Ay_dz = x(i,j,k)*Aphi_dz

           
     det=SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
     sdet = sqrt(det) 
     
     Bvecx(i,j,k) = Ay_dz/sdet
     Bvecy(i,j,k) = - Ax_dz/sdet
     Bvecz(i,j,k) = (Ax_dy-Ay_dx)/sdet

     if ( set_Avec.gt.0 ) then
        Avecx(i,j,k) = Ax
        Avecy(i,j,k) = Ay
        Avecz(i,j,k) = Az
     end if

     !Bvecx(i,j,k) = 0.0d0 
     !Bvecy(i,j,k) = 0.0d0
     !Bvecz(i,j,k) = 0.00000001/sdet

           if (CCTK_EQUALS(GRHydro_eos_type,"Polytype")) then
              call Prim2ConPolyM(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                   gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                   det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                   tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),rho(i,j,k),&
                   velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                   eps(i,j,k),press(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),&
                   w_lorentz(i,j,k))
           else
             if (evolve_temper .ne. 0) then
               call Prim2ConGenM_hot(GRHydro_eos_handle,GRHydro_reflevel,&
                    i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),gxx(i,j,k),gxy(i,j,k),&
                    gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                    det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                    tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),rho(i,j,k),&
                    velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                    eps(i,j,k),press(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),&
                    w_lorentz(i,j,k),temperature(i,j,k),y_e(i,j,k))
             else
               call Prim2ConGenM(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                    gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                    det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                    tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),rho(i,j,k),&
                    velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                    eps(i,j,k),press(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),&
                    w_lorentz(i,j,k))
             end if
           end if
           
        enddo
     enddo
  enddo
  
  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0
  Bconsrhs = 0.d0
  if (clean_divergence .ne. 0) then
    psidcrhs =0.0
  endif
  
  !Bvec = 0
  !lBvec = 0
  !Bcons = 0
  !Avec = 0
  
  return
  
end subroutine GRHydro_PoloidalMagFieldM


