 /*@@
   @file      GRHydro_RotorM.F90
   @date      Aug 15, 2011
   @author    Scott Noble, Joshua Faber, Bruno Mundim
   @desc 
   Cylindrical magnetized rotor test (see Etienne et al.).
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
   @routine    GRHydro_Rotor
   @date       Aug 11, 2011
   @author     Scott Noble, Joshua Faber, Bruno Mundim
   @desc 
   Initial data for magnetic rotor - parameters from Etienne et al. 2010
   @enddesc 
   @calls     
   @calledby   
   @history 
   Using GRHydro_ShockTubeM.F90 as a template.
   @endhistory 

@@*/

subroutine GRHydro_RotorM(CCTK_ARGUMENTS)

  use cctk

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: radius, det, rfact

  !begin EOS_Omni stuff
  CCTK_REAL :: tempEOS(1), yeEOS(1), xepsEOS(1)
  CCTK_INT :: keyerr(1), anyerr
  CCTK_INT, parameter :: have_temp = 0, n = 1
  !end EOS_Omni

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  !$OMP PARALLEL DO PRIVATE(i,j,k,radius,det,rfact,tempEOS,yeEOS,xepsEOS,keyerr,anyerr) &
  !$OMP default(none) &
  !$OMP firstprivate(nx,ny,nz,GRHydro_eos_handle, GRHydro_eos_rf_prec, grhydro_eos_type, &
  !$OMP              rotor_xc,rotor_yc, rotor_rsmooth_rel, rotor_use_smoothing, &
  !$OMP              rotor_r_rot, rotor_rhoin,rotor_pressin,rotor_rhoout,rotor_pressout, &
  !$OMP              rotor_v_max, rotor_bvcxl, rotor_bvcyl, rotor_bvczl) &
  !$OMP shared(x,y,z,rho,press,eps,vel,Bvec,w_lorentz,dens,Scon,tau,Bcons, &
  !$OMP        gxx,gxy,gxz,gyy,gyz,gzz)
  do i=1,nx
     do j=1,ny
        do k=1,nz

           radius = sqrt((x(i,j,k)-rotor_xc)**2+(y(i,j,k)-rotor_yc)**2)

           if(radius.le.rotor_r_rot) then
              rho(i,j,k) = rotor_rhoin
              press(i,j,k) = rotor_pressin
              velx(i,j,k) = -1.d0*rotor_v_max/rotor_r_rot*(y(i,j,k)-rotor_yc)
              vely(i,j,k) = +1.d0*rotor_v_max/rotor_r_rot*(x(i,j,k)-rotor_xc)
              velz(i,j,k) = 0.d0
           else if((rotor_use_smoothing.eq.1) .and. &
                ((radius.gt.rotor_r_rot) .and. &
                (radius.le.((1.0+rotor_rsmooth_rel)*rotor_r_rot)))) then
              rfact = (radius/rotor_r_rot - 1.d0) / rotor_rsmooth_rel
              rho(i,j,k) = rfact*rotor_rhoout + (1.d0-rfact)*rotor_rhoin
              press(i,j,k) = rfact*rotor_pressout + (1.d0-rfact)*rotor_pressin
              velx(i,j,k) = -1.d0*(1.d0-rfact)*rotor_v_max * (y(i,j,k)-rotor_yc) / radius
              velx(i,j,k) = +1.d0*(1.d0-rfact)*rotor_v_max * (x(i,j,k)-rotor_xc) / radius
              velz(i,j,k) = 0.d0
           else
              rho(i,j,k) = rotor_rhoout
              press(i,j,k) = rotor_pressout
              velx(i,j,k) = 0.d0
              vely(i,j,k) = 0.d0
              velz(i,j,k) = 0.d0
           endif
           keyerr = 0; anyerr = 0
           call EOS_Omni_EpsFromPress(GRHydro_eos_handle, have_temp, GRHydro_eos_rf_prec, &
                                      n, rho(i:i,j:j,k:k), eps(i:i,j:j,k:k), &
                                      tempEOS, yeEOS, press(i:i,j:j,k:k), &
                                      eps(i:i,j:j,k:k), & ! xeps argument
                                      keyerr, anyerr)
           if(anyerr .ne. 0) then
              call CCTK_WARN(0, "Error when calling EOS. After stopping swearing at me, add a decent output text.")
           end if

           Bvecx(i,j,k)=rotor_bvcxl
           Bvecy(i,j,k)=rotor_bvcyl
           Bvecz(i,j,k)=rotor_bvczl

           det=SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           if(abs(det - 1d0) .gt. 1e-8) then
              call CCTK_WARN(0, "Rotor initial data only supports flat spacetime right now")
           end if

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
  !$OMP END PARALLEL DO

end subroutine GRHydro_RotorM


