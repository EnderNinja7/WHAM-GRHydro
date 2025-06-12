 /*@@
   @file      GRHydro_ReadConformalData.F90
   @date      Fri Jul  6 12:29:27 2007
   @author    Luca Baiotti
   @desc 
   Set the missing quantities, after reading in from file initial data from conformally-flat codes (Garching)
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

 /*@@
   @routine    GRHydro_ReadConformalData
   @date       Fri Jul  6 12:31:18 2007
   @author     Luca Baiotti
   @desc 
   Set the missing quantities, after reading in from file initial data from conformally-flat codes (Garching)
   @enddesc 
   @calls     
   @calledby   
   @history 

   @endhistory 

@@*/

subroutine GRHydro_ReadConformalData(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k,handle,ierr
  CCTK_REAL :: rho_min, det

! begin EOS Omni vars
  CCTK_INT  :: n = 1
  CCTK_INT  :: poly_eoskey = 0
  CCTK_INT  :: keytemp = 0
  CCTK_INT  :: anyerr = 0
  CCTK_INT  :: keyerr(1) = 0
  CCTK_REAL :: rf_precision = 1.0d-10
  CCTK_REAL :: xpress(1) = 0.0d0
  CCTK_REAL :: xeps(1) = 0.0d0
  CCTK_REAL :: xtemp(1) = 0.0d0
  CCTK_REAL :: xye(1) = 0.0d0
  poly_eoskey = EOS_Omni_GetHandle("2D_Polytrope")
! end EOS Omni vars

  ! only gxx has been read in; copy it into gyy and gzz as well
  gyy = gxx
  gzz = gxx

  ! set the other components to zero
  gxy = 0.0
  gxz = 0.0
  gyz = 0.0

  ! set the extrinsic curvature to zero
  kxx = 0.0
  kxy = 0.0
  kxz = 0.0

  kyy = 0.0
  kyz = 0.0
  kzz = 0.0  

  ! set the shift to zero (the lapse is read in)
  betax = 0.0
  betay = 0.0
  betaz = 0.0


  ! atmosphere

  if (rho_abs_min > 0.0) then
    rho_min = rho_abs_min
  else
    call CCTK_WARN(0,"For now, with ReadConformalData you need to set rho_abs_min")
  end if

  if (initial_rho_abs_min > 0.0) then
    rho_min = initial_rho_abs_min;
  else if (initial_rho_rel_min > 0.0) then
    call CCTK_WARN(1,"Setting initial_rho_rel_min is ignored with ReadConformalData, for the time being")
  end if

  if (initial_atmosphere_factor > 0.0) rho_min = rho_min * initial_atmosphere_factor

  ierr = 0
  do i=1,cctk_lsh(1)
    do j=1,cctk_lsh(2)
      do k=1,cctk_lsh(3)

 !       write(*,*) i,j,k, rho(i,j,k)
        ! check on the read-in rho
        if (rho(i,j,k) > 0.0) then
          write(*,*) i,j,k, rho(i,j,k)
          ierr = ierr+1
        end if
        
        ! set the atmosphere
        if (rho(i,j,k) <= 0.d0) rho(i,j,k) = rho_min

!        if (rho(i,j,k) > 1.d-9)  write(*,*) i,j,k,rho(i,j,k)
      end do
    end do
  end do

!  if (ierr == 0) call CCTK_WARN(0,"rho.h5 contains only zeroes: stopping the simulation")

!  write(*,*)"rho min",  rho_min

  
!!$  where(rho < rho_min)
!!$    rho = rho_min
!!$  end where


  ! set pressure and eps from rho, using the polytropic EoS

!  handle = EOS_Handle("2D_Polytrope")

!  if (handle < 0) call CCTK_WARN(0,"For this hack you need to compile with EOS_Polytrope")

!  eos_k = EOS_Pressure(GRHydro_eos_handle, 1.0, 1.0)
!  eos_gamma = 1.0 + eos_k / EOS_SpecificIntEnergy(GRHydro_eos_handle,1.0,1.0)


!  press = eos_k * rho**eos_gamma
  
  do i=1,cctk_lsh(1)
    do j=1,cctk_lsh(2)
      do k=1,cctk_lsh(3)

        call EOS_Omni_EpsFromPress(poly_eoskey,keytemp,rf_precision,&
             n,rho(i,j,k),xeps,xtemp,xye,&
             press(i,j,k),eps(i,j,k),keyerr,anyerr)

      end do
    end do
  end do


  ! set velocities to zero

  velx(:,:,:) = 0.0
  vely(:,:,:) = 0.0
  velz(:,:,:) = 0.0
  w_lorentz = 1.0


  ! Compute conserved variables

  do i=1,cctk_lsh(1)
    do j=1,cctk_lsh(2)
      do k=1,cctk_lsh(3)
        
        call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
        
        if (CCTK_EQUALS(GRHydro_eos_type,"Polytype")) then
          call Prim2ConPoly(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
        else
          call Prim2ConGen(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
        end if
              
      end do
    end do
  end do
  

  return
  
end subroutine GRHydro_ReadConformalData
