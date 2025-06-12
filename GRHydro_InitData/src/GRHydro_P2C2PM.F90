 /*@@
   @file      GRHydro_P2C2PM.F90
   @date      Sep 25, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Luca Baiotti
   @desc 
   A test of the conservative <--> primitive variable exchange
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    p2c2pm
   @date       Sep 25, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Luca Baiotti
   @desc 
   Testing the conservative <--> primitive variable transformations.
   The values before and after should match.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine p2c2pm(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: det
  CCTK_REAL :: uxx,uxy,uxz,uyy,uyz,uzz
  CCTK_REAL :: gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send
  CCTK_REAL :: dens_send,sx_send,sy_send,sz_send,tau_send
  CCTK_REAL :: bconsx_send,bconsy_send,bconsz_send
  CCTK_REAL :: rho_send(1),velx_send,vely_send,velz_send,eps_send(1)
  CCTK_REAL :: press_send(1),w_lorentz_send,x_send,y_send,z_send,r_send
  CCTK_REAL :: bvcx_send,bvcy_send,bvcz_send,b2_send
  CCTK_REAL :: C2P_failed
  CCTK_INT :: epsnegative

! begin EOS Omni vars
  CCTK_REAL :: pmin(1), epsmin(1), local_gam(1), epsval(1)
  CCTK_INT :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress(1)=0.0d0;xtemp(1)=0.0d0;xye(1)=0.0d0;xeps(1)=0.0d0
! end EOS Omni vars

  call CCTK_WARN(1,"This test works only with Ideal_Fluid EoS")
  
  x_send = 0.0d0
  y_send = 0.0d0
  z_send = 0.0d0
  r_send = 0.0d0
  
  gxx_send = 1.0d0
  gyy_send = 1.0d0 
  gzz_send = 1.0d0
  gxy_send = 0.0d0
  gxz_send = 0.0d0
  gyz_send = 0.0d0
  
  det = 1.0d0
  
  uxx = 1.0d0
  uyy = 1.0d0 
  uzz = 1.0d0
  uxy = 0.0d0
  uxz = 0.0d0
  uyz = 0.0d0
  
  rho_send = 1.29047362d0
  velx_send   = 0.166666658d0
  vely_send   = 0.166666658d0
  velz_send   = 0.166666658d0
  eps_send  = 0.484123939d0
  
  bvcx_send   = Bx_init
  bvcy_send   = By_init
  bvcz_send   = Bz_init
  bconsx_send   = Bx_init
  bconsy_send   = By_init
  bconsz_send   = Bz_init

  w_lorentz_send = 1.d0/sqrt(1.0d0-velx_send*velx_send-vely_send*vely_send-velz_send*velz_send)

  epsnegative = 0
  
  xrho = 1.0d-10
  epsval = 1.0d0

  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho_send,eps_send,xtemp,xye,press_send,keyerr,anyerr)
  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,epsval,xtemp,xye,pmin,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  local_gam = 0.0d0
  xrho = 1.0d0
  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,epsval,xtemp,xye,local_gam,keyerr,anyerr)
  local_gam = local_gam + 1.0

  C2P_failed = 0.0d0

  write(*,*) 'P2C2PM test: the primitive variables are'
  write(*,*) '   primitive variables: '
  write(*,*) '   rho   : ',rho_send
  write(*,*) '   velx  : ',velx_send
  write(*,*) '   vely  : ',vely_send
  write(*,*) '   velz  : ',velz_send
  write(*,*) '   press : ',press_send
  write(*,*) '   eps   : ',eps_send
  write(*,*) '   W     : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send
  
  write(*,*) 'P2C2PM test: converting back to conserved variables.'
  call Prim2ConGenM(GRHydro_eos_handle,&
       gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send,det, &
       dens_send,sx_send,sy_send,sz_send,tau_send, &
       bconsx_send,bconsy_send,bconsz_send, &
       rho_send(1),velx_send,vely_send,velz_send,eps_send(1),press_send(1), &
       bvcx_send,bvcy_send,bvcz_send,w_lorentz_send) 

  write(*,*) 'P2C2PM test: initial values.'
  write(*,*) '   conservative variables: '
  write(*,*) '   dens: ',dens_send
  write(*,*) '   sx  : ',sx_send
  write(*,*) '   sy  : ',sy_send
  write(*,*) '   sz  : ',sz_send
  write(*,*) '   tau : ',tau_send
  write(*,*) '   Bconsx  : ',bconsx_send
  write(*,*) '   Bconsy  : ',bconsy_send
  write(*,*) '   Bconsz  : ',bconsz_send
  write(*,*) '   eps : ',eps_send
  write(*,*) '   W   : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send
  
  write(*,*) 'P2C2PM test: getting the associated primitive variables.'
  
  call Con2PrimGenM(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,local_gam(1),dens_send,sx_send,sy_send,sz_send, &
       tau_send,bconsx_send,bconsy_send,bconsz_send,&
       xtemp(1),xye(1),rho_send(1),velx_send,vely_send,velz_send, &
       eps_send(1),press_send(1), &      
       bvcx_send,bvcy_send,bvcz_send,b2_send,&
       w_lorentz_send, &
       gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send,&
       uxx,uxy,uxz,uyy,uyz,uzz,det,&
       epsnegative,C2P_failed)
  
  write(*,*) 'P2C2PM test: the primitive variables are'
  write(*,*) '   primitive variables: '
  write(*,*) '   rho   : ',rho_send
  write(*,*) '   velx  : ',velx_send
  write(*,*) '   vely  : ',vely_send
  write(*,*) '   velz  : ',velz_send
  write(*,*) '   press : ',press_send
  write(*,*) '   eps   : ',eps_send
  write(*,*) '   W     : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send
  
  STOP

  return

end subroutine p2c2pm
