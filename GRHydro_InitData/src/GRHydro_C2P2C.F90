 /*@@
   @file      GRHydro_C2P2C.F90
   @date      Sat Jan 26 02:44:43 2002
   @author    Luca Baiotti
   @desc 
   A test of the conservative <--> primitive variable exchange
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    c2p2c
   @date       Sat Jan 26 02:45:19 2002
   @author     Luca Baiotti
   @desc 
   Testing the conservative <--> primitive variable transformations.
   The values before and after should match.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine c2p2c(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: det
  CCTK_REAL :: uxx,uxy,uxz,uyy,uyz,uzz
  CCTK_REAL :: gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send
  CCTK_REAL :: dens_send,sx_send,sy_send,sz_send,tau_send
  CCTK_REAL :: rho_send,velx_send,vely_send,velz_send,eps_send
  CCTK_REAL :: press_send,w_lorentz_send,x_send,y_send,z_send,r_send
  CCTK_REAL :: C2P_failed
  CCTK_INT :: epsnegative

! begin EOS Omni vars
  CCTK_REAL :: pmin(1), epsmin(1), epsval(1)
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
  
  dens_send = 1.29047362d0
  sx_send   = 0.166666658d0
  sy_send   = 0.166666658d0
  sz_send   = 0.166666658d0
  tau_send  = 0.484123939d0
  
  eps_send = 1.0d-6
  press_send = 6.666666666666667d-7
  w_lorentz_send = 1.0d0

  epsnegative = 0
  
  xrho = 1.0d-10
  epsval = 1.0d0

  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,epsval,xtemp,xye,pmin,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  C2P_failed = 0.0d0

  write(*,*) 'C2P2C test: initial values.'
  write(*,*) '   conservative variables: '
  write(*,*) '   dens: ',dens_send
  write(*,*) '   sx  : ',sx_send
  write(*,*) '   sy  : ',sy_send
  write(*,*) '   sz  : ',sz_send
  write(*,*) '   tau : ',tau_send
  write(*,*) '   eps : ',eps_send
  write(*,*) '   W   : ',w_lorentz_send
  
  write(*,*) 'C2P2C test: getting the associated primitive variables.'
  call Con2PrimGen(GRHydro_eos_handle,dens_send,sx_send,sy_send,sz_send, &
       tau_send,rho_send,velx_send,vely_send,velz_send, &
       eps_send,press_send,w_lorentz_send, &
       uxx,uxy,uxz,uyy,uyz,uzz,det,x_send,y_send,z_send,r_send,&
       epsnegative,xrho(1),pmin(1),epsmin(1),GRHydro_init_data_reflevel,C2P_failed)
  
  write(*,*) 'C2P2C test: the primitive variables are'
  write(*,*) '   primitive variables: '
  write(*,*) '   rho   : ',rho_send
  write(*,*) '   velx  : ',velx_send
  write(*,*) '   vely  : ',vely_send
  write(*,*) '   velz  : ',velz_send
  write(*,*) '   press : ',press_send
  write(*,*) '   eps   : ',eps_send
  write(*,*) '   W     : ',w_lorentz_send
  
  write(*,*) 'C2P2C test: converting back to conserved variables.'
  call Prim2ConGen(GRHydro_eos_handle, &
       gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send,det, &
       dens_send,sx_send,sy_send,sz_send,tau_send,rho_send, &
       velx_send,vely_send,velz_send,eps_send,press_send,w_lorentz_send) 
  
  write(*,*) 'C2P2C test: the conserved variables are'
  write(*,*) '   conservative variables: '
  write(*,*) '   dens: ',dens_send
  write(*,*) '   sx  : ',sx_send
  write(*,*) '   sy  : ',sy_send
  write(*,*) '   sz  : ',sz_send
  write(*,*) '   tau : ',tau_send
  write(*,*) '   eps : ',eps_send
  write(*,*) '   W   : ',w_lorentz_send
  
  STOP

  return

end subroutine c2p2c
