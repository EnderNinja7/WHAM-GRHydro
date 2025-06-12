 /*@@
   @file      GRHydro_Con2Prim.F90
   @date      Sat Jan 26 02:49:32 2002
   @author    Luca Baiotti
   @desc 
   A test of the conservative to primitive exchange.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    GRHydro_con2primtest
   @date       Sat Jan 26 02:49:58 2002
   @author     Luca Baiotti
   @desc 
   A test of the conservative to primitive variable solver.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_Init_Data_RefinementLevel(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS

  GRHydro_init_data_reflevel = aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))

end subroutine GRHydro_Init_Data_RefinementLevel


subroutine GRHydro_con2primtest(CCTK_ARGUMENTS)
  
  implicit  none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer didit,i,j,k,nx,ny,nz
  CCTK_REAL det,uxx,uxy,uxz,uyy,uyz,uzz
  CCTK_REAL dens_send,sx_send,sy_send,sz_send,tau_send
  CCTK_REAL rho_send,velx_send,vely_send,velz_send,eps_send
  CCTK_REAL press_send,w_lorentz_send,x_send,y_send,z_send,r_send
  CCTK_REAL pmin, epsmin
  CCTK_REAL C2P_failed
  logical epsnegative

! begin EOS Omni vars
  integer :: n = 1
  integer :: keytemp = 0
  integer :: anyerr = 0
  integer :: keyerr(1) = 0
  real*8  :: xpress = 0.0d0
  real*8  :: xeps = 0.0d0
  real*8  :: xtemp = 0.0d0
  real*8  :: xye = 0.0d0
! end EOS Omni vars

  call CCTK_WARN(1,"For this test, remember to use a polytropic EoS and to set eos_gamma = 2.0 and eos_k = 100.0")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  x_send = 0.0d0
  y_send = 0.0d0
  z_send = 0.0d0
  r_send = 0.0d0
  
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
  
  rho_send = 1.0d0
  velx_send = 0.0d0
  vely_send = 0.0d0
  velz_send = 0.0d0
  eps_send = 1.0d-6
  press_send = 6.666666666666667d-7
  w_lorentz_send = 1.0d0

  epsnegative = .false.

  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       GRHydro_rho_min,xeps,xtemp,xye,pmin,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       GRHydro_rho_min,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  C2P_failed = 0.d0

  write(*,*) 'Con2Prim test: converting to primitive variables'
  call Con2Prim_pt(GRHydro_eos_handle,dens_send,sx_send,sy_send,sz_send,&
       tau_send,rho_send,velx_send,vely_send,velz_send,&
       eps_send,press_send,w_lorentz_send, &
       uxx,uxy,uxz,uyy,uyz,uzz,det,x_send,y_send,z_send,r_send,&
       epsnegative,GRHydro_rho_min, pmin, epsmin, GRHydro_init_data_reflevel, C2P_failed)
  write(*,*) 'Con2Prim test: the primitive variables are'
  write(*,*) '   primitive variables: '
  write(*,*) '   rho   : ',rho_send
  write(*,*) '   velx  : ',velx_send
  write(*,*) '   vely  : ',vely_send
  write(*,*) '   velz  : ',velz_send
  write(*,*) '   eps   : ',eps_send
  write(*,*) '   press : ',press_send
  write(*,*) '   w_lor : ',w_lorentz_send
  
  STOP

  return

end subroutine GRHydro_con2primtest
