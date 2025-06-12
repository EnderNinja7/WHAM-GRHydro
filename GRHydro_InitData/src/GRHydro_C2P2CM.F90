 /*@@
   @file      GRHydro_C2P2CM.F90
   @date      Sep 23, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Luca Baiotti
   @desc 
   A test of the conservative <--> primitive variable exchange
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"


 /*@@
   @routine    c2p2cM
   @date       Sep 23, 2010
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

subroutine c2p2cM(CCTK_ARGUMENTS)


  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: det, sdet, invdet
  CCTK_REAL :: uxx,uxy,uxz,uyy,uyz,uzz
  CCTK_REAL :: gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send
  CCTK_REAL :: dens_send,sx_send,sy_send,sz_send,tau_send
  CCTK_REAL :: bconsx_send, bconsy_send, bconsz_send
  CCTK_REAL :: entropy_send, entropycons_send;
  CCTK_REAL :: rho_send,velx_send,vely_send,velz_send,eps_send
  CCTK_REAL :: press_send,w_lorentz_send
  CCTK_REAL :: bvcx_send, bvcy_send, bvcz_send, b2_send
  CCTK_REAL :: C2P_failed
  CCTK_REAL :: BdotB, Bdotv, b2, beta_mag
  CCTK_INT :: epsnegative

! begin EOS Omni vars
  CCTK_REAL :: pmin(1), epsmin(1), local_gam(1), epsval(1)
  CCTK_INT :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress(1)=0.0d0;xtemp(1)=0.0d0;xye(1)=0.0d0;xeps(1)=0.0d0
! end EOS Omni vars

  call CCTK_WARN(1,"This test works only with Ideal_Fluid EoS")
  
  gxx_send = gxx_init
  gxy_send = gxy_init
  gxz_send = gxz_init
  gyy_send = gyy_init
  gyz_send = gyz_init
  gzz_send = gzz_init
  
  det = SPATIAL_DETERMINANT(gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send)
  sdet = sqrt(det)
  invdet = 1.d0 / det
  uxx = (-gyz_send**2 + gyy_send*gzz_send)*invdet
  uxy = (gxz_send*gyz_send - gxy_send*gzz_send)*invdet
  uyy = (-gxz_send**2 + gxx_send*gzz_send)*invdet
  uxz = (-gxz_send*gyy_send + gxy_send*gyz_send)*invdet
  uyz = (gxy_send*gxz_send - gxx_send*gyz_send)*invdet
  uzz = (-gxy_send**2 + gxx_send*gyy_send)*invdet

! Initialize the velocity as GRHydro_Con2PrimM_pt may use these
! values as initial guess for its Newton-Raphson procedure.
!  velx_send = 0.1d0
!  vely_send = 0.1d0 
!  velz_send = 0.1d0
!  
!  dens_send = 1.29047362d0
!  sx_send   = 0.166666658d0
!  sy_send   = 0.166666658d0
!  sz_send   = 0.166666658d0
!  tau_send  = 0.484123939d0

  velx_send = velx_init
  vely_send = vely_init
  velz_send = velz_init

  dens_send = dens_init
  sx_send   = sx_init
  sy_send   = sy_init
  sz_send   = sz_init
  tau_send  = tau_init

  bvcx_send   = Bx_init
  bvcy_send   = By_init
  bvcz_send   = Bz_init
  bconsx_send   = sdet*Bx_init
  bconsy_send   = sdet*By_init
  bconsz_send   = sdet*Bz_init

  rho_send = rho_init
  eps_send = eps_init
  press_send = press_init
  w_lorentz_send = sqrt(1.0d0-(gxx_send*velx_send**2+gyy_send*vely_send**2+ &
                               gzz_send*velz_send**2 &
                               +2.0d0*(gxy_send*velx_send*vely_send+ &
                                       gxz_send*velx_send*velz_send+ &
                                       gyz_send*vely_send*velz_send &
                                      ) &
                               ) &
                        ) 
  w_lorentz_send = 1.0d0/w_lorentz_send


  epsnegative = 0
  
  xrho = 1.0d-10
  epsval = 1.0d0

  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,epsval,xtemp,xye,pmin,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  local_gam = 0.0d0
  xrho = 1.0d0
  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,epsval,xtemp,xye,local_gam,keyerr,anyerr)
  local_gam = local_gam + 1.0

  if(use_c2p_with_entropy_eqn.eq.1)then
    entropy_send = ( local_gam(1) - 1.0d0 ) * eps_send * &
                     rho_send**(2.0d0 - local_gam(1))
    entropycons_send = sdet * w_lorentz_send * entropy_send
  endif

  C2P_failed = 0.0d0

  write(*,*) 'C2P2CM test: metric values.'
  write(*,*) '    gxx: ', gxx_send
  write(*,*) '    gxy: ', gxy_send
  write(*,*) '    gxz: ', gxz_send
  write(*,*) '    gyy: ', gyy_send
  write(*,*) '    gyz: ', gyz_send
  write(*,*) '    gzz: ', gzz_send
  write(*,*) '    uxx: ', uxx
  write(*,*) '    uxy: ', uxy
  write(*,*) '    uxz: ', uxz
  write(*,*) '    uyy: ', uyy
  write(*,*) '    uyz: ', uyz
  write(*,*) '    uzz: ', uzz
  write(*,*) '    det: ', det

  write(*,*) 'C2P2CM test: initial primitive guess values.'
  write(*,*) '   primitive variables: '
  write(*,*) '   rho   : ',rho_send
  write(*,*) '   velx  : ',velx_send
  write(*,*) '   vely  : ',vely_send
  write(*,*) '   velz  : ',velz_send
  write(*,*) '   press : ',press_send
  write(*,*) '   eps   : ',eps_send
  if(use_c2p_with_entropy_eqn.eq.1)then
    write(*,*) '   entropy  : ',entropy_send
  endif
  write(*,*) '   W     : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send
  write(*,*) '   C2P_failed : ',C2P_failed

  write(*,*) 'C2P2CM test: initial values.'
  write(*,*) '   conservative variables: '
  write(*,*) '   dens: ',dens_send
  write(*,*) '   sx  : ',sx_send
  write(*,*) '   sy  : ',sy_send
  write(*,*) '   sz  : ',sz_send
  write(*,*) '   tau : ',tau_send
  write(*,*) '   Bconsx  : ',bconsx_send
  write(*,*) '   Bconsy  : ',bconsy_send
  write(*,*) '   Bconsz  : ',bconsz_send
  if(use_c2p_with_entropy_eqn.eq.1)then
    write(*,*) '   entropycons  : ',entropycons_send
  endif
  write(*,*) '   eps : ',eps_send
  write(*,*) '   W   : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send

! Is this point magnetically dominated?
   BdotB = gxx_send*bvcx_send**2+gyy_send*bvcy_send**2+gzz_send*bvcz_send**2&
     +2.0*(gxy_send*bvcx_send*bvcy_send+gxz_send*bvcx_send*bvcy_send+ &
           gyz_send*bvcy_send*bvcz_send) 
   Bdotv = gxx_send*bvcx_send*velx_send+gyy_send*bvcy_send*vely_send  &
          +gzz_send*bvcz_send*velz_send  &
          +gxy_send*(bvcx_send*vely_send+bvcy_send*velx_send)  &
          +gxz_send*(bvcx_send*velz_send+bvcz_send*velx_send)  &
          +gyz_send*(bvcy_send*velz_send+bvcz_send*vely_send)

   b2 = BdotB/w_lorentz_send**2 + Bdotv**2

   beta_mag = 2.0*press_send/b2;

  write(*,*) '   BdotB   : ',BdotB
  write(*,*) '   Bdotv   : ',Bdotv
  write(*,*) '      b2   : ',b2
  write(*,*) 'beta_mag   : ',beta_mag
  
  if(use_c2p_with_entropy_eqn.eq.0)then
     write(*,*) 'C2P2CM test: getting the associated primitive variables.'
     call Con2PrimGenM(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,local_gam(1),dens_send,sx_send,sy_send,sz_send, &
       tau_send,bconsx_send,bconsy_send,bconsz_send,xtemp(1),xye(1),rho_send,velx_send,vely_send,velz_send, &
       eps_send,press_send,bvcx_send,bvcy_send,bvcz_send,b2_send,w_lorentz_send, &
       gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send,&
       uxx,uxy,uxz,uyy,uyz,uzz,det,&
       epsnegative,C2P_failed)
  else
     write(*,*) 'C2P2CMee test: getting the associated primitive variables.'
     call Con2PrimGenMee(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,local_gam(1),dens_send,sx_send,sy_send,sz_send, &
       tau_send,bconsx_send,bconsy_send,bconsz_send,entropycons_send,xtemp(1),xye(1),rho_send,velx_send,vely_send,velz_send, &
       eps_send,press_send,bvcx_send,bvcy_send,bvcz_send,b2_send,w_lorentz_send, &
       gxx_send,gxy_send,gxz_send,gyy_send,gyz_send,gzz_send,&
       uxx,uxy,uxz,uyy,uyz,uzz,det,&
       epsnegative,C2P_failed)
  endif

  if(use_c2p_with_entropy_eqn.eq.1)then
    entropy_send = ( local_gam(1) - 1.0d0 ) * eps_send * &
                     rho_send**(2.0d0 - local_gam(1))
    entropycons_send = sdet * w_lorentz_send * entropy_send
  endif
  
  write(*,*) 'C2P2CM test: the primitive variables are'
  write(*,*) '   primitive variables: '
  write(*,*) '   rho   : ',rho_send
  write(*,*) '   velx  : ',velx_send
  write(*,*) '   vely  : ',vely_send
  write(*,*) '   velz  : ',velz_send
  write(*,*) '   press : ',press_send
  write(*,*) '   eps   : ',eps_send
  if(use_c2p_with_entropy_eqn.eq.1)then
    write(*,*) '   entropy  : ',entropy_send
  endif
  write(*,*) '   W     : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send
  write(*,*) '   C2P_failed : ',C2P_failed
  
  write(*,*) 'C2P2CM test: converting back to conserved variables.'
  call Prim2ConGenM(GRHydro_eos_handle,gxx_send, gxy_send, gxz_send, gyy_send, gyz_send, gzz_send, det, &
       dens_send, sx_send, sy_send, sz_send, tau_send, bconsx_send, bconsy_send, bconsz_send, rho_send, &
       velx_send, vely_send, velz_send, eps_send, press_send, bvcx_send, bvcy_send, bvcz_send, w_lorentz_send) 
  
  write(*,*) 'C2P2CM test: the conserved variables are'
  write(*,*) '   conservative variables: '
  write(*,*) '   dens: ',dens_send
  write(*,*) '   sx  : ',sx_send
  write(*,*) '   sy  : ',sy_send
  write(*,*) '   sz  : ',sz_send
  write(*,*) '   tau : ',tau_send
  write(*,*) '   Bconsx  : ',bconsx_send
  write(*,*) '   Bconsy  : ',bconsy_send
  write(*,*) '   Bconsz  : ',bconsz_send
  if(use_c2p_with_entropy_eqn.eq.1)then
    write(*,*) '   entropycons  : ',entropycons_send
  endif
  write(*,*) '   eps : ',eps_send
  write(*,*) '   W   : ',w_lorentz_send
  write(*,*) '   Bvecx  : ',bvcx_send
  write(*,*) '   Bvecy  : ',bvcy_send
  write(*,*) '   Bvecz  : ',bvcz_send
  
  STOP

  return

end subroutine c2p2cM
