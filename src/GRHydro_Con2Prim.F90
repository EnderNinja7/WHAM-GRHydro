/*@@
   @file      GRHydro_RegisterVars.c
   @date      Sat Jan 26 01:06:01 2002
   @author    The GRHydro Developers
   @desc 
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "GRHydro_Macros.h"


#define ITER_TOL (1.0e-8)
#define MAXITER (50)

 /*@@
   @routine    Conservative2Primitive
   @date       Sat Jan 26 01:08:33 2002
   @author     Ian Hawke
   @desc 
   Wrapper routine that converts from conserved to primitive variables
   at every grid cell centre.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Trimmed and altered from the GR3D routines, original author Mark Miller.
   2007?: Bruno excluded the points in the atmosphere and excision region from the computation.
   Aug. 2008: Luca added a check on whether a failure at a given point may be disregarded, 
   because that point will then be restricted from a finer level. This should be completely 
   safe only if *regridding happens at times when all levels are evolved.*
   Feb. 2009: The above procedure proved to be wrong, so Luca implemented another one. 
   When a failure occurs, it is temporarily ignored, except for storing the location of where 
   it occured in a mask. At the end, after a Carpet restriction, the mask is checked and if 
   it still contains failures, the run is aborted with an error message. Only used routines 
   have been updated to use this procedure.
   @endhistory 

@@*/

subroutine Conservative2Primitive(CCTK_ARGUMENTS)

  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, sdet, pmin, epsmin, dummy1, dummy2, temp22
  logical :: epsnegative
  character*256 :: warnline

  CCTK_REAL :: oob, b2, d2, s2, bscon, bxhat, byhat, bzhat, bhatscon
  CCTK_REAL :: Wm, Wm0, Wm_resid, Wmold
  CCTK_REAL :: s2m, s2m0, s2m_resid, s2mold, s2max, taum
  CCTK_INT :: niter
  
  CCTK_REAL :: local_min_tracer
  CCTK_REAL :: local_perc_ptol

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  !Temp data holders for de-averaging
  CCTK_REAL, DIMENSION(cctk_lsh(1), cctk_lsh(2), cctk_lsh(3)) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8&
  , temp9, temp10, dens_avg, tau_avg, scon1_avg, scon2_avg, scon3_avg

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars

  if(evolve_temper.ne.0) then
     call Conservative2PrimitiveHot(CCTK_PASS_FTOF)
     return
  endif


  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  !Initialize de-average temps
  temp1 = 0
  temp22 = 0
  temp2 = 0
  temp3 = 0
  temp4 = 0
  temp5 = 0
  temp6 = 0
  temp7 = 0
  temp8 = 0
  temp9 = 0
  temp10 = 0
  dens_avg = 0
  tau_avg = 0
  scon1_avg = 0
  scon2_avg = 0
  scon3_avg = 0
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0.0d0
  end if


  !
  !call de-averaging step
  call apply(dens, nx, ny, nz, 0, temp1)
  call apply(temp1, nx, ny, nz, 1, temp2)
  call apply(temp2, nx, ny, nz, 2, dens_avg)

  call apply(tau, nx, ny, nz, 0, temp3)
  call apply(temp3, nx, ny, nz, 1, temp4)
  call apply(temp4, nx, ny, nz, 2, tau_avg)

  call apply(scon(:,:,:,1), nx, ny, nz, 0, temp5)
  call apply(temp5, nx, ny, nz, 1, temp6)
  call apply(temp6, nx, ny, nz, 2, scon1_avg)

  call apply(scon(:,:,:,2), nx, ny, nz, 0, temp7)
  call apply(temp7, nx, ny, nz, 1, temp8)
  call apply(temp8, nx, ny, nz, 2, scon2_avg)

  call apply(scon(:,:,:,3), nx, ny, nz, 0, temp9)
  call apply(temp9, nx, ny, nz, 1, temp10)
  call apply(temp10, nx, ny, nz, 2, scon3_avg)

  !For plotting dens and rho output
  open (unit=10,file="out.csv",action="write")

  ! this is a poly call
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
         GRHydro_rho_min,xeps,xtemp,xye,pmin,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       GRHydro_rho_min,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, epsnegative, anyerr, keyerr, keytemp,&
  !$OMP warnline, dummy1, dummy2)
  !do k = 1+GRHydro_stencil, nz-GRHydro_stencil
  !  do j = 1+GRHydro_stencil, ny-GRHydro_stencil
  !    do i = 1+GRHydro_stencil, nx-GRHydro_stencil
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        sdet = sdetg(i,j,k)
        !Fix as per GRHydro_Con2PrimM.F90
        if(sqrtdet_thr.gt.0d0 .and. sdetg(i,j,k).ge.sqrtdet_thr) then
          d2 = dens(i,j,k)**2
                  s2 = uxx*scon1_avg(i,j,k)**2 + uyy*scon2_avg(i,j,k)**2 &
                                            + uzz*scon3_avg(i,j,k)**2 &
                     + 2.0d0*uxy*scon1_avg(i,j,k)*scon2_avg(i,j,k) &
                     + 2.0d0*uxz*scon1_avg(i,j,k)*scon3_avg(i,j,k)&
                     + 2.0d0*uyz*scon2_avg(i,j,k)*scon3_avg(i,j,k) 
                 oob = 1.0d0/sqrt(b2)
               bxhat = 0.0d0
               byhat = 0.0d0
               bzhat = 0.0d0
               bhatscon = bxhat*scon(i,j,k,1)+byhat*scon2_avg(i,j,k) &
                                             +bzhat*scon3_avg(i,j,k)
               bscon = 0.0d0
               b2 = 0.0d0
               ! Initial guesses for iterative procedure to find Wm:
               Wm0 = sdet*sqrt(bhatscon**2+d2)
               s2m0 = (Wm0**2*s2+bhatscon**2*(b2+2.0d0*Wm0)) &
                    / (Wm0+b2)**2 
               Wm = sdet*sqrt(s2m0+d2)
               s2m = (Wm**2*s2+bscon**2*(b2+2.0d0*Wm)) &
                   / (Wm+b2)**2 
               s2m_resid = 1.0d60
               Wm_resid = 1.0d60
               niter = 0
               do while((s2m_resid.ge.ITER_TOL.and.Wm_resid.ge.ITER_TOL).and.&
                        niter.le.MAXITER)
                 Wmold = Wm
                 s2mold = s2m
                 Wm = sdet*sqrt(s2m+d2)
                 s2m = (Wm**2*s2+bscon**2*(b2+2.0d0*Wm)) &
                     / (Wm+b2)**2 
                 Wm_resid = abs(Wmold-Wm)
                 s2m_resid = abs(s2mold-s2m)
                 niter = niter + 1
               end do
               !TODO: abort execution if niter .eq. MAXITER and warn user
               taum = tau_avg(i,j,k) - 0.5d0*sdet*b2 -0.5d0*(b2*s2-bscon**2)/ &
                                                       (sdet*(Wm+b2)**2) 
               s2max = taum*(taum+2.0d0*dens(i,j,k))
               if(taum.lt.GRHydro_tau_min)then
                 tau_avg(i,j,k) = GRHydro_tau_min + 0.5d0*sdet*b2 + 0.5d0* &
                              (b2*s2-bscon**2)/(sdet*(Wm+b2)**2)
               end if
               if(s2.gt.s2max) then
                 scon(i,j,k,1) = scon(i,j,k,1)*sqrt(s2max/s2)
                 scon(i,j,k,2) = scon(i,j,k,2)*sqrt(s2max/s2)
                 scon(i,j,k,3) = scon(i,j,k,3)*sqrt(s2max/s2)
               end if
        endif

#if 0
         if (dens_avg(i,j,k).ne.dens_avg(i,j,k)) then
            !$OMP CRITICAL
            call CCTK_WARN(1,"dens NAN at entry of Con2Prim")
            write(warnline,"(A10,i5)") "reflevel: ",GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,"(1P10E15.6)") rho(i,j,k),dens_avg(i,j,k),eps(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a32,2i3)') 'hydro_excision_mask, atmosphere_mask: ', hydro_excision_mask(i,j,k), atmosphere_mask(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_ERROR("Aborting!!!")
            STOP
            !$OMP END CRITICAL
         endif
#endif


#if 0
        !$OMP CRITICAL
         if (rho(i,j,k).le.0.0d0) then
            call CCTK_WARN(1,"rho less than zero at entry of Con2Prim")
            write(warnline,"(A10,i5)") "reflevel: ",GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,"(1P10E15.6)") rho(i,j,k),dens_avg(i,j,k),eps(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_ERROR("Aborting!!!")
            STOP
         endif
         !$OMP END CRITICAL
#endif

        
        !do not compute if in atmosphere region!
        if (atmosphere_mask(i,j,k) .gt. 0) cycle
        
        epsnegative = .false.

        call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdetg(i,j,k)*sdetg(i,j,k),&
             g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),&
             g23(i,j,k),g33(i,j,k))        

        ! In excised region, set to atmosphere!    
        if (GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .gt. 0)) then
           !SET_ATMO_MIN(dens(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k)) !sqrt(det)*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
           SET_ATMO_MIN(dens_avg(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k))
           SET_ATMO_MIN(rho(i,j,k), GRHydro_rho_min, r(i,j,k)) !GRHydro_rho_min
           scon(i,j,k,:) = 0.d0
           vup(i,j,k,:) = 0.d0
           w_lorentz(i,j,k) = 1.d0

           if(evolve_temper.ne.0) then
              ! set hot atmosphere values
              temperature(i,j,k) = grhydro_hot_atmo_temp
              y_e(i,j,k) = grhydro_hot_atmo_Y_e
              y_e_con(i,j,k) = y_e(i,j,k) * dens_avg(i,j,k)
              keytemp = 1
              call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                   press(i,j,k),keyerr,anyerr)
           else
              keytemp = 0
              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)

              call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
           endif

           ! w_lorentz=1, so the expression for tau reduces to:
           tau(i,j,k)  = sdetg(i,j,k) * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens(i,j,k)
           
           cycle
        endif
        
!!$        if (det < 0.e0) then  
!!$          call CCTK_ERROR("nan produced (det negative)")
!!$          STOP
!!$        end if

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), tracer(i,j,k,itracer), &
                   dens_avg(i,j,k))

              if (use_min_tracer .ne. 0) then
                if (tracer(i,j,k,itracer) .le. local_min_tracer) then
                  tracer(i,j,k,itracer) = local_min_tracer
                end if
              end if

           enddo
           
        endif
        
        if(evolve_Y_e.ne.0) then
           Y_e(i,j,k) = max(min(Y_e_con(i,j,k) / dens_avg(i,j,k),GRHydro_Y_e_max),&
                GRHydro_Y_e_min)
        endif

         !if ( dens(i,j,k) .le. sqrt(det)*GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
         !IF_BELOW_ATMO(dens_avg(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, GRHydro_atmo_tolerance, r(i,j,k)) then
           !SET_ATMO_MIN(dens(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k)) !sqrt(det)*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
         !  SET_ATMO_MIN(dens_avg(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k))
         !  SET_ATMO_MIN(rho(i,j,k), GRHydro_rho_min, r(i,j,k)) !GRHydro_rho_min
         !  scon(i,j,k,:) = 0.d0
         !  scon1_avg(i,j,k) = 0.d0
         !  scon2_avg(i,j,k) = 0.d0
         !  scon3_avg(i,j,k) = 0.d0
         !  vup(i,j,k,:) = 0.d0
         !  w_lorentz(i,j,k) = 1.d0

         !  if(evolve_temper.ne.0) then
         !     ! set hot atmosphere values
         !     temperature(i,j,k) = grhydro_hot_atmo_temp
         !    y_e(i,j,k) = grhydro_hot_atmo_Y_e
         !    y_e_con(i,j,k) = y_e(i,j,k) * dens_avg(i,j,k)
         !     keytemp = 1
         !     call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
         !          rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
         !          press(i,j,k),keyerr,anyerr)
         !  else
         !    keytemp = 0
         !     call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
         !          rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)

         !     call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
         !          rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
         !  endif

         !  ! w_lorentz=1, so the expression for tau reduces to:
         !  tau_avg(i,j,k)  = sdetg(i,j,k) * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens_avg(i,j,k)
         !  
         !  cycle

         ! end if

         IF_BELOW_ATMO(dens(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, GRHydro_atmo_tolerance, r(i,j,k)) then
           !SET_ATMO_MIN(dens(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k)) !sqrt(det)*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
           SET_ATMO_MIN(dens(i,j,k), sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k))
           SET_ATMO_MIN(rho(i,j,k), GRHydro_rho_min, r(i,j,k)) !GRHydro_rho_min
           scon(i,j,k,:) = 0.d0
           scon1_avg(i,j,k) = 0.d0
           scon2_avg(i,j,k) = 0.d0
           scon3_avg(i,j,k) = 0.d0
           vup(i,j,k,:) = 0.d0
           w_lorentz(i,j,k) = 1.d0

           if(evolve_temper.ne.0) then
              ! set hot atmosphere values
              temperature(i,j,k) = grhydro_hot_atmo_temp
              y_e(i,j,k) = grhydro_hot_atmo_Y_e
              y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
              keytemp = 1
              call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                   press(i,j,k),keyerr,anyerr)
           else
              keytemp = 0
              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)

              call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
           endif

           ! w_lorentz=1, so the expression for tau reduces to:
           tau_avg(i,j,k)  = sdetg(i,j,k) * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens_avg(i,j,k)
           tau(i,j,k)  = sdetg(i,j,k) * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens(i,j,k)

           cycle

         end if

         if(evolve_temper.eq.0) then
          
            !write (*,*) dens_avg(i, j, k), "=?=", dens(i, j, k)
              call Con2Prim_pt(int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),&
                 GRHydro_eos_handle, dens(i,j,k),scon(i,j,k, 1),scon(i,j,k, 2), &
                 scon(i,j,k, 3),tau(i,j,k),rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2), &
                 vup(i,j,k,3),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                 uxx,uxy,uxz,uyy,uyz,uzz,sdetg(i,j,k),x(i,j,k),y(i,j,k), &
                 z(i,j,k),r(i,j,k),epsnegative,GRHydro_rho_min,pmin, epsmin, & 
                 GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
            
                 ! RePrimand call
                 
                 !temp22 = dens(i,j,k)
                 !call GRHydro_RPR_Con2Prim_pt(dens(i,j,k),scon(i,j,k, 1),scon(i,j,k, 2), &
                 !scon(i,j,k, 3),tau(i,j,k), g11(i,j,k), g12(i,j,k), g13(i,j,k), g22(i,j,k), g23(i,j,k), g33(i,j,k), rho(i,j,k), &
                 !eps(i,j,k),press(i,j,k), & 
                 !vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3),&
                 !GRHydro_C2P_failed(i,j,k), w_lorentz(i,j,k))

                 if(GRHydro_C2P_failed(i,j,k).ge.1) then
                  write(*,*) "FAILED AT X:", X(i,j,k), " Y:", y(i,j,k), " Z:", z(i,j,k)
                  write(10,*) "FAILED AT X:", X(i,j,k), " Y:", y(i,j,k), " Z:", z(i,j,k), " dens:", dens(i,j,k), "tau: ", tau(i,j,k), "scon: ", scon(i,j,k,1), ",", scon(i,j,k,2), ",", scon(i,j,k,3),&
                  " g11: ", g11(i,j,k), " g12: ", g12(i,j,k), " g13: ", g13(i,j,k), " g22: ", g22(i,j,k), " g23: ", g23(i,j,k), " g33: ", g33(i,j,k)
                endif
                 !write(10,*) temp22, ",", dens(i,j,k)
                 !foo
                 


          
        else
            call CCTK_ERROR("Should never reach this point!")
            STOP
         
        endif


        if (epsnegative) then  

#if 0
          ! cott 2010/03/30:
          ! Set point to atmosphere, but continue evolution -- this is better than using
          ! the poly EOS -- it will lead the code to crash if this happens inside a (neutron) star,
          ! but will allow the job to continue if it happens in the atmosphere or in a
          ! zone that contains garbage (i.e. boundary, buffer zones)
          ! Ultimately, we want this fixed via a new carpet mask presently under development
!           GRHydro_C2P_failed(i,j,k) = 1

          !$OMP CRITICAL
          call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0! ')
          write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
          call CCTK_WARN(GRHydro_NaN_verbose+2,warnline)
          write(warnline,'(a20,3g16.7)') 'xyz location: ',&
               x(i,j,k),y(i,j,k),z(i,j,k)
          call CCTK_WARN(GRHydro_NaN_verbose+2,warnline)
          write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
          call CCTK_WARN(GRHydro_NaN_verbose+2,warnline)
          call CCTK_WARN(GRHydro_NaN_verbose+2,"Setting the point to atmosphere")
          !$OMP END CRITICAL

          ! for safety, let's set the point to atmosphere
          dens(i,j,k) = ATMOMIN(sdetg(i,j,k)*GRHydro_rho_min, r(i,j,k)) !sqrt(det)*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
          rho(i,j,k) = ATMOMIN(GRHydro_rho_min, r(i,j,k)) !GRHydro_rho_min
          scon(i,j,k,:) = 0.d0
          vup(i,j,k,:) = 0.d0
          w_lorentz(i,j,k) = 1.d0

           call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)

           call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)

          ! w_lorentz=1, so the expression for tau reduces to:
          tau(i,j,k)  = sdetg(i,j,k) * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens(i,j,k)

#else
          ! cott 2010/03/27:      
          ! Honestly, this should never happen. We need to flag the point where
          ! this happened as having led to failing con2prim.

          !$OMP CRITICAL
          call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0, trying polytype.')
          !$OMP END CRITICAL
          call Con2Prim_ptPolytype(GRHydro_polytrope_handle, &
               dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
               tau(i,j,k), rho(i,j,k), vup(i,j,k,1), vup(i,j,k,2), &
               vup(i,j,k,3), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k), &
               uxx, uxy, uxz, uyy, uyz, uzz, sdetg(i,j,k), x(i,j,k), y(i,j,k), &
               z(i,j,k), r(i,j,k),GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
#endif          

        end if


!! Some debugging convenience output
# if 0
         !$OMP CRITICAL
         if (dens(i,j,k).ne.dens(i,j,k)) then
            call CCTK_WARN(1,"NAN at exit of Con2Prim")
            write(warnline,"(A10,i5)") "reflevel: ",GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,"(1P10E15.6)") rho(i,j,k),dens(i,j,k),eps(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_ERROR("Aborting!!!")
            STOP
         endif
         if (rho(i,j,k).ne.rho(i,j,k)) then
            call CCTK_WARN(1,"NAN in rho at exit of Con2Prim")
            write(warnline,"(A10,i5)") "reflevel: ",GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,"(1P10E15.6)") rho(i,j,k),dens(i,j,k),eps(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_ERROR("Aborting!!!")
            STOP
         endif
         !$OMP END CRITICAL
         
         !$OMP CRITICAL
         if (rho(i,j,k).le.0.0d0) then
            call CCTK_WARN(1,"rho less than zero at exit of Con2Prim")
            write(warnline,"(A10,i5)") "reflevel: ",GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,"(1P10E15.6)") rho(i,j,k),dens(i,j,k),eps(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_ERROR("Aborting!!!")
            STOP
         endif
         !$OMP END CRITICAL
#endif
         

      end do
    end do
  end do
  !$OMP END PARALLEL DO
  close(10)

  return
  
end subroutine Conservative2Primitive

/*@@
@routine    Con2Prim_pt
@date       Sat Jan 26 01:09:39 2002
@author     The GRHydro Develpoers    
@desc 
Given all the appropriate data, recover the primitive variables at
a single point.
@enddesc 
 @calls     
@calledby   
@history 
@endhistory 

@@*/

subroutine Con2Prim_pt(cctk_iteration,ii,jj,kk,&
     handle, dens, sx, sy, sz, tau, rho, velx, vely, &
     velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, sdetg, x, y, z, r, epsnegative, GRHydro_rho_min, pmin, epsmin, &
     GRHydro_reflevel, GRHydro_C2P_failed)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL dens, sx, sy, sz, tau, rho, velx, vely, velz, epsilon, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, sdetg, isdetg, w_lorentz, x, &
       y, z, r, GRHydro_rho_min
  
  CCTK_REAL s2, f, df, vlowx, vlowy, vlowz
  CCTK_INT count, i, handle, GRHydro_reflevel
  CCTK_REAL GRHydro_C2P_failed, dummy1, dummy2
  CCTK_REAL udens, usx, usy, usz, utau, pold, pnew, &
            temp1, drhobydpress, depsbydpress, dpressbydeps, dpressbydrho, pmin, epsmin

  CCTK_INT cctk_iteration,ii,jj,kk
  character(len=200) warnline
  logical epsnegative, mustbisect

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr
  CCTK_REAL :: xpress,xtemp,xye,tmp,plow
  n=1;keytemp=0;anyerr=0;keyerr=0
  xpress=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars


!!$  Undensitize the variables 

  isdetg = 1.0d0/sdetg
  udens = dens * isdetg
  usx = sx * isdetg
  usy = sy * isdetg
  usz = sz * isdetg
  utau = tau * isdetg
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.*usx*usz*uxz + 2.*usy*usz*uyz
  
  
!!$  Set initial guess for pressure:
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
                      rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)
  !pold = max(1.d-10,xpress)
  ! This is the lowest admissible pressure, otherwise we are off the physical regime
  ! Sometimes, we may end up with bad initial guesses. In that case we need to start off with something
  ! reasonable at least
  plow = max(pmin, sqrt(s2) - utau - udens)
  
  ! Start out with some reasonable initial guess
  pold = max(plow+1.d-10,xpress)

  mustbisect = .false.

!!$  Check that the variables have a chance of being physical

  if( (utau + pold + udens)**2 - s2 .le. 0.0d0) then

    if (c2p_reset_pressure .ne. 0) then
      pold = sqrt(s2 + c2p_reset_pressure_to_value) - utau - udens
    else 
      !$OMP CRITICAL
      call CCTK_WARN(GRHydro_NaN_verbose, "c2p failed and being told not to reset the pressure")
      call CCTK_ERROR("c2p failed and being told not to reset the pressure")
      !STOP
      GRHydro_C2P_failed = 1
      !$OMP END CRITICAL
    endif
  endif
  
!!$  Calculate rho and epsilon 

!define temporary variables to speed up

  rho = udens * sqrt( (utau + pold + udens)**2 - s2)/(utau + pold + udens)
  w_lorentz = (utau + pold + udens) / sqrt( (utau + pold + udens)**2 - s2)
  epsilon = (sqrt( (utau + pold + udens)**2 - s2) - pold * w_lorentz - &
       udens)/udens
  
!!$  Calculate the function

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
                      rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)

  f = pold - xpress

  if (f .ne. f) then 
     ! Ok, this yielded nonsense, let's enforce bisection!
     mustbisect = .true.
  endif
  
!!$Find the root
  
  count = 0
  pnew = pold
  do while ( ((abs(pnew - pold)/abs(pnew) .gt. GRHydro_perc_ptol) .and. &
       (abs(pnew - pold) .gt. GRHydro_del_ptol))  .or. &
       (count .lt. GRHydro_countmin))
    count = count + 1

    if (count > GRHydro_countmax) then
      GRHydro_C2P_failed = 1

      !$OMP CRITICAL
      call CCTK_WARN(1, 'count > GRHydro_countmax! ')
      write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a28,i8)') 'cctk_iteration:', cctk_iteration
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3i5,3g16.7)') 'ijk, xyz location: ',ii,jj,kk,x,y,z
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'radius: ',r
      call CCTK_WARN(1,warnline)
      call CCTK_WARN(1,"Setting the point to atmosphere")
      !$OMP END CRITICAL


      ! for safety, let's set the point to atmosphere
      SET_ATMO_MIN(rho, GRHydro_rho_min, r)
      udens = rho
      dens = sdetg * rho
      pnew = pmin
      epsilon = epsmin
      ! w_lorentz=1, so the expression for utau reduces to:
      utau  = rho + rho*epsmin - udens
      sx = 0.d0
      sy = 0.d0
      sz = 0.d0
      s2 = 0.d0
      usx = 0.d0
      usy = 0.d0
      usz = 0.d0
      w_lorentz = 1.d0
      goto 51
    end if

    call EOS_Omni_DPressByDRho(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,xtemp,xye,dpressbydrho,keyerr,anyerr)

    call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,xtemp,xye,dpressbydeps,keyerr,anyerr)

    temp1 = (utau+udens+pnew)**2 - s2
    drhobydpress = udens * s2 / (sqrt(temp1)*(udens+utau+pnew)**2)
    depsbydpress = pnew * s2 / (rho * (udens + utau + pnew) * temp1)
    df = 1.0d0 - dpressbydrho*drhobydpress - &
         dpressbydeps*depsbydpress


    pold = pnew
    
    ! Try to obtain new pressure via Newton-Raphson.
    
    pnew = pold - f/df
    
    ! Check if Newton-Raphson resulted in something reasonable!
    
    if (c2p_resort_to_bisection.ne.0) then 
    
      tmp = (utau + pnew + udens)**2 - s2
      plow = max(pmin, sqrt(s2) - utau - udens)
      
      if (pnew .lt. plow .or. tmp .le. 0.0d0 .or. mustbisect) then
      
         ! Ok, Newton-Raphson ended up finding something unphysical.
         ! Let's try to find our root via bisection (which converges slower but is more robust)
      
         pnew = (plow + pold) / 2
         tmp = (utau + pnew + udens)**2 - s2
         
         mustbisect = .false.
      end if
    
    else
      
      ! This is the standard algorithm without resorting to bisection.
    
      pnew = max(pmin, pnew)
      tmp = (utau + pnew + udens)**2 - s2
    
    endif
    
    ! Check if we are still in the physical range now.
    ! If not, we set the C2P failed mask, and set pnew = pold (which will exit the loop).
    
    if ((tmp .le. 0.0d0) .and. GRHydro_C2P_failed .eq. 0) then
      call CCTK_WARN(GRHydro_NaN_verbose, "c2p failed and we are not in physical range")
      GRHydro_C2P_failed = 1
      pnew = pold
    endif

    
!!$    Recalculate primitive variables and function
        
    rho = udens * sqrt( tmp)/(utau + pnew + udens)

!! Some debugging convenience output
#if 0
    if (rho .ne. rho) then
    !$OMP CRITICAL
      write(warnline,'(a38, i16)') 'NAN in rho! Root finding iteration: ', count
      call CCTK_WARN(1,warnline)
      write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'radius: ',r
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'udens: ', udens
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'dens: ', dens
      call CCTK_WARN(1,warnline)
      !write(warnline,'(a20,g16.7)') 'rho-old: ', rhoold
      !call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'pnew: ', pnew
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'pold: ', pold
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'xpress: ', xpress
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'f: ', f
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'df: ', df
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'f/df: ', f/df
      call CCTK_WARN(1,warnline)
      write(warnline,'(a30,g16.7)') '(utau + pnew + udens)**2 - s2: ', (utau + pnew + udens)**2 - s2
      call CCTK_WARN(1,warnline)
      write(warnline,'(a30,g16.7)') '(utau + pnew + udens): ', (utau + pnew + udens)
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'drhobydpress: ', drhobydpress
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'depsbydpress: ', depsbydpress
      call CCTK_ERROR(warnline)
      STOP
   !$OMP END CRITICAL
  endif
#endif

    w_lorentz = (utau + pnew + udens) / sqrt( tmp)
    epsilon = (sqrt( tmp) - pnew * w_lorentz - &
         udens)/udens

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)
    
    f = pnew - xpress

    if (f .ne. f) then 
       ! Ok, this yielded nonsense, let's enforce bisection!
       mustbisect = .true.
    endif

  enddo
  
!!$  Polish the root

  do i=1,GRHydro_polish

    call EOS_Omni_DPressByDRho(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,xtemp,xye,dpressbydrho,keyerr,anyerr)

    call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,xtemp,xye,dpressbydeps,keyerr,anyerr)

    temp1 = (utau+udens+pnew)**2 - s2
    drhobydpress = udens * s2 / (sqrt(temp1)*(udens+utau+pnew)**2)
    depsbydpress = pnew * s2 / (rho * (udens + utau + pnew) * temp1)
    df = 1.0d0 - dpressbydrho*drhobydpress - &
         dpressbydeps*depsbydpress
    pold = pnew
    pnew = pold - f/df
    
!!$    Recalculate primitive variables and function

    rho = udens * sqrt( (utau + pnew + udens)**2 - s2)/(utau + pnew + udens)
    w_lorentz = (utau + pnew + udens) / sqrt( (utau + pnew + udens)**2 - &
         s2)
    epsilon = (sqrt( (utau + pnew + udens)**2 - s2) - pnew * w_lorentz - &
         udens)/udens

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n, &
         rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)
    
    f = pold - xpress


  enddo

!!$  Calculate primitive variables from root

  !if (rho .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
  IF_BELOW_ATMO(rho, GRHydro_rho_min, GRHydro_atmo_tolerance, r) then
    SET_ATMO_MIN(rho, GRHydro_rho_min, r) !GRHydro_rho_min
    udens = rho
    dens = sdetg * rho
!    epsilon = (sqrt( (utau + pnew + udens)**2) - pnew -  udens)/udens
    epsilon = epsmin
    ! w_lorentz=1, so the expression for utau reduces to:
    utau  = rho + rho*epsmin - udens
    sx = 0.d0
    sy = 0.d0
    sz = 0.d0
    s2 = 0.d0
    usx = 0.d0
    usy = 0.d0
    usz = 0.d0
    w_lorentz = 1.d0
  end if

51 press = pnew
  vlowx = usx / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowy = usy / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowz = usz / ( (rho + rho*epsilon + press) * w_lorentz**2)
  velx = uxx * vlowx + uxy * vlowy + uxz * vlowz
  vely = uxy * vlowx + uyy * vlowy + uyz * vlowz
  velz = uxz * vlowx + uyz * vlowy + uzz * vlowz
  
    
!!$If all else fails, use the polytropic EoS

  if(epsilon .lt. 0.0d0) then
    epsnegative = .true.
  endif

#if 0
  if (rho .le. 0.0d0) then
    !$OMP CRITICAL
      write(warnline,'(a38, i16)') 'Epsilon negative!'
      call CCTK_WARN(1,warnline)
      write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'radius: ',r
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'udens: ', udens
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'dens: ', dens
      call CCTK_WARN(1,warnline)
      !write(warnline,'(a20,g16.7)') 'rho-old: ', rhoold
      !call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'pnew: ', pnew
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'pold: ', pold
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'xpress: ', xpress
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'f: ', f
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'df: ', df
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'f/df: ', f/df
      call CCTK_WARN(1,warnline)
      write(warnline,'(a30,g16.7)') '(utau + pnew + udens)**2 - s2: ', (utau + pnew + udens)**2 - s2
      call CCTK_WARN(1,warnline)
      write(warnline,'(a30,g16.7)') '(utau + pnew + udens): ', (utau + pnew + udens)
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'drhobydpress: ', drhobydpress
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'depsbydpress: ', depsbydpress
      call CCTK_ERROR(warnline)
      STOP
   !$OMP END CRITICAL
  endif
#endif

end subroutine Con2Prim_pt


 /*@@
   @routine    Conservative2PrimitiveBoundaries
   @date       Tue Mar 12 18:04:40 2002
   @author     The GRHydro Developers    
   @desc 
        This routine is used only if the reconstruction is performed on the conserved variables. 
        It computes the primitive variables on cell boundaries.
        Since reconstruction on conservative had not proved to be very successful, 
        some of the improvements to the C2P routines (e.g. the check about 
        whether a failure happens in a point that will be restriced anyway) 
        are not implemented here yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Conservative2PrimitiveBounds(CCTK_ARGUMENTS)

  use Con2Prim_fortran_interfaces
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS
  
  CCTK_REAL, parameter :: half = 0.5d0

  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr, pmin, epsmin
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr
  logical :: epsnegative
  character(len=100) warnline
 
  CCTK_REAL :: local_min_tracer

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr
  CCTK_REAL :: xpress,xtemp,xye,xeps
  n=1;keytemp=0;anyerr=0;keyerr=0
  xpress=0.0d0;xtemp=0.0d0;xye=0.0d0;xeps=0.0d0
! end EOS Omni vars
  
  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
  end if

  ! this is a poly call
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
         GRHydro_rho_min,1.0d0,xtemp,xye,pmin,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       GRHydro_rho_min,epsmin,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3) 
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0d0
  end if
  
  do k = GRHydro_stencil, nz - GRHydro_stencil + 1
    do j = GRHydro_stencil, ny - GRHydro_stencil + 1
      do i = GRHydro_stencil, nx - GRHydro_stencil + 1

        !do not compute if in atmosphere or in an excised region
        if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
            GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0)) cycle
         
        gxxl = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))

        epsnegative = .false.

        avg_detl = sqrt(SPATIAL_DETERMINANT(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl))
        avg_detr = sqrt(SPATIAL_DETERMINANT(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr))
        call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl*avg_detl,&
             gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
        call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr*avg_detr,&
             gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), &
                   tracer(i,j,k,itracer), dens(i,j,k))
           enddo

           if (use_min_tracer .ne. 0) then
             if (tracer(i,j,k,itracer) .le. local_min_tracer) then
               tracer(i,j,k,itracer) = local_min_tracer
             end if
           end if

        endif

        if(evolve_Y_e.ne.0) then
           Y_e(i,j,k) = Y_e_con(i,j,k) / dens(i,j,k)
        endif

        call Con2Prim_pt(int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),&
             GRHydro_eos_handle, densminus(i,j,k),&
             sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k),&
             tauminus(i,j,k),rhominus(i,j,k),velxminus(i,j,k),&
             velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),&
             pressminus(i,j,k),w_lorentzminus(i,j,k),&
             uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
             x(i,j,k)-half*CCTK_DELTA_SPACE(1),&
             y(i,j,k)-half*CCTK_DELTA_SPACE(2), &
             z(i,j,k)-half*CCTK_DELTA_SPACE(3),r(i,j,k),&
             epsnegative,GRHydro_rho_min,pmin, epsmin, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
        if (epsnegative) then
          !$OMP CRITICAL
          call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0, trying polytype!')
          !$OMP END CRITICAL
          call Con2Prim_ptPolytype(GRHydro_polytrope_handle, densminus(i,j,k),&
               sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k),&
               tauminus(i,j,k),rhominus(i,j,k),velxminus(i,j,k),&
               velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),&
               pressminus(i,j,k),w_lorentzminus(i,j,k),&
               uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
               x(i,j,k)-half*CCTK_DELTA_SPACE(1),&
               y(i,j,k)-half*CCTK_DELTA_SPACE(2), &
               z(i,j,k)-half*CCTK_DELTA_SPACE(3),r(i,j,k),GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
        end if

        if (epsminus(i,j,k) .lt. 0.0d0) then
          if (GRHydro_reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
            !$OMP CRITICAL
            call CCTK_WARN(1,'Con2Prim: stopping the code.')
            call CCTK_WARN(1, '   specific internal energy just went below 0! ')
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                 x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'velocities: ',&
                 velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_WARN(GRHydro_c2p_warnlevel, "Specific internal energy negative")
            !$OMP END CRITICAL
            exit
          endif
        endif

        epsnegative = .false.
        call Con2Prim_pt(int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),&
             GRHydro_eos_handle, densplus(i,j,k),&
             sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k),&
             tauplus(i,j,k),rhoplus(i,j,k),velxplus(i,j,k),&
             velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),&
             pressplus(i,j,k),w_lorentzplus(i,j,k),&
             uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
             x(i,j,k)+half*CCTK_DELTA_SPACE(1),&
             y(i,j,k)+half*CCTK_DELTA_SPACE(2), &
             z(i,j,k)+half*CCTK_DELTA_SPACE(3),r(i,j,k),&
             epsnegative,GRHydro_rho_min,pmin, epsmin, GRHydro_reflevel,GRHydro_C2P_failed(i,j,k))
        if (epsnegative) then
          !$OMP CRITICAL
          call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0, trying polytype!!')
          !$OMP END CRITICAL
          call Con2Prim_ptPolytype(GRHydro_polytrope_handle, densplus(i,j,k),&
               sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k),&
               tauplus(i,j,k),rhoplus(i,j,k),velxplus(i,j,k),&
               velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),&
               pressplus(i,j,k),w_lorentzplus(i,j,k),&
               uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
               x(i,j,k)+half*CCTK_DELTA_SPACE(1),&
               y(i,j,k)+half*CCTK_DELTA_SPACE(2), &
               z(i,j,k)+half*CCTK_DELTA_SPACE(3),r(i,j,k),GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
        end if

        if (epsplus(i,j,k) .lt. 0.0d0) then
          if (GRHydro_reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
            !$OMP CRITICAL
            call CCTK_WARN(1,'Con2Prim: stopping the code.')
            call CCTK_WARN(1, '   specific internal energy just went below 0! ')
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                 x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'velocities: ',&
                 velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_WARN(GRHydro_c2p_warnlevel, "Specific internal energy negative")
            write(warnline,'(a25,4g15.6)') 'coordinates: x,y,z,r:',&
                 x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
            call CCTK_WARN(1,warnline)
            !$OMP END CRITICAL
          endif
        endif

      end do
    end do
  end do
  
end subroutine Conservative2PrimitiveBounds


 /*@@
   @routine    Con2PrimPolytype
   @date       Tue Mar 19 11:43:06 2002
   @author     Ian Hawke
   @desc 
   All routines below are identical to those above, just
   specialised from polytropic type EOS.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Conservative2PrimitivePolytype(CCTK_ARGUMENTS)
  
  use Con2Prim_fortran_interfaces

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz

  CCTK_REAL :: local_min_tracer
!  character(len=400) :: warnline

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

                                                                               
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0d0
  end if

  
!!$  do k = GRHydro_stencil + 1, nz - GRHydro_stencil
!!$    do j = GRHydro_stencil + 1, ny - GRHydro_stencil
!!$      do i = GRHydro_stencil + 1, nx - GRHydro_stencil
  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx

        !do not compute if in atmosphere or in an excised region
        if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
            GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0)) cycle

        call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdetg(i,j,k)*sdetg(i,j,k),&
             g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),&
             g23(i,j,k),g33(i,j,k))        

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), & 
                   tracer(i,j,k,itracer), dens(i,j,k))
           enddo

           if (use_min_tracer .ne. 0) then
             if (tracer(i,j,k,itracer) .le. local_min_tracer) then
               tracer(i,j,k,itracer) = local_min_tracer
             end if
           end if

        endif

        if(evolve_Y_e.ne.0) then
           Y_e(i,j,k) = Y_e_con(i,j,k) / dens(i,j,k)
        endif

        call Con2Prim_ptPolytype(GRHydro_eos_handle, dens(i,j,k),&
             scon(i,j,k,1),scon(i,j,k,2), &
             scon(i,j,k,3),tau(i,j,k),rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2), &
             vup(i,j,k,3),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
             uxx,uxy,uxz,uyy,uyz,uzz,sdetg(i,j,k),x(i,j,k),y(i,j,k), &
             z(i,j,k),r(i,j,k),GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
                
      end do
    end do
  end do

  !$OMP END PARALLEL DO
  
  return

end subroutine Conservative2PrimitivePolytype


/*@@
@routine    Con2Prim_ptPolytype
@date       Sat Jan 26 01:09:39 2002
@author     The GRHydro Developers   
@desc 
Given all the appropriate data, recover the primitive variables at
a single point.
@enddesc 
 @calls     
@calledby   
@history 

@endhistory 

@@*/

subroutine Con2Prim_ptPolytype(handle, dens, sx, sy, sz, tau, rho, &
     velx, vely, velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, sdetg, x, y, z, r, GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL dens, sx, sy, sz, tau, rho, velx, vely, velz, epsilon, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, sdetg, isdetg, w_lorentz, x, &
       y, z, r, GRHydro_rho_min
  
  CCTK_REAL s2, f, df, vlowx, vlowy, vlowz
  CCTK_INT count, handle, GRHydro_reflevel
  CCTK_REAL udens, usx, usy, usz, rhoold, rhonew, &
       enthalpy, denthalpy, invfac, GRHydro_C2P_failed, dummy1, dummy2
  character(len=200) warnline

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars


!!$  Undensitize the variables 

  isdetg = 1.0d0/sdetg
  udens = dens*isdetg
  usx = sx*isdetg
  usy = sy*isdetg
  usz = sz*isdetg
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.d0*usx*usy*uxy + &
       2.d0*usx*usz*uxz + 2.d0*usy*usz*uyz

!!$  Set initial guess for rho:
  rhoold = max(GRHydro_rho_min,rho)

!  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
!       rhoold,xeps,xtemp,xye,xpress,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhoold,xeps,xtemp,xye,press,xeps,keyerr,anyerr)

  enthalpy = 1.0d0 + xeps + press / rhoold

  w_lorentz = sqrt(1.d0 + s2 / ( (udens*enthalpy)**2 ))

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhoold,epsilon,xtemp,xye,press,keyerr,anyerr)

!!$  Calculate the function

  f = rhoold * w_lorentz - udens

!!$Find the root
  
  count = 0
  rhonew = rhoold

  do while ( ((abs(rhonew - rhoold)/abs(rhonew) .gt. GRHydro_perc_ptol) .and. &
       (abs(rhonew - rhoold) .gt. GRHydro_del_ptol))  .or. &
       (count .lt. GRHydro_countmin))
    count = count + 1

    if (count > GRHydro_countmax) then
      GRHydro_C2P_failed = 1

      !$OMP CRITICAL
      call CCTK_WARN(1, 'count > GRHydro_countmax! ')
      write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'radius: ',r
      call CCTK_WARN(1,warnline)
      call CCTK_WARN(1,"Setting the point to atmosphere")
      !$OMP END CRITICAL

      ! for safety, let's set the point to atmosphere
      SET_ATMO_MIN(rhonew, GRHydro_rho_min, r)
      udens = rhonew
      dens = sdetg * rhonew

      call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
           rhonew,1.0d0,xtemp,xye,press,keyerr,anyerr)

      call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
           rhonew,xeps,xtemp,xye,press,epsilon,keyerr,anyerr)

      sx = 0.d0
      sy = 0.d0
      sz = 0.d0
      s2 = 0.d0
      usx = 0.d0
      usy = 0.d0
      usz = 0.d0
      w_lorentz = 1.d0

      goto 50
    end if

!!$    Ian has the feeling that this is an error in general. But for the
!!$    2D_Polytrope it gives the right answers.

      call EOS_Omni_DPressByDrho(handle,keytemp,GRHydro_eos_rf_prec,n,&
           rhonew,epsilon,xtemp,xye,denthalpy,keyerr,anyerr)
      denthalpy = denthalpy / rhonew

    df = w_lorentz - rhonew * s2 * denthalpy / &
         (w_lorentz*(udens**2)*(enthalpy**3))

    rhoold = rhonew
    rhonew = rhoold - f/df

    if (rhonew .lt. 0.d0) then
#if 0
      GRHydro_C2P_failed = 1
      !$OMP CRITICAL
      call CCTK_WARN(GRHydro_NaN_verbose, 'rhonew went below 0')
      !$OMP END CRITICAL
#else
      rhonew = GRHydro_rho_min/100;
#endif
    end if
    
!!$    Recalculate primitive variables and function

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rhonew,epsilon,xtemp,xye,press,keyerr,anyerr)
    
    call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rhonew,xeps,xtemp,xye,press,epsilon,keyerr,anyerr)

    enthalpy = 1.0d0 + epsilon + press / rhonew

    w_lorentz = sqrt(1.d0 + s2 / ( (udens*enthalpy)**2 ))
    tau = sdetg * ((rho * enthalpy) * w_lorentz**2 - press) - dens
    
    f = rhonew * w_lorentz - udens

  enddo

!!$  Calculate primitive variables from root


  !if (rhonew .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
  IF_BELOW_ATMO(rhonew, GRHydro_rho_min, GRHydro_atmo_tolerance, r) then
    SET_ATMO_MIN(rhonew, GRHydro_rho_min, r) !GRHydro_rho_min
    udens = rhonew
    ! before 2009/10/07 dens was not reset and was used with its negative 
    ! value below; this has since been changed, but the change produces 
    ! significant changes in the results
    dens = sdetg * rhonew

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhonew,1.0d0,xtemp,xye,press,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhonew,xeps,xtemp,xye,press,epsilon,keyerr,anyerr)

    ! w_lorentz=1, so the expression for utau reduces to:
    tau  = rhonew + rhonew*epsilon - udens
    sx = 0.d0
    sy = 0.d0
    sz = 0.d0
    s2 = 0.d0
    usx = 0.d0
    usy = 0.d0
    usz = 0.d0
    w_lorentz = 1.d0
  end if

50  rho = rhonew

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhonew,epsilon,xtemp,xye,press,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhonew,xeps,xtemp,xye,press,epsilon,keyerr,anyerr)

  enthalpy = 1.0d0 + epsilon + press / rhonew

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhonew,epsilon,xtemp,xye,press,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rhonew,xeps,xtemp,xye,press,epsilon,keyerr,anyerr)

  tau = sdetg * ((rho * enthalpy) * w_lorentz**2 - press) - dens

  if (epsilon .lt. 0.0d0) then
    GRHydro_C2P_failed = 1

    !$OMP CRITICAL
    call CCTK_WARN(1, 'epsilon < 0! ')
    write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'radius: ',r
    call CCTK_WARN(1,warnline)
    call CCTK_WARN(1,"Setting the point to atmosphere")
    !$OMP END CRITICAL
    
    rho = GRHydro_rho_min
    dens = sdetg * rho

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rhonew,1.d0,xtemp,xye,press,keyerr,anyerr)
    call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rhonew,xeps,xtemp,xye,press,epsilon,keyerr,anyerr)

    ! w_lorentz=1, so the expression for tau reduces to:
    tau  = sdetg * (rho+rho*epsilon) - dens
    usx = 0.d0
    usy = 0.d0
    usz = 0.d0
    w_lorentz = 1.d0

  end if

  invfac = 1.d0 / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowx = usx * invfac
  vlowy = usy * invfac
  vlowz = usz * invfac
  velx = uxx * vlowx + uxy * vlowy + uxz * vlowz
  vely = uxy * vlowx + uyy * vlowy + uyz * vlowz
  velz = uxz * vlowx + uyz * vlowy + uzz * vlowz
  
  return

end subroutine Con2Prim_ptPolytype


 /*@@
   @routine    Cons2PrimBoundsPolytype
   @date       Tue Mar 12 18:04:40 2002
   @author     The GRHydro Developers
   @desc 
        This routine is used only if the reconstruction is performed on the conserved variables. 
        It computes the primitive variables on cell boundaries.
        Since reconstruction on conservative had not proved to be very successful, 
        some of the improvements to the C2P routines (e.g. the check about 
        whether a failure happens in a point that will be restriced anyway) are not implemented here yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Con2PrimBoundsPolytype(CCTK_ARGUMENTS)
  
  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
  end if
                                                                               
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  do k = GRHydro_stencil, nz - GRHydro_stencil + 1
    do j = GRHydro_stencil, ny - GRHydro_stencil + 1
      do i = GRHydro_stencil, nx - GRHydro_stencil + 1
        
        !do not compute if in atmosphere or in an excised region
        if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
            GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0)) cycle
        
        gxxl = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))

        avg_detl = sqrt(SPATIAL_DETERMINANT(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl))
        avg_detr = sqrt(SPATIAL_DETERMINANT(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr))
        call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl*avg_detl,&
             gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
        call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr*avg_detr,&
             gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        

        call Con2Prim_ptPolytype(GRHydro_eos_handle, densminus(i,j,k),&
             sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k),&
             tauminus(i,j,k),rhominus(i,j,k),velxminus(i,j,k),&
             velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),&
             pressminus(i,j,k),w_lorentzminus(i,j,k),&
             uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
             x(i,j,k)-0.5d0*CCTK_DELTA_SPACE(1),&
             y(i,j,k)-0.5d0*CCTK_DELTA_SPACE(2), &
             z(i,j,k)-0.5d0*CCTK_DELTA_SPACE(3),&
             r(i,j,k),GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))
        call Con2Prim_ptPolytype(GRHydro_eos_handle, densplus(i,j,k),&
             sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k),&
             tauplus(i,j,k),rhoplus(i,j,k),velxplus(i,j,k),&
             velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),&
             pressplus(i,j,k),w_lorentzplus(i,j,k),&
             uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
             x(i,j,k)+0.5d0*CCTK_DELTA_SPACE(1),&
             y(i,j,k)+0.5d0*CCTK_DELTA_SPACE(2), &
             z(i,j,k)+0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed(i,j,k))

        if(evolve_Y_e.ne.0) then
           Y_e(i,j,k) = Y_e_con(i,j,k) / dens(i,j,k)
        endif
        
      end do
    end do
  end do
  
end subroutine Con2PrimBoundsPolytype


 /*@@
   @routine    Con2PrimBoundsTracer
   @date       Mon Mar  8 13:41:55 2004
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2PrimBoundsTracer(CCTK_ARGUMENTS)
  
  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
  end if
 
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  do k = GRHydro_stencil, nz - GRHydro_stencil + 1
    do j = GRHydro_stencil, ny - GRHydro_stencil + 1
      do i = GRHydro_stencil, nx - GRHydro_stencil + 1
        
        gxxl = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))

        avg_detl = sqrt(SPATIAL_DETERMINANT(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl))
        avg_detr = sqrt(SPATIAL_DETERMINANT(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr))
        call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl*avg_detl,&
             gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
        call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr*avg_detr,&
               gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        

        do itracer=1,number_of_tracers
           call Con2Prim_ptBoundsTracer(cons_tracerplus(i,j,k,itracer), &
                tracerplus(i,j,k,itracer), &
                rhoplus(i,j,k), sqrt(1.d0 - &
                (gxxr * velxplus(i,j,k)**2 + &
                gyyr * velyplus(i,j,k)**2 + &
                gzzr * velzplus(i,j,k)**2 + &
                2.d0 * (gxyr * velxplus(i,j,k) * velyplus(i,j,k) + &
                gxzr * velxplus(i,j,k) * velzplus(i,j,k) + &
                gyzr * velyplus(i,j,k) * velzplus(i,j,k) ) ) ), &
                avg_detr)
           call Con2Prim_ptBoundsTracer(cons_tracerminus(i,j,k,itracer), &
                tracerminus(i,j,k,itracer), &
                rhominus(i,j,k), sqrt(1.d0 - &
                (gxxr * velxminus(i,j,k)**2 + &
                gyyr * velyminus(i,j,k)**2 + &
                gzzr * velzminus(i,j,k)**2 + &
                2.d0 * (gxyr * velxminus(i,j,k) * velyminus(i,j,k) + &
                gxzr * velxminus(i,j,k) * velzminus(i,j,k) + &
                gyzr * velyminus(i,j,k) * velzminus(i,j,k) ) ) ), &
                avg_detl)
        enddo

!!$        tracerplus(i,j,k) = cons_tracerplus(i,j,k) / &
!!$             sqrt(avg_detr) / rhoplus(i,j,k) * &
!!$             sqrt(1.d0 - &
!!$                   (gxxr * velxplus(i,j,k)**2 + &
!!$                    gyyr * velyplus(i,j,k)**2 + &
!!$                    gzzr * velzplus(i,j,k)**2 + &
!!$                    2.d0 * (gxyr * velxplus(i,j,k) * velyplus(i,j,k) + &
!!$                            gxzr * velxplus(i,j,k) * velzplus(i,j,k) + &
!!$                            gyzr * velyplus(i,j,k) * velzplus(i,j,k) ) ) )
!!$        tracerminus(i,j,k) = cons_tracerminus(i,j,k) / &
!!$             sqrt(avg_detl) / rhominus(i,j,k) * &
!!$             sqrt(1.d0 - &
!!$                   (gxxl * velxminus(i,j,k)**2 + &
!!$                    gyyl * velyminus(i,j,k)**2 + &
!!$                    gzzl * velzminus(i,j,k)**2 + &
!!$                    2.d0 * (gxyl * velxminus(i,j,k) * velyminus(i,j,k) + &
!!$                            gxzl * velxminus(i,j,k) * velzminus(i,j,k) + &
!!$                            gyzl * velyminus(i,j,k) * velzminus(i,j,k) ) ) )
        
      end do
    end do
  end do
  
end subroutine Con2PrimBoundsTracer


 /*@@
   @routine    Con2Prim_ptTracer
   @date       Mon Mar  8 14:26:29 2004
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2Prim_ptTracer(cons_tracer, tracer, dens) 
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL cons_tracer, tracer, dens

  tracer = cons_tracer / dens

end subroutine Con2Prim_ptTracer


 /*@@
   @routine    Con2Prim_ptTracerBounds
   @date       Mon Mar  8 14:26:29 2004
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2Prim_ptBoundsTracer(cons_tracer, tracer, rho, one_over_w_lorentz, sdet)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL cons_tracer, tracer, rho, one_over_w_lorentz, sdet

  tracer = cons_tracer / sdet / rho * one_over_w_lorentz

end subroutine Con2Prim_ptBoundsTracer


! subroutines to manage the C2P failure mask
subroutine reset_GRHydro_C2P_failed(CCTK_ARGUMENTS)
  
  use Con2Prim_fortran_interfaces

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  GRHydro_C2P_failed = 0

  return

end subroutine reset_GRHydro_C2P_failed


subroutine sync_GRHydro_C2P_failed(CCTK_ARGUMENTS)
! a dummy routine to syncronise GRHydro_C2P_failed
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  return

end subroutine sync_GRHydro_C2P_failed


subroutine check_GRHydro_C2P_failed(CCTK_ARGUMENTS)
  
  use Con2Prim_fortran_interfaces
  use GRHydro_EvolutionMask

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, nx, ny, nz
  character(len=300) warnline
  CCTK_REAL :: dummy1, dummy2
  
  ! we skip points where evolution_mask vanishes
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: evolution_mask
  pointer (evolution_mask_ptr, evolution_mask)
  CCTK_INT :: check_evolution_mask

  call GRHydro_DeclareEvolutionMask(cctkGH, evolution_mask_ptr, &
                                    check_evolution_mask)

!  call CCTK_INFO("Checking the C2P failure mask.")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  ! do not check if told not to
  if(GRHydro_reflevel .lt. GRHydro_c2p_warn_from_reflevel) return

  !$OMP PARALLEL DO PRIVATE(i,j,k, dummy1, dummy2)
  do k = 1, nz 
    do j = 1, ny 
      do i = 1, nx

        if (GRHydro_C2P_failed(i,j,k) .ge. 1) then

           if(con2prim_oct_hack.ne.0.and.&
                (x(i,j,k) .lt. -1.0d-12 .or.&
                y(i,j,k) .lt. -1.0d-12 .or.&
                z(i,j,k) .lt. -1.0d-12)) then
              cycle
           endif

           !if ( rho(i,j,k) .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
           IF_BELOW_ATMO(rho(i,j,k), GRHydro_rho_min, GRHydro_atmo_tolerance, r(i,j,k)) then
               cycle
           end if 

          ! do not collapse conditions since Fortran does not guarantee an order
          if (check_evolution_mask.ne.0) then
            if (abs(evolution_mask(i,j,k)).le.0.5d0) then 
              cycle
            end if
          end if

          !$OMP CRITICAL
          call CCTK_WARN(1,'Con2Prim failed; stopping the code.')
          call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
          write(warnline,'(a32,i2)') 'on carpet reflevel: ',GRHydro_reflevel
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,3g16.7)') 'xyz location: ',&
               x(i,j,k),y(i,j,k),z(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,g16.7)') 'radius: ',r(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,i3,1g16.7)') 'hydro_excision_mask, C2Pfailed: ', hydro_excision_mask(i,j,k), GRHydro_C2P_failed(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,4g16.7)') 'rho eps press w_lorentz: ',&
               rho(i,j,k),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,3g16.7)') 'velocities: ',&
               vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,2g16.7)') 'dens tau: ',&
               dens(i,j,k),tau(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,3g16.7)') 'scon: ',&
               scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3)
          call CCTK_WARN(1,warnline)
          if (evolve_mhd.ne.0) then
            write(warnline,'(a32,3g16.7)') 'magnetic field: ',&
                 Bvec(i,j,k,1),Bvec(i,j,k,2),Bvec(i,j,k,3)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a32,3g16.7)') 'cons magnetic field: ',&
                 Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3)
            call CCTK_WARN(1,warnline)
            if (calculate_bcom.ne.0) then
              write(warnline,'(a32,6g16.7)') 'com magnetic field, b0, b^2: ', bcom(i,j,k,1), bcom(i,j,k,2), bcom(i,j,k,3), bcom0(i,j,k), bcom_sq(i,j,k)
              call CCTK_WARN(1,warnline)
            end if
            if (clean_divergence.ne.0) then
              write(warnline,'(a32,g16.7)') 'psidc: ',psidc(i,j,k)
              call CCTK_WARN(1,warnline)
            end if
            if (track_divB.ne.0) then
              write(warnline,'(a32,g16.7)') 'divB: ',divB(i,j,k)
              call CCTK_WARN(1,warnline)
            end if
          end if
          if (evolve_temper.ne.0) then
            write(warnline,'(a32,2g16.7)') 'temperature, entropy: ',&
                  temperature(i,j,k),entropy(i,j,k)
            call CCTK_WARN(1,warnline)
          end if
          if (evolve_entropy.ne.0) then
            write(warnline,'(a32,2g16.7)') 'entropy: ', entropy(i,j,k)
            call CCTK_WARN(1,warnline)
          end if
          if (evolve_Y_e.ne.0) then
            write(warnline,'(a32,2g16.7)') 'Y_e, Y_e_con: ',&
                  Y_e(i,j,k),Y_e_con(i,j,k)
            call CCTK_WARN(1,warnline)
          end if
          write(warnline,'(a32,6g16.7)') '3-metric: ',&
               gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
          call CCTK_WARN(1,warnline)
          write(warnline,'(a32,4g16.7)') 'lapse, shift: ',&
               alp(i,j,k),betax(i,j,k),betay(i,j,k),betaz(i,j,k)
          call CCTK_WARN(1,warnline)
          if (CCTK_EQUALS(GRHydro_c2p_failed_action, "terminate")) then
            call CCTK_TerminateNext(cctkGH)
          else if (CCTK_EQUALS(GRHydro_c2p_failed_action, "abort")) then
            call CCTK_ERROR("Aborting.")
            STOP
          else
            call CCTK_ERROR("Internal error, unknown action")
            STOP
          end if
          !$OMP END CRITICAL

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  return

end subroutine check_GRHydro_C2P_failed

/*@@
@routine    apply
@date       Fri Jul 5 14:47:00 2024
@author     Dhruv Srivastava
@desc
     Applies the de-averaging step of the WHAM process, for
     improved accuracy in Con2Prim
@enddesc
@calls
@calledby
@history

@endhistory

@@*/

!Apply function as in weno.c++
subroutine apply(data, nx, ny, nz, dirn, a_center_xyz)
  implicit none

  !Declare variables
  integer, intent(in) :: nx, ny, nz, dirn
  real*8, dimension(nx, ny, nz), intent(in) :: data
  real*8, dimension(nx, ny, nz), intent(out) :: a_center_xyz

  real*8, dimension(nx, ny, nz) :: temp_xyz

  real*8, dimension(5) :: A = 0


  real*8 :: A0, A1, A2, A3, A4, beta1, beta2, beta3, wbarplus1, wbarplus2, wbarplus3, iwbarplussum, wplus1, wplus2, wplus3, weno_eps

  integer :: i, j, k, p, q, GRHydro_stencil
  integer, dimension(5,3) :: ijk


  real*8, dimension(3,6) :: beta_shu
  real*8, dimension(3,5) ::  weno_coeffs_de_avg


  !Define beta_shu (same for all in WHAM) and WENO coeffs (for de-averaging here)
  beta_shu = reshape((/4.0/3.0, 4.0/3.0, 10.0/3.0,&
                      -19.0/3.0, -13.0/3.0, -31.0/3.0,&
                      25.0/3.0, 13.0/3.0, 25.0/3.0,&
                      11.0/3.0, 5.0/3.0, 11.0/3.0,&
                      -31.0/3.0, -13.0/3.0, -19.0/3.0,&
                      10.0/3.0, 4.0/3.0, 4.0/3.0/), shape(beta_shu))

  weno_coeffs_de_avg = reshape((/ -1.0/24.0, 0.,0.,&
                                  1.0/12.0, -1.0/24.0, 0.,&
                                  23.0/24.0, 13.0/12.0, 23.0/24.0,&
                                  0.,-1.0/24.0,1.0/12.0,&
                                  0.,0.,-1.0/24.0&
  /), shape(weno_coeffs_de_avg))

  temp_xyz = 0
  GRHydro_stencil = 3
  weno_eps = 1e-10

  i = 0
  j = 0
  k = 0

  loopz: do k = GRHydro_stencil, nz - GRHydro_stencil+1
      loopy: do j = GRHydro_stencil, ny - GRHydro_stencil+1
          loopx: do i = GRHydro_stencil, nx - GRHydro_stencil+1
              !This is in place of the definition of ijk[5] in weno.C++. F90 does not have ternary operators, so I had to do this in a complicated way
              !If there's a better way to do this, please let me know
              select case (dirn)
              case (0) !Solve by x-direction
                  ijk = reshape((/i-2, i-1, i, i+1, i+2, j, j, j, j, j, k, k, k, k, k/), shape(ijk))
              case (1) !Solve by y-direction
                  ijk = reshape((/i, i, i, i, i, j-2, j-1, j, j+1, j+2, k, k, k, k, k/), shape(ijk))
              case (2) !Solve by z-direction
                  ijk = reshape((/i, i, i, i, i, j, j, j, j, j, k-2, k-1, k, k+1, k+2/), shape(ijk))
              end select

              !Again, Cannot define functions inside a function in F90, so yet another convoluted way to do something simple
              A0 = data(ijk(1, 1), ijk(1,2), ijk(1, 3))
              A1 = data(ijk(2, 1), ijk(2,2), ijk(2, 3))
              A2 = data(ijk(3, 1), ijk(3,2), ijk(3, 3))
              A3 = data(ijk(4, 1), ijk(4,2), ijk(4, 3))
              A4 = data(ijk(5, 1), ijk(5,2), ijk(5, 3))

              A = reshape((/A0, A1, A2, A3, A4/), shape(A))



              !Calculate smoothness factors
              beta1  = beta_shu(1,1)*SQR(A(1))&
                + beta_shu(1,2)*A(1)*A(2)&
                + beta_shu(1,3)*SQR(A(2))&
                + beta_shu(1,4)*A(1)*A(3)&
                + beta_shu(1,5)*A(2)*A(3)&
                + beta_shu(1,6)*SQR(A(3))

              beta2  = beta_shu(2,1)*SQR(A(2))&
                + beta_shu(2,2)*A(2)*A(3)&
                + beta_shu(2,3)*SQR(A(3))&
                + beta_shu(2,4)*A(2)*A(4)&
                + beta_shu(2,5)*A(3)*A(4)&
                + beta_shu(2,6)*SQR(A(4))
              beta3  = beta_shu(3,1)*SQR(A(3))&
                + beta_shu(3,2)*A(3)*A(4)&
                + beta_shu(3,3)*SQR(A(4))&
                + beta_shu(3,4)*A(3)*A(5)&
                + beta_shu(3,5)*A(4)*A(5)&
                + beta_shu(3,6)*SQR(A(5))

              wbarplus1 = -9.0/80.0 / SQR(weno_eps + beta1)
              wbarplus2 = 49.0/40.0 / SQR(weno_eps + beta2)
              wbarplus3 = -9.0/80.0 / SQR(weno_eps + beta3)

              iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3)

              wplus1 = wbarplus1 * iwbarplussum
              wplus2 = wbarplus2 * iwbarplussum
              wplus3 = wbarplus3 * iwbarplussum


              !Calculate the reconstruction
              temp_xyz(ijk(3, 1), ijk(3,2), ijk(3, 3)) = 0

              looprecon: do p = 1, 5
                  temp_xyz(ijk(3, 1), ijk(3,2), ijk(3, 3)) = &
                  temp_xyz(ijk(3, 1), ijk(3,2), ijk(3, 3)) + ((wplus1 * weno_coeffs_de_avg(1,p) + &
                  wplus2 * weno_coeffs_de_avg(2,p) + wplus3 * weno_coeffs_de_avg(3,p)) * A(p))
              end do looprecon
          end do loopx
      end do loopy
  end do loopz
  a_center_xyz = temp_xyz

  contains

  !Gives square of real*8 number
  real*8 function SQR(x)
    implicit none
    real*8, intent(in)::x
    SQR = x*x
  end function SQR

end subroutine apply