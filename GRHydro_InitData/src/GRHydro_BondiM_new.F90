 /*@@
   @file      GRHydro_Bondi.F90 
   @date      Wed Jan 13 13:00:49 EST 2010
   @author    Scott C. Noble
   @desc 
   Hydro initial data for the relativistic Bondi solution about 
   a single Schwarzschild black hole. 
   @enddesc 
 @@*/

/*
   Calculates the Bondi solution, or the spherically symmetric hydrostationary 
   solution to a fluid on a static fixed background spacetime. We assume that one can
   calculate a radius "r" from the grid and that with respect to this radial coordinate, 
   the solution satisfies 

   d (\rho u^r) / dr    = 0 

   Assumes that the equation of state is  P = K \rho^\Gamma   and K  is set by 
   the location of the sonic point.     


 -- Implicitly assumes that there is no spin in the geometry as there is no Bondi 
    solution for spinning black holes.  If a spin is specified, a spherically symmetric 
    is still assumed but the 4-velocity is set consistently with the spinning spacetime. 
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

# define M_PI     3.14159265358979323846d0  /* pi */

!!$Newton-Raphson parameters:

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

subroutine GRHydro_BondiM_Iso(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz, imin, jb,N_points
  CCTK_REAL :: ONEmTINY, tiny
  PARAMETER (N_points=100000,ONEmTINY=0.999999d0,tiny=1.0d-12)
  CCTK_REAL :: M, Msq, Mdot, rs, gam, rmin_bondi, rmax_bondi, cs_sq,cs,vs_sq,vs,rhos,gmo,hs, Kval, Qdot
  CCTK_REAL :: psonic, riso_s
  CCTK_REAL :: logrmin,dlogr,rhotmp,utmp,vtmp,rspher
  CCTK_REAL :: r_bondi(N_points), logr_bondi(N_points), rho_bondi(N_points), u_bondi(N_points), v_bondi(N_points)
  CCTK_REAL :: drhodr, det, sdet, rhocheck, rhocheck2, riso, rnew, rsch, ucheck
  CCTK_REAL :: uiso, uisocheck, vcheck, ucheck2, vcheck2, xhat,yhat, zhat, xp, yp, zp
  CCTK_REAL :: f,df,ddf,a,b,c,rsm,roverm,dudr,uisocheck2,auiso,buiso
  CCTK_REAL :: bondi_bsmooth, bmag, bsonic, psonicmag

  character(400) :: debug_message

  !!$set_bondi_parameters
  M = bondi_central_mass(1)
  Msq  = M*M
  Mdot = mdot_sonicpt_bondi
  rs = r_sonicpt_bondi
  riso_s = 0.5d0*(rs-M+sqrt(rs*(rs-2.0d0*M)))
  gam = gl_gamma

  write(debug_message,'(a,2f22.14)') "Bondi_pars: M, mdot_sonicpt_bondi,",&
                                                  M,mdot_sonicpt_bondi
  call CCTK_INFO(debug_message)
  write(debug_message,'(a,2f22.14)') "Bondi_pars: r_sonicpt_bondi,gl_gamma", &
                                                    r_sonicpt_bondi,gl_gamma 
  call CCTK_INFO(debug_message)
  write(debug_message,'(a,f22.14)') "Bondi_pars: riso_s", riso_s
  call CCTK_INFO(debug_message)

  rmin_bondi = M * bondi_rmin(1)
  rmax_bondi = M * bondi_rmax(1)

  cs_sq  =  M / ( 2.*rs - 3.*M ) 

  if( cs_sq > (gam - 1.)) then 
     cs_sq = gam - 1.
     rs = 0.5 * M * ( 3. + 1./cs_sq ) 
  endif
  
  cs     =  sqrt(cs_sq)
  vs_sq  =  M / ( 2. * rs )  
  vs     =  sqrt(vs_sq)
  rhos   =  Mdot / ( 4. * M_PI * vs * rs * rs ) 
  gmo  =  gam - 1.
  hs     =  1. / ( 1. - cs_sq / (gam - 1.) )
  Kval      = hs * cs_sq * rhos**(-gmo) / gam 
  Qdot   = hs * hs * ( 1. - 3. * vs_sq ) 
  ! Get the pressure value psonic at the sonic point
  psonic = Kval * rhos**gam
  
  logrmin = log10(rmin_bondi)
  dlogr = (log10(rmax_bondi) - logrmin)/(1.*(N_points-1))

  write(debug_message,'(a,4f22.14)') "Bondi pars: cs,vs,rhos,hs",cs,vs,rhos,hs
  call CCTK_INFO(debug_message)
  write(debug_message,'(a,2f22.14)') "Bondi pars: Kval,Qdot ",Kval,Qdot
  call CCTK_INFO(debug_message)
  write(debug_message,'(a,2f22.14)') "Bondi pars: logrmin,dlogr",logrmin,dlogr
  call CCTK_INFO(debug_message)

  rhotmp=1.0d30
  imin=1

  do i=1,N_points   
     logr_bondi(i) = logrmin + dlogr*(i-1)
     r_bondi(i) = 10.**(logr_bondi(i))
     utmp = abs(r_bondi(i) - r_sonicpt_bondi) 
     if (utmp < rhotmp) then 
        rhotmp = utmp 
        imin = i 
     endif
  enddo

!!$  rhotmp = -1.  !!$ start with guess
  rhotmp=rhos      !!$ start with value at sonic point!


  do i=imin,N_points
     rspher = r_bondi(i) 
     call find_bondi_solution( rspher, rhotmp, utmp, vtmp, rs, rhos, M, Mdot, Kval, gam, Qdot )
     if(rhotmp < initial_rho_abs_min) then
        rhotmp = initial_rho_abs_min
        utmp = Kval * rhotmp**gl_gamma  / (gl_gamma - 1.)
     endif
     rho_bondi(i) = rhotmp
     u_bondi(i)   = utmp
     v_bondi(i)   = vtmp
  end do
  
!!$  rhotmp = -1.
  rhotmp=rhos      !!$ start with value at sonic point!
  
  do i=imin-1,1,-1
     rspher = r_bondi(i)
     call find_bondi_solution( rspher, rhotmp, utmp, vtmp, rs, rhos, M, Mdot, Kval, gam, Qdot )
     if(rhotmp < initial_rho_abs_min) then
        rhotmp = initial_rho_abs_min
        utmp = K * rhotmp**gl_gamma  / (gl_gamma - 1.)
     endif
     rho_bondi(i) = rhotmp
     u_bondi(i)   = utmp
     v_bondi(i)   = vtmp
  enddo

  if (CCTK_MyProc(cctkGH) == 0) then 
    open (47,file="bondi.asc",form="formatted")
    do i=1,N_points   
      write(47,'(i5,4f22.14)')i,r_bondi(i),rho_bondi(i),&
                               u_bondi(i),v_bondi(i)
    end do
    close(47)
  end if

!!$  write(debug_message,'(a,4f22.14)') "i=1:",r_bondi(1),rho_bondi(1),&
!!$                                            u_bondi(1),v_bondi(1)
!!$  call CCTK_INFO(debug_message)
!!$  write(debug_message,'(a,4f22.14)') "i=100:",r_bondi(100),rho_bondi(100),&
!!$                                              u_bondi(100),v_bondi(100)
!!$  call CCTK_INFO(debug_message)
!!$  write(debug_message,'(a,4f22.14)') "i=1000:",r_bondi(1000),rho_bondi(1000),&
!!$                                               u_bondi(1000),v_bondi(1000)
!!$  call CCTK_INFO(debug_message)
!!$  write(debug_message,'(a,4f22.14)') "i=1500:",r_bondi(1500),rho_bondi(1500),&
!!$                                               u_bondi(1500),v_bondi(1500)
!!$  call CCTK_INFO(debug_message)

!!$    // find the derivative near r=M in isotropic coords = r=9/4M in schwarzschild; 
  rnew = 2.25 * M
  j = floor ((log10(rnew) - logrmin) / dlogr + 1.0)
!!$  j = NINT((log10(rnew) - logrmin) / dlogr + 1.0)
  rhocheck = rho_bondi(j)
  call find_bondi_solution(rnew,rhocheck, ucheck, vcheck, rs, rhos, M, Mdot, Kval, gam, Qdot )
  uisocheck = 4.0*vcheck/3.0
  
!!$ the previous point was r=M in isotropic coords = r=9/4M in schwarzschild;  this one is r=1.01M in isotropic
  rnew = 0.25 * 3.02**2 * M/1.01
  j = floor((log10(rnew) - logrmin) / dlogr + 1.0)
!!$ j = NINT((log10(rnew) - logrmin) / dlogr + 1.0)
  rhocheck2 = rho_bondi(j)
  call find_bondi_solution( rnew, rhocheck2, ucheck2, vcheck2, rs, rhos, M, Mdot, Kval, gam, Qdot )
  uisocheck2 = vcheck2 / (1.0 - 1.0/2.02) / (1.0+ 1.0/2.02)

  drhodr = 100.0*(rhocheck2-rhocheck)/M 

!!$ Don't divide by M here, to simplify the math
  dudr   = 100.0*(uisocheck2-uisocheck)

  write(debug_message,'(a,3f22.14)') 'Rhocheck:',rhocheck,rhocheck2,drhodr
  call CCTK_INFO(debug_message)
  write(debug_message,'(a,3f22.14)') 'Ucheck:',uisocheck,uisocheck2,dudr
  call CCTK_INFO(debug_message)

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)



 ! Note that: B^r(t,r) = bondi_bmag M^2 / ( sqrt(\Lambda) \lambda r^2 )
 ! according to eq. 82 of PRD82 084031 (2010). Assuming Schwarzschild 
 ! coordinates B^r = bondi_bmag M^2/r^2 sqrt(1-2M/r). We can show that
 ! b^\mu b_\mu = (B^r)^2 and from the definition of plasma beta parameter
 ! we can find that bondi_bmag = sqrt(2P/\beta) r^2/M 1/sqrt(1-2M/r) 
 if(set_bondi_beta_sonicpt.ne.0) then
  bmag = sqrt(2.0d0*psonic/bondi_beta_sonicpt)*rs**2/M &
              /sqrt(1.0d0-2.0d0*M/rs)
 else
  bmag = bondi_bmag
 end if 
 bsonic = bmag*(M/rs)**2 * sqrt(1.0d0-2.0d0*M/rs)
 psonicmag = 0.5d0*bsonic**2

 write(debug_message,'(a,2f22.14)')'Bondi pars: bondi_bmag,bondi_beta_sonicpt',&
                                                bmag,bondi_beta_sonicpt
 call CCTK_INFO(debug_message)
 write(debug_message,'(a,2f22.14)')'Bondi pars: rs,bsonic',rs,bsonic
 call CCTK_INFO(debug_message)
 write(debug_message,'(a,2f22.14)')'Bondi pars: psonic,psonicmag',&
                                                psonic,psonicmag
 call CCTK_INFO(debug_message)

  do i=1,nx
     do j=1,ny
        do k=1,nz
           
           xp=x(i,j,k)
           yp=y(i,j,k)
           zp=z(i,j,k)
           
           riso = sqrt(xp*xp + yp*yp + zp*zp +1.0e-16)
           xhat = xp/riso
           yhat = yp/riso
           zhat = zp/riso
           roverm = riso/M

           bondi_bsmooth = 1.0d0
           
           if(roverm > ONEmTINY) then
              rsch = 0.25 * ( 2.*riso + M)**2 / riso
!!$              jb = floor( 0.5 +  (log10(rsch) - logrmin) / dlogr ) 
              jb = floor( (log10(rsch) - logrmin) / dlogr  + 1.0)
!!$              jb = NINT((log10(rsch) - logrmin) / dlogr + 1.0)

              if(jb >= N_points) jb = N_points-1

              rhotmp = rho_bondi(jb)+(rho_bondi(jb+1)-rho_bondi(jb))*&
                  (rsch-r_bondi(jb))/(r_bondi(jb+1)-r_bondi(jb))
!!$              call find_bondi_solution_bracket( rsch,rhotmp, utmp, vtmp, rs, rhos, M, Mdot, Kval, gam, Qdot, &
!!$                   rho_bondi(jb),rho_bondi(jb+1)) 
              call find_bondi_solution( rsch,rhotmp, utmp, vtmp, rs, rhos, M, Mdot, Kval, gam, Qdot)
              rho(i,j,k) = rhotmp
              uiso = vtmp / (1.0 - M/2.0/riso) / (1.0+ M/2.0/riso)
           else
              if(roverm > 0.5d0*ONEmTINY) then
                 rho(i,j,k) = rhocheck+drhodr*riso*(riso-M)/M
              else 
                 rho(i,j,k) = (rhocheck-drhodr*M/4.0)*(1.-cos(2.*M_PI*riso/M))/2.0
                 bondi_bsmooth = 8.0d0*riso**3 
              endif
              utmp = Kval * rho(i,j,k)**( gam ) / (gam - 1.)

!!$ match to uiso and dudr at roverm=1
!!$ a R + b R^3 --->   a+b = uisocheck;  a+3b = dudr
!!$ b = (dudr-uisocheck)/2;  a=3*uisocheck-dudr)/2

              auiso = 1.5*uisocheck-0.5*dudr
              buiso = 0.5*dudr-0.5*uisocheck
              uiso =  roverm*(auiso+buiso*roverm**2)
           endif
           eps(i,j,k) = utmp/rhotmp
           
           w_lorentz(i,j,k) = sqrt(1.0+gxx(i,j,k) * uiso**2)
           velx(i,j,k) = -1.0*uiso * xhat / w_lorentz(i,j,k)
           vely(i,j,k) = -1.0*uiso * yhat / w_lorentz(i,j,k)
           velz(i,j,k) = -1.0*uiso * zhat / w_lorentz(i,j,k)
           
           det=SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
           sdet = sqrt(det)

           Bvecx(i,j,k) = bondi_bsmooth*bmag*M**2*xhat/sdet/riso**2
           Bvecy(i,j,k) = bondi_bsmooth*bmag*M**2*yhat/sdet/riso**2
           Bvecz(i,j,k) = bondi_bsmooth*bmag*M**2*zhat/sdet/riso**2


           call Prim2ConGenM(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k), &
                          gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                          det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k), &
                    tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k), &
                          rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k), &
                          eps(i,j,k),press(i,j,k), &
                  Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),w_lorentz(i,j,k))

           if(evolve_entropy.ne.0) then
             entropy(i,j,k) = press(i,j,k) * rho(i,j,k)**(-gmo)
             entropycons(i,j,k) = sdet * entropy(i,j,k) * w_lorentz(i,j,k)
           end if

!!$        write(48,'(3f22.14)')riso,uiso,bondi_bsmooth*bondi_bmag
           
        end do
     end do
  end do

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0
  Bconsrhs = 0.d0
  if(evolve_entropy.ne.0) then
    entropyrhs = 0.d0
  end if

  return

end subroutine GRHydro_BondiM_Iso


