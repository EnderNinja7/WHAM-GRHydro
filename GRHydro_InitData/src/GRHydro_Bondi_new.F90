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

subroutine GRHydro_Bondi_Iso(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz, imin, jb,N_points
  CCTK_REAL :: ONEmTINY, tiny
  PARAMETER (N_points=100000,ONEmTINY=0.999999d0,tiny=1.0d-12)
  CCTK_REAL :: M, Msq, Mdot, rs, gam, rmin_bondi, rmax_bondi, cs_sq,cs,vs_sq,vs,rhos,gtemp,hs, Kval, Qdot
  CCTK_REAL :: logrmin,dlogr,rhotmp,utmp,vtmp,rspher
  CCTK_REAL :: r_bondi(N_points), logr_bondi(N_points), rho_bondi(N_points), u_bondi(N_points), v_bondi(N_points)
  CCTK_REAL :: drhodr, det, rhocheck, rhocheck2, riso, rnew, rsch, ucheck
  CCTK_REAL :: uiso, uisocheck, vcheck, ucheck2, vcheck2, xhat,yhat, zhat, xp, yp, zp
  CCTK_REAL :: f,df,ddf,a,b,c,rsm,roverm,dudr,uisocheck2,auiso,buiso

  character(400) :: debug_message

  !!$set_bondi_parameters
  M = bondi_central_mass(1)
  Msq  = M*M
  Mdot = mdot_sonicpt_bondi
  rs = r_sonicpt_bondi
  gam = gl_gamma

  write(debug_message,'(a,4f22.14)') "Bondi_pars: ",M,mdot_sonicpt_bondi, &
                                                    r_sonicpt_bondi,gl_gamma
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
  gtemp  =  gam - 1.
  hs     =  1. / ( 1. - cs_sq / (gam - 1.) )
  Kval      = hs * cs_sq * rhos**(-gtemp) / gam 
  Qdot   = hs * hs * ( 1. - 3. * vs_sq ) 
  
  logrmin = log10(rmin_bondi)
  dlogr = (log10(rmax_bondi) - logrmin)/(1.*(N_points-1))

  write(debug_message,'(a,8f22.14)') "More pars: ",cs,vs,rhos,hs,Kval,Qdot,&
                                                   logrmin,dlogr
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

  if(CCTK_MyProc(cctkGH) .eq. 0) then
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
!!$  call find_bondi_solution_bracket(rnew,rhocheck, ucheck, vcheck, rs, rhos, M, Mdot, Kval, gam, Qdot, &
!!$       rho_bondi(j),rho_bondi(j+1))
  call find_bondi_solution(rnew,rhocheck, ucheck, vcheck, rs, rhos, M, Mdot, Kval, gam, Qdot)
  uisocheck = 4.0*vcheck/3.0
  
!!$ the previous point was r=M in isotropic coords = r=9/4M in schwarzschild;  this one is r=1.01M in isotropic
  rnew = 0.25 * 3.02**2 * M/1.01
  j = floor((log10(rnew) - logrmin) / dlogr + 1.0)
!!$ j = NINT((log10(rnew) - logrmin) / dlogr + 1.0)
  rhocheck2 = rho_bondi(j)
!!$  call find_bondi_solution_bracket( rnew, rhocheck2, ucheck2, vcheck2, rs, rhos, M, Mdot, Kval, gam, Qdot, &
!!$       rho_bondi(j),rho_bondi(j+1))
  call find_bondi_solution( rnew, rhocheck2, ucheck2, vcheck2, rs, rhos, M, Mdot, Kval, gam, Qdot)
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


  write(debug_message,'(a,3f22.14)') 'Lower bound coordinates:',x(1,1,1),y(1,1,1),z(1,1,1)
  call CCTK_INFO(debug_message)
  write(debug_message,'(a,3f22.14)') 'Upper bound coordinates:',x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz)
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

           call Prim2ConGen(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k), &
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k), &
                tau(i,j,k),rho(i,j,k), &
                velx(i,j,k),vely(i,j,k),velz(i,j,k), &
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
           
!           if(riso.gt.1.014d0.and.riso.lt.1.015) then 
           if(z(i,j,k).ge.-1.0d0.and.z(i,j,k).le.1.0d0.and. &
              x(i,j,k).ge.6.5d0.and.x(i,j,k).le.7.5d0.and. &
              y(i,j,k).ge.1.0d0.and.y(i,j,k).le.1.5d0 ) then 
             write(debug_message,'(a,15f22.14)') 'Point to check:', &
               x(i,j,k),y(i,j,k),z(i,j,k),riso,gxx(i,j,k),dens(i,j,k),&
               tau(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),rho(i,j,k),eps(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k)
             call CCTK_INFO(debug_message)
           end if 

        end do
     end do
  end do

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0

  return

end subroutine GRHydro_Bondi_Iso

subroutine find_bondi_solution(r, rho, u, v, rs, rhos, M, Mdot, Kval, gam, Qdot )
  implicit none
  
  CCTK_REAL :: r, rho, u, v, rs, rhos, M, Mdot, Kval, gam, Qdot
  CCTK_REAL :: ur,r_sol, rho_old
  CCTK_REAL :: f, df, dx, x_old, resid, jac
  CCTK_REAL :: errx, x_orig, x
  CCTK_INT :: n_iter, i_extra, doing_extra, keep_iterating, i_increase
  CCTK_REAL :: vp, h, hp, term

  CCTK_REAL :: newt_tol_b,small_bondi
  CCTK_INT :: max_newt_iter_b, extra_newt_iter_b

  max_newt_iter_b = 30
  newt_tol_b = 1.0e-15
  extra_newt_iter_b = 2
  small_bondi = 1.0e-20

!!$  if(r>8.1043 .and. r<8.1044)write(*,*)'init guess:',r,rho
!!$  write(*,*)'init guess:',r,rho

  if (rho < 0.) then  
     if( r > 0.9*rs .and. r < 1.1*rs ) then 
        rho = rhos
     else
        if(r < rs) then
           ur = r**(-0.5)
        else     
           ur = 0.5*r**(-1.5)
        endif
        rho = Mdot / (4.*M_PI * r * r * ur)
     endif
  endif
  
!!$  if(r>8.1043 .and. r<8.1044)
!!$  write(*,*)'init guess:',r,ur,rho,rs,rhos,Mdot
!!$  if(r<1.0001e-10)write(*,*)'init guess:',r,ur,rho,rs,rhos,Mdot
!!$  write(*,*)'init guess:',r,ur,rho,rs,rhos,Mdot

  !!$ set global variables needed by residual function:
  r_sol = r 


!!$ Use Newton's method to find rho:
!!$  gnr_bondi( rho, NEWT_DIM_B, bondi_resid)     
  errx = 1.0
  df=1.0
  f=1.0
  doing_extra = 0

  rho_old=rho
  x=rho

  n_iter = 0

!!$ Start the Newton-Raphson iterations : 
  keep_iterating = 1
  do while( keep_iterating == 1 ) 
     
     hp  =  Kval * gam *  x**(gam - 2.)   !!$   //   dh/drho 
     h   =  1. +  hp * x / ( gam - 1. )
     v   =  Mdot / ( 4. * M_PI * r_sol * r_sol * x )    
     vp  =  -v / x   !!$ //  dv/drho
     term = 1. - 2.*M/r_sol + v*v
     resid  =  -Qdot  +  h * h * term
     jac =  2. * h *( hp*term + h*v*vp )
     dx = -resid / jac
     f = 0.5*resid**2
     df = -2. * f
     
     /* Save old values before calculating the new: */
     errx = 0.
     x_old = x
     
     x=x+dx

!!$     if(r>8.1043 .and. r<8.1044)write(*,*)'iter:',x,dx,resid,jac
     
!!$    /****************************************/
!!$    /* Calculate the convergence criterion */
!!$    /****************************************/
     
!!$    /* For the new criterion, always look at relative error in indep. variable: */
!!$    // METHOD specific:
     if(x==0) then
        errx = abs(dx)
        x = small_bondi
     else
        errx = abs(dx/x)
     endif
     
!!$    /*****************************************************************************/
!!$    /* If we've reached the tolerance level, then just do a few extra iterations */
!!$    /*  before stopping                                                          */
!!$    /*****************************************************************************/

!!$     if(r>8.1043 .and. r<8.1044)write(*,*)'iter3:',errx,newt_tol_b,keep_iterating, &
!!$          doing_extra,i_extra

     if((abs(errx)<=newt_tol_b) .and. (doing_extra == 0) .and. (extra_newt_iter_b > 0)) &
          doing_extra=1

    if( doing_extra == 1 ) i_extra=i_extra+1 
    
    if( ((abs(errx) <= newt_tol_b).and.(doing_extra == 0)) .or. &
         (i_extra > extra_newt_iter_b) .or. (n_iter >= (max_newt_iter_b-1)) ) &
         keep_iterating = 0

!!$     if(r>8.1043 .and. r<8.1044)write(*,*)'iter4:',errx,newt_tol_b,keep_iterating, &
!!$          doing_extra,i_extra

    
    n_iter=n_iter+1
    
 end do

 rho=x
 
!!$ Calculate other quantities:
 u = Kval * rho**( gam ) / (gam - 1.)
 v = Mdot / ( 4. * M_PI * r * r * rho )

!!$ if(r>8.1043 .and. r<8.1044)
!!$ write(*,*)'final:',r,rho,u,v

 
 return
 
end subroutine find_bondi_solution

    subroutine find_bondi_solution_bracket(r, rho, u, v, rs, rhos, M, Mdot, Kval, gam, Qdot,rho1,rho2 )
      implicit none
    
      CCTK_REAL rho1,rho2
      CCTK_REAL :: r, rho, u, v, rs, rhos, M, Mdot, Kval, gam, Qdot
      CCTK_REAL :: ur,r_sol, rho_old
      CCTK_REAL :: f, df, dx, x_old, resid, jac
      CCTK_REAL :: errx, x_orig, x
      CCTK_INT :: n_iter, i_extra, doing_extra, keep_iterating, i_increase
      CCTK_REAL :: vp, h, hp, term
      
      CCTK_REAL :: newt_tol_b,small_bondi
      CCTK_INT :: max_newt_iter_b, extra_newt_iter_b

      if(rho.gt.rho1.or.rho.lt.rho2) then
         write(6,*)'find_bondi_solution_bracket: Very bad rho! (rholow,rhoup,rho) = ',rho1,rho2,rho
         stop
      endif

         write(6,*)'find_bondi_solution_bracket: (rholow,rhoup,rho) = ',rho1,rho2,rho
    
      max_newt_iter_b = 30
      newt_tol_b = 1.0e-15
      extra_newt_iter_b = 2
      small_bondi = 1.0e-20
      
!!$  if(r>8.1043 .and. r<8.1044)write(*,*)'init guess:',r,rho
!!$  write(*,*)'init guess:',r,rho
      
!!$      if (rho < 0.) then  
!!$         if( r > 0.9*rs .and. r < 1.1*rs ) then 
!!$            rho = rhos
!!$         else
!!$            if(r < rs) then
!!$               ur = r**(-0.5)
!!$            else     
!!$               ur = 0.5*r**(-1.5)
!!$            endif
!!$            rho = Mdot / (4.*M_PI * r * r * ur)
!!$         endif
!!$      endif
      
!!$  if(r>8.1043 .and. r<8.1044)
!!$  write(*,*)'init guess:',r,ur,rho,rs,rhos,Mdot
!!$  if(r<1.0001e-10)write(*,*)'init guess:',r,ur,rho,rs,rhos,Mdot
!!$  write(*,*)'init guess:',r,ur,rho,rs,rhos,Mdot
      
!!$ set global variables needed by residual function:

      r_sol = r 
      
      
!!$ Use Newton's method to find rho:
!!$  gnr_bondi( rho, NEWT_DIM_B, bondi_resid)     
      errx = 1.0
      df=1.0
      f=1.0
      doing_extra = 0
      
      rho_old=rho
      x=rho
      
      n_iter = 0
      
!!$ Start the Newton-Raphson iterations : 
      keep_iterating = 1
      do while( keep_iterating == 1 ) 
         
         hp  =  Kval * gam *  x**(gam - 2.)   !!$   //   dh/drho 
         h   =  1. +  hp * x / ( gam - 1. )
         v   =  Mdot / ( 4. * M_PI * r_sol * r_sol * x )    
         vp  =  -v / x   !!$ //  dv/drho
         term = 1. - 2.*M/r_sol + v*v
         resid  =  -Qdot  +  h * h * term
         jac =  2. * h *( hp*term + h*v*vp )
         dx = -resid / jac
         f = 0.5*resid**2
         df = -2. * f
         
         errx = 0.
         x_old = x
         
         x=x+dx

         if(x.gt.rho1.or.x.lt.rho2) then
!!$            write(6,*)'Bad rho! ',rho1,rho2,rho
            if(x.gt.rho1)x=0.5*(x_old+rho1)
            if(x.lt.rho2)x=0.5*(x_old+rho2)
         endif

         if(x.gt.rho1.or.x.lt.rho2) then
            write(6,*)'find_bondi_solution_bracket: Bad rho, bad! (rholow,rhoup,rho) = ',rho1,rho2,rho
            stop
         endif
         
!!$     if(r>8.1043 .and. r<8.1044)write(*,*)'iter:',x,dx,resid,jac
         
!!$    /****************************************/
!!$    /* Calculate the convergence criterion */
!!$    /****************************************/
         
!!$    /* For the new criterion, always look at relative error in indep. variable: */
!!$    // METHOD specific:
         if(x==0) then
            errx = abs(dx)
            x = small_bondi
         else
            errx = abs(dx/x)
         endif
         
!!$    /*****************************************************************************/
!!$    /* If we've reached the tolerance level, then just do a few extra iterations */
!!$    /*  before stopping                                                          */
!!$    /*****************************************************************************/
         
!!$     if(r>8.1043 .and. r<8.1044)write(*,*)'iter3:',errx,newt_tol_b,keep_iterating, &
!!$          doing_extra,i_extra
         
         if((abs(errx)<=newt_tol_b) .and. (doing_extra == 0) .and. (extra_newt_iter_b > 0)) &
              doing_extra=1
         
         if( doing_extra == 1 ) i_extra=i_extra+1 
         
         if( ((abs(errx) <= newt_tol_b).and.(doing_extra == 0)) .or. &
              (i_extra > extra_newt_iter_b) .or. (n_iter >= (max_newt_iter_b-1)) ) then
            keep_iterating = 0
            if(n_iter>=(max_newt_iter_b-1))write(6,*)'find_bondi_solution_bracket: Extra iterations!'
         endif
         
!!$     if(r>8.1043 .and. r<8.1044)write(*,*)'iter4:',errx,newt_tol_b,keep_iterating, &
!!$          doing_extra,i_extra
         
         
         n_iter=n_iter+1
         
      end do
      
      rho=x
      
!!$ Calculate other quantities:
      u = Kval * rho**( gam ) / (gam - 1.)
      v = Mdot / ( 4. * M_PI * r * r * rho )
      
!!$ if(r>8.1043 .and. r<8.1044)
!!$ write(*,*)'final:',r,rho,u,v
      
      
      return
      
    end subroutine find_bondi_solution_bracket

