 /*@@
   @file      GRHydro_ShockTube.F90
   @date      Sat Jan 26 02:53:25 2002
   @author    Ian Hawke
   @desc 
   Initial data of the shock tube type.
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

 /*@@
   @routine    GRHydro_shocktube
   @date       Sat Jan 26 02:53:49 2002
   @author     Ian Hawke
   @desc 
   Initial data for shock tubes. Either diagonal or parallel to
   a coordinate axis. Either Sods problem or the standard shock tube.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Expansion and alteration of the test code from GRAstro_Hydro, 
   written by Mark Miller.
   @endhistory 

@@*/

subroutine GRHydro_shocktube(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  do i=1,nx
    do j=1,ny
      do k=1,nz
        if (CCTK_EQUALS(shocktube_type,"diagshock")) then
          direction = x(i,j,k) - shock_xpos + &
               y(i,j,k) - shock_ypos + z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"xshock")) then
          direction = x(i,j,k) - shock_xpos
        else if (CCTK_EQUALS(shocktube_type,"yshock")) then
          direction = y(i,j,k) - shock_ypos
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
          direction = z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"sphere")) then
          direction = sqrt((x(i,j,k)-shock_xpos)**2+&
                           (y(i,j,k)-shock_ypos)**2+&
                           (z(i,j,k)-shock_zpos)**2)-shock_radius
        end if
        if (CCTK_EQUALS(shock_case,"Simple")) then
          rhol = 10.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 2.d0
          epsr = 1.d-6
        else if (CCTK_EQUALS(shock_case,"Sod")) then
          rhol = 1.d0
          rhor = 0.125d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 1.5d0
          epsr = 1.2d0
        else if (CCTK_EQUALS(shock_case,"Balsaralike1")) then
          rhol = 1.0d0
          rhor = 0.125d0
          velxl = 0.0d0
          velxr = 0.0d0
          velyl = 0.0d0
          velyr = 0.0d0
          velzl = 0.0d0
          velzr = 0.0d0
          epsl = 1.0d0/rhol
          epsr = 0.1d0/rhor
!!$This line only for polytrope, k=1
!!$          epsr = 0.375d0
        else if (CCTK_EQUALS(shock_case,"Blast")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 1500.d0
          epsr = 1.5d-2
        else
          call CCTK_WARN(0,"Shock case not recognized")
        end if

        if ( ((change_shock_direction==0).and.(direction .lt. 0.0d0)).or.& 
             ((change_shock_direction==1).and.(direction .gt. 0.0d0)) ) then
          rho(i,j,k) = rhol
          velx(i,j,k) = velxl
          vely(i,j,k) = velyl
          velz(i,j,k) = velzl
          eps(i,j,k) = epsl
        else
          rho(i,j,k) = rhor
          velx(i,j,k) = velxr
          vely(i,j,k) = velyr
          velz(i,j,k) = velzr
          eps(i,j,k) = epsr
        end if

        det=SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

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
    enddo
    enddo
  enddo

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0

  return
  
end subroutine GRHydro_shocktube

subroutine GRHydro_shocktube_hot(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr, templ, tempr, yel, yer
  CCTK_REAL :: vlowx, vlowy, vlowz, w, tenthalpy

! begin EOS Omni vars                                                                                                  
  CCTK_INT :: n,keytemp,anyerr,keyerr(1)
  CCTK_INT :: handle
  character(len=256) :: warnline
  ! handle for nuclear EOS
  handle=4
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
! end EOS Omni vars      
  
  nx = cctk_ash(1)
  ny = cctk_ash(2)
  nz = cctk_ash(3)
  
  call CCTK_INFO("Setting up initial data for hot shocktube")

  if(.not.CCTK_EQUALS(Y_e_evolution_method,"GRHydro").or.&
       .not.CCTK_EQUALS(temperature_evolution_method,"GRHydro")) then
     call CCTK_WARN(0,"Must have Y_e_evolution_method and temperature_evolution_method set to GRHydro")
  endif

  if (nuceos_read_table.eq.0) then
     call CCTK_WARN(0,"You must read in a nuclear EOS table for initial data shocktube_hot to work!")
  endif

  if (.not.CCTK_EQUALS(GRHydro_eos_table,"nuc_eos").or..not.CCTK_EQUALS(GRHydro_eos_type,"General")) then
     call CCTK_WARN(0,"You must set GRHydro::GRHydro_eos_table = nuc_eos and GRHydro::GRHydro_eos_type = General!")
  endif

  do i=1,nx
    do j=1,ny
      do k=1,nz
        if (CCTK_EQUALS(shocktube_type,"diagshock")) then
          direction = x(i,j,k) - shock_xpos + &
               y(i,j,k) - shock_ypos + z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"xshock")) then
          direction = x(i,j,k) - shock_xpos
        else if (CCTK_EQUALS(shocktube_type,"yshock")) then
          direction = y(i,j,k) - shock_ypos
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
          direction = z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"sphere")) then
          direction = sqrt((x(i,j,k)-shock_xpos)**2+&
                           (y(i,j,k)-shock_ypos)**2+&
                           (z(i,j,k)-shock_zpos)**2)-shock_radius
        end if
        if (CCTK_EQUALS(shock_case,"Simple")) then
          rhol = 1.62e-8
          rhor = 1.62e-9
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          templ = 8.0d0
          tempr = 0.6d0
          yel = 0.48d0
          yer = 0.48d0
        else
          call CCTK_WARN(0,"Shock case not recognized")
        end if

        if ( ((change_shock_direction==0).and.(direction .lt. 0.0d0)).or.& 
             ((change_shock_direction==1).and.(direction .gt. 0.0d0)) ) then
          rho(i,j,k) = rhol
          velx(i,j,k) = velxl
          vely(i,j,k) = velyl
          velz(i,j,k) = velzl
          temperature(i,j,k) = templ
          y_e(i,j,k) = yel
        else
          rho(i,j,k) = rhor
          velx(i,j,k) = velxr
          vely(i,j,k) = velyr
          velz(i,j,k) = velzr
          temperature(i,j,k) = tempr
          y_e(i,j,k) = yer
        end if

        det=SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

        ! call EOS with
        keytemp = 1
        call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
             rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),press(i,j,k),keyerr,anyerr)

        if(anyerr.ne.0) then
           call CCTK_WARN(1,"Error in Initial Data EOS call!")
           write(warnline,"(A10,i8)") "keyerr= ",keyerr
           call CCTK_WARN(0,warnline)
        endif

        ! set up conserved variables
        w = 1.0d0 / &
             sqrt(1.0d0 - (gxx(i,j,k)*vel(i,j,k,1)*vel(i,j,k,1) &
                + gyy(i,j,k)*vel(i,j,k,2)*vel(i,j,k,2) &
                + gzz(i,j,k)*vel(i,j,k,3)*vel(i,j,k,3) ) )

        vlowx = gxx(i,j,k)*vel(i,j,k,1) &
             + gxy(i,j,k)*vel(i,j,k,2)  &
             + gxz(i,j,k)*vel(i,j,k,3)
        vlowy = gxy(i,j,k)*vel(i,j,k,1) &
             + gyy(i,j,k)*vel(i,j,k,2)  &
             + gyz(i,j,k)*vel(i,j,k,3)
        vlowz = gxz(i,j,k)*vel(i,j,k,1) &
             + gyz(i,j,k)*vel(i,j,k,2)  &
             + gzz(i,j,k)*vel(i,j,k,3)


        dens(i,j,k) = sqrt(det)*w*rho(i,j,k)

        tenthalpy = 1.0d0 + eps(i,j,k) + press(i,j,k) / rho(i,j,k)

        tau(i,j,k) = sqrt(det)*( (rho(i,j,k)*(1.0d0+eps(i,j,k))+press(i,j,k))*w*w - press(i,j,k)) &
             - dens(i,j,k)

        w_lorentz(i,j,k) = w

        scon(i,j,k,1) = sqrt(det)*rho(i,j,k)*tenthalpy*(w**2) &
             *vlowx
        scon(i,j,k,2) = sqrt(det)*rho(i,j,k)*tenthalpy*(w**2) &
             *vlowy
        scon(i,j,k,3) = sqrt(det)*rho(i,j,k)*tenthalpy*(w**2) &
             *vlowz
        
        Y_e_con(i,j,k) = dens(i,j,k)*Y_e(i,j,k)
           
     enddo
    enddo
  enddo

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0

  return
  
end subroutine GRHydro_shocktube_hot
