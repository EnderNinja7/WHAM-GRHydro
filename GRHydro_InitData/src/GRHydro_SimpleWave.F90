
 /*@@
   @file      GRHydro_SimpleWave.F90
   @date      Thu Aug  2 15:17:35 2007
   @author    Luca Baiotti
   @desc 
   Initial data for a simple wave with sinusoidal initial function for the velocity
     See Anile, Miller, Motta, Formation and damping of relativistic strong shocks,  
         Phys. Fluids 26, 1450 (1983)
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

 /*@@
   @routine    GRHydro_SimpleWave
   @date       Thu Aug  2 15:20:28 2007
   @author     Luca Baiotti
   @desc 
   Initial data for a simple wave with sinusoidal initial function for the velocity
     See Anile, Miller, Motta, Formation and damping of relativistic strong shocks,  
         Phys. Fluids 26, 1450 (1983)
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 

@@*/

subroutine GRHydro_SimpleWave(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: dr, k1, k2, k3, k4, in_data, old_data, source_data, new_data, c_0, det, pi

! begin EOS Omni vars
  CCTK_INT  :: n = 1
  CCTK_INT  :: keytemp = 0
  CCTK_INT  :: anyerr = 0
  CCTK_INT  :: keyerr(1) = 0
  CCTK_REAL :: xpress(1) = 0.0d0
  CCTK_REAL :: xeps(1) = 0.0d0
  CCTK_REAL :: xtemp(1) = 0.0d0
  CCTK_REAL :: xye(1) = 0.0d0
  CCTK_REAL :: rf_precision = 1.0d-10
! end EOS Omni vars
  
  call CCTK_INFO("Setting initial data for a simple wave as Anile Miller Motta")

  call CCTK_WARN(1, "The simple-wave initial-data routine works only for unigrid and on one node.")
  
  pi = 4.d0 * atan(1.0)

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  velx(:,:,:) = 0.0
  rho = 0.0
  eps = 0.0
  press = 0.0
  w_lorentz = 0.0
  simple_tmp = 0.0

  i = 0
  do while (minval(velx(:,1,1)) >= 0.0)
    i=i+1
    if (x(i,1,1) > 0.0) velx(i,1,1) = simple_wave_v_max * sin(pi * x(i,1,1))
  end do

  velx(i,1,1) = 0.0 ! set the first term that became negative to zero

  c_0 = simple_wave_constant_c_0 ! a parameter

    ! compute quantities at v=0

  simple_eps_0 = 9.d0*c_0**2 /(4.d0 * (1.d0 - 3.d0 * c_0**2))

  simple_rho_0 = simple_eps_0**3.d0

!  press_0 = eps_0**4.d0
  
!   do j = 1,i
   do j = 1,CCTK_LSH(1)
  
    simple_tmp(j,1,1) = ( 1.d0 + c_0*sqrt(3.d0) )/( 1.d0 - c_0*sqrt(3.d0) ) &
         * ( (1+velx(j,1,1))/(1-velx(j,1,1)) )**(1/(2.d0*sqrt(3.d0)))
    
    c_s(j,1,1) = ( simple_tmp(j,1,1) - 1.d0 ) / ( sqrt(3.d0)* ( simple_tmp(j,1,1) + 1.d0 ) )
    
    eps(j,1,1) = 9.d0*c_s(j,1,1)**2 /(4.d0 * (1.d0 - 3.d0 * c_s(j,1,1)**2))    
   
    rho(j,1,1) = eps(j,1,1)**3.d0
! write(*,*) j, x(j,1,1), rho(j,1,1)
!    rho(j,1,1) = rho_abs_min * rho(j,1,1)
    
    press(j,1,1) = eps(j,1,1)**(4.d0)
    
    w_lorentz(j,1,1) = 1.d0/sqrt(1.d0-velx(j,1,1)**2) ! flat spacetime
    
  end do

  !arrays
  gxx = 1.d0; gyy=1.d0; gzz=1.d0; gxy=0.d0; gxz=0.d0; gyz=0.d0
  
  

!!$  do i = 1,CCTK_LSH(1)
!!$    
!!$    write(*,*) i, x(i,1,1), velx(i,1,1), rho(i,1,1)
!!$    write(*,*) eps(i,1,1), press(i,1,1), w_lorentz(i,1,1)
!!$
!!$  end do



  do i=1,nx
    
    ! atmosphere

     if ( (rho(i,1,1) < GRHydro_rho_min).OR.(velx(i,1,1) < 0) ) then
        rho(i,1,1) = rho_abs_min
        !      rho(i,1,1) = 1.0 !the value of rho_min for the initial data
        eps(i,1,1) = rho_abs_min**(1.d0/3.d0)
        velx(i,1,1) = 0.d0
        w_lorentz(i,1,1) = 1.d0
        
        xeps = 1.0d0
        call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,rf_precision,n,&
             rho(i,1,1),xeps,xtemp,xye,press(i,1,1),keyerr,anyerr)
        
        ! polytrope only (initial data)
     end if

!    write(*,*) 'p',i, x(i,1,1), rho(i,1,1)**(4.d0/3.d0)/press(i,1,1)    
    
    call SpatialDet(gxx(i,1,1),gxy(i,1,1),gxz(i,1,1),&
         gyy(i,1,1),gyz(i,1,1),gzz(i,1,1),det)
    
!    if (CCTK_EQUALS(GRHydro_eos_type,"Polytype")) then
    ! always use polytype for initial data
      call Prim2ConPoly(GRHydro_polytrope_handle,gxx(i,1,1),gxy(i,1,1),&
           gxz(i,1,1),gyy(i,1,1),gyz(i,1,1),gzz(i,1,1),&
           det, dens(i,1,1),sx(i,1,1),sy(i,1,1),sz(i,1,1),&
           tau(i,1,1),rho(i,1,1),&
           velx(i,1,1),vely(i,1,1),velz(i,1,1),&
           eps(i,1,1),press(i,1,1),w_lorentz(i,1,1))
!!$    else
!!$      call Prim2ConGen(GRHydro_eos_handle,gxx(i,1,1),gxy(i,1,1),&
!!$           gxz(i,1,1),gyy(i,1,1),gyz(i,1,1),gzz(i,1,1),&
!!$           det, dens(i,1,1),sx(i,1,1),sy(i,1,1),sz(i,1,1),&
!!$           tau(i,1,1),rho(i,1,1),&
!!$           velx(i,1,1),vely(i,1,1),velz(i,1,1),&
!!$           eps(i,1,1),press(i,1,1),w_lorentz(i,1,1))
!!$    end if
!!$    

!      write(*,*) 'd',i, x(i,1,1), rho(i,1,1)**(4.d0/3.d0)/press(i,1,1)

  enddo
  
  
  ! planar symmetry
  do j=1,ny
    do k=1,nz
      
      velx(:,j,k)   = velx(:,1,1)
      rho(:,j,k)    = rho(:,1,1)
      eps(:,j,k)    = eps(:,1,1)
      press(:,j,k)  = press(:,1,1)
      w_lorentz(:,j,k)  = w_lorentz(:,1,1)

      sx(:,j,k)   = sx(:,1,1)
      dens(:,j,k)   = dens(:,1,1)
      tau(:,j,k)   = tau(:,1,1)
     
    enddo
  enddo
  
  
  vely(:,:,:) = 0.d0
  velz(:,:,:) = 0.d0
  
  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0

  simple_tmp = rho



!!$  do i = 1,CCTK_LSH(1)
!!$
!!$    do j = 1,CCTK_LSH(2)
!!$
!!$      do k = 1,CCTK_LSH(3)
!!$
!!$        write(*,*) i, x(i,j,k), velx(i,j,k)
!!$        write(*,*) eps(i,j,k), press(i,j,k), rho(i,j,k)
!!$
!!$      end do
!!$    end do
!!$  end do


call CCTK_INFO("Finished initial data")

  return
  
end subroutine GRHydro_SimpleWave


subroutine GRHydro_SimpleWave_Analysis(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (velx(CCTK_LSH(1),1,1) == 0.0) then
    simple_eps_0 = eps(CCTK_LSH(1),1,1)
  else
    call CCTK_WARN(1,"The wave has reached the outer boundary: the computation of simple_eps is now wrong")
  end if

  simple_rho = (rho - simple_rho_0)/simple_rho_0
  simple_eps = (eps - simple_eps_0)/simple_eps_0


  return

end subroutine GRHydro_SimpleWave_Analysis
