 /*@@
   @file      GRHydro_TOV.F90
   @date      Wed Feb 13 02:53:25 2002
   @author    Scott Hawley
   @desc 
   Initial data of the TOV type.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    GRHydro_TOV
   @date       Sat Jan 26 02:53:49 2002
   @author     Scott Hawley
   @desc 
   Initial data for TOV stars.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Alteration of the GRHydro_shockwave written by Ian Hawke
   @endhistory 

@@*/

subroutine GRHydro_TOV(CCTK_ARGUMENTS)

  implicit  none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  real*8   rho0
  common / com_tov / rho0
  
  integer   nr, n_succ
  integer i,j,k,nx,ny,nz
  integer CCTK_Equals
  CCTK_REAL rmax, dr
  CCTK_REAL xlb, xub
  CCTK_REAL ylb, yub
  CCTK_REAL zlb, zub
  CCTK_REAL tol
  parameter ( tol = 1e-10 )

  CCTK_REAL, allocatable, dimension (:) :: gtt, grr

  CCTK_REAL det
  integer   rc

#include "EOS_Base.inc"

#ifdef HAVE_ODEPACK
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

!-----------------------------------------------------------
! Set up the 1D spherical grid
!-----------------------------------------------------------
!   decide maximum radius to integrate out to:
!           Just use the maximum diagonal width of the grid
!  BUG: this may not be the proper procedure call
    call CCTK_CoordRange(cctkGH, xlb, xub, -1, "x", "cart3d")
    call CCTK_CoordRange(cctkGH, ylb, yub, -1, "y", "cart3d")
    call CCTK_CoordRange(cctkGH, zlb, zub, -1, "z", "cart3d")
    rmax = sqrt( (xub-xlb)**2 + (yub-ylb)**2 + (zub-zlb)**2 )

    nr = tov_nr
    dr = rmax / (nr - 1)   
    do j = 1, nr
      tov_r(j) = (j-1) * dr
    enddo

    
!-----------------------------------------------------------
!  Assign data at r=0
!-----------------------------------------------------------
    rho0 = GRHydro_rho_central
    tov_rho(1) = rho0
    tov_alpha(1) = 1.0
    tov_m(1)   = 0.0

!-----------------------------------------------------------
! Get the 1D TOV data
!-----------------------------------------------------------
      call tov_lsoda (tov_rho, tov_p, tov_m, tov_alpha, &
                tov_r,  nr, tol , n_succ, rc)

!-----------------------------------------------------------
! Handle any errors from integration routine
!-----------------------------------------------------------
   if (rc .ne. 0) then 
     call CCTK_WARN(0,"Fatal error in TOV integration")
   endif
   if (n_succ .lt. nr) then 
     call CCTK_WARN(0,"TOV integration did not extend as far out as requested")
   endif
   nr = n_succ

!-----------------------------------------------------------
! fill in additional physical info
!-----------------------------------------------------------
  do j=1,nr
!    BUG: the call looks something like this...
    tov_rho(j) = EOS_RestMassDens(GRHydro_eos_handle,tov_p(j))
  end do

!-----------------------------------------------------------
! Generate 3D data from 1D data
!  BUG: These Spin1D functions do not yet exist!
!-----------------------------------------------------------
  call Spin1D(tov_rho, tov_r, nr, rho, x, y, z, nx, ny, nz)
  call Spin1D(tov_p, tov_r, nr, press, x, y, z, nx, ny, nz)
  call Spin1D(tov_rho, tov_r, nr, eps, x, y, z, nx, ny, nz)

  allocate(grr(nr), gtt(nr))
  do i = 1, nr
     grr(i) = tov_r(i) / (tov_r(i) - 2*tov_m(i))
     gtt(i) = - tov_alpha(i)**2
  enddo
  call Spin1DMetric(tov_m, tov_r, nr, gxx, gxy, gxz, gyy, gyz, gzz, &
      x, y, z, nx, ny, nz)

  call Spin1DMetric(grr, gtt, gxx, gxy, gxz, gyy, gyz, gzz, nx, ny, nz, &
                          x, y, z, r)
  deallocate(gtt)
  deallocate(grr)

  do i=1,nx
    do j=1,ny
      do k=1,nz
          velx(i,j,k) = 0.0
          vely(i,j,k) = 0.0
          velz(i,j,k) = 0.0
          w_lorentz(i,j,k) = 1.0
      end do
    end do
  enddo

!-----------------------------------------------------------
! Set additional variables
!-----------------------------------------------------------
! Conserved Variables
  det = 1.0d0
  call prim2con(GRHydro_eos_handle,gxx, gxy, gxz, gyy, gyz, gzz, det, &
       dens, sx, sy, sz, tau, rho, &
       velx, vely, velz, eps, press, w_lorentz)

! Copy data to other time steps
  do i=1,nx
    do j=1,ny
      do k=1,nz

!        Primitive variables
         rho_p(i,j,k) = rho(i,j,k)
         press_p(i,j,k) = press(i,j,k)
         eps_p(i,j,k) = eps(i,j,k)
         velx_p(i,j,k) = velx(i,j,k)
         vely_p(i,j,k) = vely(i,j,k)
         velz_p(i,j,k) = velz(i,j,k)
         w_lorentz_p(i,j,k) = w_lorentz(i,j,k)

!        Conserved variables
         dens_p(i,j,k) = dens(i,j,k)
         tau_p(i,j,k) = tau(i,j,k)
         sx_p(i,j,k) = sx(i,j,k)
         sy_p(i,j,k) = sy(i,j,k)
         sz_p(i,j,k) = sz(i,j,k)

      end do
    end do
  end do
  
! BUG:  Why am I doing this?  
  densrhs = 0.d0
  sxrhs = 0.d0
  syrhs = 0.d0
  szrhs = 0.d0
  taurhs = 0.d0

  return
  
end subroutine GRHydro_TOV



! =======================================================================

      subroutine Spin1D(rad_func, r, nr, cart_func, x, y, z, &
                        nx, ny, nz)
!     ................................................................
!   Takes a "radial grid function" and creates a 3D grid function
!   BUG/FEATURE: This doesnt do any fancy interpolation at all!
!
      implicit none

      integer nr, nx, ny, nz
      CCTK_REAL   rad_func(nr), r(nr)
      CCTK_REAL   cart_func(nx,ny,nz),  x(nx,ny,nz), y(nx,ny,nz), z(nx,ny,nz)

      integer i,j,k, ri
      real    rloc

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

!      Given a point in 3D cartesian space, substitute next-further-out data 
!      value in 1D spherical space
              rloc = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
              ri = 1
100           continue
              if (r(ri) .lt. rloc) then
                 ri = ri+1
              endif

              cart_func(i,j,k) =  rad_func(ri)
          end do
        end do
      end do

      end subroutine Spin1D



! =======================================================================


       subroutine Spin1DMetric(grr,gtt,gxx,gxy, &
          gxz,gyy,gyz,gzz,nx,ny,nz,x,y,z,r)

!     ................................................................
!     this subroutine converts spherical metric components to cartesian
!
!    This code was totally ripped off of Francisco Guzmans spheretocart
!    code.  Used with his persmission. 

      implicit NONE
      integer nx,ny,nz
      integer i,j,k

      CCTK_REAL x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz),r(nx,ny,nz)
      CCTK_REAL grr(nx,ny,nz),gtt(nx,ny,nz)
      CCTK_REAL gxx(nx,ny,nz),gxy(nx,ny,nz),gxz(nx,ny,nz)
      CCTK_REAL gyy(nx,ny,nz),gyz(nx,ny,nz),gzz(nx,ny,nz)

      CCTK_REAL xp,yp,zp,rp

!     define derivatives drx = (dr/dx)
      CCTK_REAL drx,dry,drz
      CCTK_REAL dtx,dty,dtz
      CCTK_REAL dpx,dpy,dpz

! this equation for gxx gyy gzz assumes that the metric is diagonal
! and that grr = 1 at the origin, but should compute for all points
! including the origin and the line x=y=0.


      do k = 1, nz
         do j = 1, ny
            do i = 1, nx

               xp = x(i,j,k)
               yp = y(i,j,k)
               zp = z(i,j,k)
               rp = r(i,j,k)

               if (rp.eq.0) then

                  gxx(i,j,k) = 1.
                  gyy(i,j,k) = 1.
                  gzz(i,j,k) = 1.
                  gxy(i,j,k) = 0.
                  gxz(i,j,k) = 0.
                  gyz(i,j,k) = 0.

               else

                  gxx(i,j,k) = (xp/rp)**2 * grr(i,j,k) + &
                       (zp**2 + yp**2)/(rp**2)
                  gyy(i,j,k) = (yp/rp)**2 * grr(i,j,k) + &
                       (xp**2 + zp**2)/(rp**2)
                  gzz(i,j,k) = (zp/rp)**2 * grr(i,j,k) + &
                       (xp**2 + yp**2)/(rp**2)

                  if ((xp.eq.0) .and. (yp.eq.0)) then

                     gxy(i,j,k) = 0.
                     gyz(i,j,k) = 0.
                     gxz(i,j,k) = 0.

                  else

                     gxy(i,j,k) = (xp * yp/(rp**2) * grr(i,j,k) &
                         + (xp*yp*zp**2)/(rp**2*(xp**2+yp**2)) &
                         - (yp*xp)/(xp**2+yp**2))

                     gxz(i,j,k) = (xp * zp/(rp**2) &
                         * grr(i,j,k) - (xp * zp/(rp**2))) 

                     gyz(i,j,k) = (yp * zp/(rp**2) &
                         * grr(i,j,k) - (yp * zp/(rp**2)))

                  endif

               endif

            enddo
          enddo
        enddo
#else
        write(6,*) 'GRHydro_TOV: Does nothing.  Recompile with odepack'
#endif

      return
      end



