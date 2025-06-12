! $Header$

! Extrapolate the difference between the grid function var and oldvar to
! var where indicated by mask.  Use the normal direction given by dirx,
! diry, and dirz for extrapolation.

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "maskvalues.h"

subroutine excision_extrapolate (ierr, var, oldvar, &
     mask, dirx, diry, dirz, ni, nj, nk, var0)
  
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  ! Arguments
  
  ! out: zero for success, nonzero for error
  integer      :: ierr

  ! in: value of var inside mask.
  CCTK_REAL    :: var0

  ! in: array sizes for grid functions
  !     (you can pass in cctk_lsh(:) for these)
  integer      :: ni,nj,nk
  
  ! inout: grid function that should be interpolated
  CCTK_REAL    :: var(ni,nj,nk)
  
  ! in: other grid function for interpolation
  CCTK_REAL    :: oldvar(ni,nj,nk)
  
  ! in: mask
  CCTK_REAL    :: mask(ni,nj,nk)
  
  ! in: normal directions to use for interpolation
  CCTK_REAL    :: dirx(ni,nj,nk), diry(ni,nj,nk), dirz(ni,nj,nk)

  integer i,j,k
  integer ii,jj,kk
  
  if (CCTK_IsThornActive(CCTK_THORNSTRING) == 0) then
     call CCTK_WARN (0, "The routine excision_extrapolate was called, but thorn " // CCTK_THORNSTRING // " is not active")
  end if
  
  do k=2,nk-1
     do j=2,nj-1
        do i=2,ni-1
           if (abs(mask(i,j,k)-MASK_BOUNDARY)<0.01) then
              ii = i + int(dirx(i,j,k))
              jj = j + int(diry(i,j,k))
              kk = k + int(dirz(i,j,k))
              var(i,j,k) = oldvar(i,j,k) + var(ii,jj,kk) - oldvar(ii,jj,kk)
           else if (abs(mask(i,j,k)-MASK_EXCISED)<0.01) then
              var(i,j,k) = var0
           end if
        end do
     end do
  end do

  ierr = 0

end subroutine excision_extrapolate
