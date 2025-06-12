! $Header$

! Take a mask that contains only the values 0 and 1.
! Return a mask where the outermost 0s have been replaced by 0.5s.

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "maskvalues.h"

subroutine excision_findboundary (ierr, mask, ni, nj, nk)
  
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Arguments
  
  ! out: zero for success, nonzero for error
  integer      :: ierr

  ! in: array sizes for grid functions
  !     (you can pass in cctk_lsh(:) for these)
  integer      :: ni,nj,nk
  
  ! inout: mask
  CCTK_REAL    :: mask(ni,nj,nk)

! Internal variables.

  integer i,j,k
  integer ii,jj,kk
  logical bnd
  
  if (CCTK_IsThornActive(CCTK_THORNSTRING) == 0) then
     call CCTK_WARN (0, "The routine excision_findboundary was called, but thorn " // CCTK_THORNSTRING // " is not active")
  end if

! Loop over grid points.

  do k=2,nk-1
     do j=2,nj-1
        do i=2,ni-1

!          Check if we are in an excised point.

           if (abs(mask(i,j,k)-MASK_EXCISED)<0.01) then

              bnd = .false.

!             If any neighbour is active, mask the point
!             as a boundary point.

              do kk=k-1,k+1
                 do jj=j-1,j+1
                    do ii=i-1,i+1
                       bnd = bnd.or.abs(mask(ii,jj,kk)-MASK_ACTIVE)<0.01
                    end do
                 end do
              end do

              if (bnd) mask(i,j,k) = MASK_BOUNDARY

           end if

        end do
     end do
  end do
  
  ierr = 0
  
end subroutine excision_findboundary
