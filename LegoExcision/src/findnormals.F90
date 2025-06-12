! $Header$

! Take a mask where the boundary has been been marked.
! Return information about the normals in these locations.

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "maskvalues.h"

subroutine excision_findnormals (ierr, mask, dirx, diry, dirz, ni, nj, nk)

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Arguments
  
  ! out: zero for success, nonzero for error
  integer      :: ierr
  
  ! in: array sizes for grid functions
  !     (you can pass in cctk_lsh(:) for these)
  integer      :: ni,nj,nk

  ! in: mask
  CCTK_REAL    :: mask(ni,nj,nk)
  
  ! out: normal directions to use for interpolation
  CCTK_REAL    :: dirx(ni,nj,nk), diry(ni,nj,nk), dirz(ni,nj,nk)

  ! Internal variables.

  integer i,j,k,ii,jj,kk

  CCTK_REAL sx,sy,sz,smag
  CCTK_REAL vx,vy,vz
  CCTK_REAL dot,dotmax
  
  if (CCTK_IsThornActive(CCTK_THORNSTRING) == 0) then
     call CCTK_WARN (0, "The routine excision_findnormals was called, but thorn " // CCTK_THORNSTRING // " is not active")
  end if

  ! Initialise direction arrays to zero.

  dirx = 0
  diry = 0
  dirz = 0

  ! Loop over all grid points.

  do k=2,nk-1
     do j=2,nj-1
        do i=2,ni-1

           ! Check if the point is on the boundary, if it
           ! is not forget about it.

           if (abs(mask(i,j,k)-MASK_BOUNDARY)<0.01) then

              ! Initialise the vector (sx,sy,sz) to zero.

              sx = 0.0D0
              sy = 0.0D0
              sz = 0.0D0

              ! Loop around nearest neighbours.

              do kk=-1,1
                 do jj=-1,1
                    do ii=-1,1

                       ! If a given neighbour is active, then add its
                       ! relative coordinates to (sx,sy,sz).

                       if (abs(mask(i+ii,j+jj,k+kk)-MASK_ACTIVE)<0.01) then
                          sx = sx + dble(ii)
                          sy = sy + dble(jj)
                          sz = sz + dble(kk)
                       end if

                    end do
                 end do
              end do

              ! Find magnitude of (sx,sy,sz).

              smag = sqrt(sx**2 + sy**2 + sz**2)

              ! Normalize (sx,sy,sz).

              if (smag /= 0.0D0) then
                 sx = sx/smag
                 sy = sy/smag
                 sz = sz/smag
              else
                 ! There is a tie.  Break it.
                 sx = 1
                 sy = 0
                 sz = 0
              end if

              ! Now we have a unit vector pointing in the average
              ! direction on the unmasked neighbours.  This is
              ! as close to the normal direction as we can get.

              ! What we need to do now is find that unmasked
              ! neighbour that is closest to this direction.
              ! For this I loop again over unmasked neighbours,
              ! and find the normalized dot product between
              ! (sx,sy,sz) and the direction of a given neighbour.
              ! The largest dot product will give us our best
              ! normal direction.

              dotmax = 0.0D0

              dirx(i,j,k) = 0.0D0
              diry(i,j,k) = 0.0D0
              dirz(i,j,k) = 0.0D0

              do kk=-1,1
                 do jj=-1,1
                    do ii=-1,1

                       ! If the point is not active, forget about it.

                       if (abs(mask(i+ii,j+jj,k+kk)-MASK_ACTIVE)<0.01) then

                          ! Find normalized dot product.

                          vx = dble(ii)
                          vy = dble(jj)
                          vz = dble(kk)

                          dot = (sx*vx + sy*vy + sz*vz) &
                              / sqrt(vx**2 + vy**2 + vz**2)

!                         If the new product is larger than the
!                         largest so far, redefine the normal.

                          if (dot > dotmax) then

                             dotmax = dot

                             dirx(i,j,k) = vx
                             diry(i,j,k) = vy
                             dirz(i,j,k) = vz

                          end if

                       end if

                    end do
                 end do
              end do

              ! Sanity check.  If the check fails then
              ! something is very wrong.

              ii = int(dirx(i,j,k))
              jj = int(diry(i,j,k))
              kk = int(dirz(i,j,k))

              if (ii==0 .and. jj==0 .and. kk==0) then
                 call CCTK_WARN (0, "Could not decide for a normal direction")
              end if

              if (abs(mask(i+ii,j+jj,k+kk)-MASK_ACTIVE)>0.01) then
                 call CCTK_WARN (0, "Mask boundary layer is too thick")
              end if

           end if
           
        end do
     end do
  end do
  
  ! Success

  ierr = 0
  
end subroutine excision_findnormals
