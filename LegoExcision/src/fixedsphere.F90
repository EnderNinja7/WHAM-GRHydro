! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine Lego_FixedSphere(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: ione = 1
  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  
  integer ierr

  if (CCTK_IsThornActive(CCTK_THORNSTRING) == 0) then
     call CCTK_WARN (0, "The routine Lego_FixedSphere was called, but thorn " // CCTK_THORNSTRING // " is not active")
  end if
  
  if (CCTK_EQUALS(fixed_excision,"none")) return
  
  ! Default: no excision
  emask = 1

  ! First excision region
  if (num_fixed_regions >= 1) then

     ! Only excise if there would be more than a single grid point excised
     if (fixed_size > minval(CCTK_DELTA_SPACE(:))) then

        if (CCTK_EQUALS(fixed_excision,"sphere")) then
           where ((x - fixed_origin_x)**2 + (y - fixed_origin_y)**2 + (z - fixed_origin_z)**2 < fixed_size**2)
              emask = 0
           end where
        else if (CCTK_EQUALS(fixed_excision,"cube")) then
           where (max(abs(x - fixed_origin_x), abs(y - fixed_origin_y), abs(z - fixed_origin_z)) < fixed_size)
              emask = 0
           end where
        end if

     end if

  end if

  ! Second excision region
  if (num_fixed_regions >= 2) then

     ! Only excise if there would be more than a single grid point excised
     if (fixed2_size > minval(CCTK_DELTA_SPACE(:))) then

        if (CCTK_EQUALS(fixed_excision,"sphere")) then
           where ((x - fixed2_origin_x)**2 + (y - fixed2_origin_y)**2 + (z - fixed2_origin_z)**2 < fixed2_size**2)
              emask = 0
           end where
        else if (CCTK_EQUALS(fixed_excision,"cube")) then
           where (max(abs(x - fixed2_origin_x), abs(y - fixed2_origin_y), abs(z - fixed2_origin_z)) < fixed2_size)
              emask = 0
           end where
        end if

     end if

  end if
  
  call excision_findboundary (ierr, emask, cctk_lsh(1), cctk_lsh(2), cctk_lsh(3))
  ierr = Boundary_SelectGroupForBC (cctkGH, faces, ione, -ione, "spacemask::mask", "None")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  
end subroutine Lego_FixedSphere
