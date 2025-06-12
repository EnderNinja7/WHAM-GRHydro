#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine GRHydro_Refluxing_CaptureFluxes(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: poison = 42.0d+42

  ! Shift the fluxes, so that each cell contains both its
  ! cell-centered value and its left face flux. (GRHydro uses the
  ! opposite convention, where each cell contains its right face
  ! flux.) In other words, here we define that flux_i is the flux to
  ! the left of density_i, whereas GRHydro defines that flux_(i-1) is
  ! the flux to the left of density_i.
  integer, parameter :: grhydro_offset = -1

  integer :: imin(3), imax(3), ioff(3)

  ! Region with valid data
  imin(:) = 1           + cctk_nghostzones(:)
  imax(:) = cctk_lsh(:) - cctk_nghostzones(:)
  ! There is one more flux value in this direction
  imax(flux_direction) = imax(flux_direction) + 1

  ! Offset between source and destination
  ioff(:) = 0
  ioff(flux_direction) = grhydro_offset

  ! Capture fluxes from GRHydro
  call capture(flux(:,:,:,3*index_dens+flux_direction), densflux)
  call capture(flux(:,:,:,3*index_sx  +flux_direction), sxflux  )
  call capture(flux(:,:,:,3*index_sy  +flux_direction), syflux  )
  call capture(flux(:,:,:,3*index_sz  +flux_direction), szflux  )
  call capture(flux(:,:,:,3*index_tau +flux_direction), tauflux )
  if (evolve_Y_e.ne.0) then
     call capture(flux(:,:,:,3*index_ye+flux_direction), Y_e_con_flux)
  end if
  if (evolve_MHD.ne.0) then
     call capture(flux(:,:,:,3*index_Bconsx+flux_direction), Bconsxflux)
     call capture(flux(:,:,:,3*index_Bconsy+flux_direction), Bconsyflux)
     call capture(flux(:,:,:,3*index_Bconsz+flux_direction), Bconszflux)
  end if

contains

  subroutine capture(refluxing, grhydro)
    CCTK_REAL, intent(out) :: refluxing(:,:,:)
    CCTK_REAL, intent(in)  :: grhydro(:,:,:)
    integer, parameter :: rk = kind(refluxing)
    integer   :: di,dj,dk
    integer   :: i,j,k
    CCTK_REAL :: avg_alp

    ! Poison Refluxing flux variable

    !$omp parallel do private(i,j,k)
    do k=1,size(refluxing,3)
       do j=1,size(refluxing,2)
          do i=1,size(refluxing,1)
             refluxing(i,j,k) = poison
          end do
       end do
    end do

    ! Copy interior of GRHydro flux variable to Refluxing flux variable

    di = ioff(1)
    dj = ioff(2)
    dk = ioff(3)
    !$omp parallel do private(i,j,k, avg_alp)
    do k=imin(3),imax(3)
       do j=imin(2),imax(2)
          do i=imin(1),imax(1)
             avg_alp = 0.5_rk * (alp(i,j,k) + alp(i+di,j+dj,k+dk))
             refluxing(i,j,k) = avg_alp * grhydro(i+di,j+dj,k+dk)
          end do
       end do
    end do
  end subroutine capture
end subroutine GRHydro_Refluxing_CaptureFluxes
