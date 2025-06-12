 /*@@
   @file      GRHydro_ReconstructTest.F90
   @date      Sat Jan 26 02:51:49 2002
   @author    Luca Baiotti
   @desc 
   A test of the reconstruction algorithm.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    GRHydro_reconstruction_test
   @date       Sat Jan 26 02:52:17 2002
   @author     Luca Baiotti
   @desc 
   Tests the reconstruction method by giving smoothly varying initial
   data with a shock along the x axis.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_reconstruction_test(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer i, j, k, nx, ny, nz, ierr

  logical, dimension(:,:,:), allocatable :: trivial_rp

  allocate(trivial_rp(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Allocation problems with trivial_rp")
  end if
  trivial_rp = .false.

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  xoffset = 1
  yoffset = 0
  zoffset = 0
  
  do k=1, cctk_lsh(3)
    do j=1, cctk_lsh(2)
      do i=1, cctk_lsh(1)
        if (i < 8) then 
          rho(i,j,k) = 5+sin(real(i+j+k)/2.)
        else 
          rho(i,j,k) = 1+sin(real(i+j+k)/2.)
        end if
      end do
    end do
  end do
  
  call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
       rho, rhoplus, rhominus, trivial_rp, hydro_excision_mask)
  deallocate(trivial_rp)
  
  open(unit=11,file='original.dat',status='unknown')
  do  i=1, CCTK_LSH(1)
    write(11,*) i, rho(i,4,7)
    write(*,*) i, rho(i,4,7)
  end do
  close(11)
  
  open(unit=12,file='extensions.dat',status='unknown')
  do  i=2, CCTK_LSH(1)-1
    write(12,*) i-0.5, rhominus(i,4,7)
    write(12,*) i+0.5, rhoplus(i,4,7)
  end do
  close(12)
  
  
  write(*,*) 'test_reconstruction: done tvdreconstruct.'
  write(*,*) 'data written in files original.dat and extension.dat'
    
  stop
  
  return
  
end subroutine GRHydro_reconstruction_test
