#include <cassert>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
void
GRHydro_Refluxing_SetupVars(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ierr = 0;

  ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::dens"),
      use_MoL_slow_multirate_sector, index_dens);
  ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::scon[0]"),
      use_MoL_slow_multirate_sector, index_sx);
  ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::scon[1]"),
      use_MoL_slow_multirate_sector, index_sy);
  ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::scon[2]"),
      use_MoL_slow_multirate_sector, index_sz);
  ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::tau"),
      use_MoL_slow_multirate_sector, index_tau);
  if(CCTK_Equals(Y_e_evolution_method, "GRHydro")) {
    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::Y_e_con"),
        use_MoL_slow_multirate_sector, index_ye);
  }
  if(CCTK_Equals(Bvec_evolution_method, "GRHydro")) {
    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::Bcons[0]"),
        use_MoL_slow_multirate_sector, index_Bconsx);
    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::Bcons[1]"),
        use_MoL_slow_multirate_sector, index_Bconsy);
    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("GRHydro::Bcons[2]"),
        use_MoL_slow_multirate_sector, index_Bconsz);
  }

  assert(!ierr);
}
