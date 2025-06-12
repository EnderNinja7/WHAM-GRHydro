#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
void GRHydro_Refluxing_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS

  bool const evol_ye  = CCTK_Equals(Y_e_evolution_method, "GRHydro");
  // NOTE: with the vector potential there is no point in refluxing the B field
  bool const evol_mhd = CCTK_Equals(Bvec_evolution_method, "GRHydro");

  int const expected_nvars = 5 + evol_ye + 3*evol_mhd;
  if(expected_nvars != nvars) {
    char msg[512];
    CCTK_VParamWarn("GRHydro_Refluxing", "GRHydro_Refluxing::nvars should "
        "be %d, but is %d", expected_nvars, nvars);
  }
}
