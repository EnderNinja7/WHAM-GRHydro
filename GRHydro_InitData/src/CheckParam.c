#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void GRHydro_InitData_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (timelevels < 2)
  {
      CCTK_PARAMWARN("You have to set 'HydroBase::timelevels to at least 2");
  }

  if(CCTK_Equals(Bvec_evolution_method,"GRHydro") && 
     ((CCTK_Equals(initial_hydro,"ony_atmo")) ||
      (CCTK_Equals(initial_hydro,"read_conformal")) ||
      (CCTK_Equals(initial_hydro,"simple_wave")) ||
      (CCTK_Equals(initial_data,"con2primtest")) ||
      (CCTK_Equals(initial_data,"reconstruction_test")) ||
      (CCTK_Equals(shocktube_type,"sphere")))) 
    {
      CCTK_PARAMWARN("That test not yet implemented in MHD!");
    }

  if(!CCTK_Equals(Bvec_evolution_method,"GRHydro") && 
     ((CCTK_Equals(shock_case,"Balsara0"   )) ||
      (CCTK_Equals(shock_case,"Balsara1"   )) ||
      (CCTK_Equals(shock_case,"Balsara2"   )) ||
      (CCTK_Equals(shock_case,"Balsara3"   )) ||
      (CCTK_Equals(shock_case,"Balsara4"   )) ||
      (CCTK_Equals(shock_case,"Balsara5"   )) ||
      (CCTK_Equals(shock_case,"Alfven"     )) ||
      (CCTK_Equals(shock_case,"Komissarov1")) ||
      (CCTK_Equals(shock_case,"Komissarov2")) ||
      (CCTK_Equals(shock_case,"Komissarov3")) ||
      (CCTK_Equals(shock_case,"Komissarov4")) ||
      (CCTK_Equals(shock_case,"Komissarov5")) ||
      (CCTK_Equals(shock_case,"Komissarov6")) ||
      (CCTK_Equals(shock_case,"Komissarov7")) ||
      (CCTK_Equals(shock_case,"Komissarov8")) ||
      (CCTK_Equals(shock_case,"Komissarov9")) ||
      (CCTK_Equals(initial_hydro,"cylexp")) 
      ))
    {
      CCTK_PARAMWARN("That test requires MHD!  Set Bvec_evolution_method=GRHYDRO!");
    }

  /* Checks for Bondi solution initial data : */
  if(CCTK_Equals(Bvec_evolution_method,"GRHydro") && 
     CCTK_Equals(initial_hydro,"hydro_bondi_solution") )
    {
      CCTK_PARAMWARN("Please set initial_hydro='magnetized_bondi_solution' instead to initialize the magnetic field correctly!");
    }

  if(!CCTK_Equals(Bvec_evolution_method,"GRHydro") && 
     CCTK_Equals(initial_hydro,"magnetized_bondi_solution") )
    {
      CCTK_PARAMWARN("Please set initial_hydro='hydro_bondi_solution' instead to NOT initialize the magnetic field!");
    }

  if((CCTK_Equals(initial_hydro,"magnetized_bondi_solution") || CCTK_Equals(initial_hydro,"magnetized_bondi_solution_iso")) &&
     !CCTK_Equals(initial_Bvec, "magnetized Bondi") )
    {
      CCTK_PARAMWARN("Please set initial_Bvec='magnetized Bondi' to properly initialize the magnetic field for Bondi solutions!");
    }

  if( CCTK_Equals(initial_hydro,"magnetized_bondi_solution") || CCTK_Equals(initial_hydro,"hydro_bondi_solution")  ) { 
    if( num_bondi_sols > 1 ) { 
      CCTK_PARAMWARN("Currently only one Bondi solution is supported, please change [num_bondi_sols] ");
    }
  }

  if(CCTK_Equals(Bvec_evolution_method,"GRHydro") &&
     CCTK_Equals(initial_hydro,"magnetized_bondi_solution") &&
     CCTK_Equals(entropy_evolution_method,"GRHydro") &&
     !CCTK_Equals(initial_entropy,"magnetized Bondi") ){
      CCTK_PARAMWARN("Please set initial_entropy='magnetized Bondi' in order to initialize the entropy correctly!");
  }

  if(CCTK_Equals(entropy_evolution_method,"GRHydro")){
    *evolve_entropy = 1;
  }else{
    *evolve_entropy = 0;
  }
}
