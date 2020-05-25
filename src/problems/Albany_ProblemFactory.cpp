//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Albany_ProblemFactory.hpp"

#include "Albany_HeatProblem.hpp"
#include "Albany_ThermalProblem.hpp"
#include "Albany_ThermalProblemWithSensitivities.hpp"
#include "Albany_PopulateMesh.hpp"
#include "Albany_SideLaplacianProblem.hpp"

namespace Albany
{

std::string getName(std::string const& key)
{
  if (key.size() < 3) return key;
  return key.substr(0, key.size() - 3);
}
// In "Thermal 3D", extract 3.
int getNumDim(std::string const& key)
{
  if (key.size() < 3) return -1;
  return static_cast<int>(key[key.size() - 2] - '0');
} 



bool BasicProblemFactory::provides (const std::string& key) const 
{
  return key == "Heat 1D" ||
         key == "Heat 2D" ||
         key == "Heat 3D" ||
         key == "Thermal 1D" ||
         key == "Thermal 2D" ||
         key == "Thermal 3D" ||
         key == "Thermal With Sensitivities 1D" ||
         key == "Thermal With Sensitivities 2D" ||
         key == "Thermal With Sensitivities 3D" ||
         key == "Populate Mesh" ||
         key == "Side Laplacian 3D";
}

BasicProblemFactory::obj_ptr_type
BasicProblemFactory::
create (const std::string& key,
        const Teuchos::RCP<const Teuchos_Comm>&     comm,
        const Teuchos::RCP<Teuchos::ParameterList>& topLevelParams,
        const Teuchos::RCP<ParamLib>&               paramLib) const
{
  obj_ptr_type problem;

  auto problemParams = Teuchos::sublist(topLevelParams, "Problem", true);
  auto discParams = Teuchos::sublist(topLevelParams, "Discretization");

  if (key == "Heat 1D") {
    problem = Teuchos::rcp(new HeatProblem(problemParams, paramLib, 1, comm));
  } else if (key == "Heat 2D") {
    problem = Teuchos::rcp(new HeatProblem(problemParams, paramLib, 2, comm));
  } else if (key == "Heat 3D") {
    problem = Teuchos::rcp(new HeatProblem(problemParams, paramLib, 3, comm));
  } else if (getName(key) == "Thermal") {
    problem =
        Teuchos::rcp(new ThermalProblem(problemParams, paramLib, getNumDim(key), comm));
  } else if (getName(key) == "Thermal With Sensitivities") {
    problem =
        Teuchos::rcp(new ThermalProblemWithSensitivities(problemParams, paramLib, getNumDim(key), comm));
  } else if (key == "Populate Mesh") {
    problem = Teuchos::rcp(new PopulateMesh(problemParams, discParams, paramLib));
  } else if (key == "Side Laplacian 3D") {
    problem = Teuchos::rcp(new SideLaplacian(problemParams, paramLib, 1));
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error,
      "Error! Unrecognized key in BasicProblemFactory. Did you forget to check with 'provides(key)' first?\n");
  }

  return problem;
}

} // namespace Albany
