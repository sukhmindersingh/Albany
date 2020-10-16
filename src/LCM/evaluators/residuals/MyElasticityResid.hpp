//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef MYELASTICITYRESID_HPP
#define MYELASTICITYRESID_HPP

#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_config.hpp"

#include <deal.II/lac/affine_constraints.h>

namespace LCM {
/** \brief Finite Element Interpolation Evaluator

    This evaluator interpolates nodal DOF values to quad points.

*/

template <typename EvalT, typename Traits>
class MyElasticityResid : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>
{
 public:
  MyElasticityResid(Teuchos::ParameterList& p);

  void
  postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& vm);

  void
  evaluateFields(typename Traits::EvalData d);

 private:
  using ScalarT     = typename EvalT::ScalarT;
  using MeshScalarT = typename EvalT::MeshScalarT;

  // Input:
  //PHX::MDField<const ScalarT, Cell, QuadPoint, Dim, Dim>      Stress;
  //PHX::MDField<const MeshScalarT, Cell, Node, QuadPoint, Dim> wGradBF;

  //PHX::MDField<const ScalarT, Cell>                      density;
  //PHX::MDField<const MeshScalarT, Cell, Node, QuadPoint> wBF;
  //PHX::MDField<const ScalarT, Cell, QuadPoint, Dim>      uDotDot;

  // Output:
  PHX::MDField<ScalarT, Cell, Node, Dim> ExResidual;

  PHX::MDField<ScalarT, Cell, Dim> MyResidual;

  dealii::AffineConstraints<double> constraints;

  int  numNodes;
  int  numQPs;
  int  numDims;
  bool enableTransient;
  bool hasDensity;
};
}  // namespace LCM

#endif
