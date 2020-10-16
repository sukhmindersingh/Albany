//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Phalanx_DataLayout.hpp"
#include "Teuchos_TestForException.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

//#define HARD_CODED_BODY_FORCE_ELASTICITY_RESID

namespace LCM {

//**********************************************************************
template <typename EvalT, typename Traits>
MyElasticityResid<EvalT, Traits>::MyElasticityResid(Teuchos::ParameterList& p)
    : ExResidual(
          p.get<std::string>("Residual Name"),
          p.get<Teuchos::RCP<PHX::DataLayout>>("Node Vector Data Layout")),
      MyResidual(
          p.get<std::string>("My Residual Name"),
          p.get<Teuchos::RCP<PHX::DataLayout>>("Cell Vector Data Layout"))
{
  this->addEvaluatedField(ExResidual);
  this->addEvaluatedField(MyResidual);

  this->setName("MyElasticityResid" + PHX::print<EvalT>());
}

//**********************************************************************
template <typename EvalT, typename Traits>
void
MyElasticityResid<EvalT, Traits>::postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{

  this->utils.setFieldData(ExResidual, fm);
  this->utils.setFieldData(MyResidual, fm);

  std::vector<PHX::DataLayout::size_type> dims;
  numNodes = dims[1];
  numQPs   = dims[2];
  numDims  = dims[3];
}

//**********************************************************************
template <typename EvalT, typename Traits>
void
MyElasticityResid<EvalT, Traits>::evaluateFields(
    typename Traits::EvalData workset)
{

//fe_values.get_dof_values_from_vector(values, vector);

  for (int cell = 0; cell < workset.numCells; ++cell) {
    for (int node = 0; node < numNodes; ++node) {
      for (int dim = 0; dim < numDims; dim++) {
        ExResidual(cell, node, dim) = 0.0;
      }
      for (int qp = 0; qp < numQPs; ++qp) {
        for (int i = 0; i < numDims; i++) {
          for (int dim = 0; dim < numDims; dim++) {
            ExResidual(cell, node, i) += 1.;
                //Stress(cell, qp, i, dim) * wGradBF(cell, node, qp, dim);
          }
        }
      }
    }
  }

  //Teuchos::RCP<Thyra_Vector> f = workset.f;
  //constraints.distribute_local_to_global(cell_residual, dof_indices, f);
}

//**********************************************************************
}  // namespace LCM
