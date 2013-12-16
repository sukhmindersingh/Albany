//*****************************************************************//
//    Albany 2.0:  Copyright 2012 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_RCP.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Sacado_ParameterRegistration.hpp"

#include "Intrepid_FunctionSpaceTools.hpp"

namespace Aeras {

//**********************************************************************
template<typename EvalT, typename Traits>
ShallowWaterResid<EvalT, Traits>::
ShallowWaterResid(const Teuchos::ParameterList& p,
              const Teuchos::RCP<Albany::Layouts>& dl) :
  wBF      (p.get<std::string> ("Weighted BF Name"), dl->node_qp_scalar),
  wGradBF  (p.get<std::string> ("Weighted Gradient BF Name"),dl->node_qp_gradient),
  U        (p.get<std::string> ("QP Variable Name"), dl->qp_vector),
  Ugrad    (p.get<std::string> ("Gradient QP Variable Name"), dl->qp_vecgradient),
  UDot     (p.get<std::string> ("QP Time Derivative Variable Name"), dl->qp_vector),
  Residual (p.get<std::string> ("Residual Name"), dl->node_vector)
{

  Teuchos::ParameterList* eulerList = p.get<Teuchos::ParameterList*>("Shallow Water Problem");
  gravity = eulerList->get<double>("Gravity", 1.0); //Default: Re=1

  this->addDependentField(U);
  this->addDependentField(Ugrad);
  this->addDependentField(UDot);
  this->addDependentField(wBF);
  this->addDependentField(wGradBF);

  this->addEvaluatedField(Residual);


  this->setName("Aeras::ShallowWaterResid"+PHX::TypeString<EvalT>::value);

  std::vector<PHX::DataLayout::size_type> dims;
  wGradBF.fieldTag().dataLayout().dimensions(dims);
  numNodes = dims[1];
  numQPs   = dims[2];
  numDims  = dims[3];

  U.fieldTag().dataLayout().dimensions(dims);
  vecDim  = dims[2];

//*out << " vecDim = " << vecDim << std::endl;
//*out << " numDims = " << numDims << std::endl;
//*out << " numQPs = " << numQPs << std::endl; 
//*out << " numNodes = " << numNodes << std::endl; 

  // Register Reynolds number as Sacado-ized Parameter
  Teuchos::RCP<ParamLib> paramLib = p.get<Teuchos::RCP<ParamLib> >("Parameter Library");
  new Sacado::ParameterRegistration<EvalT, SPL_Traits>("Gravity", this, paramLib);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ShallowWaterResid<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(U,fm);
  this->utils.setFieldData(Ugrad,fm);
  this->utils.setFieldData(UDot,fm);
  this->utils.setFieldData(wBF,fm);
  this->utils.setFieldData(wGradBF,fm);

  this->utils.setFieldData(Residual,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ShallowWaterResid<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  for (std::size_t i=0; i < Residual.size(); ++i) Residual(i)=0.0;

  // Depth Equation (Eq# 0)
  for (std::size_t cell=0; cell < workset.numCells; ++cell) {
    for (std::size_t qp=0; qp < numQPs; ++qp) {
      for (std::size_t node=0; node < numNodes; ++node) {
          Residual(cell,node,0) += UDot(cell,qp,0)*wBF(cell,node,qp)
                                 - U(cell,qp,0)*U(cell,qp,1)*wGradBF(cell,node,qp,0)
                                 - U(cell,qp,0)*U(cell,qp,2)*wGradBF(cell,node,qp,1);
      }
    }
  }

  // Velocity Equations (Eq# 1,2) -- u_dot = 0 only for now.
  for (std::size_t cell=0; cell < workset.numCells; ++cell) {
    for (std::size_t qp=0; qp < numQPs; ++qp) {
      for (std::size_t node=0; node < numNodes; ++node) {
                  Residual(cell,node,1) += ( UDot(cell,qp,1))*wBF(cell,node,qp);
                  Residual(cell,node,2) += ( UDot(cell,qp,2))*wBF(cell,node,qp);
      }
    }
  }
  // Velocity Equations (Eq# 1,2)
//  for (std::size_t cell=0; cell < workset.numCells; ++cell) {
//    for (std::size_t qp=0; qp < numQPs; ++qp) {
//      for (std::size_t node=0; node < numNodes; ++node) {
//                  Residual(cell,node,1) += ( UDot(cell,qp,1)
//                                             + U(cell,qp,1)*Ugrad(cell,qp,1,0) + U(cell,qp,2)*Ugrad(cell,qp,1,1)
//                                             + gravity*Ugrad(cell,qp,0,0)
//                                            )*wBF(cell,node,qp);
//                  Residual(cell,node,2) += ( UDot(cell,qp,2)
//                                             + U(cell,qp,1)*Ugrad(cell,qp,2,0) + U(cell,qp,2)*Ugrad(cell,qp,2,1)
//                                             + gravity*Ugrad(cell,qp,0,1)
//                                            )*wBF(cell,node,qp);
//      }
//    }
//  }

}

//**********************************************************************
// Provide Access to Parameter for sensitivity/optimization/UQ
template<typename EvalT,typename Traits>
typename ShallowWaterResid<EvalT,Traits>::ScalarT&
ShallowWaterResid<EvalT,Traits>::getValue(const std::string &n)
{
  return gravity;
}
//**********************************************************************
}
