//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef MYELASTICITYPROBLEM_HPP
#define MYELASTICITYPROBLEM_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Albany_AbstractProblem.hpp"

#include "PHAL_AlbanyTraits.hpp"
#include "PHAL_Dimension.hpp"
#include "PHAL_Workset.hpp"

#include "AAdapt_RC_Manager.hpp"


namespace Albany {

/*!
 * \brief Abstract interface for representing a 2-D finite element
 * problem.
 */
class MyElasticityProblem : public Albany::AbstractProblem
{
 public:
  //! Default constructor
  MyElasticityProblem(
      const Teuchos::RCP<Teuchos::ParameterList>& params_,
      const Teuchos::RCP<ParamLib>&               paramLib_,
      const int                                   numDim_,
      const Teuchos::RCP<AAdapt::rc::Manager>&    rc_mgr);

  //! Destructor
  virtual ~MyElasticityProblem();

  //! Return number of spatial dimensions
  virtual int
  spatialDimension() const
  {
    return numDim;
  }

  //! Get boolean telling code if SDBCs are utilized
  virtual bool
  useSDBCs() const
  {
    return use_sdbcs_;
  }

  //! Build the PDE instantiations, boundary conditions, and initial solution
  virtual void
  buildProblem(
      Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> meshSpecs,
      StateManager&                                            stateMgr);

  // Build evaluators
  virtual Teuchos::Array<Teuchos::RCP<const PHX::FieldTag>>
  buildEvaluators(
      PHX::FieldManager<PHAL::AlbanyTraits>&      fm0,
      const Albany::MeshSpecsStruct&              meshSpecs,
      Albany::StateManager&                       stateMgr,
      Albany::FieldManagerChoice                  fmchoice,
      const Teuchos::RCP<Teuchos::ParameterList>& responseList);

  //! Each problem must generate it's list of valid parameters
  Teuchos::RCP<const Teuchos::ParameterList>
  getValidProblemParameters() const;

  void
  getAllocatedStates(
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<
          Teuchos::RCP<Kokkos::DynRankView<RealType, PHX::Device>>>> oldState_,
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<
          Teuchos::RCP<Kokkos::DynRankView<RealType, PHX::Device>>>> newState_)
      const;

 private:
  //! Private to prohibit copying
  MyElasticityProblem(const MyElasticityProblem&);

  //! Private to prohibit copying
  MyElasticityProblem&
  operator=(const MyElasticityProblem&);

 public:
  //! Main problem setup routine. Not directly called, but indirectly by
  //! following functions
  template <typename EvalT>
  Teuchos::RCP<const PHX::FieldTag>
  constructEvaluators(
      PHX::FieldManager<PHAL::AlbanyTraits>&      fm0,
      const Albany::MeshSpecsStruct&              meshSpecs,
      Albany::StateManager&                       stateMgr,
      Albany::FieldManagerChoice                  fmchoice,
      const Teuchos::RCP<Teuchos::ParameterList>& responseList);

  void
  constructDirichletEvaluators(const Albany::MeshSpecsStruct& meshSpecs);

 protected:
  ///
  /// Boolean marking whether SDBCs are used
  bool use_sdbcs_;

  //! Boundary conditions on source term
  bool haveSource;
  int  numDim;

  //! Compute exact error in displacement solution
  bool computeError;

  std::string                   matModel;
  Teuchos::RCP<Albany::Layouts> dl;

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<
      Teuchos::RCP<Kokkos::DynRankView<RealType, PHX::Device>>>>
      oldState;
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<
      Teuchos::RCP<Kokkos::DynRankView<RealType, PHX::Device>>>>
      newState;

  Teuchos::RCP<AAdapt::rc::Manager> rc_mgr;

  /// Problem PL 
  const Teuchos::RCP<Teuchos::ParameterList> params;

};

}  // namespace Albany

#include "Albany_EvaluatorUtils.hpp"
#include "Albany_ProblemUtils.hpp"
#include "Albany_ResponseUtilities.hpp"
#include "Albany_SolutionAverageResponseFunction.hpp"
#include "Albany_SolutionMaxValueResponseFunction.hpp"
#include "Albany_SolutionTwoNormResponseFunction.hpp"
#include "Albany_Utils.hpp"

#include "DefGrad.hpp"
#include "Density.hpp"
#include "ElasticModulus.hpp"
#include "MyElasticityResid.hpp"
#include "PHAL_SaveStateField.hpp"
#include "PHAL_Source.hpp"
#include "PoissonsRatio.hpp"
#include "Strain.hpp"
#include "Stress.hpp"
//#include "ElasticityDispErrResid.hpp"

#include "Time.hpp"
//#include "CapExplicit.hpp"
//#include "CapImplicit.hpp"

template <typename EvalT>
Teuchos::RCP<const PHX::FieldTag>
Albany::MyElasticityProblem::constructEvaluators(
    PHX::FieldManager<PHAL::AlbanyTraits>&      fm0,
    const Albany::MeshSpecsStruct&              meshSpecs,
    Albany::StateManager&                       stateMgr,
    Albany::FieldManagerChoice                  fieldManagerChoice,
    const Teuchos::RCP<Teuchos::ParameterList>& responseList)
{

  using PHAL::AlbanyTraits;
  using PHX::DataLayout;
  using PHX::MDALayout;
  using std::vector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // get the name of the current element block
  std::string elementBlockName = meshSpecs.ebName;

  const int numNodes    = 4; //intrepidBasis->getCardinality();
  const int worksetSize = meshSpecs.worksetSize;

  const int numDim      = 3;//cubature->getDimension();
  const int numQPts     = 4;// cubature->getNumPoints();
  const int numVertices = 4;// cellType->getNodeCount();
     //const int numVertices = cellType->getVertexCount();

  *out << "Field Dimensions: Workset=" << worksetSize
       << ", Vertices= " << numVertices << ", Nodes= " << numNodes
       << ", QuadPts= " << numQPts << ", Dim= " << numDim << std::endl;

  // Construct standard FEM evaluators with standard field names
  dl = rcp(
      new Albany::Layouts(worksetSize, numVertices, numNodes, numQPts, numDim));
  TEUCHOS_TEST_FOR_EXCEPTION(
      dl->vectorAndGradientLayoutsAreEquivalent == false,
      std::logic_error,
      "Data Layout Usage in Mechanics problems assume vecDim = numDim");
  Albany::EvaluatorUtils<EvalT, PHAL::AlbanyTraits> evalUtils(dl);

  // Displacement Fields

  Teuchos::ArrayRCP<std::string> dof_names(1);
  dof_names[0] = "Displacement";
  Teuchos::ArrayRCP<std::string> resid_names(1);
  resid_names[0] = dof_names[0] + " Residual";

  fm0.template registerEvaluator<EvalT>(
      evalUtils.constructGatherSolutionEvaluator_noTransient(true, dof_names));

  fm0.template registerEvaluator<EvalT>(
      evalUtils.constructScatterResidualEvaluator(true, resid_names));


  // Temporary variable used numerous times below
  Teuchos::RCP<PHX::Evaluator<AlbanyTraits>> ev;

  {  // Displacement Resid
    RCP<ParameterList> p = rcp(new ParameterList("Displacement Resid"));

    // Output
    p->set<std::string>("Residual Name", "Displacement Residual");
    p->set<RCP<DataLayout>>("Node Vector Data Layout", dl->node_vector);
    p->set<RCP<DataLayout>>("Cell Vector Data Layout", dl->cell_vector);

    ev = rcp(new LCM::MyElasticityResid<EvalT, AlbanyTraits>(*p));
    fm0.template registerEvaluator<EvalT>(ev);
  }

  std::cout << "Evaluators constructed." << std::endl;

  //if (Teuchos::nonnull(rc_mgr)) rc_mgr->createEvaluators<EvalT>(fm0, dl);

  //if (fieldManagerChoice == Albany::BUILD_RESID_FM) {
    //PHX::Tag<typename EvalT::ScalarT> res_tag("Scatter", dl->dummy);
    //fm0.requireField<EvalT>(res_tag);

    //return res_tag.clone();
  //}
  return Teuchos::null;
}

#endif  // ALBANY_MYELASTICITYPROBLEM_HPP
