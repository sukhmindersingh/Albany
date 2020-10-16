//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//
#include "MyElasticityProblem.hpp"
#include "Albany_BCUtils.hpp"
#include "Albany_ProblemUtils.hpp"
#include "Albany_Utils.hpp"

Albany::MyElasticityProblem::MyElasticityProblem(
    const Teuchos::RCP<Teuchos::ParameterList>& params_,
    const Teuchos::RCP<ParamLib>&               paramLib_,
    const int                                   numDim_,
    const Teuchos::RCP<AAdapt::rc::Manager>&    rc_mgr_)
    : Albany::AbstractProblem(params_, paramLib_, numDim_),
      params(params_), 
      haveSource(false),
      numDim(numDim_),
      use_sdbcs_(false),
      rc_mgr(rc_mgr_)
{
  std::string& method = params->get("Name", "MyElasticity ");
  *out << "Problem Name = " << method << std::endl;

  haveSource = params->isSublist("Source Functions");

  matModel =
      params->sublist("Material Model").get("Model Name", "LinearElasticity");

  computeError = params->get<bool>("Compute Error", false);
}

Albany::MyElasticityProblem::~MyElasticityProblem() {}

void
Albany::MyElasticityProblem::buildProblem(
    Teuchos::ArrayRCP<Teuchos::RCP<Albany::MeshSpecsStruct>> meshSpecs,
    Albany::StateManager&                                    stateMgr)
{
  /* Construct All Phalanx Evaluators */
  TEUCHOS_TEST_FOR_EXCEPTION(
      meshSpecs.size() != 1,
      std::logic_error,
      "Problem supports one Material Block");

  fm.resize(1);

  fm[0] = Teuchos::rcp(new PHX::FieldManager<PHAL::AlbanyTraits>);
  buildEvaluators(
      *fm[0], *meshSpecs[0], stateMgr, BUILD_RESID_FM, Teuchos::null);

  if (meshSpecs[0]->nsNames.size() >
      0) {  // Build a nodeset evaluator if nodesets are present
    constructDirichletEvaluators(*meshSpecs[0]);
  }
}

Teuchos::Array<Teuchos::RCP<const PHX::FieldTag>>
Albany::MyElasticityProblem::buildEvaluators(
    PHX::FieldManager<PHAL::AlbanyTraits>&      fm0,
    const Albany::MeshSpecsStruct&              meshSpecs,
    Albany::StateManager&                       stateMgr,
    Albany::FieldManagerChoice                  fmchoice,
    const Teuchos::RCP<Teuchos::ParameterList>& responseList)
{
  // Call constructeEvaluators<EvalT>(*rfm[0], *meshSpecs[0], stateMgr);
  // for each EvalT in PHAL::AlbanyTraits::BEvalTypes
  ConstructEvaluatorsOp<MyElasticityProblem> op(
      *this, fm0, meshSpecs, stateMgr, fmchoice, responseList);
  Sacado::mpl::for_each<PHAL::AlbanyTraits::BEvalTypes> fe(op);
  return *op.tags;
}

// Dirichlet BCs
void
Albany::MyElasticityProblem::constructDirichletEvaluators(
    const Albany::MeshSpecsStruct& meshSpecs)
{
  // Construct Dirichlet evaluators for all nodesets and names
  std::vector<std::string> dirichletNames(neq);
  dirichletNames[0] = "X";
  if (neq > 1) dirichletNames[1] = "Y";
  if (neq > 2) dirichletNames[2] = "Z";
  Albany::BCUtils<Albany::DirichletTraits> dirUtils;
  dfm = dirUtils.constructBCEvaluators(
      meshSpecs.nsNames, dirichletNames, this->params, this->paramLib);
  use_sdbcs_  = dirUtils.useSDBCs();
  offsets_    = dirUtils.getOffsets();
  nodeSetIDs_ = dirUtils.getNodeSetIDs();
}


Teuchos::RCP<const Teuchos::ParameterList>
Albany::MyElasticityProblem::getValidProblemParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
      this->getGenericProblemParams("ValidElasticityProblemParams");

  validPL->sublist("Density", false, "");
  validPL->sublist("Elastic Modulus", false, "");
  validPL->sublist("Poissons Ratio", false, "");
  validPL->sublist("Material Model", false, "");

  validPL->set<bool>("Compute Error", false, "");

#ifdef ALBANY_ATO
  // Add additional parameters now for Topological Optimization.
  // ... these in an evaluator rather that in getValidProblemParameters()
  // ... as (apparently) they arn't parsable but used later
  validPL->set<bool>("avgJ", false, "");
  validPL->set<bool>("volavgJ", false, "");
  validPL->set<bool>("weighted_Volume_Averaged_J", false, "");
#endif  // ALBANY_ATO

  if (matModel == "CapExplicit" || matModel == "CapImplicit") {
    validPL->set<double>("A", false, "");
    validPL->set<double>("B", false, "");
    validPL->set<double>("C", false, "");
    validPL->set<double>("theta", false, "");
    validPL->set<double>("R", false, "");
    validPL->set<double>("kappa0", false, "");
    validPL->set<double>("W", false, "");
    validPL->set<double>("D1", false, "");
    validPL->set<double>("D2", false, "");
    validPL->set<double>("calpha", false, "");
    validPL->set<double>("psi", false, "");
    validPL->set<double>("N", false, "");
    validPL->set<double>("L", false, "");
    validPL->set<double>("phi", false, "");
    validPL->set<double>("Q", false, "");
  }

  if (matModel == "GursonSD") {
    validPL->set<double>("f0", false, "");
    validPL->set<double>("Y0", false, "");
    validPL->set<double>("kw", false, "");
    validPL->set<double>("N", false, "");
    validPL->set<double>("q1", false, "");
    validPL->set<double>("q2", false, "");
    validPL->set<double>("q3", false, "");
    validPL->set<double>("eN", false, "");
    validPL->set<double>("sN", false, "");
    validPL->set<double>("fN", false, "");
    validPL->set<double>("fc", false, "");
    validPL->set<double>("ff", false, "");
    validPL->set<double>("flag", false, "");
  }

  return validPL;
}

void
Albany::MyElasticityProblem::getAllocatedStates(
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<
        Teuchos::RCP<Kokkos::DynRankView<RealType, PHX::Device>>>> oldState_,
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<
        Teuchos::RCP<Kokkos::DynRankView<RealType, PHX::Device>>>> newState_)
    const
{
  oldState_ = oldState;
  newState_ = newState;
}
